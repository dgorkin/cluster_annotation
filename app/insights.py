"""AI insights: per-cluster cell-identity hypotheses from Claude.

Two-step pipeline per cluster:
  1. RESEARCH  — claude-fable-5 (primary) with adaptive thinking + the web_search tool, asked to
     identify the likely cell type from the dataset's biological context + the cluster's top marker
     genes/motifs, and to find 3-5 supporting references via web search. Falls back to
     claude-opus-4-8 if Fable 5 refuses or is rerouted by safety classifiers, or on an API error.
  2. STRUCTURE — claude-haiku-4-5 reformats the free-text analysis into the ClusterInsight schema
     (mechanical, low refusal risk; no tools).

Results are cached to .cache/<name>/ai_insights/cluster_<id>.json and invalidated when the marker
inputs change (inputs_hash) or the prompt template changes (PROMPT_VERSION).

The API key is read from a locked-down file (see load_api_key); it is never logged or echoed.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import stat
import tempfile
from concurrent.futures import CancelledError, ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Optional

import pandas as pd
from pydantic import BaseModel, Field

import anthropic

import data as D

# Bump whenever a prompt changes, so cached artifacts produced by an older prompt are not
# silently mixed with newer ones.
# v2 (2026-08-17): references must carry a resolvable identifier (PMID/DOI/URL) or be omitted —
#   v1 returned five unanchored citations on 2 of 3 test clusters.
# v3 (2026-08-17): primer mode, and cite sparingly — with the whole library in the prompt behind
#   convenient keys, v2 over-cited (15 refs on one cluster, including chromVAR/Signac methods
#   papers and oligodendrocyte work on an astrocyte cluster). Fabrication is impossible in primer
#   mode; padding was the remaining failure mode.
PROMPT_VERSION = 3
WEB_SEARCH_TOOL = {"type": "web_search_20260209", "name": "web_search"}
MAX_PAUSE_CONTINUATIONS = 6  # cap on web_search pause_turn loops per call


# ----------------------------------------------------------------- cost / usage accounting
# Anthropic first-party list prices, USD per million tokens (input, output), for the models
# offered in D.KNOWN_MODELS. These drive an ESTIMATE computed from the token counts the API
# reports back — it is not a billing record, and it ignores any negotiated discount. Check the
# figures against https://platform.claude.com/docs/en/pricing if they look off.
MODEL_PRICES = {
    "claude-opus-5":     (5.0, 25.0),
    "claude-opus-4-8":   (5.0, 25.0),
    "claude-opus-4-7":   (5.0, 25.0),
    "claude-opus-4-6":   (5.0, 25.0),
    "claude-fable-5":   (10.0, 50.0),
    "claude-sonnet-5":   (3.0, 15.0),
    "claude-sonnet-4-6": (3.0, 15.0),
    "claude-haiku-4-5":  (1.0,  5.0),
}
CACHE_READ_MULTIPLIER = 0.1     # cached input bills at ~0.1x the model's base input rate
CACHE_WRITE_MULTIPLIER = 1.25   # writing a 5-minute cache entry bills at ~1.25x
WEB_SEARCH_COST = 0.01          # $10 per 1,000 web searches
USAGE_TOKEN_KEYS = ("input_tokens", "cache_read_input_tokens",
                    "cache_creation_input_tokens", "output_tokens", "web_search_requests")


def _usage_of(resp) -> dict:
    """Token + server-tool counts for one API response, as a plain dict. Never raises.

    Read defensively: a refusal or an error can leave fields unset, and we would rather log a
    zero than lose the whole (already paid for) result to an AttributeError.
    """
    u = getattr(resp, "usage", None)
    stu = getattr(u, "server_tool_use", None)

    def _n(obj, field: str) -> int:
        return int(getattr(obj, field, 0) or 0)

    return {
        "model": getattr(resp, "model", None) or "?",
        # stop_reason belongs in the record: a truncated research call used to be indistinguishable
        # from a good one after the fact, which is how a "no marker data provided" annotation got
        # written for a cluster with 1,925 markers.
        "stop_reason": getattr(resp, "stop_reason", None) or "?",
        "input_tokens": _n(u, "input_tokens"),
        "cache_read_input_tokens": _n(u, "cache_read_input_tokens"),
        "cache_creation_input_tokens": _n(u, "cache_creation_input_tokens"),
        "output_tokens": _n(u, "output_tokens"),
        "web_search_requests": _n(stu, "web_search_requests"),
    }


_DATE_SUFFIX = re.compile(r"-\d{8}$")

# Adaptive thinking exists on the 4.6-and-later families only. Haiku 4.5 rejects it outright with
# `400 adaptive thinking is not supported on this model` — which is how the first primer build
# died, after its (paid) web-search research had already succeeded. Anything not listed here gets
# no `thinking` parameter at all, which every model accepts.
ADAPTIVE_THINKING_MODELS = frozenset({
    "claude-opus-5", "claude-opus-4-8", "claude-opus-4-7", "claude-opus-4-6",
    "claude-fable-5", "claude-mythos-5", "claude-sonnet-5", "claude-sonnet-4-6",
})


def supports_adaptive_thinking(model: Optional[str]) -> bool:
    if not model:
        return False
    return (model in ADAPTIVE_THINKING_MODELS
            or _DATE_SUFFIX.sub("", model) in ADAPTIVE_THINKING_MODELS)


def price_for(model: Optional[str]) -> Optional[tuple[float, float]]:
    """List price for a *served* model id, or None if we have no price for it.

    We send an alias ("claude-haiku-4-5") but the API reports the resolved snapshot back in
    `response.model` ("claude-haiku-4-5-20251001"). An exact-match lookup therefore prices those
    calls at $0 — measured live on 2026-08-17, where the Haiku structuring step silently
    contributed nothing to the total. Fall back to the alias with the date suffix stripped.
    """
    if not model:
        return None
    if model in MODEL_PRICES:
        return MODEL_PRICES[model]
    return MODEL_PRICES.get(_DATE_SUFFIX.sub("", model))


def call_cost(rec: dict) -> float:
    """Estimated USD for one call, from its token counts and the served model's list price."""
    price_in, price_out = price_for(rec.get("model")) or (0.0, 0.0)
    return (rec.get("input_tokens", 0) / 1e6 * price_in
            + rec.get("cache_read_input_tokens", 0) / 1e6 * price_in * CACHE_READ_MULTIPLIER
            + rec.get("cache_creation_input_tokens", 0) / 1e6 * price_in * CACHE_WRITE_MULTIPLIER
            + rec.get("output_tokens", 0) / 1e6 * price_out
            + rec.get("web_search_requests", 0) * WEB_SEARCH_COST)


def summarize_usage(calls: list[dict]) -> dict:
    """Roll per-call usage records up into totals plus an estimated cost.

    Cost is summed per call rather than from the totals, because a fallback (or the Haiku
    structuring step) can serve part of one cluster at a different price. `unpriced_models`
    is non-empty when a served model is missing from MODEL_PRICES, so a stale price table
    shows up as a visible gap instead of a silently low number.
    """
    calls = [c for c in calls if c]
    total = {k: sum(c.get(k, 0) for c in calls) for k in USAGE_TOKEN_KEYS}
    total["api_calls"] = len(calls)
    total["models"] = sorted({c.get("model", "?") for c in calls})
    total["unpriced_models"] = sorted({m for m in total["models"] if price_for(m) is None})
    total["est_cost_usd"] = round(sum(call_cost(c) for c in calls), 4)
    total["calls"] = calls  # per-call breakdown, kept for debugging cache behaviour
    return total


def record_cost(rec: Optional[dict]) -> float:
    """Estimated cost recorded on a cached insight/review, or 0.0 if it predates telemetry."""
    if not rec:
        return 0.0
    return float(((rec.get("_meta") or {}).get("usage") or {}).get("est_cost_usd") or 0.0)


# ----------------------------------------------------------------- structured output schema
class Reference(BaseModel):
    citation: str = Field(description="Human-readable citation: authors, title/journal, year.")
    pmid: Optional[str] = Field(default=None, description="PubMed ID if known, else null.")
    url: Optional[str] = Field(default=None, description="A working URL to the source if known.")
    supports: str = Field(description="What claim about this cluster this reference supports.")


class ClusterInsight(BaseModel):
    primary_identity: str = Field(description="Most likely cell type/state for this cluster.")
    alternative_identities: list[str] = Field(
        default_factory=list, description="Other plausible identities, most to least likely.")
    confidence: str = Field(description="One of: high, medium, low.")
    reasoning: str = Field(description="Why this identity, grounded in the markers and context.")
    key_genes: list[str] = Field(default_factory=list, description="Marker genes driving the call.")
    key_motifs: list[str] = Field(default_factory=list, description="Marker motifs/TFs driving it.")
    caveats: str = Field(description="Ambiguities, contaminating signals, or what would confirm it.")
    references: list[Reference] = Field(default_factory=list, description="3-5 supporting sources.")
    revision_note: str = Field(
        default="", description="Empty for an original annotation. On a reannotation: what changed "
        "vs the original annotation and why (or that it stands, and why).")


# ----------------------------------------------------------------- reference primer schema
class PrimerReference(BaseModel):
    key: str = Field(description="Short handle for citing this source, e.g. R1, R2.")
    citation: str = Field(description="Authors, title/journal, year.")
    pmid: Optional[str] = Field(default=None, description="PubMed ID if found, else null.")
    url: Optional[str] = Field(default=None, description="DOI url or working URL if found.")
    covers: str = Field(description="Which cell types or claims this source supports.")


class PrimerCellType(BaseModel):
    name: str = Field(description="Cell type / state expected in this sample.")
    aliases: list[str] = Field(default_factory=list, description="Other names for it.")
    canonical_markers: list[str] = Field(
        default_factory=list, description="Marker genes that identify it.")
    canonical_motifs: list[str] = Field(
        default_factory=list, description="TF motifs/families that identify it.")
    distinguishing_notes: str = Field(
        description="How to tell it apart from the types it is most often confused with.")
    reference_keys: list[str] = Field(
        default_factory=list, description="Keys of supporting references, from the library.")


class ReferencePrimer(BaseModel):
    """A vetted, dataset-level reference sheet, built once with web search.

    Per-cluster calls then run with NO web search against this, and cite BY KEY into its
    reference library — so a citation is resolved from verified data in code rather than
    produced afresh by the model for every cluster.
    """
    expected_cell_types: list[PrimerCellType] = Field(
        default_factory=list,
        description="Cell types/states plausibly present, broad enough to cover the whole sample.")
    references: list[PrimerReference] = Field(
        default_factory=list, description="The citation library, each with a resolvable id.")
    coverage_notes: str = Field(
        description="What this primer covers well, and what it may miss.")


class PrimerClusterInsight(BaseModel):
    """Per-cluster output in primer mode. Cites by key instead of emitting free-text references."""
    primary_identity: str = Field(description="Most likely cell type/state for this cluster.")
    matched_primer_type: str = Field(
        description="Name of the primer cell type this matches, or empty string if none fits.")
    alternative_identities: list[str] = Field(default_factory=list)
    confidence: str = Field(description="One of: high, medium, low.")
    reasoning: str = Field(description="Why this identity, grounded in the markers.")
    key_genes: list[str] = Field(default_factory=list)
    key_motifs: list[str] = Field(default_factory=list)
    caveats: str = Field(description="Ambiguities, contamination, what would confirm it.")
    reference_keys: list[str] = Field(
        default_factory=list,
        description="Keys of primer references supporting this call. Only keys from the library.")
    revision_note: str = Field(default="")


class CohortFlag(BaseModel):
    category: str = Field(description="One of: redundant, inconsistent, ambiguous, other.")
    clusters: list[int] = Field(default_factory=list, description="Cluster id(s) involved.")
    issue: str = Field(description="What looks off across these clusters.")
    suggestion: str = Field(description="Concrete suggestion (merge, re-examine, relabel, …).")
    severity: str = Field(description="One of: high, medium, low.")


class CohortReview(BaseModel):
    """A 'head' review over all per-cluster annotations, checking global coherence."""
    overall: str = Field(description="Overall read on whether the annotation set is coherent.")
    coverage: str = Field(description="Which expected cell types/lineages are well represented.")
    missing_expected: list[str] = Field(
        default_factory=list, description="Cell types expected in this sample but not seen.")
    flags: list[CohortFlag] = Field(
        default_factory=list, description="Specific cross-cluster issues to review.")


# ----------------------------------------------------------------- API key (locked-down file)
def load_api_key(cfg: Optional[dict] = None) -> Optional[str]:
    """Resolve the Anthropic key: env var first, else the configured secrets file. Never logged."""
    key = os.environ.get("ANTHROPIC_API_KEY")
    if key and key.strip():
        return key.strip()
    aic = D.ai_insights_cfg(cfg or {})
    path = Path(os.path.expanduser(aic["api_key_file"]))
    if not path.exists():
        return None
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        k, _, v = line.partition("=")
        if k.strip() == "ANTHROPIC_API_KEY":
            v = v.strip().strip('"').strip("'")
            return v or None
    return None


def api_key_available(cfg: Optional[dict] = None) -> bool:
    return bool(load_api_key(cfg))


_CLIENT: Optional[anthropic.Anthropic] = None

# The SDK default is a 10-minute timeout with our old max_retries=8, i.e. up to ~90 minutes of
# silent retrying on a single wedged request — during a bulk run that reads as a hang with no
# diagnostic. Bound it: a research call that has not finished in 10 minutes is not going to.
REQUEST_TIMEOUT_S = 600.0
MAX_RETRIES = 3


def get_client(cfg: Optional[dict] = None) -> anthropic.Anthropic:
    global _CLIENT
    if _CLIENT is None:
        key = load_api_key(cfg)
        if not key:
            raise RuntimeError("No Anthropic API key found (env ANTHROPIC_API_KEY or secrets file).")
        _CLIENT = anthropic.Anthropic(api_key=key, max_retries=MAX_RETRIES,
                                      timeout=REQUEST_TIMEOUT_S)
    return _CLIENT


# ----------------------------------------------------------------- prompt construction
def _marker_table(df: pd.DataFrame, cluster: int, kind: str, n: int) -> str:
    eff = D.EFFECT_COL[kind]
    top = D.top_markers(df, cluster, n, eff, ascending=False)
    if top.empty:
        return "(none)"
    cols = ["feature", eff, "pct.1", "pct.2", "delta_pct", "p_val_adj"]  # dotted names -> label access
    lines = ["\t".join(cols)]
    for _, row in top[cols].iterrows():
        lines.append("\t".join([
            str(row["feature"]), f"{row[eff]:.3f}", f"{row['pct.1']:.3f}",
            f"{row['pct.2']:.3f}", f"{row['delta_pct']:.3f}", f"{row['p_val_adj']:.2e}",
        ]))
    return "\n".join(lines)


def build_research_prompt(cfg: dict, cluster: int, genes: pd.DataFrame,
                          motifs: pd.DataFrame) -> tuple[str, str]:
    aic = D.ai_insights_cfg(cfg)
    n_ref = aic["num_references"]
    system = (
        "You are an expert developmental biologist and single-cell genomics analyst helping to "
        "annotate clusters from a single-nucleus multiome (RNA + ATAC) experiment. For the given "
        "cluster, infer the most likely cell type/state in the stated biological context, reasoning "
        "from the marker genes and marker motifs (transcription-factor activity). Then use the "
        f"web_search tool to find up to {n_ref} high-quality SUPPORTING references — prefer "
        "primary literature and well-known embryo cell atlases. EVERY reference you list must "
        "carry a resolvable identifier that you actually retrieved: a PubMed ID (preferred), "
        "else a DOI, else a working URL. Never invent, guess or reconstruct an identifier. If you "
        f"cannot anchor a source to a real identifier, leave it out — {n_ref} is a target, not a "
        "quota, and three anchored references are worth far more than five unverifiable ones. For "
        "each one, say which specific claim about this cluster it supports. State your confidence, "
        "plausible alternatives, the key markers behind the call, and caveats. Write a clear, "
        "well-organized analysis in prose; a separate step will convert it to a structured record.")
    bio = D.biological_context(cfg) or "(no biological context provided)"
    user = (
        f"BIOLOGICAL CONTEXT:\n{bio}\n\n"
        f"CLUSTER: {cluster}\n\n"
        f"TOP MARKER GENES (by {D.EFFECT_COL['gene']}; delta_pct = pct.1 - pct.2 = specificity):\n"
        f"{_marker_table(genes, cluster, 'gene', aic['top_n_genes'])}\n\n"
        f"TOP MARKER MOTIFS (by {D.EFFECT_COL['motif']}):\n"
        f"{_marker_table(motifs, cluster, 'motif', aic['top_n_motifs'])}\n\n"
        f"Identify the cell type(s) in this cluster and find {n_ref} supporting references.")
    return system, user


# ----------------------------------------------------------------- step 1: research (+ fallback)
def _text_of(resp) -> str:
    return "\n".join(b.text for b in resp.content if getattr(b, "type", None) == "text").strip()


def _used_web(resp) -> bool:
    return any(getattr(b, "type", "").startswith(("server_tool_use", "web_search"))
               for b in resp.content)


def _web_search_tool(aic: dict) -> dict:
    tool = dict(WEB_SEARCH_TOOL)
    mu = aic.get("web_search_max_uses")
    if mu:
        tool["max_uses"] = int(mu)
    return tool


def _run_research(client, model, system, user, max_tokens, effort, tool, cache: bool = True):
    """One model attempt with the web_search pause_turn loop.

    Returns (resp, text|None, calls) where `calls` holds one _usage_of() record per HTTP request
    made — the pause_turn loop can make several, and each is billed separately.
    """
    messages = [{"role": "user", "content": user}]
    # Cache the prompt prefix. Each pause_turn continuation below re-sends the whole accumulated
    # conversation, so without this every continuation re-bills all the web-search results
    # already sitting in `messages` at the full input rate. The top-level form lets the API place
    # the breakpoint on the last cacheable block, so continuation N reads everything through
    # N-1 at ~0.1x. Caveat: a breakpoint only looks back 20 content blocks, so a single turn that
    # appends more than that (a burst of searches) can still miss — watch
    # usage.cache_read_input_tokens, which is recorded per call below.
    cache_kw = {"cache_control": {"type": "ephemeral"}} if cache else {}
    calls: list[dict] = []
    resp = None
    for _ in range(MAX_PAUSE_CONTINUATIONS):
        # Stream, don't create(). The web_search_20260209 tool runs dynamic filtering (code
        # execution server-side), so a single research turn regularly outruns the HTTP timeout
        # when sent non-streaming — measured 2026-08-17: the same prompt returns in 48s without
        # the tool and times out past 300s with it. Streaming keeps the connection fed with
        # events, and get_final_message() reassembles exactly the Message create() would have
        # returned, so the pause_turn handling and usage accounting below are unchanged.
        with client.messages.stream(
            model=model,
            max_tokens=max_tokens,
            system=system,
            thinking={"type": "adaptive"},
            output_config={"effort": effort},
            tools=[tool],
            messages=messages,
            **cache_kw,
        ) as stream:
            resp = stream.get_final_message()
        calls.append(_usage_of(resp))
        if resp.stop_reason == "refusal":
            return resp, None, calls
        if resp.stop_reason == "pause_turn":
            messages.append({"role": "assistant", "content": resp.content})
            continue
        if resp.stop_reason == "max_tokens":
            # The model ran out of output budget mid-turn — with search uncapped it can spend the
            # whole allowance on thinking and searching and never write the analysis. Downstream,
            # the structuring step honestly reports "no data provided" and we cache a confident-
            # looking but empty annotation for a cluster that has thousands of markers. Treat it
            # as a failed attempt so the caller falls back or reports the cluster as failed.
            return resp, None, calls
        return resp, _text_of(resp), calls
    return resp, _text_of(resp), calls  # ran out of continuations; use what we have


def _research_with_fallback(cfg: dict, system: str, user: str,
                            label: str) -> tuple[str, str, bool, list[dict]]:
    """Run a web-search research call, primary -> fallback.

    Returns (text, model_served, web_used, calls). `calls` accumulates usage across every
    attempt, including one that refused or came back empty — a refusal after partial output is
    still billed, so it belongs in the cost record.
    """
    aic = D.ai_insights_cfg(cfg)
    client = get_client(cfg)
    tool = _web_search_tool(aic)
    cache = bool(aic["prompt_caching"])
    calls: list[dict] = []
    last = "unknown error"
    for model in (aic["primary_model"], aic["fallback_model"]):
        try:
            resp, text, attempt_calls = _run_research(
                client, model, system, user, aic["max_tokens"], aic["effort"], tool, cache)
            calls.extend(attempt_calls)
            if text is None:
                if resp is not None and resp.stop_reason == "max_tokens":
                    last = (f"truncated at max_tokens={aic['max_tokens']} from {resp.model} "
                            "before writing the analysis — raise ai_insights.max_tokens or cap "
                            "web_search_max_uses")
                else:
                    cat = (getattr(resp.stop_details, "category", None)
                           if resp is not None and resp.stop_details else None)
                    last = f"refusal (category={cat}) from {resp.model}"
                continue
            if not text:
                last = f"empty response (stop_reason={resp.stop_reason}) from {resp.model}"
                continue
            return text, resp.model, _used_web(resp), calls
        except anthropic.APIStatusError as e:
            last = f"{type(e).__name__}: {getattr(e, 'message', e)}"
            continue
    raise RuntimeError(f"research failed for {label}: {last}")


def research_cluster(cfg: dict, cluster: int, genes: pd.DataFrame,
                     motifs: pd.DataFrame) -> tuple[str, str, bool, list[dict]]:
    """Returns (analysis_text, model_that_served, web_used, usage_calls)."""
    system, user = build_research_prompt(cfg, cluster, genes, motifs)
    return _research_with_fallback(cfg, system, user, f"cluster {cluster}")


# ----------------------------------------------------------------- step 2: structure
def structure_insight(cfg: dict, analysis_text: str) -> tuple[ClusterInsight, dict]:
    """Reformat the research prose into the schema. Returns (insight, usage_record)."""
    aic = D.ai_insights_cfg(cfg)
    client = get_client(cfg)
    resp = client.messages.parse(
        model=aic["structuring_model"],
        max_tokens=aic["max_tokens"],
        system=("Convert the analysis into the structured schema. Do not invent or add facts; only "
                "restructure what is present. Copy references (citation, PMID, URL) exactly as "
                "given. If a reference carries no identifier, copy it with pmid and url left null "
                "— never fabricate, guess or complete an identifier that the analysis did not "
                "state. Put a DOI in the url field as https://doi.org/<doi>."),
        messages=[{"role": "user", "content": analysis_text}],
        output_format=ClusterInsight,
    )
    if resp.parsed_output is None:
        raise RuntimeError(f"structuring returned no parsed output (stop_reason={resp.stop_reason})")
    return resp.parsed_output, _usage_of(resp)


# ----------------------------------------------------------------- reference primer
def primer_path(cfg: dict) -> Path:
    return insight_dir(cfg) / "_reference_primer.json"


def load_primer(cfg: dict) -> Optional[dict]:
    p = primer_path(cfg)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def build_primer_prompt(cfg: dict) -> tuple[str, str]:
    bio = D.biological_context(cfg) or "(no biological context provided)"
    system = (
        "You are an expert developmental biologist preparing a reference sheet that will be used "
        "to annotate every cluster in one single-nucleus multiome (RNA + ATAC) experiment. Use the "
        "web_search tool thoroughly — this is the ONE chance to consult the literature for this "
        "dataset, and everything downstream depends on it.\n\n"
        "Produce two things:\n"
        "1. An inventory of the cell types and states plausibly present in the stated sample. Be "
        "generous and systematic: include the abundant populations, the minor ones that are easy "
        "to miss, and non-neural/other-lineage populations (vasculature, microglia, meninges, "
        "blood) where relevant. For each, give the marker genes and TF motif families that "
        "identify it, and say how to distinguish it from the types it is most often confused "
        "with — that discriminating detail is the most valuable part of this sheet.\n"
        "2. A citation library covering those types. EVERY reference must carry an identifier you "
        "actually retrieved: a PubMed ID (preferred), else a DOI, else a working URL. Never "
        "invent, guess or reconstruct an identifier; omit a source you cannot anchor. Prefer "
        "primary literature and well-known embryo/organ cell atlases. Give each reference a short "
        "key (R1, R2, ...) — per-cluster annotation will cite these keys, so the library must "
        "stand on its own.\n\n"
        "Write the sheet as clear prose organised under those two headings; a later step converts "
        "it to a structured record.")
    user = (f"BIOLOGICAL CONTEXT:\n{bio}\n\n"
            "Build the reference sheet for this sample: the expected cell types with their "
            "discriminating markers, and the anchored citation library that supports them.")
    return system, user


def build_primer(cfg: dict, force: bool = False) -> dict:
    """Build (or reuse) the dataset's reference primer. One web-grounded call per dataset."""
    if not force:
        existing = load_primer(cfg)
        if existing:
            return existing
    aic = D.ai_insights_cfg(cfg)
    system, user = build_primer_prompt(cfg)
    text, model_used, web_used, calls = _research_with_fallback(cfg, system, user, "primer")
    # Structure the primer with the PRIMARY model, not the cheap structuring model: this is one
    # call per dataset, and every cluster's annotation and citations are read out of the result, so
    # losing a reference or a discriminating note here costs far more than the price difference.
    primer, struct_model, struct_calls = _structured_call(
        cfg, "Convert the reference sheet into the structured schema. Do not invent or add facts, "
             "and do not drop any cell type or reference that the sheet lists. Copy each "
             "reference's identifier exactly as given; if a reference has no identifier, copy it "
             "with pmid and url null rather than completing one.",
        text, ReferencePrimer, [aic["primary_model"], aic["structuring_model"]])
    out = primer.model_dump()
    out["_meta"] = {
        "kind": "reference_primer",
        "model_used": model_used,
        "web_used": web_used,
        "context": D.biological_context(cfg),
        "prompt_version": PROMPT_VERSION,
        "generated_at": pd.Timestamp.now().isoformat(timespec="seconds"),
        "usage": summarize_usage(calls + struct_calls),
    }
    _write_versioned(primer_path(cfg), json.dumps(out, indent=2))
    return out


def primer_stamp(primer: Optional[dict]) -> str:
    """Short fingerprint of a primer, for the per-cluster cache key.

    Rebuilding the primer changes the evidence every cluster was annotated against, so cached
    per-cluster insights must go stale with it.
    """
    if not primer:
        return "none"
    payload = json.dumps({
        "types": sorted(t.get("name", "") for t in primer.get("expected_cell_types") or []),
        "refs": sorted(r.get("key", "") for r in primer.get("references") or []),
        "at": (primer.get("_meta") or {}).get("generated_at"),
    }, sort_keys=True)
    return hashlib.md5(payload.encode()).hexdigest()[:12]


def _primer_digest(primer: dict) -> str:
    """The primer rendered for the prompt. Stable text so it caches as a shared prefix."""
    lines = ["EXPECTED CELL TYPES IN THIS SAMPLE:"]
    for t in primer.get("expected_cell_types") or []:
        lines.append(f"- {t.get('name', '?')}"
                     + (f" (aka {', '.join(t['aliases'])})" if t.get("aliases") else ""))
        if t.get("canonical_markers"):
            lines.append(f"    markers: {', '.join(t['canonical_markers'])}")
        if t.get("canonical_motifs"):
            lines.append(f"    motifs: {', '.join(t['canonical_motifs'])}")
        if t.get("distinguishing_notes"):
            lines.append(f"    distinguish: {t['distinguishing_notes']}")
        if t.get("reference_keys"):
            lines.append(f"    refs: {', '.join(t['reference_keys'])}")
    lines.append("\nREFERENCE LIBRARY (cite these keys; do not cite anything else):")
    for r in primer.get("references") or []:
        ident = f"PMID {r['pmid']}" if r.get("pmid") else (r.get("url") or "no identifier")
        lines.append(f"- [{r.get('key', '?')}] {r.get('citation', '?')}  ({ident})"
                     + (f" — covers: {r['covers']}" if r.get("covers") else ""))
    if primer.get("coverage_notes"):
        lines.append(f"\nCOVERAGE NOTES: {primer['coverage_notes']}")
    return "\n".join(lines)


def build_primer_cluster_prompt(cfg: dict, cluster: int, genes: pd.DataFrame,
                                motifs: pd.DataFrame, primer: dict) -> tuple[str, str]:
    aic = D.ai_insights_cfg(cfg)
    n_ref = aic["num_references"]
    system = (
        "You are annotating one cluster from a single-nucleus multiome (RNA + ATAC) experiment, "
        "using a reference sheet that was prepared for this exact dataset and whose citations have "
        "already been verified. You have no web access — the sheet is your literature.\n\n"
        "Identify the cluster from its marker genes and marker motifs. Prefer a type from the "
        "sheet when one genuinely fits, and name it in matched_primer_type. If nothing on the "
        "sheet fits, say so plainly: leave matched_primer_type empty, give your best independent "
        "call in primary_identity, and set confidence low — a forced fit to a listed type is worse "
        "than an honest mismatch, because it hides a population the sheet missed.\n\n"
        f"Cite support by reference KEY only, from the library in the sheet. Never write a "
        f"citation of your own, and never cite a key that is not listed — an unsupported call with "
        f"no keys is fine and expected.\n\n"
        f"Cite SPARINGLY and specifically: at most {n_ref} keys, and only those that directly "
        f"support THIS cluster's identity. The library covers the whole dataset, so most of it is "
        f"about other cell types — do not cite a reference merely because it is available. In "
        f"particular do not cite methods or tooling papers (sequencing protocols, analysis "
        f"software) as evidence for a cell type. Three precisely relevant references are worth far "
        f"more than a dozen loosely related ones, and a padded list makes the annotation harder to "
        f"check, not easier.\n\n"
        f"State confidence, plausible alternatives, the markers behind the call, and caveats.")
    bio = D.biological_context(cfg) or "(no biological context provided)"
    # The primer goes in `system`, not the user turn. System renders before messages, so a cache
    # breakpoint at the end of it is byte-identical for every cluster and all 34 share one cached
    # prefix. Put the primer after the per-cluster markers instead and each cluster gets its own
    # prefix, which caches nothing across clusters — the whole point of building it once.
    system = f"{system}\n\nBIOLOGICAL CONTEXT:\n{bio}\n\n{_primer_digest(primer)}"
    user = (
        f"----- CLUSTER TO ANNOTATE -----\n"
        f"CLUSTER: {cluster}\n\n"
        f"TOP MARKER GENES (by {D.EFFECT_COL['gene']}; delta_pct = pct.1 - pct.2 = specificity):\n"
        f"{_marker_table(genes, cluster, 'gene', aic['top_n_genes'])}\n\n"
        f"TOP MARKER MOTIFS (by {D.EFFECT_COL['motif']}):\n"
        f"{_marker_table(motifs, cluster, 'motif', aic['top_n_motifs'])}\n\n"
        f"Identify cluster {cluster} and cite the supporting reference keys.")
    return system, user


def _resolve_reference_keys(primer: dict, keys) -> tuple[list[dict], list[str]]:
    """Turn cited keys into full references from the verified library.

    This is the point of primer mode: citations are looked up, not generated, so a cluster cannot
    emit a plausible-looking reference that does not exist. Returns (references, unknown_keys).
    """
    library = {r.get("key"): r for r in (primer.get("references") or []) if r.get("key")}
    refs, unknown = [], []
    for k in dict.fromkeys(keys or []):
        src = library.get(k) or library.get(str(k).strip())
        if not src:
            unknown.append(k)
            continue
        refs.append({"citation": src.get("citation", ""), "pmid": src.get("pmid"),
                     "url": src.get("url"), "supports": src.get("covers", "")})
    return refs, unknown


def generate_one_primer(cfg: dict, cluster: int, genes: pd.DataFrame, motifs: pd.DataFrame,
                        primer: dict) -> dict:
    """Annotate one cluster against the primer: no tools, one call, structured output directly.

    No web search means no server-side tool loop, so this needs neither streaming nor the separate
    Haiku structuring hop — the schema comes straight off the call.
    """
    aic = D.ai_insights_cfg(cfg)
    system, user = build_primer_cluster_prompt(cfg, cluster, genes, motifs, primer)
    insight, model_used, calls = _structured_call(
        cfg, system, user, PrimerClusterInsight,
        [aic["primary_model"], aic["fallback_model"]], cache=bool(aic["prompt_caching"]))
    raw = insight.model_dump()
    refs, unknown = _resolve_reference_keys(primer, raw.pop("reference_keys", []))
    matched = raw.pop("matched_primer_type", "")
    record = {**raw, "references": refs}
    record["_meta"] = {
        "cluster": cluster,
        "mode": "primer",
        "matched_primer_type": matched,
        "model_used": model_used,
        "web_used": False,
        "primer_stamp": primer_stamp(primer),
        "unresolved_reference_keys": unknown,
        "inputs_hash": inputs_hash(cfg, cluster, genes, motifs, primer=primer),
        "prompt_version": PROMPT_VERSION,
        "generated_at": pd.Timestamp.now().isoformat(timespec="seconds"),
        "usage": summarize_usage(calls),
    }
    return record


# ----------------------------------------------------------------- cache + orchestration
def _ensure_writable(p: Path) -> None:
    """Add the owner-write bit to `p` if it lacks one, so we can create/overwrite it.

    The AI-insight cache may be chmod'd read-only (dir 0555, files 0444) to protect paid
    artifacts from stray writes. We always back a file up before overwriting it, so once that
    backup exists it's safe to re-enable writes here rather than fail the (also paid) new run.
    """
    try:
        mode = p.stat().st_mode
        if not mode & stat.S_IWUSR:
            p.chmod(mode | stat.S_IWUSR)
    except OSError:
        pass  # best-effort; the caller's write will surface a real failure


def _write_versioned(path: Path, text: str) -> None:
    """Write `text`, preserving any prior version as a timestamped backup first.

    Unlike the rest of `.cache/`, AI insight artifacts are paid API runs ($-cost), not cheaply
    regenerable. The only way the app overwrites one is a deliberate **Regenerate**; before it does,
    we copy the existing file into `.../_backups/<name>.<ts>.json` so a costly run is never lost to
    an accidental click. Backups accumulate (they're small JSON) and can be pruned by hand.

    The cache may be locked read-only between runs to guard those paid artifacts; we re-enable the
    owner-write bit on the directory and target *after* backing up, so a deliberate regenerate is
    never silently discarded just because the prior run was protected.
    """
    _ensure_writable(path.parent)
    if path.exists():
        bak_dir = path.parent / "_backups"
        bak_dir.mkdir(exist_ok=True)
        ts = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        shutil.copy2(path, bak_dir / f"{path.stem}.{ts}{path.suffix}")
        _ensure_writable(path)  # prior file may be 0444; allow the overwrite now that it's backed up
    path.write_text(text)


def insight_dir(cfg: dict) -> Path:
    d = D.cache_dir(cfg) / "ai_insights"
    d.mkdir(parents=True, exist_ok=True)
    return d


def insight_path(cfg: dict, cluster: int) -> Path:
    return insight_dir(cfg) / f"cluster_{cluster}.json"


def _hash_payload(cfg: dict, cluster: int, genes: pd.DataFrame, motifs: pd.DataFrame,
                  models: Optional[tuple[str, str]] = None,
                  primer: Optional[dict] = None, mode: Optional[str] = None) -> str:
    aic = D.ai_insights_cfg(cfg)
    g = D.top_markers(genes, cluster, aic["top_n_genes"], D.EFFECT_COL["gene"], False)["feature"].tolist()
    m = D.top_markers(motifs, cluster, aic["top_n_motifs"], D.EFFECT_COL["motif"], False)["feature"].tolist()
    payload = {"v": PROMPT_VERSION, "ctx": D.biological_context(cfg), "genes": g, "motifs": m}
    if models is not None:
        payload["primary"], payload["fallback"] = models
    # Mode and primer DO change the answer, so unlike model choice they belong in the key: a
    # primer-mode annotation is not interchangeable with a per-cluster one, and rebuilding the
    # primer changes the evidence every cluster was judged against.
    if mode is not None:
        payload["mode"] = mode
        if mode == "primer":
            payload["primer"] = primer_stamp(primer)
    return hashlib.md5(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def inputs_hash(cfg: dict, cluster: int, genes: pd.DataFrame, motifs: pd.DataFrame,
                primer: Optional[dict] = None) -> str:
    """Cache key over the things that change the ANSWER: markers, context, prompt version.

    Model choice is deliberately NOT part of this. It used to be, which meant flipping the
    sidebar model dropdown marked every cached artifact stale (every cluster) even though
    not one marker had changed — a ~$19 regeneration sitting behind a one-click control. The
    model that actually served is recorded in _meta.model_used, which is where provenance
    belongs. Artifacts written under the old key are converted by migrate_inputs_hashes().
    """
    mode = D.ai_insights_cfg(cfg)["research_mode"]
    # per_cluster keeps the historical payload (no mode key) so migrated artifacts stay current.
    if mode == "per_cluster":
        return _hash_payload(cfg, cluster, genes, motifs)
    return _hash_payload(cfg, cluster, genes, motifs, primer=primer, mode=mode)


def legacy_inputs_hash(cfg: dict, cluster: int, genes: pd.DataFrame, motifs: pd.DataFrame) -> str:
    """The pre-migration key, which folded the configured primary/fallback model names in.

    Reproduced byte-for-byte (sort_keys puts the two model fields back in their old positions)
    so migrate_inputs_hashes() can tell "generated under the old scheme, markers unchanged"
    apart from "genuinely stale".
    """
    aic = D.ai_insights_cfg(cfg)
    return _hash_payload(cfg, cluster, genes, motifs,
                         (aic["primary_model"], aic["fallback_model"]))


def load_insight(cfg: dict, cluster: int) -> Optional[dict]:
    p = insight_path(cfg, cluster)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def _is_current(cfg: dict, cluster: int, genes, motifs, primer: Optional[dict] = None) -> bool:
    rec = load_insight(cfg, cluster)
    return bool(rec) and rec.get("_meta", {}).get("inputs_hash") == inputs_hash(
        cfg, cluster, genes, motifs, primer=primer)


def clusters_needing_insight(cfg: dict, clusters, genes, motifs, force: bool = False,
                             primer: Optional[dict] = None) -> list[int]:
    if force:
        return list(clusters)
    if primer is None and D.ai_insights_cfg(cfg)["research_mode"] == "primer":
        primer = load_primer(cfg)
    return [c for c in clusters if not _is_current(cfg, c, genes, motifs, primer=primer)]


def migrate_inputs_hashes(cfg: dict, clusters, genes: pd.DataFrame, motifs: pd.DataFrame,
                          apply: bool = False) -> dict:
    """Rewrite legacy model-dependent cache keys to the model-independent form, in place.

    An artifact is only rewritten when its stored key provably equals legacy_inputs_hash() for
    this config's CURRENT models — "written under the old scheme, markers unchanged". Anything
    else is left alone, so a genuinely stale artifact stays stale rather than being laundered
    into looking current.

    Run this BEFORE changing primary_model. The legacy key embeds the configured model names,
    so once the config moves the match can no longer be made and the artifacts read as stale.

    Reannotations carry no inputs_hash, so they need no migration. Returns the cluster ids in
    each bucket; pass apply=True to actually write (each rewrite is backed up first).
    """
    out: dict[str, list[int]] = {"migrated": [], "already_current": [], "stale": [], "missing": []}
    for c in clusters:
        rec = load_insight(cfg, c)
        if rec is None:
            out["missing"].append(c)
            continue
        meta = rec.get("_meta") or {}
        stored = meta.get("inputs_hash")
        if stored == inputs_hash(cfg, c, genes, motifs):
            out["already_current"].append(c)
        elif stored == legacy_inputs_hash(cfg, c, genes, motifs):
            out["migrated"].append(c)
            if apply:
                meta["inputs_hash"] = inputs_hash(cfg, c, genes, motifs)
                meta["inputs_hash_migrated_from"] = stored
                rec["_meta"] = meta
                _write_versioned(insight_path(cfg, c), json.dumps(rec, indent=2))
        else:
            out["stale"].append(c)
    return out


def generate_one(cfg: dict, cluster: int, genes: pd.DataFrame, motifs: pd.DataFrame,
                 force: bool = False, primer: Optional[dict] = None) -> dict:
    """Annotate one cluster in whichever mode the config selects.

    primer mode  — one no-tool call against the dataset's vetted reference sheet. Cheap, fast,
                   citations resolved from the library in code.
    per_cluster  — the original: streamed web-search research, then a Haiku structuring pass.
                   Slower and dearer, but does its own literature work; use it for a cluster the
                   primer cannot place.
    """
    aic = D.ai_insights_cfg(cfg)
    mode = aic["research_mode"]
    if mode == "primer" and primer is None:
        primer = load_primer(cfg)
        if primer is None:
            raise RuntimeError(
                "research_mode is 'primer' but no reference primer exists for this dataset — "
                "build it first (Build reference primer in the AI insights tab, or "
                "scripts/generate_all.py, which builds it automatically).")
    if not force and _is_current(cfg, cluster, genes, motifs, primer=primer):
        return load_insight(cfg, cluster)

    if mode == "primer":
        record = generate_one_primer(cfg, cluster, genes, motifs, primer)
    else:
        text, model_used, web_used, calls = research_cluster(cfg, cluster, genes, motifs)
        insight, struct_call = structure_insight(cfg, text)
        record = insight.model_dump()
        record["_meta"] = {
            "cluster": cluster,
            "mode": "per_cluster",
            "model_used": model_used,
            "web_used": web_used,
            "inputs_hash": inputs_hash(cfg, cluster, genes, motifs),
            "prompt_version": PROMPT_VERSION,
            "generated_at": pd.Timestamp.now().isoformat(timespec="seconds"),
            "usage": summarize_usage(calls + [struct_call]),
        }
    payload = json.dumps(record, indent=2)
    try:
        _write_versioned(insight_path(cfg, cluster), payload)
    except OSError as e:
        # The API call already cost money; never throw the result away because the cache write
        # failed. Stash it somewhere guaranteed-writable and point the caller at it to recover.
        salvage = Path(tempfile.gettempdir()) / (
            f"cluster_annotation_unsaved_{cfg['name']}_cluster_{cluster}_"
            f"{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.json")
        try:
            salvage.write_text(payload)
            hint = f"salvaged paid result to {salvage}"
        except OSError:
            hint = "could not salvage result to a temp file either"
        raise RuntimeError(f"failed to write insight for cluster {cluster}: {e} ({hint})") from e
    return record


def generate_all(cfg: dict, clusters, genes: pd.DataFrame, motifs: pd.DataFrame,
                 progress_cb: Optional[Callable[[int, Optional[str]], None]] = None,
                 force: bool = True,
                 should_stop: Optional[Callable[[], bool]] = None) -> list[tuple[int, str]]:
    """Generate insights for `clusters` concurrently. Returns a list of (cluster, error) failures.

    progress_cb(cluster, error) is invoked from the calling thread as each cluster completes
    (error is None on success) — safe to update Streamlit widgets from it.

    should_stop() is polled as clusters complete: when it turns true, clusters that have not
    started yet are cancelled. Calls already in flight are paid for and are left to finish, so
    stopping bounds future spend rather than pretending to undo the current one.
    """
    aic = D.ai_insights_cfg(cfg)
    # Load the primer once and pass it down: 34 workers each re-reading and re-hashing it is
    # wasted work, and it guarantees every cluster in a run is judged against the same sheet.
    primer = load_primer(cfg) if aic["research_mode"] == "primer" else None
    errors: list[tuple[int, str]] = []
    with ThreadPoolExecutor(max_workers=aic["max_workers"]) as ex:
        futs = {ex.submit(generate_one, cfg, c, genes, motifs, force, primer): c
                for c in clusters}
        for fut in as_completed(futs):
            c = futs[fut]
            try:
                fut.result()
                err = None
            except CancelledError:
                continue          # never started; not a failure worth reporting
            except Exception as e:  # one cluster failing must not kill the batch
                err = str(e)
                errors.append((c, err))
            if progress_cb:
                progress_cb(c, err)
            if should_stop and should_stop():
                for pending in futs:
                    pending.cancel()   # no-op on anything already running
    return errors


# ----------------------------------------------------------------- "head" cohort review
def _structured_call(cfg: dict, system: str, user: str, schema, models, cache: bool = False):
    """Single structured (no-tool) call with primary->fallback on refusal/error.

    Returns (obj, model, calls) — `calls` carries usage for every attempt, refusals included.
    `cache=True` puts a breakpoint at the end of the system prompt, which is what makes the
    primer a shared prefix across every cluster in a run.
    """
    client = get_client(cfg)
    aic = D.ai_insights_cfg(cfg)
    system_param = ([{"type": "text", "text": system, "cache_control": {"type": "ephemeral"}}]
                    if cache else system)
    calls: list[dict] = []
    last = "unknown error"
    for model in models:
        try:
            # Only send `thinking` where the model supports it: this helper is used with the
            # primary model (adaptive) and with the structuring model (Haiku 4.5, which 400s on it).
            thinking_kw = ({"thinking": {"type": "adaptive"}}
                           if supports_adaptive_thinking(model) else {})
            resp = client.messages.parse(
                model=model,
                max_tokens=aic["max_tokens"],
                system=system_param,
                messages=[{"role": "user", "content": user}],
                output_format=schema,
                **thinking_kw,
            )
            calls.append(_usage_of(resp))
            if resp.stop_reason == "refusal":
                cat = getattr(resp.stop_details, "category", None) if resp.stop_details else None
                last = f"refusal (category={cat}) from {resp.model}"
                continue
            if resp.parsed_output is None:
                last = f"no parsed output (stop_reason={resp.stop_reason}) from {resp.model}"
                continue
            return resp.parsed_output, resp.model, calls
        except anthropic.APIStatusError as e:
            last = f"{type(e).__name__}: {getattr(e, 'message', e)}"
            continue
    raise RuntimeError(f"structured call failed: {last}")


def cohort_review_path(cfg: dict) -> Path:
    return insight_dir(cfg) / "_cohort_review.json"


def load_cohort_review(cfg: dict) -> Optional[dict]:
    p = cohort_review_path(cfg)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def _cohort_summary_line(c: int, rec: dict) -> str:
    alts = ", ".join(rec.get("alternative_identities") or [])
    genes = ", ".join((rec.get("key_genes") or [])[:12])
    motifs = ", ".join((rec.get("key_motifs") or [])[:8])
    return (f"Cluster {c}: {rec.get('primary_identity', '?')} "
            f"[confidence {rec.get('confidence', '?')}]"
            + (f"; alts: {alts}" if alts else "")
            + (f"; genes: {genes}" if genes else "")
            + (f"; motifs: {motifs}" if motifs else ""))


def build_cohort_prompt(cfg: dict, present: dict[int, dict]) -> tuple[str, str]:
    system = (
        "You are a senior single-cell biologist doing a QC review of an entire set of cluster cell-"
        "type annotations from one experiment. Judge them AS A WHOLE: flag clusters that are likely "
        "the same type and over-split (redundant), annotations that are mutually inconsistent or "
        "biologically implausible together, and ambiguous calls that need another look. Assess "
        "whether the cell types/lineages expected in the stated sample are represented, and list any "
        "conspicuously missing ones. Be specific and cite cluster ids. Do not invent markers; reason "
        "from the provided annotations and the biological context.")
    bio = D.biological_context(cfg) or "(no biological context provided)"
    lines = "\n".join(_cohort_summary_line(c, present[c]) for c in sorted(present))
    user = (f"BIOLOGICAL CONTEXT:\n{bio}\n\n"
            f"PER-CLUSTER ANNOTATIONS ({len(present)} clusters):\n{lines}\n\n"
            "Review the set for global coherence, coverage, and cross-cluster issues.")
    return system, user


def cohort_review(cfg: dict, clusters, force: bool = False) -> dict:
    """Run the head review over all cached per-cluster insights. Caches a single JSON."""
    present = {c: load_insight(cfg, c) for c in clusters}
    present = {c: r for c, r in present.items() if r}
    if not present:
        raise RuntimeError("No per-cluster insights yet — generate clusters first.")
    aic = D.ai_insights_cfg(cfg)
    system, user = build_cohort_prompt(cfg, present)
    review, model_used, calls = _structured_call(
        cfg, system, user, CohortReview, [aic["primary_model"], aic["fallback_model"]])
    out = review.model_dump()
    out["_meta"] = {
        "model_used": model_used,
        "covered_clusters": sorted(present),
        "n_clusters": len(present),
        "generated_at": pd.Timestamp.now().isoformat(timespec="seconds"),
        "usage": summarize_usage(calls),
    }
    _write_versioned(cohort_review_path(cfg), json.dumps(out, indent=2))
    return out


def cluster_flags(review: Optional[dict], cluster: int) -> list[dict]:
    """Cohort-review flags that reference a given cluster."""
    if not review:
        return []
    return [f for f in (review.get("flags") or []) if cluster in (f.get("clusters") or [])]


def flagged_clusters(review: Optional[dict]) -> list[int]:
    if not review:
        return []
    out: set[int] = set()
    for f in review.get("flags") or []:
        out.update(f.get("clusters") or [])
    return sorted(out)


# ----------------------------------------------------------------- reannotation (keeps original)
def reannotation_path(cfg: dict, cluster: int) -> Path:
    return insight_dir(cfg) / f"cluster_{cluster}.reannotation.json"


def load_reannotation(cfg: dict, cluster: int) -> Optional[dict]:
    p = reannotation_path(cfg, cluster)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text())
    except (json.JSONDecodeError, OSError):
        return None


def reannotation_is_stale(cfg: dict, cluster: int) -> bool:
    """True when a reannotation describes an original that has since been regenerated.

    A reannotation is written against one specific original annotation and one specific cohort
    critique. Regenerate the original and the reannotation stays on disk describing something that
    no longer exists — still readable, still plausible, and now wrong. It records no inputs_hash of
    its own, so compare write times: an original written after its reannotation means the
    reannotation is history.
    """
    reann = load_reannotation(cfg, cluster)
    if not reann:
        return False
    orig = load_insight(cfg, cluster)
    if not orig:
        return True
    r_at = (reann.get("_meta") or {}).get("generated_at") or ""
    o_at = (orig.get("_meta") or {}).get("generated_at") or ""
    return bool(r_at) and bool(o_at) and o_at > r_at   # ISO-8601 sorts lexicographically


def build_reannotation_prompt(cfg: dict, cluster: int, genes: pd.DataFrame, motifs: pd.DataFrame,
                              original: dict, flags: list[dict]) -> tuple[str, str]:
    aic = D.ai_insights_cfg(cfg)
    n_ref = aic["num_references"]
    system = (
        "You are critically REANNOTATING one cluster in light of a cohort-level QC review of the "
        "whole annotation set. Reconsider the cell identity: if the original call stands, strengthen "
        "it; if the critique shows it is mislabeled, over-split, or ambiguous, revise it. Reason from "
        "the markers and use the web_search tool to find current supporting references (up to "
        f"{n_ref}). Every reference must carry a resolvable identifier you actually retrieved — a "
        "PubMed ID (preferred), else a DOI, else a working URL; never invent one, and omit any "
        "source you cannot anchor. In a clearly labeled 'what changed' section, state "
        "plainly how your reannotation differs from the original annotation and why — or that it is "
        "unchanged and why it stands. Write a clear prose analysis; a later step structures it.")
    bio = D.biological_context(cfg) or "(no biological context provided)"
    orig = (f"primary_identity: {original.get('primary_identity', '?')}\n"
            f"confidence: {original.get('confidence', '?')}\n"
            f"alternatives: {', '.join(original.get('alternative_identities') or []) or '—'}\n"
            f"key_genes: {', '.join(original.get('key_genes') or []) or '—'}\n"
            f"key_motifs: {', '.join(original.get('key_motifs') or []) or '—'}\n"
            f"reasoning: {original.get('reasoning', '')}")
    if flags:
        crit = "\n".join(
            f"- [{f.get('severity', '?')}/{f.get('category', '?')}] {f.get('issue', '')} "
            f"(suggestion: {f.get('suggestion', '')})" for f in flags)
    else:
        crit = "(no specific flag references this cluster; re-examine it for overall coherence)"
    user = (
        f"BIOLOGICAL CONTEXT:\n{bio}\n\n"
        f"CLUSTER: {cluster}\n\n"
        f"TOP MARKER GENES (by {D.EFFECT_COL['gene']}):\n"
        f"{_marker_table(genes, cluster, 'gene', aic['top_n_genes'])}\n\n"
        f"TOP MARKER MOTIFS (by {D.EFFECT_COL['motif']}):\n"
        f"{_marker_table(motifs, cluster, 'motif', aic['top_n_motifs'])}\n\n"
        f"ORIGINAL ANNOTATION:\n{orig}\n\n"
        f"COHORT REVIEW CRITIQUE FOR THIS CLUSTER:\n{crit}\n\n"
        "Reannotate this cluster, addressing the critique, and explain what changed vs the original.")
    return system, user


def generate_reannotation(cfg: dict, cluster: int, genes: pd.DataFrame, motifs: pd.DataFrame,
                          review: dict) -> dict:
    """Reannotate one cluster using the cohort critique. Writes a SEPARATE file; original untouched."""
    original = load_insight(cfg, cluster)
    if original is None:
        raise RuntimeError(f"no original insight for cluster {cluster} — generate it first")
    flags = cluster_flags(review, cluster)
    system, user = build_reannotation_prompt(cfg, cluster, genes, motifs, original, flags)
    text, model_used, web_used, calls = _research_with_fallback(
        cfg, system, user, f"reannotation {cluster}")
    insight, struct_call = structure_insight(cfg, text)
    rec = insight.model_dump()
    rec["_meta"] = {
        "cluster": cluster,
        "kind": "reannotation",
        "model_used": model_used,
        "web_used": web_used,
        "addressed_flags": flags,
        "cohort_generated_at": (review.get("_meta") or {}).get("generated_at"),
        "generated_at": pd.Timestamp.now().isoformat(timespec="seconds"),
        "usage": summarize_usage(calls + [struct_call]),
    }
    _write_versioned(reannotation_path(cfg, cluster), json.dumps(rec, indent=2))
    return rec


def dataset_cost(cfg: dict, clusters) -> dict:
    """Estimated spend already recorded on disk for this dataset, by artifact kind.

    Reads the cached JSONs rather than tracking a running total, so it stays correct across
    restarts and reports whatever is actually on disk. Artifacts generated before usage
    telemetry existed contribute 0 and are counted in `untracked` so the total reads as a
    lower bound instead of looking complete.
    """
    out = {"originals": 0.0, "reannotations": 0.0, "cohort_review": 0.0,
           "n_priced": 0, "untracked": 0}
    for c in clusters:
        for key, rec in (("originals", load_insight(cfg, c)),
                         ("reannotations", load_reannotation(cfg, c))):
            if rec is None:
                continue
            cost = record_cost(rec)
            out[key] += cost
            out["n_priced" if cost else "untracked"] += 1
    review = load_cohort_review(cfg)
    if review is not None:
        cost = record_cost(review)
        out["cohort_review"] = cost
        out["n_priced" if cost else "untracked"] += 1
    out["total"] = round(out["originals"] + out["reannotations"] + out["cohort_review"], 4)
    for k in ("originals", "reannotations", "cohort_review"):
        out[k] = round(out[k], 4)
    return out


def format_dataset_cost(cfg: dict, clusters) -> str:
    """One-line summary of dataset_cost(), for the runner scripts and the UI."""
    t = dataset_cost(cfg, clusters)
    line = (f"recorded spend for {cfg['name']}: originals ~${t['originals']:.2f} + "
            f"reannotations ~${t['reannotations']:.2f} + cohort ~${t['cohort_review']:.2f} "
            f"= ~${t['total']:.2f} over {t['n_priced']} artifact(s)")
    if t["untracked"]:
        line += (f"; {t['untracked']} artifact(s) predate cost telemetry and count as $0, "
                 "so this is a lower bound")
    return line


def reannotate_flagged(cfg: dict, clusters, genes: pd.DataFrame, motifs: pd.DataFrame, review: dict,
                       progress_cb: Optional[Callable[[int, Optional[str]], None]] = None,
                       should_stop: Optional[Callable[[], bool]] = None
                       ) -> list[tuple[int, str]]:
    """Reannotate each cluster in `clusters` concurrently. Returns (cluster, error) failures.

    `should_stop` behaves as in generate_all: not-yet-started clusters are cancelled.
    """
    aic = D.ai_insights_cfg(cfg)
    errors: list[tuple[int, str]] = []
    with ThreadPoolExecutor(max_workers=aic["max_workers"]) as ex:
        futs = {ex.submit(generate_reannotation, cfg, c, genes, motifs, review): c for c in clusters}
        for fut in as_completed(futs):
            c = futs[fut]
            try:
                fut.result()
                err = None
            except CancelledError:
                continue
            except Exception as e:
                err = str(e)
                errors.append((c, err))
            if progress_cb:
                progress_cb(c, err)
            if should_stop and should_stop():
                for pending in futs:
                    pending.cancel()
    return errors
