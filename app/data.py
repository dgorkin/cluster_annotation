"""Data layer: config loading, preprocessing, marker tables, and PDF page indices."""
from __future__ import annotations
import json
import os
import re
import subprocess
from pathlib import Path

import pandas as pd
import yaml

from pdf import page_texts

REQUIRED_OUTPUTS = ("gene_markers.tsv.gz", "motif_markers.tsv.gz", "feature_plot_index.tsv")

# ---------------------------------------------------------------- UI preferences
# Project-level, not per-dataset: the render DPI and image width are properties of the screen you
# are looking at, and the last-used config obviously cannot live inside a dataset's own store.
# Kept in .run/ (gitignored) beside the launcher's pidfile.
PREFS_FILE = Path(__file__).resolve().parent.parent / ".run" / "ui_prefs.json"


def load_prefs() -> dict:
    try:
        return json.loads(PREFS_FILE.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def save_prefs(values: dict) -> None:
    """Merge `values` into the stored prefs. Never raises — prefs are a convenience, and failing
    to remember a slider position must not take the app down."""
    try:
        PREFS_FILE.parent.mkdir(parents=True, exist_ok=True)
        merged = {**load_prefs(), **values}
        PREFS_FILE.write_text(json.dumps(merged, indent=2, sort_keys=True))
    except OSError:
        pass


def available_configs() -> list[Path]:
    """Real dataset configs in config/, newest-modified first.

    Excludes *.template.yaml: the shipped template points at placeholder paths, so offering it as
    something to load would just auto-load into an error on a fresh clone. It is there to copy.
    """
    d = Path(__file__).resolve().parent.parent / "config"
    return sorted((p for p in d.glob("*.yaml") if not p.name.endswith(".template.yaml")),
                  key=lambda p: -p.stat().st_mtime)


# ---------------------------------------------------------------- config
def load_config(config_path: str | Path) -> dict:
    config_path = Path(config_path).resolve()
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    # project root = parent of the config dir (default config/mydataset.yaml -> project root)
    root = config_path.parent.parent if config_path.parent.name == "config" else config_path.parent
    cfg["_root"] = str(root)
    cfg["_config_path"] = str(config_path)
    return cfg


def resolve(cfg: dict, path: str) -> Path:
    """Resolve a config path: absolute as-is; else relative to project root, then cwd."""
    p = Path(path)
    if p.is_absolute():
        return p
    root = Path(cfg["_root"])
    for base in (root, Path.cwd()):
        cand = base / p
        if cand.exists():
            return cand
    return root / p


def cache_dir(cfg: dict) -> Path:
    d = resolve(cfg, cfg.get("cache_dir", ".cache")) / cfg["name"]
    d.mkdir(parents=True, exist_ok=True)
    return d


def annotations_db(cfg: dict) -> Path:
    """Path to the per-dataset annotation SQLite DB (user data — NOT under .cache)."""
    d = resolve(cfg, cfg.get("annotations_dir", "annotations"))
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{cfg['name']}.sqlite"


# ---------------------------------------------------------------- preprocessing
def _needs_preprocess(cfg: dict) -> bool:
    cd = cache_dir(cfg)
    outs = [cd / f for f in REQUIRED_OUTPUTS]
    if not all(o.exists() for o in outs):
        return True
    oldest_out = min(os.path.getmtime(o) for o in outs)
    srcs = [resolve(cfg, cfg["markers_rds"]), resolve(cfg, cfg["motif_markers_rds"])]
    newest_src = max(os.path.getmtime(s) for s in srcs)
    return newest_src > oldest_out


def ensure_preprocessed(cfg: dict, force: bool = False) -> Path:
    """Run the R export if outputs are missing/stale. Returns the cache dir."""
    cd = cache_dir(cfg)
    if not force and not _needs_preprocess(cfg):
        return cd
    fp = cfg.get("featureplot", {})
    gene_topn = str(fp.get("gene", {}).get("topn", 20))
    motif_topn = str(fp.get("motif", {}).get("topn", 30))
    gene_sep = str(fp.get("gene", {}).get("separator", "Gapdh"))
    script = Path(cfg["_root"]) / "preprocess" / "export_markers.R"
    cmd = ["Rscript", str(script),
           str(resolve(cfg, cfg["markers_rds"])),
           str(resolve(cfg, cfg["motif_markers_rds"])),
           str(cd), gene_topn, motif_topn, gene_sep]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"export_markers.R failed:\n{res.stderr}\n{res.stdout}")
    return cd


# ---------------------------------------------------------------- marker tables
EFFECT_COL = {"gene": "avg_log2FC", "motif": "avg_diff"}


def load_markers(cfg: dict, kind: str) -> pd.DataFrame:
    cd = cache_dir(cfg)
    fn = "gene_markers.tsv.gz" if kind == "gene" else "motif_markers.tsv.gz"
    return pd.read_csv(cd / fn, sep="\t")


def top_markers(df: pd.DataFrame, cluster: int, n: int, sort_by: str,
                ascending: bool) -> pd.DataFrame:
    sub = df[df["cluster"] == cluster].copy()
    sub = sub.sort_values(sort_by, ascending=ascending, kind="mergesort")
    return sub.head(n).reset_index(drop=True)


# ---------------------------------------------------------------- feature-plot index
def load_feature_index(cfg: dict) -> pd.DataFrame:
    return pd.read_csv(cache_dir(cfg) / "feature_plot_index.tsv", sep="\t")


def clusters(cfg: dict) -> list[int]:
    idx = load_feature_index(cfg)
    return sorted(idx["cluster"].unique().tolist())


def featureplot_pdf(cfg: dict, kind: str, cluster: int) -> Path:
    key = "gene_featureplot_glob" if kind == "gene" else "motif_featureplot_glob"
    return resolve(cfg, cfg[key].replace("{cluster}", str(cluster)))


# ---------------------------------------------------------------- AI insights config
# Models offered in the sidebar selector (current + mythos-class). The configured primary/fallback
# are always added to the list too, so a YAML value outside this set still shows up.
KNOWN_MODELS = [
    "claude-opus-5", "claude-opus-4-8", "claude-fable-5", "claude-opus-4-7",
    "claude-opus-4-6", "claude-sonnet-5", "claude-sonnet-4-6", "claude-haiku-4-5",
]

AI_INSIGHTS_DEFAULTS = {
    "enabled": True,
    "auto_generate_on_load": False,  # never auto-spend on load; generation is button-only + confirmed
    "primary_model": "claude-opus-4-8",     # Fable 5 reroutes on "bio" content for this dataset
    "fallback_model": "claude-opus-4-7",    # distinct fallback if the primary refuses/reroutes
    "structuring_model": "claude-haiku-4-5",  # mechanical reformat into the schema
    "effort": "high",
    "top_n_genes": 30,
    "top_n_motifs": 30,
    "num_references": 5,
    # primer      — one web-grounded reference sheet per dataset, then cheap no-tool per-cluster
    #               calls that cite BY KEY into its verified library. Default: ~6x cheaper, much
    #               faster, and citations are resolved in code so they cannot be fabricated.
    # per_cluster — the original per-cluster web search. Slower and dearer, but does its own
    #               literature work; use it for a cluster the primer cannot place.
    "research_mode": "primer",
    "web_search_max_uses": 4,               # cap searches/cluster; None = unlimited (costlier)
    "prompt_caching": True,                 # cache the re-sent prefix in the pause_turn loop
    "max_workers": 5,
    # Covers thinking + tool use + the written analysis in one call. Too low and a search-heavy
    # cluster truncates before writing anything usable (see the note in config/mydataset.yaml).
    "max_tokens": 32000,
    "api_key_file": "~/.config/cluster_annotation/secrets.env",
}


def ai_insights_cfg(cfg: dict) -> dict:
    """AI-insights settings merged over the defaults."""
    merged = dict(AI_INSIGHTS_DEFAULTS)
    merged.update(cfg.get("ai_insights") or {})
    return merged


def biological_context(cfg: dict) -> str:
    return str(cfg.get("biological_context") or "").strip()


# ---------------------------------------------------------------- other annotations
def other_annotations(cfg: dict) -> list[dict]:
    """Configured extra per-cluster PDFs (one page per cluster, in sorted-cluster order).

    Each entry: {name: str, pdf: <path>}. `{cluster}` in pdf is substituted if present,
    otherwise the PDF is treated as a single file paged by cluster position.
    """
    return cfg.get("other_annotations") or []


def other_pdf(cfg: dict, entry: dict, cluster: int) -> Path:
    return resolve(cfg, str(entry["pdf"]).replace("{cluster}", str(cluster)))


def feature_pages(cfg: dict, kind: str, cluster: int) -> pd.DataFrame:
    """Rows (page, feature) for a cluster's feature-plot PDF; page 1 feature == __UMAP__."""
    idx = load_feature_index(cfg)
    sub = idx[(idx["kind"] == kind) & (idx["cluster"] == cluster)]
    return sub.sort_values("page").reset_index(drop=True)


# ---------------------------------------------------------------- tangram page index
def tangram_index(cfg: dict) -> dict[int, int]:
    """Map cluster -> 1-based page in the Tangram PDF (pages are ordered as strings)."""
    cd = cache_dir(cfg)
    cache = cd / "tangram_index.tsv"
    pdf = resolve(cfg, cfg["tangram_pdf"])
    if cache.exists() and os.path.getmtime(cache) >= os.path.getmtime(pdf):
        t = pd.read_csv(cache, sep="\t")
        return dict(zip(t["cluster"], t["page"]))
    mapping: dict[int, int] = {}
    for i, txt in enumerate(page_texts(pdf), start=1):
        m = re.search(r"Cluster\s+(\d+)", txt)
        if m:
            mapping[int(m.group(1))] = i
    pd.DataFrame({"cluster": list(mapping), "page": list(mapping.values())}).to_csv(
        cache, sep="\t", index=False)
    return mapping
