"""Cluster Annotation Browser — v1 (browse-only).

Run:  streamlit run app/app.py --server.port 8501 --server.address 127.0.0.1
Then: ssh -L 8501:localhost:8501 <server>   and open http://localhost:8501
"""
from __future__ import annotations
import json
import os
import sys
from pathlib import Path

import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).parent))  # allow `import data`/`import pdf`
import data as D
import pdf as P
import insights as I
import store as S
import export as E
import jobs as J

st.set_page_config(page_title="Cluster Annotation Browser", layout="wide")

EFFECT = {"gene": "avg_log2FC", "motif": "avg_diff"}


def inject_js(script: str) -> None:
    """Run a <script> in the app document.

    Uses st.html rather than the deprecated st.components.v1.html (whose stated removal date has
    already passed, so it will vanish on some future Streamlit upgrade). st.html injects inline
    instead of inside a zero-height iframe, so there is no wasted layout box. The scripts below
    still address the document via `window.parent` — correct either way, because at top level
    window.parent is window itself.

    **The IIFE wrapper is load-bearing.** Injecting inline means every script shares one global
    lexical scope, so two scripts that both declare `const doc` — or the same script re-injected on
    the next rerun — is a *parse-time* SyntaxError that silently kills the whole script, listeners
    and all. That is exactly what happened when this moved off the iframe (which gave each script
    its own document): the keyboard shortcuts never registered, with nothing in the UI to say so.
    Wrapping here rather than in each caller means a new script cannot reintroduce the bug.
    """
    st.html("<script>(function(){\n" + script + "\n})();</script>",
            unsafe_allow_javascript=True)


def block_clear_cache_shortcut():
    """Neutralize Streamlit's built-in 'C' = Clear cache hotkey.

    Clearing the cache would force costly AI insights to regenerate, so we swallow the key in the
    capture phase (before Streamlit's own handler) whenever focus isn't in a text field. A real
    Ctrl/Cmd+C copy (modifier + a non-empty selection) is left alone.
    """
    inject_js("""
const doc = window.parent.document, w = window.parent;
if (w.__noClearCache) doc.removeEventListener('keydown', w.__noClearCache, true);
w.__noClearCache = function(e){
  if (e.key !== 'c' && e.key !== 'C') return;
  const t = e.target, tag = t && t.tagName ? t.tagName.toUpperCase() : '';
  if (tag === 'INPUT' || tag === 'TEXTAREA' || (t && t.isContentEditable)) return;
  const hasSel = w.getSelection && String(w.getSelection()).length > 0;
  if ((e.ctrlKey || e.metaKey) && hasSel) return;   // preserve genuine copy
  e.stopImmediatePropagation();
  e.preventDefault();
};
doc.addEventListener('keydown', w.__noClearCache, true);  // capture: runs before Streamlit
""")


block_clear_cache_shortcut()


# --------------------------------------------------------------- loading
def load_dataset(config_path: str, force: bool = False):
    cfg = D.load_config(config_path)
    with st.spinner("Preprocessing markers (Rscript)…" if force else "Loading dataset…"):
        D.ensure_preprocessed(cfg, force=force)
    st.session_state.cfg = cfg
    st.session_state.genes = D.load_markers(cfg, "gene")
    st.session_state.motifs = D.load_markers(cfg, "motif")
    st.session_state.feat_idx = D.load_feature_index(cfg)
    st.session_state.tangram = D.tangram_index(cfg)
    st.session_state.clusters = D.clusters(cfg)
    load_cell_counts(cfg, force=force)
    aic = D.ai_insights_cfg(cfg)
    if aic["enabled"] and aic["auto_generate_on_load"] and I.api_key_available(cfg):
        start_ai_generation(cfg, force=False)  # background job; status in the AI insights tab


def load_cell_counts(cfg, force: bool = False):
    """Put ncells / %cells per cluster in session state, computing them first if needed.

    Part of loading rather than a button, because the numbers are a property of the dataset: once
    computed they are cached until the object changes. The read of the object is the expensive bit
    (tens of seconds and several GB for a whole-embryo multiome), hence the explicit spinner and
    the never-raise contract — a dataset whose object is missing or ambiguous must still open.
    """
    # `auto: false` opts out of ever reading the object on load — the CLI (./run_app.sh counts)
    # and the Recount button still work, they pass force.
    if force or (D.cell_counts_cfg(cfg)["auto"] and D.cell_counts_needed(cfg)):
        try:
            with st.spinner("Counting cells per cluster from the Seurat object — reads a large "
                            "file, usually under a minute…"):
                D.ensure_cell_counts(cfg, force=force, expected=st.session_state.clusters)
        except Exception as exc:  # noqa: BLE001 - counts are a convenience, not a prerequisite
            st.warning(f"Could not count cells per cluster: {exc}")
    st.session_state.cell_counts = D.load_cell_counts(cfg)


# --------------------------------------------------------------- paid AI work, as durable jobs
# Every paid action runs through jobs.start, on a background thread with its state on disk. The
# alternative (doing the work inline, inside st.spinner) tied a minutes-long, money-spending call
# to one script run: switching section or reloading the page hid it completely, with no way to
# tell whether it was still going.
AI_JOBS = ("primer", "generate", "cohort", "reannotate")


def start_primer_build(cfg) -> bool:
    def work(_prog, _stop):
        got = I.build_primer(cfg, force=True)
        return (f"Primer built: {len(got.get('expected_cell_types') or [])} cell types, "
                f"{len(got.get('references') or [])} references (~${I.record_cost(got):.2f}).")

    return J.start(cfg, "primer", "Build reference primer", work)


def start_ai_generation(cfg, force: bool, only: list[int] | None = None) -> bool:
    """Annotate clusters that need it (or exactly `only`), as a background job."""
    clusters = st.session_state.clusters
    genes, motifs = st.session_state.genes, st.session_state.motifs
    if only is None:
        try:
            todo = I.clusters_needing_insight(cfg, clusters, genes, motifs, force=force)
        except Exception as e:  # noqa: BLE001
            st.warning(f"AI insights: could not determine work ({e}).")
            return False
    else:
        todo = list(only)
    if not todo:
        st.info("Every cluster already has a current annotation — nothing to generate.")
        return False

    label = (f"Annotate cluster {todo[0]}" if len(todo) == 1
             else f"Annotate {len(todo)} clusters")

    def work(prog, stop):
        state = {"done": 0}

        def cb(cluster, err):
            state["done"] += 1
            prog.update(done=state["done"],
                        message=f"cluster {cluster}" + (f" failed: {err}" if err else " done"))

        errors = I.generate_all(cfg, todo, genes, motifs, progress_cb=cb, force=True,
                               should_stop=stop)
        spent = sum(I.record_cost(I.load_insight(cfg, c)) for c in todo)
        done = state["done"] - len(errors)
        msg = f"Annotated {done} of {len(todo)} cluster(s) · ~${spent:.2f} estimated."
        if errors:
            msg += (f" {len(errors)} failed ({', '.join(str(c) for c, _ in errors[:6])}"
                    f"{'…' if len(errors) > 6 else ''}) — retry from this tab.")
        return msg

    return J.start(cfg, "generate", label, work, total=len(todo), stoppable=len(todo) > 1)


def start_reannotation(cfg, targets, review) -> bool:
    if not targets:
        return False
    genes, motifs = st.session_state.genes, st.session_state.motifs
    targets = list(targets)
    label = (f"Reannotate cluster {targets[0]}" if len(targets) == 1
             else f"Reannotate {len(targets)} flagged clusters")

    def work(prog, stop):
        state = {"done": 0}

        def cb(cluster, err):
            state["done"] += 1
            prog.update(done=state["done"],
                        message=f"cluster {cluster}" + (f" failed: {err}" if err else " done"))

        errors = I.reannotate_flagged(cfg, targets, genes, motifs, review, progress_cb=cb,
                                     should_stop=stop)
        spent = sum(I.record_cost(I.load_reannotation(cfg, c)) for c in targets)
        msg = (f"Reannotated {state['done'] - len(errors)} of {len(targets)} cluster(s) "
               f"· ~${spent:.2f} estimated.")
        if errors:
            msg += f" {len(errors)} failed."
        return msg

    return J.start(cfg, "reannotate", label, work, total=len(targets),
                   stoppable=len(targets) > 1)


def start_cohort_review(cfg, clusters) -> bool:
    def work(_prog, _stop):
        rev = I.cohort_review(cfg, clusters, force=True)
        flags = rev.get("flags") or []
        return (f"Reviewed {len(clusters)} clusters: {len(flags)} flag(s) over "
                f"{len(I.flagged_clusters(rev))} cluster(s) (~${I.record_cost(rev):.2f}).")

    return J.start(cfg, "cohort", "Cohort review", work)


@st.fragment(run_every=3)
def sidebar_job_indicator(cfg):
    """One line in the sidebar while paid work is in flight, visible from every section."""
    kind = J.any_running(cfg, AI_JOBS)
    if not kind:
        return
    rec = J.read(cfg, kind) or {}
    done, total = rec.get("done") or 0, rec.get("total")
    counts = f" · {done}/{total}" if total else ""
    # Plain st.* — a fragment cannot address st.sidebar directly; the caller supplies the sidebar
    # as the active container instead.
    st.warning(f"⏳ **{rec.get('label', kind)}** running{counts}\n\n"
               "Details in **AI insights**. Safe to keep browsing.")


@st.fragment(run_every=3)
def job_banner(cfg):
    """Status for whatever paid work is in flight (or last finished), refreshed every 3s.

    In a fragment so the poll redraws just this box rather than re-running the whole page — which
    would re-render the Tangram PDF every three seconds. A finished job triggers one full rerun so
    the rest of the page (the primer summary, the overview counts) catches up.
    """
    for kind in AI_JOBS:
        rec = J.read(cfg, kind)
        if not rec:
            continue
        state, label = rec["state"], rec.get("label", kind)
        if state == J.RUNNING:
            done, total = rec.get("done") or 0, rec.get("total")
            head = f"⏳ **{label}** — running since {rec.get('started_at', '?')}"
            if total:
                st.progress(min(done / total, 1.0), text=f"{head} · {done} / {total}")
            else:
                st.info(f"{head}. This makes one long call; there is nothing to show until it "
                        "returns.")
            st.caption(f"{rec.get('message', '')} · {rec.get('primary_model', '?')} · "
                       "keeps running if you switch section or reload the page.")
            if rec.get("stoppable") and not rec.get("stopping"):
                if st.button("■ Stop after the clusters already in flight", key=f"stop_{kind}",
                             help="Cancels the clusters that have not started. Calls already in "
                                  "flight are paid for and are left to finish."):
                    J.request_stop(cfg, kind)
                    st.rerun(scope="fragment")
            elif rec.get("stopping"):
                st.caption("Stopping — waiting for the calls already in flight to return.")
        else:
            box = {J.DONE: st.success, J.FAILED: st.error}.get(state, st.warning)
            box(f"**{label}** — {state}. {rec.get('message', '')}")
            if st.button("Dismiss", key=f"dismiss_{kind}"):
                J.clear(cfg, kind)
                st.rerun()          # full rerun: the page behind this box is now stale
            # One full rerun when a job lands, so the rest of the page reflects the new files.
            seen = st.session_state.setdefault("_jobs_seen", set())
            stamp = (kind, rec.get("finished_at"))
            if stamp not in seen:
                seen.add(stamp)
                st.rerun()


def confirm_action(label: str, prompt: str, key: str, *, disabled: bool = False,
                   help: str | None = None) -> bool:
    """A button that needs a second 'Yes, I'm sure' click before it fires.

    Guards paid AI generation against accidental clicks. Three states:

      idle    the plain button
      armed   the `prompt` as an 'are you sure?', with Yes / Cancel both live
      firing  Yes has been clicked: both buttons render **disabled** and this returns True once

    The firing state exists because the confirming click used to leave a live, primary-styled
    "Yes, I'm sure" on screen while the paid work ran — inviting a second click on an action that
    costs money. Clicking Yes now reruns immediately into a cold row, and only then is the work
    launched, so there is no live button to double-click and no doubt about whether the first
    click registered.
    """
    armed, firing = f"_armed_{key}", f"_firing_{key}"
    if st.session_state.pop(firing, False):
        st.warning(prompt)
        c1, c2 = st.columns(2)
        c1.button("⏳ Starting…", key=f"{key}__yes_cold", disabled=True, width="stretch")
        c2.button("Cancel", key=f"{key}__no_cold", disabled=True, width="stretch")
        return True
    if st.session_state.get(armed):
        st.warning(prompt)
        cy, cn = st.columns(2)
        if cy.button("✅ Yes, I'm sure", key=f"{key}__yes", type="primary", width="stretch"):
            st.session_state[armed] = False
            st.session_state[firing] = True
            st.rerun()
        if cn.button("Cancel", key=f"{key}__no", width="stretch"):
            st.session_state[armed] = False
            st.rerun()
        return False
    if st.button(label, key=f"{key}__start", disabled=disabled, help=help,
                 width="stretch"):
        st.session_state[armed] = True
        st.rerun()
    return False


@st.cache_data(show_spinner=False)
def _file_bytes(path: str, _mtime: float) -> bytes:
    """File bytes, cached on (path, mtime), so download buttons don't re-read multi-MB PNGs every
    rerun. `_mtime` busts the cache when the page is re-rendered at a new DPI."""
    with open(path, "rb") as fh:
        return fh.read()


@st.cache_data(show_spinner=False)
def _page_count(path: str, _mtime: float) -> int:
    """Page count, cached on (path, mtime), so we don't reopen the PDF with fitz every rerun."""
    return P.page_count(path)


def _state_stamp(cfg, clusters) -> tuple[float, float, float]:
    """(annotation-db mtime, newest insight mtime, cell-counts mtime) — the overview's cache key.

    The insight directory's own mtime is not enough: regenerating rewrites existing files in
    place, which does not touch the directory. Stat the files instead.
    """
    try:
        db = D.annotations_db(cfg).stat().st_mtime
    except OSError:
        db = 0.0
    try:
        ins = max((p.stat().st_mtime for p in I.insight_dir(cfg).glob("cluster_*.json")),
                  default=0.0)
    except OSError:
        ins = 0.0
    try:
        counts = D.cell_counts_path(cfg).stat().st_mtime
    except OSError:
        counts = 0.0
    return db, ins, counts


@st.cache_data(show_spinner=False)
def cluster_overview(_cfg, _clusters, _name: str, _db_mtime: float, _ins_mtime: float,
                     _counts_mtime: float):
    """One row per cluster: what has been decided, what the AI said, what you've marked.

    Underscore-prefixed args are excluded from Streamlit's cache key; the mtimes are what
    invalidate it, so saving an annotation or regenerating an insight refreshes this while a mere
    slider move does not re-read 34 JSON files and ~140 SQLite rows.
    """
    review = I.load_cohort_review(_cfg)
    flagged = set(I.flagged_clusters(review))
    annotations = S.all_annotations(_cfg)
    counts = D.load_cell_counts(_cfg)
    rows = []
    for c in _clusters:
        a = annotations.get(c, {})
        rec = I.load_insight(_cfg, c)
        reann = I.load_reannotation(_cfg, c)
        rows.append({
            "cluster": c,
            "ncells": a.get("ncells") if a.get("ncells") is not None
            else counts.get(c, {}).get("ncells"),
            "reviewed": bool(a.get("reviewed")),
            "label": a.get("annot_type") or "",
            "abbrev": a.get("annot_abbrev") or "",
            "order": a.get("annot_order"),
            "ai_identity": (rec or {}).get("primary_identity") or "",
            "conf": (rec or {}).get("confidence") or "",
            "flagged": c in flagged,
            "reannotated": reann is not None,
            "reann_stale": I.reannotation_is_stale(_cfg, c),
            "stars_gene": len(S.selected_markers(_cfg, c, "gene")),
            "stars_motif": len(S.selected_markers(_cfg, c, "motif")),
            "notes": S.comment_count(_cfg, c),
            "refs": len((rec or {}).get("references") or []),
            "cost": I.record_cost(rec) + I.record_cost(reann),
        })
    # Identifier-safe column names throughout: itertuples() renames anything that isn't a valid
    # Python identifier to a positional _N, so "AI identity" / "★genes" silently break attribute
    # access. Display names live in the dataframe's column_config instead.
    return pd.DataFrame(rows)


def cluster_label(row) -> str:
    """Selector label carrying state, so the dropdown isn't 34 bare integers."""
    bits = [f"{row.cluster:>2}",
            "✓" if row.reviewed else ("·" if row.label else " ")]
    name = row.label or row.ai_identity or "—"
    bits.append(name[:44] + ("…" if len(name) > 44 else ""))
    tags = []
    if row.flagged:
        tags.append("🚩")
    if row.stars_gene or row.stars_motif:
        tags.append(f"★{row.stars_gene + row.stars_motif}")
    if row.notes:
        tags.append(f"💬{row.notes}")
    return " ".join(bits) + ("  " + " ".join(tags) if tags else "")


def img_with_download(png_path: str, caption: str, key: str, width: int | None = None):
    # width=None -> fill container (old behaviour); int -> cap display width so big PDF
    # renders fit a normal browser window. Full-res stays available via the download button.
    if width:
        st.image(png_path, caption=caption, width=width)  # hover ⤢ = fullscreen zoom
    else:
        st.image(png_path, caption=caption, width="stretch")
    st.download_button("⤓ Download full-resolution PNG", _file_bytes(png_path, os.path.getmtime(png_path)),
                       file_name=Path(png_path).name, mime="image/png", key=key)


# Every shortcut works by clicking a real button, so the key and the click go through exactly the
# same code path and a shortcut cannot drift from what the button does. (key, label, what it does) —
# this list is both the JS keymap and the on-screen legend, so they cannot disagree.
SHORTCUTS = [
    ("ArrowLeft",  "←",     "fp_prev",     "Previous feature plot"),
    ("ArrowRight", "→",     "fp_next",     "Next feature plot"),
    ("[",          "[",     "cl_prev",     "Previous cluster"),
    ("]",          "]",     "cl_next",     "Next cluster"),
    ("u",          "u",     "next_unrev",  "Jump to the next unreviewed cluster"),
]


def keyboard_nav(sections: list[str]):
    """Inject one keydown listener that clicks the button behind each shortcut.

    Digits 1-9 and 0 select the 1st-10th section via the hidden buttons in the shortcuts expander.
    Ignores key events while a text field is focused, so typing a label or a config path isn't
    hijacked. Re-injected each run; replaces the prior listener rather than stacking another.
    """
    keymap = {k: f".st-key-{key} button" for k, _, key, _ in SHORTCUTS}
    for i, _ in enumerate(sections[:10]):
        keymap[str((i + 1) % 10)] = f".st-key-sec_{i} button"
    inject_js(f"""
const doc = window.parent.document, w = window.parent;
const keymap = {json.dumps(keymap)};
if (w.__kbNav) doc.removeEventListener('keydown', w.__kbNav);
w.__kbNav = function(e){{
  if (e.ctrlKey || e.metaKey || e.altKey) return;
  const t = e.target, tag = t && t.tagName ? t.tagName.toUpperCase() : '';
  if (tag === 'INPUT' || tag === 'TEXTAREA' || (t && t.isContentEditable)) return;
  const sel = keymap[e.key];
  if (!sel) return;
  const el = doc.querySelector(sel);
  if (el && !el.disabled) {{ e.preventDefault(); el.click(); }}
}};
doc.addEventListener('keydown', w.__kbNav);
const mark = function(){{
  const el = doc.getElementById('kb-status');
  if (el) el.textContent = '⌨️ keyboard listener: active';
}};
mark();
w.requestAnimationFrame(mark);   // the badge may mount just after this script runs
""")


def shortcut_key(sections: list[str], section: str):
    """The on-screen legend, plus the hidden buttons the digit shortcuts click.

    A key press has to land on a real widget, so each digit gets a button. They are hidden with CSS
    rather than tucked inside the collapsed expander: Streamlit marks a closed expander's contents
    `inert`, and an inert element may ignore a synthetic click — which would make the shortcuts work
    only while the legend happened to be open.
    """
    with st.sidebar.expander("⌨️ Keyboard shortcuts"):
        st.markdown("\n".join(f"- **{shown}** — {what}" for _, shown, _, what in SHORTCUTS))
        st.markdown("- **1**…**9**, **0** — jump to a section:\n"
                    + "\n".join(f"    - **{(i + 1) % 10}** {name}"
                                + ("  ← here" if name == section else "")
                                for i, name in enumerate(sections[:10])))
        st.caption("Shortcuts are ignored while you're typing in a field. Streamlit's own "
                   "**C** = clear cache is disabled here, since that would discard paid AI work.")
        # Rewritten to "active" by the injected listener. If the injection ever breaks again, this
        # says so instead of the keys just quietly doing nothing.
        st.html('<span id="kb-status" style="font-size:0.8em;opacity:0.7">'
                '⌨️ keyboard listener: not detected</span>')
    for i, name in enumerate(sections[:10]):
        st.sidebar.button(name, key=f"sec_{i}", on_click=_go_to_section, args=(name,))
    hidden = ", ".join(f".st-key-sec_{i}" for i in range(len(sections[:10])))
    st.html(f"<style>{hidden} {{ display: none; }}</style>")


def _go_to_section(name: str):
    st.session_state["section_pick"] = name
    st.session_state["section_last"] = name


# --------------------------------------------------------------- annotation of record
def tsv_pref_keys(cfg) -> tuple[str, str]:
    """Pref keys for the TSV path and the auto-write flag — per dataset, since the path is."""
    return f"tsv_path:{cfg['name']}", f"tsv_auto:{cfg['name']}"


def current_tsv_path(cfg) -> str:
    """The configured TSV target: what you last typed in the export panel, else `output_tsv:` from
    the YAML, else a dated default beside the annotation DB."""
    pk, _ = tsv_pref_keys(cfg)
    return str(D.load_prefs().get(pk) or E.tsv_path(cfg))


def tsv_auto_on(cfg) -> bool:
    """Whether saving an annotation also rewrites the TSV. Defaults on when the dataset config
    names an `output_tsv:` — asking for the file is asking for it to be kept current."""
    _, ak = tsv_pref_keys(cfg)
    return bool(D.load_prefs().get(ak, bool(cfg.get("output_tsv"))))


def write_tsv_now(cfg, only_annotated: bool | None = None, quiet: bool = False):
    """Write the TSV to the configured path. Never raises: a save must not be lost to a bad path.

    only_annotated=None follows the Export panel's "Include unlabelled clusters" toggle, so the
    file a save rewrites has the same shape as one written by the button.
    """
    if only_annotated is None:
        only_annotated = not st.session_state.get("exp_all", False)
    try:
        path, n = E.write_tsv(cfg, st.session_state.clusters, st.session_state.genes,
                              st.session_state.motifs, out_path=current_tsv_path(cfg),
                              only_annotated=only_annotated)
        if not quiet:
            st.toast(f"Wrote {n} row(s) to {Path(path).name}")
        return path, n
    except Exception as exc:  # noqa: BLE001
        st.warning(f"Could not write the TSV: {exc}")
        return None, 0


def sidebar_annotation(cfg, cluster: int, counts: dict):
    """The annotation form itself — in the sidebar, so it stays reachable from every section.

    Everything in the main pane is evidence; this is where a decision gets written down. It lives
    beside the evidence rather than in a section of its own so you can type the label while looking
    at the marker table or the spatial projection that justifies it. The AI annotation is offered as
    a prefill: a first draft to accept or edit, not something to read in one pane and retype here.
    """
    a = S.get_annotation(cfg, cluster)
    sugg = E.suggest_from_insight(cfg, cluster)
    n = counts.get(cluster, {})

    with st.sidebar.expander(f"✍️ Annotation — cluster {cluster}", expanded=True):
        if sugg:
            # Buttons inside st.form only submit, so the prefill lives outside it: it writes the
            # widgets' session_state and reruns, and session_state takes precedence over `value=`.
            if st.button("↙ Prefill from AI", width="stretch", key="an_prefill",
                         help=f"AI proposes: {sugg.get('annot_type', '—')}. Copies it and its "
                              "PMIDs into the form. Nothing is saved until you press Save."):
                for k, v in sugg.items():
                    st.session_state[f"an_{k}_{cluster}"] = v
                st.rerun()

        with st.form(f"annot_form_{cluster}"):
            vals = {
                "annot_type": st.text_input(
                    "annot_type — full label", value=a.get("annot_type") or "",
                    key=f"an_annot_type_{cluster}", placeholder="Forebrain progenitors"),
                "annot_abbrev": st.text_input(
                    "annot_abbrev — short label", value=a.get("annot_abbrev") or "",
                    key=f"an_annot_abbrev_{cluster}", placeholder="Fb"),
                "annot_origin": st.text_input(
                    "annot_origin — germ layer / lineage", value=a.get("annot_origin") or "",
                    key=f"an_annot_origin_{cluster}", placeholder="Ectoderm"),
                "reference_omg_type": st.text_input(
                    "reference_omg_type — matching atlas type",
                    value=a.get("reference_omg_type") or "",
                    key=f"an_reference_omg_type_{cluster}", placeholder="Telencephalon"),
                "annot_trajectory": st.text_input(
                    "annot_trajectory", value=a.get("annot_trajectory") or "",
                    key=f"an_annot_trajectory_{cluster}", placeholder="-"),
                "annot_order": st.number_input(
                    "annot_order — figure order", min_value=0, step=1,
                    value=int(a["annot_order"]) if a.get("annot_order") is not None else 0,
                    key=f"an_annot_order_{cluster}",
                    help="0 = unset; unnumbered clusters sort to the end of the export."),
                "refs": st.text_input(
                    "refs", value=a.get("refs") or "", key=f"an_refs_{cluster}",
                    placeholder="PMIDs:26371318,25820448"),
                # Prefilled from the counts computed off the Seurat object; a saved value wins.
                "ncells": st.number_input(
                    "ncells", min_value=0, step=1,
                    value=int(a.get("ncells") if a.get("ncells") is not None
                              else n.get("ncells") or 0),
                    key=f"an_ncells_{cluster}",
                    help="0 = unset. Prefilled from the Seurat object when available."),
                "pct_cells": st.number_input(
                    "%cells", min_value=0.0, max_value=1.0, step=0.0001, format="%.6f",
                    value=float(a.get("pct_cells") if a.get("pct_cells") is not None
                                else n.get("pct_cells") or 0.0),
                    key=f"an_pct_cells_{cluster}",
                    help="Fraction, not percent (0.0508 = 5.08%). Prefilled from the object."),
            }
            reviewed = st.checkbox("Mark this cluster reviewed", value=bool(a.get("reviewed")),
                                   key=f"an_reviewed_{cluster}")
            if st.form_submit_button("💾 Save annotation", type="primary", width="stretch"):
                # 0 is the "unset" sentinel for the numeric fields — store NULL so the export
                # leaves the cell blank rather than asserting a real zero.
                clean = dict(vals)
                for k in ("annot_order", "ncells"):
                    clean[k] = int(clean[k]) or None
                clean["pct_cells"] = float(clean["pct_cells"]) or None
                S.set_annotation(cfg, cluster, clean, reviewed=reviewed)
                st.toast(f"Annotation saved for cluster {cluster}")
                if tsv_auto_on(cfg):
                    write_tsv_now(cfg)
        if n:
            st.caption(f"Object says cluster {cluster} has **{n['ncells']:,}** cells "
                       f"({n['pct_cells'] * 100:.2f}%).")
        stars_g = len(S.selected_markers(cfg, cluster, "gene"))
        stars_m = len(S.selected_markers(cfg, cluster, "motif"))
        st.caption(f"★ {stars_g} gene(s), ★ {stars_m} motif(s), "
                   f"💬 {S.comment_count(cfg, cluster)} note(s) on this cluster. "
                   "The ★ markers and these notes fill the remaining export columns.")


def cell_counts_panel(cfg, clusters):
    """What ncells / %cells were computed from, and a way to redo it."""
    st.markdown("### Cells per cluster")
    counts = st.session_state.get("cell_counts") or {}
    meta = D.cell_counts_meta(cfg)
    cc = D.cell_counts_cfg(cfg)
    if counts:
        total = meta.get("total_cells") or sum(v["ncells"] for v in counts.values())
        src = (f"`{Path(str(meta['object'])).name}` · column `{meta['cluster_column']}`"
               if meta.get("object") else f"`{D.cell_counts_path(cfg)}`")
        st.caption(f"**{total:,}** cells over {len(counts)} clusters, from {src}"
                   + (f" · computed {meta['generated_at']}" if meta.get("generated_at") else "")
                   + ". These prefill `ncells` / `%cells`; anything you type in the form wins.")
        missing = [c for c in clusters if c not in counts]
        if missing:
            st.warning(f"No count for cluster(s) {missing} — the cluster ids in the object and in "
                       "the marker tables don't line up completely.")
    elif cc["seurat_object"]:
        st.caption(f"Not counted yet, from `{cc['seurat_object']}`."
                   + ("" if cc["auto"] else " `cell_counts.auto` is off, so loading won't do it — "
                                            "use the button below or `./run_app.sh counts`."))
    else:
        st.caption("No Seurat object configured, so `ncells` / `%cells` are yours to type. Set "
                   "`cell_counts.seurat_object` (or `seurat_object`) in the dataset config to have "
                   "them counted, or point `cell_counts.tsv` at a table you already have.")
    if cc["seurat_object"] and not cc["tsv"] and st.button(
            "↻ Recount from the Seurat object",
            help="Reads the object again — tens of seconds and several GB. Needed only if the "
                 "object or its clustering changed."):
        load_cell_counts(cfg, force=True)
        st.rerun()


def export_panel(cfg, clusters, genes, motifs):
    """xlsx and TSV export. Sits with the all-clusters view: exporting is a whole-dataset action."""
    st.markdown("### Export")
    done = S.annotated_clusters(cfg)
    st.caption(f"{len(done)} of {len(clusters)} clusters have a label. Both formats carry the same "
               "13 columns as the annotation template.")
    only_annotated = not st.toggle(
        "Include unlabelled clusters", key="exp_all",
        help="Off: only clusters with a label. On: every cluster, so you can see what is "
             "still outstanding.")

    xl, tsv = st.columns(2)
    with xl:
        if st.button("⤓ Build xlsx export", width="stretch", type="primary"):
            try:
                path, n = E.write_xlsx(cfg, clusters, genes, motifs, only_annotated=only_annotated)
                st.session_state["export_path"] = str(path)
                st.success(f"Wrote {n} row(s) to `{path}`")
            except Exception as exc:
                st.error(f"Export failed: {exc}")
        saved = st.session_state.get("export_path")
        if saved and Path(saved).exists():
            st.download_button("⤓ Download the xlsx", _file_bytes(saved, os.path.getmtime(saved)),
                               file_name=Path(saved).name, key="dl_xlsx", width="stretch",
                               mime="application/vnd.openxmlformats-officedocument."
                                    "spreadsheetml.sheet")
    with tsv:
        pk, ak = tsv_pref_keys(cfg)
        path_str = st.text_input("Output TSV", value=current_tsv_path(cfg), key="tsv_path_in",
                                 help="Written in place each time — a fixed path a downstream "
                                      "script can keep reading. Set `output_tsv:` in the dataset "
                                      "config to make it the default.")
        auto = st.toggle("Rewrite it on every save", value=tsv_auto_on(cfg), key="tsv_auto_in",
                         help="Keeps the file in step with the annotation DB, so you never have "
                              "to remember to export.")
        if (path_str, auto) != (current_tsv_path(cfg), tsv_auto_on(cfg)):
            D.save_prefs({pk: path_str, ak: auto})
        if st.button("⤓ Write TSV now", width="stretch"):
            path, n = write_tsv_now(cfg, only_annotated=only_annotated, quiet=True)
            if path:
                st.success(f"Wrote {n} row(s) to `{path}`")
        existing = Path(current_tsv_path(cfg))
        if existing.exists():
            st.caption("Last written "
                       f"{pd.Timestamp.fromtimestamp(existing.stat().st_mtime):%Y-%m-%d %H:%M}")
            st.download_button("⤓ Download the TSV",
                               _file_bytes(str(existing), os.path.getmtime(existing)),
                               file_name=existing.name, key="dl_tsv", width="stretch",
                               mime="text/tab-separated-values")

    with st.expander("Preview the export"):
        st.dataframe(E.build_table(cfg, clusters, genes, motifs, only_annotated=only_annotated),
                     hide_index=True, width="stretch")
        st.caption("`key_marker_genes` / `select_marker_motifs` come from the **★** stars on the "
                   "marker tables, `comments` from your notes, and `ncells` / `%cells` from the "
                   "Seurat object unless you typed your own.")


# --------------------------------------------------------------- sidebar
st.sidebar.title("🧬 Cluster Annotation")
_prefs = D.load_prefs()
_configs = D.available_configs()
if not _configs:
    # First run after a fresh clone: only the template exists, and it points at placeholders.
    st.info(
        "**No dataset configured yet.** Copy the template and edit the paths in it:\n\n"
        "```bash\ncp config/dataset.template.yaml config/mydataset.yaml\n"
        "./run_app.sh doctor config/mydataset.yaml\n```\n\n"
        "Then reload this page. See the README for what each field means.")
    st.stop()
_names = [str(p) for p in _configs]
_last = _prefs.get("config_path")
_default_idx = _names.index(_last) if _last in _names else 0
# A dropdown of the configs that exist, not a path to retype. Free-text remains available for a
# config living outside config/.
cfg_path = st.sidebar.selectbox("Dataset config", _names, index=_default_idx,
                                format_func=lambda p: Path(p).name, key="cfg_pick")
col_a, col_b = st.sidebar.columns(2)
if col_a.button("Load dataset", type="primary"):
    load_dataset(cfg_path)
    D.save_prefs({"config_path": cfg_path})
if col_b.button("Re-preprocess", help="Force re-export markers from RDS"):
    load_dataset(cfg_path, force=True)
    D.save_prefs({"config_path": cfg_path})

# Auto-load on first visit so the app is usable without a click, and so a browser refresh (which
# clears session_state) doesn't dump you back to a "click Load dataset" screen. Guarded: a broken
# config falls back to the manual prompt with the error shown rather than a bare traceback.
if "cfg" not in st.session_state and not st.session_state.get("autoload_failed"):
    try:
        load_dataset(cfg_path)
    except Exception as exc:  # noqa: BLE001 - surface it, don't crash the page
        st.session_state["autoload_failed"] = True
        st.error(f"Could not auto-load `{Path(cfg_path).name}`: {exc}")

if "cfg" not in st.session_state:
    st.info("Pick a dataset config in the sidebar and click **Load dataset**.")
    st.stop()

cfg = st.session_state.cfg
clusters = st.session_state.clusters
st.sidebar.success(f"**{cfg['name']}** · res {cfg.get('resolution','?')} · {len(clusters)} clusters")
with st.sidebar:                  # a fragment writes to the active container, and
    sidebar_job_indicator(cfg)    # cannot address st.sidebar itself

# Cluster picker, labelled with state (decision / AI call / flags / stars / notes) so progress
# across 34 clusters is visible instead of being 34 bare integers.
_ov = cluster_overview(cfg, clusters, cfg["name"], *_state_stamp(cfg, clusters))
_labels = {int(r.cluster): cluster_label(r) for r in _ov.itertuples()}
cluster = st.sidebar.selectbox("Cluster", clusters, index=0,
                              format_func=lambda c: _labels.get(c, str(c)), key="cluster_pick")


def _step_cluster(delta: int):
    """Walk the cluster list, wrapping. A callback so the selectbox follows before the rerun."""
    cur = clusters.index(st.session_state.get("cluster_pick", clusters[0]))
    st.session_state["cluster_pick"] = clusters[(cur + delta) % len(clusters)]


_cp, _cn = st.sidebar.columns(2)
_cp.button("◀ Cluster", key="cl_prev", on_click=_step_cluster, args=(-1,), width="stretch")
_cn.button("Cluster ▶", key="cl_next", on_click=_step_cluster, args=(1,), width="stretch")

_todo = [int(r.cluster) for r in _ov.itertuples() if not r.reviewed]
if _todo:
    def _next_unreviewed():
        st.session_state["cluster_pick"] = next((c for c in _todo if c > cluster), _todo[0])

    st.sidebar.button(f"⏭ Next unreviewed ({len(_todo)} left)", width="stretch", key="next_unrev",
                      on_click=_next_unreviewed,
                      help="Jump to the lowest-numbered cluster not yet marked reviewed.")

# The annotation of record — in the sidebar so it can be filled in while browsing any section.
sidebar_annotation(cfg, cluster, st.session_state.get("cell_counts") or {})

# Per-cluster notes (stored in the SQLite annotation DB; shown in the Comments tab).
with st.sidebar.form("note_form", clear_on_submit=True):
    _note = st.text_area("📝 Note for this cluster", height=90, placeholder="Add a note…",
                         label_visibility="visible")
    if st.form_submit_button("Add note") and _note.strip():
        S.add_comment(cfg, cluster, _note)
        st.toast(f"Note added to cluster {cluster}")
st.sidebar.caption(f"💬 {S.comment_count(cfg, cluster)} note(s) for cluster {cluster} "
                   "— see the **Comments** tab.")

render_dpi = st.sidebar.slider("Image render DPI", 72, 300, int(_prefs.get("render_dpi", 150)), 6,
                               help="Higher = sharper but slower. Use fullscreen (⤢) or download to zoom.")
disp_w = st.sidebar.slider("Max image width (px)", 400, 1600, int(_prefs.get("disp_w", 800)), 50,
                           help="On-screen size for PDF renders (Tangram / UMAP / feature / other). "
                                "Full resolution is always available via the download button.")
# Remember the slider positions across sessions — they describe your screen, not the dataset, and
# re-tuning them on every launch was pure friction. Written only when they actually change.
if (render_dpi, disp_w) != (_prefs.get("render_dpi"), _prefs.get("disp_w")):
    D.save_prefs({"render_dpi": render_dpi, "disp_w": disp_w})

# Model choice: remembered per dataset in .run/ui_prefs.json, so it survives a browser refresh,
# a second tab, and an app restart. It used to live in session_state only, which meant a reload
# silently reverted to the YAML default — you would pick Opus 5, refresh, and spend on the old
# model without anything on screen disagreeing with you.
_aic = D.ai_insights_cfg(cfg)             # YAML values, before any override
_model_prefs = D.load_prefs()
_pm_key, _fm_key = f"ai_primary:{cfg['name']}", f"ai_fallback:{cfg['name']}"
_opts = list(dict.fromkeys(
    D.KNOWN_MODELS + [_aic["primary_model"], _aic["fallback_model"],
                      _model_prefs.get(_pm_key), _model_prefs.get(_fm_key)]))
_opts = [m for m in _opts if m]
# Seed the widgets from prefs once per session and let `key=` carry them afterwards. Passing both
# `index=` and `key=` makes Streamlit warn that one of them is being ignored — the ambiguity that
# made this control look broken.
for _k, _stored, _yaml_val in ((("ai_primary", _model_prefs.get(_pm_key), _aic["primary_model"])),
                               (("ai_fallback", _model_prefs.get(_fm_key), _aic["fallback_model"]))):
    if st.session_state.get(_k) not in _opts:
        st.session_state[_k] = _stored if _stored in _opts else _yaml_val

with st.sidebar.expander("AI insights — models"):
    _pm = st.selectbox("Primary model", _opts, key="ai_primary")
    _fm = st.selectbox("Fallback model", _opts, key="ai_fallback")
    _workers = st.slider("Parallel clusters (workers)", 1, 16, int(_aic["max_workers"]), 1,
                         help="Clusters generated concurrently. Higher = faster up to your API "
                              "rate limits (the SDK auto-retries throttling); ~8-10 is a good ceiling.")
    if (_pm, _fm) != (_model_prefs.get(_pm_key), _model_prefs.get(_fm_key)):
        D.save_prefs({_pm_key: _pm, _fm_key: _fm})
    st.caption(
        f"Remembered for **{cfg['name']}** (the YAML default is `{_aic['primary_model']}`). The "
        "**AI insights** tab shows what the next run will use. Changing a model marks cached "
        "annotations stale — they refresh only when you regenerate, never on load.")

# The override the rest of the app reads. Written onto the loaded config so every insights call
# picks it up, since they all resolve models through D.ai_insights_cfg(cfg).
cfg.setdefault("ai_insights", {})
cfg["ai_insights"]["primary_model"] = _pm
cfg["ai_insights"]["fallback_model"] = _fm
cfg["ai_insights"]["max_workers"] = _workers

cache_imgs = D.cache_dir(cfg) / "img"


# --------------------------------------------------------------- marker table panel
def marker_panel(kind: str):
    df = st.session_state.genes if kind == "gene" else st.session_state.motifs
    eff = EFFECT[kind]
    c1, c2 = st.columns([1, 2])
    n = c1.slider(f"Top N {kind}s", 5, 100, 20, 5, key=f"n_{kind}")
    sort_opts = {
        "p_val_adj (most significant)": ("p_val_adj", True),
        f"{eff} (effect size ↓)": (eff, False),
        "delta_pct = pct.1 − pct.2 (specificity ↓)": ("delta_pct", False),
    }
    sort_label = c2.radio("Sort by", list(sort_opts), horizontal=True, key=f"s_{kind}")
    col, asc = sort_opts[sort_label]

    selected = S.selected_markers(cfg, cluster, kind)
    only_sel = st.checkbox(f"★ selected only ({len(selected)})", key=f"mk_only_{kind}",
                           help="Show only the markers you've starred for this cluster.")
    if only_sel:
        sub = (df[(df["cluster"] == cluster) & (df["feature"].isin(selected))]
               .sort_values(col, ascending=asc, kind="mergesort").reset_index(drop=True))
    else:
        sub = D.top_markers(df, cluster, n, col, asc)

    show = ["feature", eff, "pct.1", "pct.2", "delta_pct", "p_val_adj"]
    disp = sub[show].copy()
    disp.insert(0, "★", disp["feature"].isin(selected))
    edited = st.data_editor(
        disp, hide_index=True, width="stretch",
        key=f"mk_ed_{kind}_{cluster}_{only_sel}", disabled=show,  # only ★ is editable
        column_config={
            "★": st.column_config.CheckboxColumn("★", help="Mark a maximally-informative marker",
                                                 default=False),
            eff: st.column_config.NumberColumn(format="%.3f"),
            "delta_pct": st.column_config.NumberColumn("Δpct", format="%.3f"),
            "p_val_adj": st.column_config.NumberColumn(format="%.2e"),
        })
    S.set_markers(cfg, cluster, kind, set(disp["feature"]), set(edited.loc[edited["★"], "feature"]))
    st.caption(f"Cluster {cluster}: {len(df[df.cluster==cluster])} total {kind} markers · "
               f"{len(selected)} starred. Tick ★ to keep a marker; **★ selected only** filters to them.")


# --------------------------------------------------------------- feature-plot stepper
def feature_stepper(kind: str, width: int):
    pdf = D.featureplot_pdf(cfg, kind, cluster)
    if not pdf.exists():
        st.warning(f"Feature-plot PDF not found: {pdf}")
        return
    feats = D.feature_pages(cfg, kind, cluster)
    feats = feats[feats.feature != "__UMAP__"].reset_index(drop=True)
    if feats.empty:
        st.info(f"No {kind} feature plots for this cluster (no markers).")
        return

    labels = [f"p{r.page}: {r.feature}" for r in feats.itertuples()]
    sel_key = f"fpsel_{kind}_{cluster}"
    if st.session_state.get(sel_key) not in labels:
        st.session_state[sel_key] = labels[0]

    def _step(delta: int):  # runs as a callback (before rerun) so the selectbox follows
        cur = labels.index(st.session_state[sel_key])
        st.session_state[sel_key] = labels[(cur + delta) % len(labels)]

    c1, c2, c3 = st.columns([1, 1, 4])
    c1.button("◀ Prev", key="fp_prev", on_click=_step, args=(-1,), width="stretch")
    c2.button("Next ▶", key="fp_next", on_click=_step, args=(1,), width="stretch")
    choice = c3.selectbox("Jump to feature", labels, key=sel_key, label_visibility="collapsed")

    i = labels.index(choice)
    row = feats.iloc[i]
    page, feat = int(row.page), row.feature
    st.caption(f"{kind.capitalize()} plot {i + 1} / {len(labels)} · "
               "step with the ◀ ▶ buttons or the ← → arrow keys")
    png = P.render_page(pdf, page, render_dpi, cache_imgs)
    img_with_download(png, f"Cluster {cluster} · {kind} · {feat} (page {page})",
                      key=f"dl_{kind}", width=width)


# --------------------------------------------------------------- main tabs
# --------------------------------------------------------------- AI insights panel
def render_flag_list(flags: list):
    order = {"high": 0, "medium": 1, "low": 2}
    dot = {"high": "🔴", "medium": "🟡", "low": "🟢"}
    for f in sorted(flags, key=lambda x: order.get((x.get("severity") or "").lower(), 3)):
        sev = (f.get("severity") or "?").lower()
        cl = ", ".join(str(c) for c in (f.get("clusters") or []))
        st.markdown(f"{dot.get(sev, '⚪')} **{f.get('category', '?')}** · clusters {cl} "
                    f"· _{f.get('severity', '?')}_")
        st.write(f.get("issue", ""))
        if f.get("suggestion"):
            st.caption(f"→ {f['suggestion']}")


def render_insight(rec: dict):
    meta = rec.get("_meta", {})
    conf = (rec.get("confidence") or "").lower()
    badge = {"high": "🟢", "medium": "🟡", "low": "🔴"}.get(conf, "⚪")
    st.subheader(rec.get("primary_identity", "—"))
    if rec.get("revision_note"):
        st.info(f"**What changed vs original:** {rec['revision_note']}")
    st.markdown(f"{badge} **Confidence:** {rec.get('confidence', '?')}")
    alts = rec.get("alternative_identities") or []
    if alts:
        st.markdown("**Alternatives:** " + ", ".join(alts))
    st.markdown("**Reasoning**")
    st.write(rec.get("reasoning", ""))
    c1, c2 = st.columns(2)
    if rec.get("key_genes"):
        c1.markdown("**Key genes**\n\n" + ", ".join(f"`{g}`" for g in rec["key_genes"]))
    if rec.get("key_motifs"):
        c2.markdown("**Key motifs**\n\n" + ", ".join(f"`{m}`" for m in rec["key_motifs"]))
    if rec.get("caveats"):
        st.markdown("**Caveats**")
        st.write(rec["caveats"])
    refs = rec.get("references") or []
    if refs:
        st.markdown("**References** — web-grounded but AI-selected; *verify before citing*.")
        for i, r in enumerate(refs, 1):
            cite = r.get("citation", "")
            pmid = r.get("pmid")
            url = r.get("url") or (f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else None)
            head = f"[{cite}]({url})" if url else cite
            if pmid:
                head += f"  ·  PMID {pmid}"
            st.markdown(f"{i}. {head}")
            if r.get("supports"):
                st.caption(f"    ↳ {r['supports']}")
    served = meta.get("model_used", "?")
    when = meta.get("generated_at", "?")
    web = "web-searched" if meta.get("web_used") else "no web results"
    st.caption(f"Generated by **{served}** ({web}) · {when}. "
               "Footer shows the model that actually served — a reroute to the fallback is visible here.")
    st.caption(usage_line(meta))


def usage_line(meta: dict) -> str:
    """Token/cost line for an insight footer. Empty-ish for artifacts predating telemetry."""
    u = meta.get("usage") or {}
    if not u:
        return "_Cost not recorded — this annotation predates usage telemetry._"
    bits = [f"**~${u.get('est_cost_usd', 0):.3f}** est.",
            f"{u.get('api_calls', 0)} API call(s)",
            f"{u.get('input_tokens', 0):,} in / {u.get('output_tokens', 0):,} out tok"]
    if u.get("cache_read_input_tokens"):
        bits.append(f"{u['cache_read_input_tokens']:,} tok read from cache (billed ~0.1×)")
    if u.get("cache_creation_input_tokens"):
        bits.append(f"{u['cache_creation_input_tokens']:,} tok cache-write")
    if u.get("web_search_requests"):
        bits.append(f"{u['web_search_requests']} web search(es)")
    line = " · ".join(bits)
    if u.get("unpriced_models"):
        line += f"  ⚠️ no list price on file for {', '.join(u['unpriced_models'])} — cost understated"
    return line


def render_cohort(rev: dict):
    meta = rev.get("_meta", {})
    st.markdown("**Overall**")
    st.write(rev.get("overall", ""))
    st.markdown("**Coverage**")
    st.write(rev.get("coverage", ""))
    missing = rev.get("missing_expected") or []
    if missing:
        st.markdown("**Conspicuously missing cell types**")
        st.markdown("\n".join(f"- {m}" for m in missing))
    flags = rev.get("flags") or []
    if flags:
        st.markdown(f"**Cross-cluster flags ({len(flags)})**")
        render_flag_list(flags)
    else:
        st.success("No cross-cluster issues flagged.")
    st.caption(f"Reviewed {meta.get('n_clusters', '?')} clusters by **{meta.get('model_used', '?')}** "
               f"· {meta.get('generated_at', '?')}")
    st.caption(usage_line(meta))


# --------------------------------------------------------------- all-clusters overview
def overview_panel(cfg, clusters, overview):
    """The cross-cluster view: progress at a glance, plus search over the AI annotations.

    Everything else in the app is one cluster at a time, which makes it hard to notice that three
    clusters got three different mesothelial labels, or to find which cluster mentions notochord.
    """
    done = int(overview["reviewed"].sum())
    labelled = int((overview["label"] != "").sum())
    st.markdown(f"**{labelled}** of {len(clusters)} labelled · **{done}** marked reviewed · "
                f"**{int(overview['flagged'].sum())}** flagged by the cohort review · "
                f"~${overview['cost'].sum():.2f} of AI annotation recorded")

    q = st.text_input("🔎 Search the AI annotations", key="ov_q",
                      placeholder="e.g. notochord, microglia, Sim1 — matches identity, "
                                  "reasoning, caveats and key markers")
    view = overview
    if q.strip():
        needle = q.strip().lower()

        def hit(c: int) -> bool:
            for rec in (I.load_insight(cfg, c), I.load_reannotation(cfg, c)):
                if not rec:
                    continue
                blob = " ".join([
                    str(rec.get("primary_identity") or ""), str(rec.get("reasoning") or ""),
                    str(rec.get("caveats") or ""), " ".join(rec.get("alternative_identities") or []),
                    " ".join(rec.get("key_genes") or []), " ".join(rec.get("key_motifs") or []),
                ]).lower()
                if needle in blob:
                    return True
            return False

        view = overview[[hit(int(c)) for c in overview["cluster"]]]
        st.caption(f"{len(view)} cluster(s) mention “{q.strip()}”"
                   + (" — nothing matched" if view.empty else ""))

    show = ["cluster", "ncells", "reviewed", "label", "abbrev", "order", "ai_identity", "conf",
            "flagged", "reannotated", "reann_stale", "stars_gene", "stars_motif", "notes",
            "refs", "cost"]
    st.dataframe(
        view[show], hide_index=True, width="stretch",
        column_config={
            "ncells": st.column_config.NumberColumn(
                "cells", format="%d",
                help="Cells in the cluster, from the Seurat object (or what you typed)"),
            "reviewed": st.column_config.CheckboxColumn("✓", help="Marked reviewed"),
            "flagged": st.column_config.CheckboxColumn("🚩", help="Cohort review flagged it"),
            "reannotated": st.column_config.CheckboxColumn("↻", help="Has a reannotation"),
            "reann_stale": st.column_config.CheckboxColumn(
                "↻stale", help="Reannotation predates the current original — it describes a "
                               "version that has since been regenerated"),
            "label": st.column_config.TextColumn("your label", width="medium"),
            "abbrev": st.column_config.TextColumn("abbrev"),
            "ai_identity": st.column_config.TextColumn("AI identity", width="large"),
            "stars_gene": st.column_config.NumberColumn("★genes"),
            "stars_motif": st.column_config.NumberColumn("★motifs"),
            "refs": st.column_config.NumberColumn("refs", help="Anchored references on the AI call"),
            "cost": st.column_config.NumberColumn("$", format="%.3f",
                                                  help="Recorded AI spend for this cluster"),
        })
    st.caption("Read-only. Edit a cluster's label in the sidebar's **✍️ Annotation** form; the ★ "
               "counts come from the marker tables.")


st.header(f"Cluster {cluster}")
# Ordered the way a cluster is actually worked through: where is it (UMAP, spatial), what does it
# express (feature plots, markers), what else is known (other annotations), what does the AI think,
# then the cross-cluster views. The annotation form itself is in the sidebar, reachable from all of
# them.
SECTIONS = ["UMAP highlight", "Spatial (Tangram)", "Feature plots", "Marker genes",
            "Marker motifs", "Other annotations", "AI insights", "Cohort review", "Comments",
            "All clusters"]
# Lazy sections: only the selected section's code runs per rerun. st.tabs would execute ALL ten
# bodies every rerun (re-rendering the 26 MB Tangram, feature plots, etc. on every interaction) —
# the main cause of the freezing. segmented_control looks/behaves like a tab bar.
# Seeded through session_state rather than `default=`, because the shortcut buttons write that
# key — passing both makes Streamlit warn that the default is being ignored.
st.session_state.setdefault("section_pick", SECTIONS[0])
section = st.segmented_control("Section", SECTIONS, key="section_pick",
                               label_visibility="collapsed")
if section is None:  # clicking the active chip deselects it; keep the last section, like tabs
    section = st.session_state.get("section_last", SECTIONS[0])
st.session_state["section_last"] = section
shortcut_key(SECTIONS, section)   # sidebar legend + the buttons the digit keys click

if section == "All clusters":
    overview_panel(cfg, clusters, _ov)
    st.divider()
    cell_counts_panel(cfg, clusters)
    st.divider()
    export_panel(cfg, clusters, st.session_state.genes, st.session_state.motifs)

elif section == "Marker genes":
    marker_panel("gene")

elif section == "Marker motifs":
    marker_panel("motif")

elif section == "Spatial (Tangram)":
    tmap = st.session_state.tangram
    if cluster in tmap:
        png = P.render_page(D.resolve(cfg, cfg["tangram_pdf"]), tmap[cluster], render_dpi, cache_imgs)
        img_with_download(png, f"Tangram spatial projection · cluster {cluster}",
                          key="dl_tangram", width=disp_w)
    else:
        st.warning(f"No Tangram page found for cluster {cluster}.")

elif section == "UMAP highlight":
    pdf = D.featureplot_pdf(cfg, "gene", cluster)  # page 1 = UMAP highlight
    if pdf.exists():
        png = P.render_page(pdf, 1, render_dpi, cache_imgs)
        img_with_download(png, f"UMAP — cluster {cluster} highlighted", key="dl_umap", width=disp_w)
    else:
        st.warning(f"Feature-plot PDF not found: {pdf}")

elif section == "Feature plots":
    kind = st.radio("Feature type", ["gene", "motif"], horizontal=True,
                    format_func=str.capitalize, key="fp_kind")
    feature_stepper(kind, disp_w)

elif section == "Other annotations":
    items = D.other_annotations(cfg)
    if not items:
        st.info("No extra annotation PDFs configured. Add them under `other_annotations:` "
                "in the dataset YAML (a `name` and a `pdf` per entry, one page per cluster).")
    else:
        names = [it.get("name") or Path(str(it["pdf"])).stem for it in items]
        sel = st.selectbox("Annotation set", names, key="other_sel")
        entry = items[names.index(sel)]
        pdf = D.other_pdf(cfg, entry, cluster)
        if not pdf.exists():
            st.warning(f"PDF not found: {pdf}")
        else:
            npages = _page_count(str(pdf), os.path.getmtime(pdf))
            page = clusters.index(cluster) + 1  # positional: sorted-cluster order -> 1..N
            if page > npages:
                st.warning(f"'{sel}' has {npages} page(s); no page for cluster {cluster} "
                           f"(position {page} in sorted order).")
            else:
                png = P.render_page(pdf, page, render_dpi, cache_imgs)
                img_with_download(png, f"{sel} · cluster {cluster} (page {page})",
                                  key="dl_other", width=disp_w)

elif section == "AI insights":
    aic = D.ai_insights_cfg(cfg)
    if not aic["enabled"]:
        st.info("AI insights are disabled for this dataset (`ai_insights.enabled: false`).")
    elif not I.api_key_available(cfg):
        st.warning(
            "No Anthropic API key found, so AI insights are off. Add a key to "
            f"`{aic['api_key_file']}` (one line: `ANTHROPIC_API_KEY=sk-ant-…`, file `chmod 600`). "
            "Generate a dedicated key at https://console.claude.com to track this app's usage.")
    else:
        review = I.load_cohort_review(cfg)
        flags = I.cluster_flags(review, cluster)
        genes_df, motifs_df = st.session_state.genes, st.session_state.motifs
        st.caption("Generation makes **paid** API calls — each action below asks you to confirm, and "
                   "any current annotation is backed up before it's overwritten.")
        # State the models and the output ceiling here, at the point of spending: a run that used a
        # model you didn't intend, or truncated against a low max_tokens, is only obvious afterwards.
        _pricing = "" if I.price_for(aic["primary_model"]) else "  ⚠️ no list price on file — cost will be understated"
        st.caption(
            f"🤖 Next run uses **{aic['primary_model']}** (fallback `{aic['fallback_model']}`, "
            f"structuring `{aic['structuring_model']}`) · effort `{aic['effort']}` · "
            f"max_tokens `{aic['max_tokens']:,}`. Change the models in the sidebar.{_pricing}")
        if aic["max_tokens"] < 16000:
            st.warning(
                f"`max_tokens` is {aic['max_tokens']:,} for this dataset. A search-heavy call can "
                "spend that entire budget on thinking and web search and get truncated before it "
                "writes any analysis — which is what a *'truncated at max_tokens'* failure means. "
                "Raise `ai_insights.max_tokens` to 32000 in the dataset config.")
        st.caption(f"💰 {I.format_dataset_cost(cfg, clusters)}")

        # Anything in flight is shown here first, and blocks the other paid buttons: two runs on
        # one dataset would race on the same cache files, and it is never what you meant to buy.
        job_banner(cfg)
        _busy = J.any_running(cfg, AI_JOBS)
        _busy_help = (f"Waiting for the running job ({_busy}) to finish." if _busy else None)

        # In primer mode the dataset's reference sheet is the prerequisite for every per-cluster
        # annotation, so it needs a way in from the UI — not just the command-line runner.
        if aic["research_mode"] == "primer":
            primer = I.load_primer(cfg)
            if primer is None:
                st.info(
                    "**No reference primer yet.** In `primer` mode every cluster is annotated "
                    "against a dataset-level reference sheet built by one literature-search call. "
                    "Build it first — per-cluster annotation cannot run without it.")
            else:
                pm = primer.get("_meta") or {}
                nref = len(primer.get("references") or [])
                anchored = sum(1 for r in (primer.get("references") or [])
                               if r.get("pmid") or r.get("url"))
                st.success(
                    f"**Reference primer** · {len(primer.get('expected_cell_types') or [])} cell "
                    f"types · {nref} references ({anchored} with a PMID/DOI) · built "
                    f"{pm.get('generated_at', '?')} by {pm.get('model_used', '?')} · "
                    f"~${I.record_cost(primer):.2f}")
            if confirm_action(
                    "🔨 Rebuild reference primer" if primer else "🔨 Build reference primer",
                    "Build the dataset reference sheet? This makes one paid literature-search "
                    "call — the largest single call in the pipeline. Rebuilding it changes the "
                    "evidence every cluster is judged against, so all per-cluster annotations "
                    "will read as needing regeneration.",
                    key="primer_build", disabled=bool(_busy), help=_busy_help):
                if start_primer_build(cfg):
                    st.rerun()      # straight into the running banner
            with st.expander("What the primer contains"):
                if primer is None:
                    st.caption("Nothing yet.")
                else:
                    st.caption(primer.get("coverage_notes") or "")
                    st.dataframe(pd.DataFrame([
                        {"cell type": t.get("name", ""),
                         "markers": ", ".join(t.get("canonical_markers") or []),
                         "distinguish from": t.get("distinguishing_notes", "")}
                        for t in (primer.get("expected_cell_types") or [])],
                    ), hide_index=True, width="stretch")
                    st.dataframe(pd.DataFrame([
                        {"key": r.get("key", ""), "citation": r.get("citation", ""),
                         "id": r.get("pmid") or r.get("url") or "—",
                         "covers": r.get("covers", "")}
                        for r in (primer.get("references") or [])],
                    ), hide_index=True, width="stretch")
            st.divider()
        if confirm_action(
                f"Regenerate original — cluster {cluster}",
                f"Regenerate the AI annotation for cluster {cluster}? This makes a paid API call and "
                "overwrites the current annotation (a timestamped backup is kept).",
                key="ai_regen_one", disabled=bool(_busy), help=_busy_help):
            if start_ai_generation(cfg, force=True, only=[cluster]):
                st.rerun()
        reannot_help = None if review is not None else "Run a cohort review first (Cohort review tab)."
        if confirm_action(
                "Reannotate (use cohort review)",
                f"Reannotate cluster {cluster} using the cohort review? This makes a paid API call "
                "and writes a separate revised annotation (the original is kept).",
                key="ai_reannot", disabled=(review is None or bool(_busy)),
                help=reannot_help or _busy_help):
            if start_reannotation(cfg, [cluster], review):
                st.rerun()
        if confirm_action(
                "Regenerate ALL originals",
                f"Regenerate ALL {len(clusters)} cluster annotations? This re-runs the entire paid "
                "AI pass and overwrites every current annotation (each backed up first).",
                key="ai_regen_all", disabled=bool(_busy), help=_busy_help):
            if start_ai_generation(cfg, force=True):
                st.rerun()

        st.markdown("### 1 · Original annotation")
        rec = I.load_insight(cfg, cluster)
        if rec is None:
            st.info("No insight yet for this cluster. Use **Regenerate original** above "
                    "(or **Regenerate ALL originals** to do every cluster). Generation never runs "
                    "automatically.")
        else:
            render_insight(rec)

        st.divider()
        st.markdown("### 2 · Cohort review — flags for this cluster")
        if review is None:
            st.caption("No cohort review yet — run it in the **Cohort review** tab.")
        elif not flags:
            st.success("No cohort flags reference this cluster.")
        else:
            render_flag_list(flags)

        st.divider()
        st.markdown("### 3 · Reannotation")
        re = I.load_reannotation(cfg, cluster)
        if re is None:
            st.caption("None yet. Use **Reannotate (use cohort review)** above — it keeps the "
                       "original above and writes a separate revised annotation here.")
        else:
            if I.reannotation_is_stale(cfg, cluster):
                st.warning(
                    "⚠️ This reannotation is older than the original above — the original has been "
                    "regenerated since, so this describes a version that no longer exists. Treat "
                    "it as history, or reannotate again to refresh it.")
            render_insight(re)
            # Without this a reannotation is inert: you can read it but never say "this is the one
            # I accept", so the revised call never reaches the export.
            st.divider()
            pc1, pc2 = st.columns([1, 2])
            if pc1.button("✓ Accept as my annotation", key="promote_reann", width="stretch",
                          help="Copies this reannotation's cell type and PMIDs into the "
                               "Annotation form for this cluster."):
                vals = E.suggest_from_insight(cfg, cluster, prefer_reannotation=True)
                if vals:
                    S.set_annotation(cfg, cluster, vals)
                    for k, v in vals.items():          # keep the open form in sync
                        st.session_state[f"an_{k}_{cluster}"] = v
                    st.success(f"Saved as the annotation for cluster {cluster}: "
                               f"**{vals.get('annot_type', '')}**. Refine it in the sidebar's "
                               "**✍️ Annotation** form.")
                else:
                    st.warning("Nothing to copy from this reannotation.")
            _cur = S.get_annotation(cfg, cluster).get("annot_type")
            pc2.caption(f"Your current label: **{_cur}**" if _cur
                        else "You have no label recorded for this cluster yet.")

elif section == "Cohort review":
    aic = D.ai_insights_cfg(cfg)
    if not aic["enabled"]:
        st.info("AI insights are disabled for this dataset (`ai_insights.enabled: false`).")
    elif not I.api_key_available(cfg):
        st.warning("No Anthropic API key found — see the **AI insights** tab to set one up.")
    else:
        n_done = sum(1 for c in clusters if I.load_insight(cfg, c) is not None)
        st.caption(f"A 'head' pass that reviews all {n_done}/{len(clusters)} generated cluster "
                   "annotations together — flags over-split/redundant clusters, inconsistencies, "
                   "and missing expected cell types. One call, no web search.")
        job_banner(cfg)
        _busy = J.any_running(cfg, AI_JOBS)
        _busy_help = (f"Waiting for the running job ({_busy}) to finish." if _busy else None)
        if confirm_action(
                "Run cohort review",
                f"Run the cohort review over all {n_done} generated annotations? This makes a paid "
                "API call.",
                key="cohort_run", disabled=(n_done == 0 or bool(_busy)),
                help=_busy_help):
            if start_cohort_review(cfg, clusters):
                st.rerun()
        if n_done == 0:
            st.info("Generate per-cluster insights first (AI insights tab), then run the review.")
        rev = I.load_cohort_review(cfg)
        if rev is not None:
            cov = rev.get("_meta", {}).get("n_clusters")
            if cov is not None and cov < n_done:
                st.warning(f"Cached review covered {cov} clusters; {n_done} now exist — "
                           "re-run to include the rest.")
            flagged = I.flagged_clusters(rev)
            if flagged:
                st.markdown(f"**{len(flagged)} flagged cluster(s):** {flagged}")
                if confirm_action(
                        "Reannotate all flagged clusters",
                        f"Reannotate all {len(flagged)} flagged clusters? This makes a paid API call "
                        "per cluster; each original is kept and a separate reannotation written.",
                        key="reannot_all",
                        help="Keeps each original; writes a separate reannotation per cluster. "
                             "Review all three (original · flags · reannotation) in the AI insights tab.",
                        disabled=bool(_busy)):
                    if start_reannotation(cfg, flagged, rev):
                        st.rerun()
            render_cohort(rev)
        elif n_done:
            st.info("No cohort review yet. Click **Run cohort review**.")

elif section == "Comments":
    show_all = st.toggle("Show all clusters", key="cmt_all")
    st.subheader("Comments — all clusters" if show_all else f"Comments — cluster {cluster}")
    st.caption("Add notes from the **📝 Note for this cluster** box in the sidebar. Newest at the bottom.")
    rows = S.list_comments(cfg, None if show_all else cluster)
    if not rows:
        st.info("No notes yet. Add one from the sidebar.")
    else:
        for r in rows:  # chronological: oldest first
            meta = r["created_at"] + (f" · cluster {r['cluster']}" if show_all else "")
            cc1, cc2 = st.columns([12, 1])
            with cc1:
                st.caption(meta)
                st.write(r["text"])
            if cc2.button("🗑", key=f"cmt_del_{r['id']}", help="Delete this note"):
                S.delete_comment(cfg, r["id"])
                st.rerun()
            st.divider()


# One keydown listener for the whole app, injected last so every button it targets exists.
keyboard_nav(SECTIONS)
