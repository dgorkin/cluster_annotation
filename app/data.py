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
CELL_COUNTS_FILE = "cluster_cells.tsv"

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


def rscript_bin(cfg: dict) -> str:
    """The Rscript to run the preprocessing scripts with. Both are base-R only, so the system one
    is normally right; `rscript:` in the config points at a specific install if needed."""
    return str(cfg.get("rscript") or "Rscript")


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
    cmd = [rscript_bin(cfg), str(script),
           str(resolve(cfg, cfg["markers_rds"])),
           str(resolve(cfg, cfg["motif_markers_rds"])),
           str(cd), gene_topn, motif_topn, gene_sep]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"export_markers.R failed:\n{res.stderr}\n{res.stdout}")
    return cd


# ---------------------------------------------------------------- cells per cluster
# ncells and %cells are two of the export columns, and typing 42 pairs of them by hand from a
# Seurat object is exactly the kind of transcription a tool should do. They are computed once per
# dataset, in the same preprocessing step as the marker tables, and cached.
CELL_COUNTS_DEFAULTS = {
    "seurat_object": None,    # .rds / .rdata holding the object the clusters came from
    "cluster_column": "auto",  # metadata column with the cluster ids; auto = match the markers
    "object_var": "auto",     # variable name inside a save()-format .rdata
    "auto": True,             # compute during preprocessing when the cache is missing/stale
    "tsv": None,              # a precomputed cluster,ncells[,pct_cells] TSV — skips R entirely
}


def cell_counts_cfg(cfg: dict) -> dict:
    """Cell-count settings merged over the defaults.

    `seurat_object` falls back to the top-level `seurat_object` / `seurat_rdata` keys, so a config
    that already names the object doesn't have to repeat the path.
    """
    merged = dict(CELL_COUNTS_DEFAULTS)
    merged.update(cfg.get("cell_counts") or {})
    if not merged["seurat_object"]:
        merged["seurat_object"] = cfg.get("seurat_object") or cfg.get("seurat_rdata")
    return merged


def cell_counts_path(cfg: dict) -> Path:
    """Where the computed counts live — the override TSV if one is configured, else the cache."""
    cc = cell_counts_cfg(cfg)
    return resolve(cfg, cc["tsv"]) if cc["tsv"] else cache_dir(cfg) / CELL_COUNTS_FILE


def cell_counts_needed(cfg: dict) -> bool:
    """True when the counts could be computed but haven't been (or are older than the object)."""
    cc = cell_counts_cfg(cfg)
    if cc["tsv"] or not cc["seurat_object"]:
        return False            # supplied directly, or nothing to compute from
    obj = resolve(cfg, cc["seurat_object"])
    if not obj.exists():
        return False
    out = cache_dir(cfg) / CELL_COUNTS_FILE
    return not out.exists() or os.path.getmtime(out) < os.path.getmtime(obj)


def ensure_cell_counts(cfg: dict, force: bool = False,
                       expected: list[int] | None = None) -> Path | None:
    """Count cells per cluster from the Seurat object, unless that is already done.

    Returns the counts file, or None when the dataset has no object configured. Raises with the R
    script's own message on failure — a wrong cluster column is worth stopping for, since silently
    exporting counts from the wrong grouping would be worse than no counts at all.
    """
    cc = cell_counts_cfg(cfg)
    if cc["tsv"]:
        p = resolve(cfg, cc["tsv"])
        if not p.exists():
            raise RuntimeError(f"cell_counts.tsv does not exist: {p}")
        return p
    if not cc["seurat_object"]:
        return None
    obj = resolve(cfg, cc["seurat_object"])
    if not obj.exists():
        raise RuntimeError(f"Seurat object not found: {obj}")
    cd = cache_dir(cfg)
    out = cd / CELL_COUNTS_FILE
    if not force and not cell_counts_needed(cfg):
        return out if out.exists() else None
    script = Path(cfg["_root"]) / "preprocess" / "export_cell_counts.R"
    cmd = [rscript_bin(cfg), str(script), str(obj), str(cd),
           str(cc["cluster_column"] or "auto"), str(cc["object_var"] or "auto"),
           ",".join(str(c) for c in (expected or []))]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"export_cell_counts.R failed:\n{res.stderr}\n{res.stdout}")
    return out


def load_cell_counts(cfg: dict) -> dict[int, dict]:
    """{cluster: {'ncells': int, 'pct_cells': float}}, or {} if counts aren't available.

    Never raises: counts are a convenience that prefills two export columns, so a malformed file
    must not stop the app from opening. `pct_cells` is derived when the file omits it.
    """
    path = cell_counts_path(cfg)
    try:
        df = pd.read_csv(path, sep="\t")
    except (OSError, ValueError, pd.errors.ParserError):
        return {}
    if "cluster" not in df.columns or "ncells" not in df.columns:
        return {}
    df = df.dropna(subset=["cluster", "ncells"])
    if "pct_cells" not in df.columns:
        total = df["ncells"].sum()
        df["pct_cells"] = df["ncells"] / total if total else None
    out = {}
    for r in df.itertuples():
        try:
            out[int(r.cluster)] = {"ncells": int(r.ncells), "pct_cells": float(r.pct_cells)}
        except (TypeError, ValueError):
            continue
    return out


def cell_counts_meta(cfg: dict) -> dict:
    """What the counts were computed from (object, column, total) — for the UI caption. {} if
    unknown, e.g. when the counts came from a hand-supplied TSV."""
    try:
        return json.loads((cache_dir(cfg) / "cluster_cells.meta.json").read_text())
    except (OSError, json.JSONDecodeError):
        return {}


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
    "primary_model": "claude-opus-5",       # current Opus; thinking is on by default here
    "fallback_model": "claude-opus-4-8",    # used automatically if the primary declines/reroutes
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
