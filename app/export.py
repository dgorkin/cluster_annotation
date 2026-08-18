"""Export the annotation of record to an xlsx matching the E9_annotations template.

The 13 columns come from three places, deliberately not duplicated in storage:

  annot_* / reference_omg_type / refs / ncells / %cells   the `annotations` table (user's decision)
  key_marker_genes / select_marker_motifs                 starred markers (the ★ checkboxes)
  comments                                                the `comments` table, joined

Anything the user hasn't filled in is left blank rather than guessed — the AI insight is offered
as a *prefill* in the UI, so what lands here is what the user accepted, not what a model proposed.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

import data as D
import insights as I
import store as S

# Column order and spelling of the template (the annotation template spreadsheet).
# '%cells' keeps its awkward name because downstream scripts key on the template's header.
EXPORT_COLUMNS = ["cluster", "annot_origin", "annot_trajectory", "annot_order", "annot_abbrev",
                  "annot_type", "reference_omg_type", "key_marker_genes", "select_marker_motifs",
                  "refs", "ncells", "%cells", "comments"]


def _sorted_stars(cfg: dict, cluster: int, kind: str, markers: pd.DataFrame) -> str:
    """Starred features for a cluster, ordered by effect size rather than by when they were
    clicked, so the exported list reads like the marker table the user was looking at."""
    stars = S.selected_markers(cfg, cluster, kind)
    if not stars:
        return ""
    eff = D.EFFECT_COL[kind]
    sub = markers[(markers["cluster"] == cluster) & (markers["feature"].isin(stars))]
    ordered = sub.sort_values(eff, ascending=False)["feature"].tolist()
    # Anything starred but absent from the current marker table (markers re-exported since) still
    # belongs in the output — dropping the user's pick silently would be worse than odd ordering.
    ordered += sorted(stars - set(ordered))
    return ",".join(ordered)


def _comments_text(cfg: dict, cluster: int) -> str:
    rows = S.list_comments(cfg, cluster)
    return " | ".join(r["text"].replace("\n", " ").strip() for r in rows if r["text"].strip())


def build_table(cfg: dict, clusters, genes: pd.DataFrame, motifs: pd.DataFrame,
                only_annotated: bool = True) -> pd.DataFrame:
    """Assemble the export table.

    only_annotated=True emits just the clusters that have a label, which is the normal case: a
    spreadsheet padded with 20 blank rows is harder to review than one with 14 real ones. Pass
    False to get every cluster, e.g. to see what is still outstanding.
    """
    annotations = S.all_annotations(cfg)
    labelled = S.annotated_clusters(cfg)
    rows = []
    for c in clusters:
        if only_annotated and c not in labelled:
            continue
        a = annotations.get(c, {})
        rows.append({
            "cluster": c,
            "annot_origin": a.get("annot_origin") or "",
            "annot_trajectory": a.get("annot_trajectory") or "",
            "annot_order": a.get("annot_order"),
            "annot_abbrev": a.get("annot_abbrev") or "",
            "annot_type": a.get("annot_type") or "",
            "reference_omg_type": a.get("reference_omg_type") or "",
            "key_marker_genes": _sorted_stars(cfg, c, "gene", genes),
            "select_marker_motifs": _sorted_stars(cfg, c, "motif", motifs),
            "refs": a.get("refs") or "",
            "ncells": a.get("ncells"),
            "%cells": a.get("pct_cells"),
            "comments": _comments_text(cfg, c),
        })
    df = pd.DataFrame(rows, columns=EXPORT_COLUMNS)
    # Order by the figure order the user assigned, falling back to cluster id for unnumbered rows,
    # so the sheet comes out in the order the figures will use.
    if not df.empty:
        df = (df.assign(_ord=df["annot_order"].fillna(1e9))
                .sort_values(["_ord", "cluster"], kind="mergesort")
                .drop(columns="_ord").reset_index(drop=True))
    return df


def suggest_from_insight(cfg: dict, cluster: int, prefer_reannotation: bool = True) -> dict:
    """Prefill values proposed by the AI annotation, for the user to accept or edit.

    Returns only the fields the insight can speak to; the rest stay for the user. PMIDs are
    formatted as the template does ("PMIDs:123,456") and only real identifiers are included —
    a reference without one contributes nothing rather than an empty placeholder.
    """
    rec = (I.load_reannotation(cfg, cluster) if prefer_reannotation else None) \
        or I.load_insight(cfg, cluster)
    if not rec:
        return {}
    pmids = [str(r.get("pmid")).strip() for r in (rec.get("references") or []) if r.get("pmid")]
    out = {"annot_type": (rec.get("primary_identity") or "").strip()}
    if pmids:
        out["refs"] = "PMIDs:" + ",".join(dict.fromkeys(pmids))
    return out


def export_path(cfg: dict, out_dir: Optional[str] = None) -> Path:
    # Default beside the annotation DB, honouring the same `annotations_dir` key the store uses —
    # a hardcoded "annotations" here silently ignored the config override.
    d = D.resolve(cfg, out_dir or cfg.get("annotations_dir", "annotations"))
    d.mkdir(parents=True, exist_ok=True)
    stamp = pd.Timestamp.now().strftime("%Y%m%d")
    return d / f"{cfg['name']}_annotations_{stamp}.xlsx"


def write_xlsx(cfg: dict, clusters, genes: pd.DataFrame, motifs: pd.DataFrame,
               out_dir: Optional[str] = None, only_annotated: bool = True) -> tuple[Path, int]:
    """Write the export and return (path, row_count). Sheet name matches the template's."""
    df = build_table(cfg, clusters, genes, motifs, only_annotated=only_annotated)
    path = export_path(cfg, out_dir)
    sheet = (cfg.get("name") or "annotations")[:31]  # xlsx sheet names cap at 31 chars
    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        df.to_excel(xl, sheet_name=sheet, index=False)
        ws = xl.sheets[sheet]
        for i, col in enumerate(df.columns, 1):
            widest = max([len(str(col))] + [len(str(v)) for v in df[col].head(200) if v is not None])
            ws.column_dimensions[ws.cell(row=1, column=i).column_letter].width = min(widest + 2, 60)
        ws.freeze_panes = "A2"
    return path, len(df)
