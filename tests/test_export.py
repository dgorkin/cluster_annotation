#!/usr/bin/env python
"""Offline tests for the annotation store and the xlsx export. No API key, no network.

    ~/.conda/envs/cluster_annotation/bin/python tests/test_export.py
"""
import json
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import export as E    # noqa: E402
import insights as I  # noqa: E402
import store as S     # noqa: E402


def _fixture(tmp):
    cfg = {"name": "t", "cache_dir": tmp, "_root": tmp, "annotations_dir": tmp}
    genes = pd.DataFrame({"feature": ["Six3", "Lhx2", "Emx2"], "cluster": [0, 0, 0],
                          "avg_log2FC": [3.0, 2.0, 1.0]})
    motifs = pd.DataFrame({"feature": ["LHX", "EMX"], "cluster": [0, 0],
                           "avg_diff": [2.0, 1.0]})
    return cfg, genes, motifs


def test_export_columns_match_the_template_exactly():
    """Downstream scripts key on the template's header, including the odd '%cells' spelling."""
    wb_path = ROOT / "ref_materials" / "E9_annotations_06092026.xlsx"
    if not wb_path.exists():
        print("      (template not present; skipping)")
        return
    import openpyxl
    ws = openpyxl.load_workbook(wb_path, data_only=True).worksheets[0]
    template = [c.value for c in next(ws.iter_rows(max_row=1))]
    assert E.EXPORT_COLUMNS == template, f"\n  ours: {E.EXPORT_COLUMNS}\n  tmpl: {template}"


def test_annotation_roundtrip_and_field_validation():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, _, _ = _fixture(tmp)
        assert S.get_annotation(cfg, 0) == {}
        S.set_annotation(cfg, 0, {"annot_type": "Forebrain progenitors", "annot_abbrev": "Fb",
                                  "annot_order": 1, "pct_cells": 0.0508})
        a = S.get_annotation(cfg, 0)
        assert a["annot_type"] == "Forebrain progenitors" and a["annot_order"] == 1
        assert a["reviewed"] == 0, "saving a field must not mark a cluster reviewed"

        # Upsert must update in place, not insert a second row, and must not clear other columns.
        S.set_annotation(cfg, 0, {"annot_origin": "Ectoderm"}, reviewed=True)
        a = S.get_annotation(cfg, 0)
        assert a["annot_origin"] == "Ectoderm" and a["annot_type"] == "Forebrain progenitors"
        assert a["reviewed"] == 1
        assert len(S.all_annotations(cfg)) == 1

        # A typo'd column must fail loudly rather than vanish.
        try:
            S.set_annotation(cfg, 0, {"annot_typ": "oops"})
        except ValueError as e:
            assert "annot_typ" in str(e)
        else:
            raise AssertionError("unknown field was accepted")


def test_starred_markers_export_in_effect_size_order():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, genes, motifs = _fixture(tmp)
        # Star them in an order that is NOT the effect-size order.
        S.set_markers(cfg, 0, "gene", {"Six3", "Lhx2", "Emx2"}, {"Emx2", "Six3"})
        S.set_annotation(cfg, 0, {"annot_type": "Forebrain progenitors"})
        df = E.build_table(cfg, [0], genes, motifs)
        assert df.loc[0, "key_marker_genes"] == "Six3,Emx2", df.loc[0, "key_marker_genes"]


def test_starred_marker_missing_from_current_table_is_still_exported():
    """Markers can be re-exported from the RDS; a prior pick must not silently disappear."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg, genes, motifs = _fixture(tmp)
        S.set_markers(cfg, 0, "gene", {"Gone"}, {"Gone"})   # not in `genes`
        S.set_markers(cfg, 0, "gene", {"Six3"}, {"Six3"})
        S.set_annotation(cfg, 0, {"annot_type": "x"})
        got = E.build_table(cfg, [0], genes, motifs).loc[0, "key_marker_genes"]
        assert got == "Six3,Gone", got


def test_only_annotated_clusters_are_exported_by_default():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, genes, motifs = _fixture(tmp)
        S.set_annotation(cfg, 1, {"annot_type": "Midbrain progenitors"})
        S.set_annotation(cfg, 2, {"annot_trajectory": "-"})   # no label -> not exported
        assert list(E.build_table(cfg, [0, 1, 2], genes, motifs)["cluster"]) == [1]
        assert list(E.build_table(cfg, [0, 1, 2], genes, motifs,
                                 only_annotated=False)["cluster"]) == [0, 1, 2]


def test_rows_come_out_in_figure_order():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, genes, motifs = _fixture(tmp)
        S.set_annotation(cfg, 5, {"annot_type": "c", "annot_order": 1})
        S.set_annotation(cfg, 1, {"annot_type": "a", "annot_order": 3})
        S.set_annotation(cfg, 3, {"annot_type": "b"})          # unnumbered -> last
        df = E.build_table(cfg, [1, 3, 5], genes, motifs)
        assert list(df["cluster"]) == [5, 1, 3], list(df["cluster"])


def test_comments_are_joined_into_one_cell():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, genes, motifs = _fixture(tmp)
        S.add_comment(cfg, 0, "first note")
        S.add_comment(cfg, 0, "second\nnote with newline")
        S.set_annotation(cfg, 0, {"annot_type": "x"})
        got = E.build_table(cfg, [0], genes, motifs).loc[0, "comments"]
        assert got == "first note | second note with newline", got


def test_prefill_from_insight_formats_pmids_and_skips_unanchored():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, _, _ = _fixture(tmp)
        I.insight_dir(cfg)
        I.insight_path(cfg, 0).write_text(json.dumps({
            "primary_identity": "Embryonic microglia",
            "references": [{"citation": "a", "pmid": "24316888"},
                           {"citation": "b", "pmid": None, "url": "https://x"},
                           {"citation": "c", "pmid": "24316888"}],   # duplicate
            "_meta": {}}))
        got = E.suggest_from_insight(cfg, 0)
        assert got["annot_type"] == "Embryonic microglia"
        assert got["refs"] == "PMIDs:24316888", got
        assert E.suggest_from_insight(cfg, 99) == {}, "absent insight must yield no prefill"


def test_prefill_prefers_the_reannotation():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, _, _ = _fixture(tmp)
        I.insight_dir(cfg)
        I.insight_path(cfg, 0).write_text(json.dumps(
            {"primary_identity": "original call", "_meta": {}}))
        I.reannotation_path(cfg, 0).write_text(json.dumps(
            {"primary_identity": "revised call", "_meta": {}}))
        assert E.suggest_from_insight(cfg, 0)["annot_type"] == "revised call"
        assert E.suggest_from_insight(cfg, 0, prefer_reannotation=False)["annot_type"] \
            == "original call"


def test_xlsx_writes_and_reads_back_with_the_template_header():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, genes, motifs = _fixture(tmp)
        S.set_annotation(cfg, 0, {"annot_type": "Forebrain progenitors", "annot_abbrev": "Fb",
                                  "annot_origin": "Ectoderm", "annot_order": 1,
                                  "refs": "PMIDs:26371318", "ncells": 3095,
                                  "pct_cells": 0.0508393836854036})
        S.set_markers(cfg, 0, "gene", {"Six3", "Lhx2"}, {"Six3", "Lhx2"})
        path, n = E.write_xlsx(cfg, [0], genes, motifs, out_dir=tmp)
        assert path.exists() and n == 1
        back = pd.read_excel(path)
        assert list(back.columns) == E.EXPORT_COLUMNS
        assert back.loc[0, "annot_abbrev"] == "Fb"
        assert back.loc[0, "key_marker_genes"] == "Six3,Lhx2"
        assert abs(back.loc[0, "%cells"] - 0.0508393836854036) < 1e-12
        assert back.loc[0, "ncells"] == 3095


def test_export_honours_annotations_dir_instead_of_a_hardcoded_path():
    """Regression: export_path hardcoded "annotations", so a config override was ignored and the
    file landed in the real project dir. Caught by an AppTest run against a throwaway store."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg, genes, motifs = _fixture(tmp)
        cfg["annotations_dir"] = tmp
        cfg.pop("_root", None)
        cfg["_root"] = tmp
        assert E.export_path(cfg).parent == Path(tmp), E.export_path(cfg)
        S.set_annotation(cfg, 0, {"annot_type": "x"})
        path, _ = E.write_xlsx(cfg, [0], genes, motifs)
        assert path.parent == Path(tmp), path


def test_export_of_an_empty_store_is_a_valid_empty_sheet():
    """Exporting before annotating anything should produce headers, not crash."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg, genes, motifs = _fixture(tmp)
        path, n = E.write_xlsx(cfg, [0, 1], genes, motifs, out_dir=tmp)
        assert n == 0
        assert list(pd.read_excel(path).columns) == E.EXPORT_COLUMNS


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL  {t.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"  ERROR {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
