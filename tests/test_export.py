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

import data as D     # noqa: E402
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


def test_cell_counts_prefill_ncells_and_pct_when_the_user_typed_none():
    """The counts are measurements: exporting blanks when we know them is just lost information."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg, genes, motifs = _fixture(tmp)
        (Path(tmp) / "t").mkdir(parents=True, exist_ok=True)
        D.cell_counts_path(cfg).write_text("cluster\tncells\tpct_cells\n0\t2834\t0.0536\n")
        S.set_annotation(cfg, 0, {"annot_type": "x"})
        row = E.build_table(cfg, [0], genes, motifs).loc[0]
        assert row["ncells"] == 2834 and abs(row["%cells"] - 0.0536) < 1e-12

        # A typed value is a decision and must win over the computed one.
        S.set_annotation(cfg, 0, {"ncells": 99, "pct_cells": 0.5})
        row = E.build_table(cfg, [0], genes, motifs).loc[0]
        assert row["ncells"] == 99 and row["%cells"] == 0.5


def test_cell_counts_derives_pct_when_the_file_only_has_counts():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, _, _ = _fixture(tmp)
        D.cell_counts_path(cfg).write_text("cluster\tncells\n0\t30\n1\t70\n")
        got = D.load_cell_counts(cfg)
        assert got[0]["ncells"] == 30 and abs(got[1]["pct_cells"] - 0.7) < 1e-12


def test_missing_or_broken_cell_counts_are_not_fatal():
    """Counts are a convenience; a malformed file must not stop the app or the export."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg, genes, motifs = _fixture(tmp)
        assert D.load_cell_counts(cfg) == {}          # nothing written yet
        D.cell_counts_path(cfg).write_text("not\ta\tcounts\tfile\n")
        assert D.load_cell_counts(cfg) == {}
        S.set_annotation(cfg, 0, {"annot_type": "x"})
        assert pd.isna(E.build_table(cfg, [0], genes, motifs).loc[0, "ncells"])


def test_tsv_export_matches_the_xlsx_content_and_header():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, genes, motifs = _fixture(tmp)
        S.set_annotation(cfg, 0, {"annot_type": "Forebrain progenitors", "annot_abbrev": "Fb",
                                  "ncells": 3095, "pct_cells": 0.0508393836854036})
        S.set_markers(cfg, 0, "gene", {"Six3", "Lhx2"}, {"Six3", "Lhx2"})
        path, n = E.write_tsv(cfg, [0], genes, motifs, out_path=str(Path(tmp) / "out.tsv"))
        assert path.exists() and n == 1
        back = pd.read_csv(path, sep="\t")
        assert list(back.columns) == E.EXPORT_COLUMNS
        assert back.loc[0, "key_marker_genes"] == "Six3,Lhx2"
        assert abs(back.loc[0, "%cells"] - 0.0508393836854036) < 1e-12
        # Blank, not the string "nan" — read.delim in R has to see NA.
        assert "nan" not in path.read_text()


def test_tsv_path_prefers_explicit_then_config_then_a_default():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, _, _ = _fixture(tmp)
        assert E.tsv_path(cfg).suffix == ".tsv"
        assert E.tsv_path(cfg).parent == Path(tmp)
        cfg["output_tsv"] = "{name}_annotations.tsv"
        assert E.tsv_path(cfg).name == "t_annotations.tsv"
        assert E.tsv_path(cfg, str(Path(tmp) / "explicit.tsv")).name == "explicit.tsv"


def test_tsv_rewrites_the_same_path_rather_than_dating_a_new_one():
    """The point of output_tsv is one file a downstream script keeps reading."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg, genes, motifs = _fixture(tmp)
        cfg["output_tsv"] = "{name}_annotations.tsv"
        S.set_annotation(cfg, 0, {"annot_type": "first"})
        p1, _ = E.write_tsv(cfg, [0], genes, motifs)
        S.set_annotation(cfg, 0, {"annot_type": "second"})
        p2, _ = E.write_tsv(cfg, [0], genes, motifs)
        assert p1 == p2
        assert pd.read_csv(p2, sep="\t").loc[0, "annot_type"] == "second"


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
