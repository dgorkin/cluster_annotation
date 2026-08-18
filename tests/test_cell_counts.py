#!/usr/bin/env python
"""Offline tests for counting cells per cluster. No API key, no network, no real Seurat object.

    ~/.conda/envs/cluster_annotation/bin/python tests/test_cell_counts.py

The R script is exercised for real against a tiny stand-in object: a plain list carrying the
`meta.data` / `active.ident` attributes that the script reads, which is all it looks at. Tests
needing R skip themselves when `Rscript` isn't on PATH.
"""
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import data as D  # noqa: E402

HAVE_R = shutil.which("Rscript") is not None

# 6 cells: clusters 0,0,0,1,1,2. `barcode` is unique per cell and `nCount` continuous, so both must
# be rejected as cluster columns; `stage` has too few values to match; idents mirror the clustering.
MAKE_OBJECT = r"""
args <- commandArgs(trailingOnly = TRUE)
cl <- factor(c(0, 0, 0, 1, 1, 2))
md <- data.frame(
  barcode         = paste0("cell", 1:6),
  nCount          = c(10, 20, 30, 40, 50, 60),
  stage           = factor(c("a", "a", "a", "b", "b", "b")),
  seurat_clusters = cl,
  clust_r1.2      = cl,
  stringsAsFactors = FALSE)
obj <- structure(list(), class = "Seurat", meta.data = md, active.ident = cl)
if (length(args) >= 2 && args[2] == "save") {
  second <- obj
  save(obj, second, file = args[1])       # save()-format archive with TWO candidate objects
} else {
  saveRDS(obj, args[1])
}
"""


def _make_object(path: Path, fmt: str = "rds") -> None:
    script = path.parent / "make_object.R"
    script.write_text(MAKE_OBJECT)
    res = subprocess.run(["Rscript", str(script), str(path), fmt],
                         capture_output=True, text=True)
    assert res.returncode == 0, res.stderr


def _cfg(tmp: str, **cell_counts) -> dict:
    return {"name": "t", "cache_dir": tmp, "_root": str(ROOT), "annotations_dir": tmp,
            "cell_counts": cell_counts}


def test_auto_detection_picks_the_column_that_matches_the_marker_clusters():
    """Auto mode verifies itself: the chosen column's values ARE the cluster ids. Name-guessing
    would happily pick `stage` on an object where the clustering column is called something else."""
    if not HAVE_R:
        print("      (no Rscript; skipping)")
        return
    with tempfile.TemporaryDirectory() as tmp:
        obj = Path(tmp) / "obj.rds"
        _make_object(obj)
        cfg = _cfg(tmp, seurat_object=str(obj))
        D.ensure_cell_counts(cfg, expected=[0, 1, 2])
        counts = D.load_cell_counts(cfg)
        assert counts[0] == {"ncells": 3, "pct_cells": 0.5}, counts
        assert counts[2]["ncells"] == 1 and abs(counts[2]["pct_cells"] - 1 / 6) < 1e-12
        meta = D.cell_counts_meta(cfg)
        assert meta["cluster_column"] == "seurat_clusters", meta
        assert meta["total_cells"] == 6


def test_a_named_column_is_used_as_given():
    if not HAVE_R:
        print("      (no Rscript; skipping)")
        return
    with tempfile.TemporaryDirectory() as tmp:
        obj = Path(tmp) / "obj.rds"
        _make_object(obj)
        cfg = _cfg(tmp, seurat_object=str(obj), cluster_column="clust_r1.2")
        D.ensure_cell_counts(cfg, expected=[0, 1, 2])
        assert D.cell_counts_meta(cfg)["cluster_column"] == "clust_r1.2"


def test_no_matching_column_fails_loudly_and_names_the_candidates():
    """Counts attached to the wrong clusters is the one outcome worth crashing over."""
    if not HAVE_R:
        print("      (no Rscript; skipping)")
        return
    with tempfile.TemporaryDirectory() as tmp:
        obj = Path(tmp) / "obj.rds"
        _make_object(obj)
        cfg = _cfg(tmp, seurat_object=str(obj))
        try:
            D.ensure_cell_counts(cfg, expected=[0, 1, 2, 3, 4])   # not this object's clustering
        except RuntimeError as e:
            assert "could not identify the cluster column" in str(e), str(e)
            assert "seurat_clusters" in str(e), "the error should list what it looked at"
        else:
            raise AssertionError("mismatched clusters were accepted")
        assert D.load_cell_counts(cfg) == {}, "a failed run must not leave counts behind"


def test_an_ambiguous_rdata_archive_asks_for_object_var():
    if not HAVE_R:
        print("      (no Rscript; skipping)")
        return
    with tempfile.TemporaryDirectory() as tmp:
        obj = Path(tmp) / "obj.rdata"
        _make_object(obj, fmt="save")
        cfg = _cfg(tmp, seurat_object=str(obj))
        try:
            D.ensure_cell_counts(cfg, expected=[0, 1, 2])
        except RuntimeError as e:
            assert "object_var" in str(e), str(e)
        else:
            raise AssertionError("an ambiguous archive was silently resolved")
        # Naming one resolves it — and proves the save()-format path works at all.
        cfg = _cfg(tmp, seurat_object=str(obj), object_var="second")
        D.ensure_cell_counts(cfg, expected=[0, 1, 2])
        assert D.load_cell_counts(cfg)[0]["ncells"] == 3
        assert D.cell_counts_meta(cfg)["object_var"] == "second"


def test_counts_are_reused_until_the_object_changes():
    if not HAVE_R:
        print("      (no Rscript; skipping)")
        return
    with tempfile.TemporaryDirectory() as tmp:
        obj = Path(tmp) / "obj.rds"
        _make_object(obj)
        cfg = _cfg(tmp, seurat_object=str(obj))
        assert D.cell_counts_needed(cfg)
        D.ensure_cell_counts(cfg, expected=[0, 1, 2])
        assert not D.cell_counts_needed(cfg), "a fresh cache must not trigger another read"
        obj.touch()
        assert D.cell_counts_needed(cfg), "an object newer than the cache must trigger a re-read"


def test_a_supplied_tsv_is_used_and_the_object_never_read():
    with tempfile.TemporaryDirectory() as tmp:
        supplied = Path(tmp) / "mine.tsv"
        supplied.write_text("cluster\tncells\n0\t10\n1\t30\n")
        cfg = _cfg(tmp, seurat_object="/does/not/exist.rds", tsv=str(supplied))
        assert not D.cell_counts_needed(cfg)
        assert D.ensure_cell_counts(cfg) == supplied
        counts = D.load_cell_counts(cfg)
        assert counts[1]["ncells"] == 30 and abs(counts[1]["pct_cells"] - 0.75) < 1e-12


def test_no_object_configured_is_simply_no_counts():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _cfg(tmp)
        assert D.ensure_cell_counts(cfg) is None
        assert D.load_cell_counts(cfg) == {}
        assert not D.cell_counts_needed(cfg)


def test_legacy_top_level_object_keys_are_honoured():
    """Existing configs name the object as `seurat_rdata:`; they must keep working untouched."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg = {"name": "t", "cache_dir": tmp, "_root": str(ROOT),
               "seurat_rdata": "inputs/object/whatever.rdata"}
        assert D.cell_counts_cfg(cfg)["seurat_object"] == "inputs/object/whatever.rdata"


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
