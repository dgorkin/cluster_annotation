#!/usr/bin/env python
"""Offline tests for primer mode. No API key, no network.

The point of primer mode is that citations are RESOLVED from a verified library rather than
generated per cluster, and that the primer is a cache-shareable prefix. Both are structural
properties, so they can be tested without spending anything.

    ~/.conda/envs/cluster_annotation/bin/python tests/test_primer.py
"""
import copy
import json
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import data as D      # noqa: E402
import insights as I  # noqa: E402

PRIMER = {
    "expected_cell_types": [
        {"name": "Microglia", "aliases": ["brain macrophages"],
         "canonical_markers": ["Cx3cr1", "P2ry12"], "canonical_motifs": ["SPI1"],
         "distinguishing_notes": "vs border-associated macrophages: lacks Mrc1",
         "reference_keys": ["R1"]},
        {"name": "Gliogenic radial glia", "aliases": [], "canonical_markers": ["Nfia", "Slc1a3"],
         "canonical_motifs": ["NFI"], "distinguishing_notes": "vs neurogenic RG: Nfia high",
         "reference_keys": ["R2"]},
    ],
    "references": [
        {"key": "R1", "citation": "Kierdorf et al. 2013", "pmid": "23334579", "url": None,
         "covers": "microglial origin"},
        {"key": "R2", "citation": "Piper et al. 2010", "pmid": None,
         "url": "https://doi.org/10.1523/x", "covers": "NFIA in gliogenesis"},
    ],
    "coverage_notes": "forebrain focus",
    "_meta": {"generated_at": "2026-08-17T00:00:00"},
}


def _fixture(tmp, mode="primer"):
    cfg = {"name": "t", "cache_dir": tmp, "_root": tmp,
           "biological_context": "E15.5 mouse forebrain",
           "ai_insights": {"research_mode": mode, "primary_model": "claude-opus-5",
                           "fallback_model": "claude-opus-4-8", "top_n_genes": 5,
                           "top_n_motifs": 5}}
    genes = pd.DataFrame({"feature": ["Cx3cr1", "P2ry12"], "cluster": [0, 0],
                          "avg_log2FC": [3.0, 2.0], "pct.1": [.9, .8], "pct.2": [.1, .1],
                          "delta_pct": [.8, .7], "p_val_adj": [1e-9, 1e-8]})
    motifs = pd.DataFrame({"feature": ["SPI1"], "cluster": [0], "avg_diff": [2.0],
                           "pct.1": [.9], "pct.2": [.1], "delta_pct": [.8], "p_val_adj": [1e-9]})
    return cfg, genes, motifs


def test_primer_is_in_system_so_clusters_share_one_cached_prefix():
    """If the primer sat after the per-cluster markers, every cluster would get a distinct prefix
    and the cache breakpoint would share nothing — defeating the point of building it once."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg, g, m = _fixture(tmp)
        s0, u0 = I.build_primer_cluster_prompt(cfg, 0, g, m, PRIMER)
        assert "REFERENCE LIBRARY" in s0 and "REFERENCE LIBRARY" not in u0
        assert "TOP MARKER GENES" in u0 and "TOP MARKER GENES" not in s0
        # The system prompt must be byte-identical across clusters for the prefix to be shared.
        s1, _ = I.build_primer_cluster_prompt(cfg, 1, g, m, PRIMER)
        assert s0 == s1, "system prompt varies by cluster; the shared prefix would not cache"


def test_cited_keys_resolve_to_verified_references():
    refs, unknown = I._resolve_reference_keys(PRIMER, ["R2", "R1"])
    assert [r["pmid"] for r in refs] == [None, "23334579"]
    assert refs[0]["url"] == "https://doi.org/10.1523/x"
    assert refs[0]["supports"] == "NFIA in gliogenesis"
    assert unknown == []


def test_invented_keys_cannot_become_citations():
    """The structural guarantee: a key not in the library yields no reference at all."""
    refs, unknown = I._resolve_reference_keys(PRIMER, ["R1", "R99", "Smith et al. 2020"])
    assert len(refs) == 1 and refs[0]["pmid"] == "23334579"
    assert unknown == ["R99", "Smith et al. 2020"]


def test_duplicate_keys_are_collapsed():
    refs, _ = I._resolve_reference_keys(PRIMER, ["R1", "R1", "R1"])
    assert len(refs) == 1


def test_primer_stamp_changes_when_the_primer_changes():
    assert I.primer_stamp(None) == "none"
    base = I.primer_stamp(PRIMER)
    rebuilt = copy.deepcopy(PRIMER)
    rebuilt["_meta"]["generated_at"] = "2026-08-18T00:00:00"
    assert I.primer_stamp(rebuilt) != base
    extra = copy.deepcopy(PRIMER)
    extra["references"].append({"key": "R3", "citation": "x", "covers": "y"})
    assert I.primer_stamp(extra) != base


def test_mode_and_primer_are_part_of_the_cache_key():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, g, m = _fixture(tmp)
        pc = copy.deepcopy(cfg); pc["ai_insights"]["research_mode"] = "per_cluster"
        h_primer = I.inputs_hash(cfg, 0, g, m, primer=PRIMER)
        h_percluster = I.inputs_hash(pc, 0, g, m)
        assert h_primer != h_percluster, "a primer annotation is not interchangeable with a per-cluster one"
        rebuilt = copy.deepcopy(PRIMER)
        rebuilt["_meta"]["generated_at"] = "2026-08-18T00:00:00"
        assert I.inputs_hash(cfg, 0, g, m, primer=rebuilt) != h_primer


def test_per_cluster_hash_is_unchanged_so_migrated_artifacts_stay_current():
    """Adding primer mode must not invalidate the artifacts already migrated to the new key."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg, g, m = _fixture(tmp, mode="per_cluster")
        assert I.inputs_hash(cfg, 0, g, m) == I._hash_payload(cfg, 0, g, m)


def test_primer_mode_without_a_primer_fails_loudly():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, g, m = _fixture(tmp)
        try:
            I.generate_one(cfg, 0, g, m, force=True)
        except RuntimeError as e:
            assert "primer" in str(e).lower(), e
        else:
            raise AssertionError("should refuse to annotate with no primer")


def test_primer_roundtrips_to_disk():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, _, _ = _fixture(tmp)
        assert I.load_primer(cfg) is None
        I.primer_path(cfg).write_text(json.dumps(PRIMER))
        got = I.load_primer(cfg)
        assert got and len(got["references"]) == 2
        assert I.primer_stamp(got) == I.primer_stamp(PRIMER)


def test_digest_lists_every_reference_key_the_model_may_cite():
    digest = I._primer_digest(PRIMER)
    for r in PRIMER["references"]:
        assert f"[{r['key']}]" in digest
    assert "PMID 23334579" in digest
    assert "https://doi.org/10.1523/x" in digest
    assert "do not cite anything else" in digest.lower()


def test_adaptive_thinking_capability_is_checked_per_model():
    """Regression: _structured_call sent thinking={"type":"adaptive"} unconditionally, so using it
    with the Haiku structuring model returned 400 'adaptive thinking is not supported on this
    model' — after the primer's paid web-search research had already succeeded."""
    assert I.supports_adaptive_thinking("claude-opus-5")
    assert I.supports_adaptive_thinking("claude-opus-4-8-20260115"), "dated snapshots too"
    assert not I.supports_adaptive_thinking("claude-haiku-4-5")
    assert not I.supports_adaptive_thinking("claude-haiku-4-5-20251001")
    assert not I.supports_adaptive_thinking(None)
    # Every model the app can be configured with must be classified, not guessed at.
    for m in D.KNOWN_MODELS:
        assert isinstance(I.supports_adaptive_thinking(m), bool)


def test_stale_reannotation_is_detected():
    """A reannotation is written against one original. Regenerate that original and the
    reannotation still sits on disk describing a version that no longer exists."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg, _, _ = _fixture(tmp)
        I.insight_dir(cfg)

        def write(path, at):
            path.write_text(json.dumps({"primary_identity": "x", "_meta": {"generated_at": at}}))

        write(I.insight_path(cfg, 0), "2026-08-17T10:00:00")
        assert I.reannotation_is_stale(cfg, 0) is False, "no reannotation -> not stale"

        write(I.reannotation_path(cfg, 0), "2026-08-17T11:00:00")
        assert I.reannotation_is_stale(cfg, 0) is False, "written after its original -> current"

        write(I.insight_path(cfg, 0), "2026-08-17T12:00:00")
        assert I.reannotation_is_stale(cfg, 0) is True, "original regenerated since -> history"

        write(I.reannotation_path(cfg, 1), "2026-08-17T11:00:00")
        assert I.reannotation_is_stale(cfg, 1) is True, "no surviving original -> stale"


def test_research_mode_default_is_primer():
    assert D.AI_INSIGHTS_DEFAULTS["research_mode"] == "primer"


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
