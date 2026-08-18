#!/usr/bin/env python
"""Generate AI insights for every cluster, then run the cohort review.

Run from the project root with the cluster_annotation conda env:
    ~/.conda/envs/cluster_annotation/bin/python scripts/generate_all.py [config] [workers]

Defaults: config = config/mydataset.yaml, workers = ai_insights.max_workers from the config.
Idempotent: clusters whose cached insight is already current (same markers + model) are skipped.
Prints one line per completed cluster so it can be tailed when run in the background.
"""
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import data as D
import insights as I


def _resolve_config(path: str) -> str:
    """Fail with instructions rather than a stack trace when no config is given/found."""
    if Path(path).exists():
        return path
    available = sorted((ROOT / "config").glob("*.yaml"))
    available = [p for p in available if not p.name.endswith(".template.yaml")]
    print(f"ERROR: no such config: {path}", flush=True)
    if available:
        print("  available: " + ", ".join(str(p) for p in available), flush=True)
    else:
        print("  no dataset configs yet. Create one:\n"
              "    cp config/dataset.template.yaml config/mydataset.yaml", flush=True)
    sys.exit(2)


def main():
    config = sys.argv[1] if len(sys.argv) > 1 else str(ROOT / "config" / "mydataset.yaml")
    workers = int(sys.argv[2]) if len(sys.argv) > 2 else None

    cfg = D.load_config(_resolve_config(config))
    if workers:
        cfg.setdefault("ai_insights", {})["max_workers"] = workers
    if not I.api_key_available(cfg):
        print("ERROR: no Anthropic API key found — cannot generate.", flush=True)
        sys.exit(1)

    # Build the marker cache from the source RDS if it's missing/stale (the app does this on load;
    # the script must too, so it works on a fresh or cleared cache rather than failing in load_markers).
    print("ensuring preprocessed marker cache…", flush=True)
    D.ensure_preprocessed(cfg)

    aic = D.ai_insights_cfg(cfg)
    genes = D.load_markers(cfg, "gene")
    motifs = D.load_markers(cfg, "motif")
    clusters = D.clusters(cfg)

    # In primer mode the dataset's reference sheet is the prerequisite for everything else: build
    # it first (one web-grounded call) so all clusters are judged against the same vetted library.
    primer = None
    if aic["research_mode"] == "primer":
        primer = I.load_primer(cfg)
        if primer is None:
            print("building the reference primer (one web-grounded call)…", flush=True)
            t = time.time()
            primer = I.build_primer(cfg)
            print(f"primer ok in {time.time()-t:.0f}s: "
                  f"{len(primer.get('expected_cell_types') or [])} cell types, "
                  f"{len(primer.get('references') or [])} references "
                  f"(~${I.record_cost(primer):.3f})", flush=True)
        else:
            print(f"reusing the reference primer from "
                  f"{(primer.get('_meta') or {}).get('generated_at')} "
                  f"({len(primer.get('expected_cell_types') or [])} cell types, "
                  f"{len(primer.get('references') or [])} references)", flush=True)

    todo = I.clusters_needing_insight(cfg, clusters, genes, motifs, force=False, primer=primer)
    print(f"dataset={cfg['name']} clusters={len(clusters)} to_generate={len(todo)} "
          f"mode={aic['research_mode']} primary={aic['primary_model']} "
          f"fallback={aic['fallback_model']} workers={aic['max_workers']}", flush=True)

    t0 = time.time()
    done = [0]
    spent = [0.0]

    def cb(c, err):
        done[0] += 1
        tag = "ok" if not err else f"FAIL: {err[:120]}"
        # Only price a success: on failure nothing was written, and reading the cluster's
        # previous insight would bill this run for an earlier one.
        cost = I.record_cost(I.load_insight(cfg, c)) if not err else 0.0
        spent[0] += cost
        print(f"[{done[0]}/{len(todo)}] cluster {c}: {tag}  "
              f"(~${cost:.3f}, run total ~${spent[0]:.2f}, +{time.time()-t0:.0f}s)", flush=True)

    errors = I.generate_all(cfg, todo, genes, motifs, progress_cb=cb, force=True) if todo else []
    print(f"generation finished in {time.time()-t0:.0f}s — {len(errors)} failure(s)"
          + (f": {errors}" if errors else ""), flush=True)
    if done[0]:
        ok = done[0] - len(errors)
        print(f"this run: ~${spent[0]:.2f} for {ok} cluster(s)"
              + (f", ~${spent[0]/ok:.3f}/cluster" if ok else "")
              + (" (failed clusters' partial spend is not recorded)" if errors else ""), flush=True)

    print("running cohort review over all cached clusters…", flush=True)
    try:
        rev = I.cohort_review(cfg, clusters, force=True)
        meta = rev["_meta"]
        print(f"cohort review ok: reviewed {meta['n_clusters']} clusters by {meta['model_used']}; "
              f"{len(rev.get('flags') or [])} flag(s); flagged clusters = {I.flagged_clusters(rev)} "
              f"(~${I.record_cost(rev):.3f})", flush=True)
        print(f"missing_expected = {rev.get('missing_expected')}", flush=True)
    except Exception as e:
        print(f"cohort review FAILED: {e}", flush=True)

    print(I.format_dataset_cost(cfg, clusters), flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
