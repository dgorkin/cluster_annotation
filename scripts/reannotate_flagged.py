#!/usr/bin/env python
"""Reannotate clusters the cohort review flagged, keeping each original.

Run from the project root with the cluster_annotation conda env:
    python scripts/reannotate_flagged.py [config] [workers] [--clusters 1,2,3]

--clusters restricts the run to an explicit list instead of every flagged cluster. That is usually
what you want: reannotation escalates a cluster to full web-search research, which helps where the
IDENTITY is uncertain but cannot resolve a "these clusters look like the same population" flag —
that is a re-clustering decision, and re-clustering is out of scope here. Spending on the
redundancy-flagged clusters buys nothing.

Requires a cohort review to exist (run scripts/generate_all.py first, or the Cohort review tab).
Writes cluster_<id>.reannotation.json per cluster — originals (cluster_<id>.json) are untouched.
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
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    config = args[0] if args else str(ROOT / "config" / "mydataset.yaml")
    workers = int(args[1]) if len(args) > 1 else None
    explicit = None
    if "--clusters" in sys.argv:
        raw = sys.argv[sys.argv.index("--clusters") + 1]
        explicit = [int(x) for x in raw.replace(" ", "").split(",") if x != ""]

    cfg = D.load_config(_resolve_config(config))
    if workers:
        cfg.setdefault("ai_insights", {})["max_workers"] = workers
    if not I.api_key_available(cfg):
        print("ERROR: no Anthropic API key found — cannot reannotate.", flush=True)
        sys.exit(1)

    review = I.load_cohort_review(cfg)
    if review is None:
        print("ERROR: no cohort review found — run scripts/generate_all.py first.", flush=True)
        sys.exit(1)

    flagged = I.flagged_clusters(review)
    if explicit is not None:
        unknown = sorted(set(explicit) - set(D.clusters(cfg)))
        if unknown:
            print(f"ERROR: not clusters in this dataset: {unknown}", flush=True)
            sys.exit(1)
        not_flagged = sorted(set(explicit) - set(flagged))
        targets = explicit
        print(f"reannotating {len(targets)} explicitly requested cluster(s): {targets}"
              + (f"  (note: {not_flagged} carry no cohort flag)" if not_flagged else ""),
              flush=True)
    else:
        targets = flagged
    if not targets:
        print("No flagged clusters — nothing to reannotate.", flush=True)
        return

    # Build the marker cache from the source RDS if it's missing/stale (the app does this on load;
    # the script must too, so it works on a fresh or cleared cache rather than failing in load_markers).
    print("ensuring preprocessed marker cache…", flush=True)
    D.ensure_preprocessed(cfg)

    genes = D.load_markers(cfg, "gene")
    motifs = D.load_markers(cfg, "motif")
    aic = D.ai_insights_cfg(cfg)
    print(f"dataset={cfg['name']} flagged={len(targets)} clusters={targets} "
          f"primary={aic['primary_model']} fallback={aic['fallback_model']} "
          f"workers={aic['max_workers']}", flush=True)

    t0 = time.time()
    done = [0]
    spent = [0.0]

    def cb(c, err):
        done[0] += 1
        tag = "ok" if not err else f"FAIL: {err[:120]}"
        # Only price a success: on failure the just-written file doesn't exist, and reading the
        # cluster's previous reannotation would bill this run for an earlier one.
        cost = I.record_cost(I.load_reannotation(cfg, c)) if not err else 0.0
        spent[0] += cost
        print(f"[{done[0]}/{len(targets)}] cluster {c}: {tag}  "
              f"(~${cost:.3f}, run total ~${spent[0]:.2f}, +{time.time()-t0:.0f}s)", flush=True)

    errors = I.reannotate_flagged(cfg, targets, genes, motifs, review, progress_cb=cb)
    print(f"reannotation finished in {time.time()-t0:.0f}s — {len(errors)} failure(s)"
          + (f": {errors}" if errors else ""), flush=True)
    print(f"this run: ~${spent[0]:.2f} for {done[0] - len(errors)} reannotation(s)"
          + (" (failed clusters' partial spend is not recorded)" if errors else ""), flush=True)
    print(I.format_dataset_cost(cfg, D.clusters(cfg)), flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
