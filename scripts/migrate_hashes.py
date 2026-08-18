#!/usr/bin/env python
"""One-time migration: make cached AI-insight cache keys independent of model choice.

The original inputs_hash folded the configured primary/fallback model names in, so switching
model marked every cached annotation stale even though no marker had changed. This rewrites
those keys to the marker-only form, leaving provenance to _meta.model_used.

Dry run (default), then apply:
    ~/.conda/envs/cluster_annotation/bin/python scripts/migrate_hashes.py config/mydataset.yaml
    ~/.conda/envs/cluster_annotation/bin/python scripts/migrate_hashes.py config/mydataset.yaml --apply

RUN THIS BEFORE CHANGING primary_model — the legacy key embeds the model names, so after the
config moves the artifacts can no longer be recognised and will read as stale.
Each rewritten file is backed up to .cache/<name>/ai_insights/_backups/ first.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import data as D      # noqa: E402
import insights as I  # noqa: E402


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
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    apply = "--apply" in sys.argv
    config = args[0] if args else str(ROOT / "config" / "mydataset.yaml")

    cfg = D.load_config(_resolve_config(config))
    D.ensure_preprocessed(cfg)
    genes, motifs = D.load_markers(cfg, "gene"), D.load_markers(cfg, "motif")
    clusters = D.clusters(cfg)
    aic = D.ai_insights_cfg(cfg)

    print(f"dataset={cfg['name']} clusters={len(clusters)} "
          f"models={aic['primary_model']}/{aic['fallback_model']} "
          f"mode={'APPLY' if apply else 'dry run'}")

    res = I.migrate_inputs_hashes(cfg, clusters, genes, motifs, apply=apply)
    for bucket in ("migrated", "already_current", "stale", "missing"):
        ids = res[bucket]
        if ids:
            print(f"  {bucket:<16} {len(ids):>3}  {ids}")

    if res["stale"]:
        print("\n  NOTE: 'stale' artifacts match neither key — their markers or biological_context "
              "changed, or they were generated under different models. They are left untouched "
              "and will regenerate on the next run.")
    if not apply and res["migrated"]:
        print(f"\n  dry run — re-run with --apply to rewrite {len(res['migrated'])} file(s)")
    elif apply and res["migrated"]:
        left = I.clusters_needing_insight(cfg, clusters, genes, motifs)
        print(f"\n  migrated {len(res['migrated'])} file(s); "
              f"{len(left)} of {len(clusters)} cluster(s) now need regeneration")
    print("DONE")


if __name__ == "__main__":
    main()
