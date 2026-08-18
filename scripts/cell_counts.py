#!/usr/bin/env python
"""Count cells per cluster from the dataset's Seurat object, ahead of opening the app.

    ./run_app.sh counts config/mydataset.yaml [--force]

The app does this itself the first time a dataset loads, but reading a whole-embryo multiome object
takes tens of seconds and several GB — so it is worth doing from the command line if you would
rather not wait on a page load, or if you want to see which metadata column was matched.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import data as D  # noqa: E402


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    force = "--force" in sys.argv[1:]
    configs = [Path(a) for a in args] or sorted(D.available_configs())
    if not configs:
        print("no dataset config found — pass one: ./run_app.sh counts config/mydataset.yaml")
        return 1
    failed = 0
    for path in configs:
        cfg = D.load_config(str(path))
        cc = D.cell_counts_cfg(cfg)
        print(f"\n-- {path.name} ({cfg.get('name')})")
        if cc["tsv"]:
            print(f"   counts come from {cc['tsv']} — nothing to compute")
            continue
        if not cc["seurat_object"]:
            print("   no Seurat object configured (cell_counts.seurat_object) — skipping")
            continue
        try:
            # The marker tables define the cluster ids, and the column in the object is chosen by
            # matching them, so preprocessing has to have run first.
            D.ensure_preprocessed(cfg)
            out = D.ensure_cell_counts(cfg, force=force, expected=D.clusters(cfg))
        except Exception as exc:  # noqa: BLE001 - report and carry on to the next config
            failed += 1
            print(f"   FAILED: {exc}")
            continue
        counts = D.load_cell_counts(cfg)
        meta = D.cell_counts_meta(cfg)
        print(f"   {sum(v['ncells'] for v in counts.values()):,} cells over {len(counts)} clusters"
              f" -> {out}")
        if meta:
            print(f"   object var '{meta.get('object_var')}', cluster column "
                  f"'{meta.get('cluster_column')}'")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
