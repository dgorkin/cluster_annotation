#!/usr/bin/env python
"""Preflight checks: catch a broken setup as one clear report instead of a mid-run traceback.

    ./run_app.sh doctor              # all configs
    ./run_app.sh doctor config/mydataset.yaml

Checks the interpreter and imports, Rscript, the secrets file's existence and permissions (never
its contents), and every input path each dataset config points at. Exits non-zero if anything is
FAIL, so it is usable as a gate in a script.
"""
import importlib
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

OK, WARN, FAIL = "ok  ", "WARN", "FAIL"
counts = {OK: 0, WARN: 0, FAIL: 0}


def report(level: str, what: str, detail: str = "") -> None:
    counts[level] += 1
    mark = {OK: "  ✓", WARN: "  !", FAIL: "  ✗"}[level]
    print(f"{mark} {what}" + (f"\n      {detail}" if detail else ""))


def check_python() -> None:
    print("\n-- interpreter and packages")
    report(OK, f"python {sys.version.split()[0]} at {sys.executable}")
    for mod, why in [("streamlit", "the UI"), ("fitz", "PDF rendering (pymupdf)"),
                     ("pandas", "marker tables"), ("yaml", "config parsing"),
                     ("openpyxl", "xlsx export"), ("anthropic", "AI insights")]:
        try:
            m = importlib.import_module(mod)
            report(OK, f"import {mod} {getattr(m, '__version__', '')}".rstrip())
        except Exception as e:
            report(FAIL, f"import {mod} — needed for {why}", f"{type(e).__name__}: {e}")


def check_rscript() -> None:
    print("\n-- Rscript (marker export from the source RDS)")
    path = shutil.which("Rscript")
    if not path:
        report(FAIL, "Rscript not on PATH",
               "preprocess/export_markers.R cannot run, so a dataset with no cached marker "
               "tables cannot be loaded at all.")
        return
    try:
        v = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=30)
        report(OK, f"Rscript at {path}", (v.stdout or v.stderr).strip().splitlines()[0])
    except Exception as e:
        report(WARN, f"Rscript at {path} but --version failed", str(e))


def check_secrets(cfg=None) -> None:
    print("\n-- Anthropic API key")
    import data as D
    import insights as I
    aic = D.ai_insights_cfg(cfg or {})
    path = Path(os.path.expanduser(aic["api_key_file"]))
    if os.environ.get("ANTHROPIC_API_KEY"):
        report(WARN, "ANTHROPIC_API_KEY is set in the environment",
               "It takes precedence over the secrets file, and an exported key is visible to "
               "other processes of this user. Prefer the chmod-600 file.")
    if not path.exists():
        report(WARN, f"no secrets file at {path}",
               "AI insights will be disabled; the rest of the app works. Create it with:\n"
               f"      printf 'ANTHROPIC_API_KEY=sk-ant-...\\n' > {path} && chmod 600 {path}")
        return
    mode = path.stat().st_mode
    if mode & (stat.S_IRGRP | stat.S_IROTH | stat.S_IWGRP | stat.S_IWOTH):
        report(FAIL, f"{path} is readable or writable beyond its owner "
                     f"({stat.filemode(mode)})",
               f"On a shared server that exposes the key to other users. Fix: chmod 600 {path}")
    else:
        report(OK, f"secrets file {path} ({stat.filemode(mode)})")
    # Confirm a key can be parsed out, without printing any part of it.
    report(OK if I.load_api_key(cfg or {}) else FAIL,
           "key parses from the secrets file" if I.load_api_key(cfg or {})
           else f"no ANTHROPIC_API_KEY= line found in {path}")


def check_config(path: Path) -> None:
    import data as D
    print(f"\n-- config {path.name}")
    try:
        cfg = D.load_config(str(path))
    except Exception as e:
        report(FAIL, f"{path.name} does not parse", f"{type(e).__name__}: {e}")
        return
    name = cfg.get("name", "?")
    report(OK, f"parses; dataset name '{name}', cache -> {D.cache_dir(cfg)}")

    for key in ("markers_rds", "motif_markers_rds", "tangram_pdf"):
        if key not in cfg:
            report(WARN, f"{key} not set")
            continue
        p = D.resolve(cfg, cfg[key])
        report(OK if p.exists() else FAIL, f"{key}", str(p) if not p.exists() else "")

    # Per-cluster feature plots: check the first cluster only, since a missing pattern shows up
    # there and checking 34x2 files on every doctor run is noise.
    cached = (D.cache_dir(cfg) / "feature_plot_index.tsv").exists()
    if cached:
        clusters = D.clusters(cfg)
        for kind in ("gene", "motif"):
            p = D.featureplot_pdf(cfg, kind, clusters[0])
            report(OK if p.exists() else FAIL, f"{kind} feature plots (cluster {clusters[0]})",
                   str(p) if not p.exists() else "")
        report(OK, f"marker cache present ({len(clusters)} clusters)")
    else:
        report(WARN, "no marker cache yet",
               "First load will run preprocess/export_markers.R via Rscript.")

    for entry in D.other_annotations(cfg):
        p = D.other_pdf(cfg, entry, 0)
        report(OK if p.exists() else WARN,
               f"other_annotations '{entry.get('name', '?')}'", str(p) if not p.exists() else "")

    aic = D.ai_insights_cfg(cfg)
    report(OK, f"models: primary {aic['primary_model']}, fallback {aic['fallback_model']}, "
               f"structuring {aic['structuring_model']}")
    report(OK, f"effort {aic['effort']}, max_uses {aic['web_search_max_uses']}, "
               f"caching {aic['prompt_caching']}, workers {aic['max_workers']}")

    import insights as I
    unpriced = [m for m in (aic["primary_model"], aic["fallback_model"], aic["structuring_model"])
                if I.price_for(m) is None]
    if unpriced:
        report(WARN, f"no list price on file for {', '.join(unpriced)}",
               "Cost telemetry will understate spend for these. Add them to "
               "insights.MODEL_PRICES.")


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    if args:
        configs = [Path(a) for a in args]
    else:
        # Skip *.template.yaml: it ships with placeholder paths on purpose, and reporting those as
        # failures is the first thing a new user would see after cloning.
        configs = sorted(p for p in (ROOT / "config").glob("*.yaml")
                         if not p.name.endswith(".template.yaml"))
        if not configs:
            print("\n-- dataset configs")
            report(WARN, "no dataset config yet",
                   "Create one from the template, then re-run:\n"
                   "      cp config/dataset.template.yaml config/mydataset.yaml\n"
                   "      ./run_app.sh doctor config/mydataset.yaml")

    print(f"cluster-annotation doctor — {ROOT}")
    check_python()
    check_rscript()
    first = None
    for c in configs:
        check_config(c)
        if first is None:
            try:
                import data as D
                first = D.load_config(str(c))
            except Exception:
                pass
    check_secrets(first)

    print(f"\n{counts[OK]} ok, {counts[WARN]} warning(s), {counts[FAIL]} failure(s)")
    if counts[FAIL]:
        print("Fix the failures above before running the app or a generation pass.")
        return 1
    print("Ready." if not counts[WARN] else "Usable, with warnings noted above.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
