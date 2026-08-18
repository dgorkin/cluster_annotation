#!/usr/bin/env python
"""Side-by-side comparison of two models on the same clusters, same pipeline.

Both arms run with identical settings apart from primary_model, so differences are the model's
rather than the harness's. Results go to scratch cache dirs — the dataset's real artifacts under
.cache/<name>/ai_insights/ are never touched, so this is safe to run against a finished dataset.

    python scripts/compare_models.py config/mydataset.yaml 28,2,10 claude-opus-4-8 claude-opus-5

Costs one full research pass per cluster per arm (~$0.45 each at Opus rates), so the default of
3 clusters x 2 arms is roughly $2.70. Wall time ~8 min per cluster, run `workers` at a time.
"""
import json
import sys
import time
import copy
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import data as D      # noqa: E402
import insights as I  # noqa: E402

WORKERS = 3


def run_arm(cfg_base: dict, clusters, genes, motifs, primary: str, fallback: str,
            scratch: str) -> dict:
    cfg = copy.deepcopy(cfg_base)
    cfg["cache_dir"] = scratch
    cfg["ai_insights"]["primary_model"] = primary
    cfg["ai_insights"]["fallback_model"] = fallback
    cfg["ai_insights"]["max_workers"] = WORKERS

    print(f"\n=== arm: {primary} (fallback {fallback}) ===", flush=True)
    t0 = time.time()
    done = [0]

    def cb(c, err):
        done[0] += 1
        rec = I.load_insight(cfg, c)
        cost = I.record_cost(rec) if not err else 0.0
        tag = "ok" if not err else f"FAIL: {err[:100]}"
        print(f"  [{done[0]}/{len(clusters)}] cluster {c}: {tag} "
              f"(~${cost:.3f}, +{time.time()-t0:.0f}s)", flush=True)

    errors = I.generate_all(cfg, clusters, genes, motifs, progress_cb=cb, force=True)
    out = {"primary": primary, "fallback": fallback, "wall_s": round(time.time() - t0),
           "errors": dict(errors), "clusters": {}}
    for c in clusters:
        rec = I.load_insight(cfg, c)
        if rec is None:
            continue
        u = rec["_meta"].get("usage") or {}
        out["clusters"][c] = {
            "identity": rec.get("primary_identity"),
            "confidence": rec.get("confidence"),
            "alternatives": rec.get("alternative_identities") or [],
            "key_genes": rec.get("key_genes") or [],
            "caveats": rec.get("caveats", ""),
            "reasoning": rec.get("reasoning", ""),
            "n_refs": len(rec.get("references") or []),
            "pmids": [r.get("pmid") for r in (rec.get("references") or []) if r.get("pmid")],
            # Keep the reference objects in full — citation quality is the thing a reviewer has
            # to judge by eye, and a pmid count alone hides "5 plausible-looking citations with
            # no resolvable identifier", which is exactly what turned up on 2026-08-17.
            "references": rec.get("references") or [],
            "served": rec["_meta"].get("model_used"),
            "cost": u.get("est_cost_usd"),
            "in_tok": u.get("input_tokens"),
            "cache_read": u.get("cache_read_input_tokens"),
            "out_tok": u.get("output_tokens"),
            "searches": u.get("web_search_requests"),
        }
    return out


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
    config = args[0] if args else str(ROOT / "config" / "mydataset.yaml")
    clusters = [int(x) for x in args[1].split(",")] if len(args) > 1 else [28, 2, 10]
    model_a = args[2] if len(args) > 2 else "claude-opus-4-8"
    model_b = args[3] if len(args) > 3 else "claude-opus-5"

    cfg = D.load_config(_resolve_config(config))
    if not I.api_key_available(cfg):
        print("ERROR: no Anthropic API key found.")
        sys.exit(1)
    D.ensure_preprocessed(cfg)
    genes, motifs = D.load_markers(cfg, "gene"), D.load_markers(cfg, "motif")
    aic = D.ai_insights_cfg(cfg)

    print(f"dataset={cfg['name']} clusters={clusters} arms=[{model_a}, {model_b}] "
          f"effort={aic['effort']} max_uses={aic['web_search_max_uses']} "
          f"caching={aic['prompt_caching']}")
    print(f"estimated spend ~${0.45 * len(clusters) * 2:.2f}", flush=True)

    # The dataset's existing annotations, for reference. Generated under the OLD pipeline
    # (uncapped searches, no caching) so treat content differences as indicative, not causal.
    existing = {c: I.load_insight(cfg, c) for c in clusters}

    arms = []
    for primary, fallback in ((model_a, model_b), (model_b, model_a)):
        with tempfile.TemporaryDirectory() as scratch:
            arms.append(run_arm(cfg, clusters, genes, motifs, primary, fallback, scratch))

    print(f"\n{'='*100}\nCOMPARISON — {cfg['name']}\n{'='*100}")
    for c in clusters:
        print(f"\n--- cluster {c} " + "-" * 84)
        old = existing.get(c) or {}
        if old:
            print(f"  on disk (old pipeline, {(old.get('_meta') or {}).get('model_used')}):")
            print(f"      {old.get('primary_identity')}  [{old.get('confidence')}]")
        for arm in arms:
            d = arm["clusters"].get(c)
            if not d:
                print(f"  {arm['primary']:<18} FAILED: {arm['errors'].get(c, 'no record')[:80]}")
                continue
            print(f"  {arm['primary']:<18} {d['identity']}")
            print(f"  {'':<18} conf={d['confidence']}  refs={d['n_refs']}  "
                  f"pmids={len(d['pmids'])}  served={d['served']}")
            print(f"  {'':<18} ~${d['cost']:.3f}  in={d['in_tok']:,} cache_rd={d['cache_read']:,} "
                  f"out={d['out_tok']:,} srch={d['searches']}")
            print(f"  {'':<18} alts: {', '.join(d['alternatives'][:3]) or '—'}")

    print(f"\n{'='*100}\nTOTALS\n{'='*100}")
    for arm in arms:
        costs = [d['cost'] for d in arm['clusters'].values() if d.get('cost')]
        tot = sum(costs)
        print(f"  {arm['primary']:<18} ~${tot:.2f} over {len(costs)} cluster(s) "
              f"(~${tot/len(costs):.3f}/cluster)  wall {arm['wall_s']}s  "
              f"failures={len(arm['errors'])}")

    out = ROOT / "logs" / f"model_comparison_{cfg['name']}.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps({"dataset": cfg["name"], "clusters": clusters,
                               "existing": {c: (existing[c] or {}) for c in clusters},
                               "arms": arms}, indent=2, default=str))
    print(f"\nfull records (reasoning, caveats, references) written to {out}")
    print("DONE")


if __name__ == "__main__":
    main()
