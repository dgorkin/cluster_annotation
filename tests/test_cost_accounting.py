#!/usr/bin/env python
"""Offline tests for usage/cost telemetry and prompt-cache request shape.

No API key and no network needed — the research loop is driven with a fake client, so this
pins the parts that are otherwise only observable by spending money.

    ~/.conda/envs/cluster_annotation/bin/python tests/test_cost_accounting.py
"""
import json
import sys
import tempfile
import types
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import insights as I  # noqa: E402


# ----------------------------------------------------------------- fakes
def _usage(inp=0, cache_read=0, cache_write=0, out=0, searches=0):
    return types.SimpleNamespace(
        input_tokens=inp,
        cache_read_input_tokens=cache_read,
        cache_creation_input_tokens=cache_write,
        output_tokens=out,
        server_tool_use=types.SimpleNamespace(web_search_requests=searches, web_fetch_requests=0),
    )


def _resp(stop_reason, model="claude-opus-4-8", usage=None, text="analysis text"):
    return types.SimpleNamespace(
        stop_reason=stop_reason, model=model, usage=usage or _usage(), stop_details=None,
        content=[types.SimpleNamespace(type="text", text=text)],
    )


class _FakeClient:
    """Returns scripted responses and records the kwargs of every request.

    Mirrors the streaming interface the research loop uses: messages.stream(...) as a context
    manager, with get_final_message() reassembling the Message.
    """

    def __init__(self, script):
        self.calls = []
        outer = self

        class _Stream:
            def __init__(self, resp):
                self._resp = resp

            def __enter__(self):
                return self

            def __exit__(self, *exc):
                return False

            def get_final_message(self):
                return self._resp

        class _Messages:
            def stream(self, **kw):
                outer.calls.append(kw)
                return _Stream(script[len(outer.calls) - 1])

        self.messages = _Messages()


# ----------------------------------------------------------------- tests
def test_call_cost_matches_hand_computation():
    rec = {"model": "claude-opus-4-8", "input_tokens": 10_000,
           "cache_read_input_tokens": 100_000, "cache_creation_input_tokens": 20_000,
           "output_tokens": 4_000, "web_search_requests": 3}
    # 0.05 base in + 0.05 cache read (0.1x) + 0.125 cache write (1.25x) + 0.10 out + 0.03 search
    assert abs(I.call_cost(rec) - 0.355) < 1e-9, I.call_cost(rec)


def test_cache_read_is_ten_times_cheaper_than_fresh_input():
    fresh = {"model": "claude-opus-4-8", "input_tokens": 1_000_000}
    cached = {"model": "claude-opus-4-8", "cache_read_input_tokens": 1_000_000}
    assert abs(I.call_cost(fresh) - 5.0) < 1e-9
    assert abs(I.call_cost(cached) - 0.5) < 1e-9


def test_summarize_usage_totals_and_prices_per_call():
    calls = [
        {"model": "claude-opus-4-8", "input_tokens": 1_000_000, "output_tokens": 0},
        {"model": "claude-haiku-4-5", "input_tokens": 1_000_000, "output_tokens": 0},
        None,  # dropped
    ]
    s = I.summarize_usage(calls)
    assert s["api_calls"] == 2
    assert s["input_tokens"] == 2_000_000
    assert s["models"] == ["claude-haiku-4-5", "claude-opus-4-8"]
    assert s["unpriced_models"] == []
    # Priced per call ($5 + $1), NOT by applying one model's rate to the total.
    assert abs(s["est_cost_usd"] - 6.0) < 1e-6, s["est_cost_usd"]


def test_dated_snapshot_ids_are_priced_like_their_alias():
    """The API resolves the alias we send into a dated snapshot on the response, e.g.
    claude-haiku-4-5 -> claude-haiku-4-5-20251001. Observed live 2026-08-17 pricing at $0."""
    assert I.price_for("claude-haiku-4-5-20251001") == (1.0, 5.0)
    assert I.price_for("claude-opus-4-8-20260115") == (5.0, 25.0)
    assert I.price_for("claude-haiku-4-5") == (1.0, 5.0)
    assert I.price_for("totally-unknown-model") is None
    assert I.price_for(None) is None

    dated = {"model": "claude-haiku-4-5-20251001", "input_tokens": 3_646, "output_tokens": 1_437}
    alias = {**dated, "model": "claude-haiku-4-5"}
    assert abs(I.call_cost(dated) - I.call_cost(alias)) < 1e-12
    assert abs(I.call_cost(dated) - 0.010831) < 1e-6, I.call_cost(dated)
    # ...and it must not be reported as a pricing gap.
    assert I.summarize_usage([dated])["unpriced_models"] == []


def test_unknown_model_is_flagged_rather_than_silently_free():
    s = I.summarize_usage([{"model": "claude-nextgen-9", "input_tokens": 1_000_000,
                            "web_search_requests": 2}])
    assert s["unpriced_models"] == ["claude-nextgen-9"]
    # Token cost unknown (0), but the per-search charge is model-independent and still counted.
    assert abs(s["est_cost_usd"] - 0.02) < 1e-9, s["est_cost_usd"]


def test_usage_of_survives_a_response_with_no_usage():
    rec = I._usage_of(types.SimpleNamespace())
    assert rec["model"] == "?"
    assert all(rec[k] == 0 for k in I.USAGE_TOKEN_KEYS)


def test_pause_turn_loop_records_one_usage_record_per_http_call():
    script = [_resp("pause_turn", usage=_usage(inp=40_000, searches=3)),
              _resp("pause_turn", usage=_usage(inp=5_000, cache_read=40_000, searches=2)),
              _resp("end_turn", usage=_usage(inp=2_000, cache_read=45_000, out=4_000))]
    client = _FakeClient(script)
    resp, text, calls = I._run_research(client, "claude-opus-4-8", "sys", "user",
                                        8000, "high", I.WEB_SEARCH_TOOL, cache=True)
    assert text == "analysis text"
    assert len(calls) == 3, "each continuation is a separate billed call"
    assert len(client.calls) == 3
    s = I.summarize_usage(calls)
    assert s["web_search_requests"] == 5
    assert s["cache_read_input_tokens"] == 85_000


def test_cache_control_is_sent_when_enabled_and_omitted_when_disabled():
    on = _FakeClient([_resp("end_turn")])
    I._run_research(on, "claude-opus-4-8", "sys", "user", 8000, "high",
                    I.WEB_SEARCH_TOOL, cache=True)
    assert on.calls[0]["cache_control"] == {"type": "ephemeral"}

    off = _FakeClient([_resp("end_turn")])
    I._run_research(off, "claude-opus-4-8", "sys", "user", 8000, "high",
                    I.WEB_SEARCH_TOOL, cache=False)
    assert "cache_control" not in off.calls[0], "prompt_caching: false must not send the param"


def test_refusal_returns_no_text_but_still_records_its_usage():
    client = _FakeClient([_resp("refusal", usage=_usage(inp=1_200, out=5))])
    resp, text, calls = I._run_research(client, "claude-opus-4-8", "sys", "user",
                                        8000, "high", I.WEB_SEARCH_TOOL, cache=True)
    assert text is None
    assert len(calls) == 1, "a refusal can still be billed, so it belongs in the record"


def test_research_streams_rather_than_blocking_on_create():
    """Non-streaming + web_search outruns the HTTP timeout (measured 2026-08-17), so the
    research loop must use messages.stream(). A fake without create() proves it never falls back."""
    client = _FakeClient([_resp("end_turn")])
    assert not hasattr(client.messages, "create"), "fake must not offer the non-streaming path"
    _, text, calls = I._run_research(client, "claude-opus-4-8", "sys", "user", 8000, "high",
                                     I.WEB_SEARCH_TOOL, cache=True)
    assert text == "analysis text" and len(calls) == 1


def test_client_bounds_its_retry_budget():
    """max_retries=8 at the 10-minute default meant ~90 min of silent retrying per wedged call."""
    assert I.MAX_RETRIES <= 3, I.MAX_RETRIES
    assert I.REQUEST_TIMEOUT_S <= 600, I.REQUEST_TIMEOUT_S


def test_truncated_research_is_a_failure_not_a_silent_empty_annotation():
    """Regression: a cluster with ~1,900 gene markers was cached as "Unable to determine - no
    marker data or analysis provided" for $0.52. The research call had hit max_tokens mid-turn
    after 16 searches; the truncated text flowed into structuring, which honestly reported no
    data. A truncated attempt must not become an annotation."""
    client = _FakeClient([_resp("max_tokens", usage=_usage(inp=1_102, out=6_616, searches=16))])
    resp, text, calls = I._run_research(client, "claude-opus-4-8", "sys", "user", 8000, "high",
                                        I.WEB_SEARCH_TOOL, cache=True)
    assert text is None, "truncated output must not be passed on as an analysis"
    assert len(calls) == 1, "the truncated call was still billed, so it stays in the record"
    assert calls[0]["stop_reason"] == "max_tokens", "stop_reason must be recorded to diagnose this"


def test_usage_records_stop_reason():
    assert I._usage_of(_resp("end_turn"))["stop_reason"] == "end_turn"
    assert I._usage_of(types.SimpleNamespace())["stop_reason"] == "?"


def test_max_tokens_default_leaves_room_for_search_heavy_clusters():
    import data as D
    assert D.AI_INSIGHTS_DEFAULTS["max_tokens"] >= 16000, D.AI_INSIGHTS_DEFAULTS["max_tokens"]


def test_continuations_are_capped():
    script = [_resp("pause_turn") for _ in range(I.MAX_PAUSE_CONTINUATIONS + 3)]
    client = _FakeClient(script)
    I._run_research(client, "claude-opus-4-8", "sys", "user", 8000, "high",
                    I.WEB_SEARCH_TOOL, cache=True)
    assert len(client.calls) == I.MAX_PAUSE_CONTINUATIONS


def test_dataset_cost_sums_disk_and_counts_pre_telemetry_artifacts():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = {"name": "t", "cache_dir": tmp, "_root": tmp}
        d = I.insight_dir(cfg)

        def write(path, cost):
            meta = {"usage": {"est_cost_usd": cost}} if cost is not None else {}
            path.write_text(json.dumps({"primary_identity": "x", "_meta": meta}))

        write(I.insight_path(cfg, 0), 0.30)
        write(I.insight_path(cfg, 1), 0.20)
        write(I.insight_path(cfg, 2), None)          # predates telemetry
        write(I.reannotation_path(cfg, 1), 0.40)
        write(I.cohort_review_path(cfg), 0.10)

        t = I.dataset_cost(cfg, [0, 1, 2])
        assert abs(t["originals"] - 0.50) < 1e-9, t
        assert abs(t["reannotations"] - 0.40) < 1e-9, t
        assert abs(t["cohort_review"] - 0.10) < 1e-9, t
        assert abs(t["total"] - 1.00) < 1e-9, t
        assert t["n_priced"] == 4 and t["untracked"] == 1, t
        assert "lower bound" in I.format_dataset_cost(cfg, [0, 1, 2])
        assert d.exists()


def _hash_fixture(tmp):
    """Minimal cfg + marker frames for exercising the cache-key functions."""
    import pandas as pd
    cfg = {"name": "t", "cache_dir": tmp, "_root": tmp, "biological_context": "E9.5 mouse",
           "ai_insights": {"primary_model": "claude-opus-4-8",
                           "fallback_model": "claude-opus-4-7"}}
    frame = pd.DataFrame({"feature": ["A", "B", "C"], "cluster": [0, 0, 0],
                          "avg_log2FC": [3.0, 2.0, 1.0], "avg_diff": [3.0, 2.0, 1.0]})
    return cfg, frame, frame


def test_inputs_hash_is_independent_of_model_choice():
    """Switching model used to mark every cached artifact in a dataset stale."""
    import copy
    with tempfile.TemporaryDirectory() as tmp:
        cfg, g, m = _hash_fixture(tmp)
        other = copy.deepcopy(cfg)
        other["ai_insights"]["primary_model"] = "claude-opus-5"
        assert I.inputs_hash(cfg, 0, g, m) == I.inputs_hash(other, 0, g, m)
        # ...whereas the legacy key moved with the model, which is the bug being migrated away.
        assert I.legacy_inputs_hash(cfg, 0, g, m) != I.legacy_inputs_hash(other, 0, g, m)
        # Markers and context must still invalidate.
        ctx = copy.deepcopy(cfg)
        ctx["biological_context"] = "E15.5 mouse forebrain"
        assert I.inputs_hash(ctx, 0, g, m) != I.inputs_hash(cfg, 0, g, m)


def test_migration_converts_legacy_keys_and_leaves_stale_ones_alone():
    with tempfile.TemporaryDirectory() as tmp:
        cfg, g, m = _hash_fixture(tmp)
        I.insight_dir(cfg)

        def write(cluster, stored_hash):
            I.insight_path(cfg, cluster).write_text(json.dumps(
                {"primary_identity": "x", "_meta": {"inputs_hash": stored_hash}}))

        write(0, I.legacy_inputs_hash(cfg, 0, g, m))   # old scheme -> migrate
        write(1, I.inputs_hash(cfg, 1, g, m))          # already new -> leave
        write(2, "0" * 32)                             # genuinely stale -> leave

        dry = I.migrate_inputs_hashes(cfg, [0, 1, 2, 3], g, m, apply=False)
        assert dry == {"migrated": [0], "already_current": [1], "stale": [2], "missing": [3]}, dry
        # A dry run must not touch disk.
        assert json.loads(I.insight_path(cfg, 0).read_text())["_meta"]["inputs_hash"] \
            == I.legacy_inputs_hash(cfg, 0, g, m)

        I.migrate_inputs_hashes(cfg, [0, 1, 2, 3], g, m, apply=True)
        meta0 = json.loads(I.insight_path(cfg, 0).read_text())["_meta"]
        assert meta0["inputs_hash"] == I.inputs_hash(cfg, 0, g, m)
        assert meta0["inputs_hash_migrated_from"]      # provenance of the rewrite is kept
        assert json.loads(I.insight_path(cfg, 2).read_text())["_meta"]["inputs_hash"] == "0" * 32

        # After migration cluster 0 survives a model switch; cluster 2 stays stale.
        import copy
        moved = copy.deepcopy(cfg)
        moved["ai_insights"]["primary_model"] = "claude-opus-5"
        assert I.clusters_needing_insight(moved, [0, 1, 2], g, m) == [2]


def test_record_cost_tolerates_old_and_malformed_records():
    assert I.record_cost(None) == 0.0
    assert I.record_cost({}) == 0.0
    assert I.record_cost({"_meta": {}}) == 0.0
    assert I.record_cost({"_meta": {"usage": {}}}) == 0.0
    assert I.record_cost({"_meta": {"usage": {"est_cost_usd": 0.42}}}) == 0.42


def test_every_selectable_model_has_a_list_price():
    """The sidebar lets the user switch models; an unpriced one silently under-reports cost."""
    import data as D
    missing = [m for m in D.KNOWN_MODELS if m not in I.MODEL_PRICES]
    assert not missing, f"add list prices for {missing}"


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
        except Exception as e:  # noqa: BLE001 - report, don't abort the suite
            failed += 1
            print(f"  ERROR {t.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
