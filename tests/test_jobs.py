#!/usr/bin/env python
"""Offline tests for the durable background-job layer. No API key, no network, no Streamlit.

    ~/.conda/envs/cluster_annotation/bin/python tests/test_jobs.py

The point of jobs.py is that a paid, minutes-long call stays visible across script re-runs, page
reloads and other tabs — so what matters here is that the state file always tells the truth,
including when the job crashes or the process that owned it is gone.
"""
import sys
import tempfile
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "app"))

import jobs as J  # noqa: E402


def _cfg(tmp: str, name: str = "t") -> dict:
    return {"name": name, "cache_dir": tmp, "_root": tmp}


def _wait(pred, timeout: float = 5.0) -> bool:
    """Poll until pred() or timeout — the work runs on its own thread."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if pred():
            return True
        time.sleep(0.02)
    return False


def test_no_record_before_anything_runs():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _cfg(tmp)
        assert J.read(cfg, "primer") is None
        assert not J.is_running(cfg, "primer")
        assert J.any_running(cfg, ("primer", "generate")) is None


def test_a_job_reports_running_then_done_with_its_summary():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _cfg(tmp)
        release = threading.Event()

        def work(prog, _stop):
            prog.update(message="searching")
            release.wait(5)
            return "Primer built: 12 cell types."

        assert J.start(cfg, "primer", "Build reference primer", work)
        assert _wait(lambda: J.is_running(cfg, "primer")), "should report running while in flight"
        # Wait for the job's own first progress write, rather than racing it.
        assert _wait(lambda: (J.read(cfg, "primer") or {}).get("message") == "searching")
        rec = J.read(cfg, "primer")
        # The record has to stand on its own: another session reads only this file.
        assert rec["label"] == "Build reference primer" and rec["message"] == "searching"
        assert rec["primary_model"] and rec["started_at"]

        release.set()
        assert _wait(lambda: J.read(cfg, "primer")["state"] == J.DONE)
        rec = J.read(cfg, "primer")
        assert rec["message"] == "Primer built: 12 cell types."
        assert rec["finished_at"] and not J.is_running(cfg, "primer")


def test_a_crashing_job_records_failed_rather_than_staying_running():
    """A job that raises must not leave 'running' on disk — that would be a lie forever."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _cfg(tmp)

        def work(_prog, _stop):
            raise RuntimeError("truncated at max_tokens=8000")

        assert J.start(cfg, "primer", "Build reference primer", work)
        assert _wait(lambda: J.read(cfg, "primer")["state"] == J.FAILED)
        assert "truncated at max_tokens" in J.read(cfg, "primer")["message"]


def test_starting_twice_is_refused_so_a_paid_call_cannot_double_run():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _cfg(tmp)
        release, ran = threading.Event(), []

        def work(_prog, _stop):
            ran.append(1)
            release.wait(5)
            return "ok"

        assert J.start(cfg, "primer", "first", work)
        assert _wait(lambda: J.is_running(cfg, "primer"))
        assert J.start(cfg, "primer", "second", work) is False, "second launch must be refused"
        release.set()
        assert _wait(lambda: J.read(cfg, "primer")["state"] == J.DONE)
        assert len(ran) == 1, ran


def test_progress_is_readable_from_outside_the_thread():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _cfg(tmp)
        step, release = threading.Event(), threading.Event()

        def work(prog, _stop):
            prog.update(done=7, message="cluster 12 done")
            step.set()
            release.wait(5)
            return "done"

        J.start(cfg, "generate", "Annotate 45 clusters", work, total=45)
        assert step.wait(5)
        assert _wait(lambda: (J.read(cfg, "generate") or {}).get("done") == 7)
        rec = J.read(cfg, "generate")
        assert rec["total"] == 45 and rec["message"] == "cluster 12 done"
        release.set()
        assert _wait(lambda: J.read(cfg, "generate")["state"] == J.DONE)


def test_stop_is_visible_to_the_job_and_noted_in_the_record():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _cfg(tmp)
        started, saw_stop = threading.Event(), threading.Event()

        def work(_prog, should_stop):
            started.set()
            for _ in range(250):
                if should_stop():
                    saw_stop.set()
                    return "3 of 45 clusters annotated."
                time.sleep(0.02)
            return "ran to completion"

        J.start(cfg, "generate", "Annotate 45 clusters", work, total=45, stoppable=True)
        assert started.wait(5)
        J.request_stop(cfg, "generate")
        assert saw_stop.wait(5), "the job never observed the stop request"
        assert _wait(lambda: J.read(cfg, "generate")["state"] == J.DONE)
        msg = J.read(cfg, "generate")["message"]
        assert msg.startswith("Stopped early."), msg
        # The flag must not leak into the next run.
        assert not J.stop_requested(cfg, "generate")


def test_a_record_owned_by_a_dead_process_reads_as_interrupted():
    """The app restarting kills the thread; 'running' on disk must not be believed."""
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _cfg(tmp)
        J._write(cfg, "primer", {"kind": "primer", "label": "Build reference primer",
                                 "state": J.RUNNING, "pid": 2 ** 22, "started_at": "x",
                                 "done": 0, "total": None, "message": "…"})
        rec = J.read(cfg, "primer")
        assert rec["state"] == J.INTERRUPTED, rec
        assert "restarted" in rec["message"]
        assert not J.is_running(cfg, "primer")
        # And a fresh run is allowed, since nothing is actually holding it.
        assert J.start(cfg, "primer", "retry", lambda _p, _s: "ok")
        assert _wait(lambda: J.read(cfg, "primer")["state"] == J.DONE)


def test_clear_dismisses_a_finished_job_but_never_a_running_one():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _cfg(tmp)
        release = threading.Event()
        J.start(cfg, "primer", "Build reference primer", lambda _p, _s: release.wait(5) or "ok")
        assert _wait(lambda: J.is_running(cfg, "primer"))
        J.clear(cfg, "primer")
        assert J.is_running(cfg, "primer"), "a running job must not be dismissable"
        release.set()
        assert _wait(lambda: J.read(cfg, "primer")["state"] == J.DONE)
        J.clear(cfg, "primer")
        assert J.read(cfg, "primer") is None


def test_jobs_are_scoped_per_dataset():
    with tempfile.TemporaryDirectory() as tmp:
        one, two = _cfg(tmp, "ds_one"), _cfg(tmp, "ds_two")
        release = threading.Event()
        J.start(one, "primer", "one", lambda _p, _s: release.wait(5) or "ok")
        assert _wait(lambda: J.is_running(one, "primer"))
        assert not J.is_running(two, "primer"), "another dataset must be unaffected"
        assert J.start(two, "primer", "two", lambda _p, _s: "ok")
        release.set()
        assert _wait(lambda: J.read(one, "primer")["state"] == J.DONE)


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
