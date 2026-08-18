"""Durable background jobs for the paid AI calls.

Streamlit re-runs the whole script on every interaction, so work started *inside* a run only has
status for as long as that run lasts. A primer build takes minutes: switch section or reload the
page and the spinner disappears, the confirm button resets, and there is no way to tell whether a
paid call is still in flight — which is how the same literature search gets paid for twice.

So a job runs on its own thread and records its state in a file under the dataset's cache. The
thread outlives the script run that started it (the app is one long-lived server process), and the
state file is readable from any session — a different tab, or the same tab after a reload.

Two rules for anything running here:
  * the job function must not touch `st.*` — it has no script context. It returns a summary string;
    the UI reads state from the file.
  * every exit path writes a terminal state, so "running" can never be a lie told by a crash.

A server restart does kill the thread. That is detected rather than believed: the record carries
the pid that owns it, and a record claiming to run under a pid that is gone reads as `interrupted`.
"""
from __future__ import annotations

import json
import os
import threading
from datetime import datetime
from pathlib import Path
from typing import Callable, Optional

import data as D

RUNNING, DONE, FAILED, INTERRUPTED = "running", "done", "failed", "interrupted"
_LOCK = threading.Lock()          # serializes read-modify-write of the state files
_STOP: set[tuple[str, str]] = set()   # (dataset, kind) asked to stop; in-process, like the threads


def _dir(cfg: dict) -> Path:
    d = D.cache_dir(cfg) / "jobs"
    d.mkdir(parents=True, exist_ok=True)
    return d


def path(cfg: dict, kind: str) -> Path:
    return _dir(cfg) / f"{kind}.json"


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _write(cfg: dict, kind: str, rec: dict) -> None:
    """Replace the record. Never raises: losing a status line must not kill a paid job."""
    try:
        with _LOCK:
            tmp = path(cfg, kind).with_suffix(".tmp")
            tmp.write_text(json.dumps(rec, indent=2))
            tmp.replace(path(cfg, kind))       # atomic: a reader never sees half a record
    except OSError:
        pass


def _pid_alive(pid: Optional[int]) -> bool:
    if not pid:
        return False
    try:
        os.kill(int(pid), 0)
    except (OSError, ValueError):
        return False
    return True


def read(cfg: dict, kind: str) -> Optional[dict]:
    """The job record, or None if this dataset has never run this job.

    A `running` record whose owning process is gone is reported as `interrupted` — the state on
    disk is what the last live writer managed to say, not necessarily what is true now.
    """
    try:
        rec = json.loads(path(cfg, kind).read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if rec.get("state") == RUNNING and not _pid_alive(rec.get("pid")):
        rec["state"] = INTERRUPTED
        rec["message"] = ("The app was restarted while this was running, so its outcome is "
                          "unknown. Anything it finished writing was kept; re-run to complete "
                          "the rest.")
    return rec


def is_running(cfg: dict, kind: str) -> bool:
    rec = read(cfg, kind)
    return bool(rec and rec["state"] == RUNNING)


def any_running(cfg: dict, kinds) -> Optional[str]:
    """The first of `kinds` currently running — paid AI work is serialized in the UI, because two
    overlapping runs on one dataset would race on the same cache files."""
    return next((k for k in kinds if is_running(cfg, k)), None)


def clear(cfg: dict, kind: str) -> None:
    """Dismiss a finished record. A running one is left alone."""
    if is_running(cfg, kind):
        return
    try:
        path(cfg, kind).unlink()
    except OSError:
        pass


# ---------------------------------------------------------------- stop requests
def request_stop(cfg: dict, kind: str) -> None:
    """Ask a running job to stop. Only meaningful for jobs that do per-cluster work: the call
    already in flight is paid for and cannot be pulled back, but the ones not yet started can."""
    _STOP.add((cfg["name"], kind))
    rec = read(cfg, kind)
    if rec and rec["state"] == RUNNING:
        rec["stopping"] = True
        _write(cfg, kind, rec)


def stop_requested(cfg: dict, kind: str) -> bool:
    return (cfg["name"], kind) in _STOP


# ---------------------------------------------------------------- running a job
class Progress:
    """Handed to the job function so it can report where it is. Writes through to the state file."""

    def __init__(self, cfg: dict, kind: str):
        self._cfg, self._kind = cfg, kind

    def update(self, *, done: Optional[int] = None, total: Optional[int] = None,
               message: Optional[str] = None) -> None:
        rec = read(self._cfg, self._kind)
        if not rec or rec["state"] != RUNNING:
            return
        if done is not None:
            rec["done"] = int(done)
        if total is not None:
            rec["total"] = int(total)
        if message is not None:
            rec["message"] = message
        _write(self._cfg, self._kind, rec)


def start(cfg: dict, kind: str, label: str,
          fn: Callable[[Progress, Callable[[], bool]], Optional[str]],
          *, total: Optional[int] = None, stoppable: bool = False) -> bool:
    """Run `fn` on a background thread, tracked in `<cache>/<name>/jobs/<kind>.json`.

    Returns False if this job is already running (the caller should not double-spend). `fn` is
    called as fn(progress, should_stop) and may return a summary line for the UI.

    The thread is a daemon: `./run_app.sh stop` should stop the app, not hang on a 20-minute
    generation. The cost is that a stopped server leaves an `interrupted` record, which `read()`
    reports honestly.
    """
    if is_running(cfg, kind):
        return False
    _STOP.discard((cfg["name"], kind))
    aic = D.ai_insights_cfg(cfg)
    _write(cfg, kind, {
        "kind": kind, "label": label, "state": RUNNING, "pid": os.getpid(),
        "started_at": _now(), "finished_at": None, "done": 0, "total": total,
        "message": "Starting…", "stopping": False,
        # What it is actually running with, recorded at launch: the sidebar can be changed while a
        # job is in flight, so the config at read time is not necessarily the one being billed.
        "primary_model": aic["primary_model"], "fallback_model": aic["fallback_model"],
        "effort": aic["effort"], "max_tokens": aic["max_tokens"], "stoppable": bool(stoppable),
    })

    def _run() -> None:
        prog = Progress(cfg, kind)
        try:
            summary = fn(prog, lambda: stop_requested(cfg, kind))
            state, message = DONE, (summary or "Finished.")
            if stop_requested(cfg, kind):
                message = f"Stopped early. {message}"
        except Exception as exc:  # noqa: BLE001 - a failed job must still report, not vanish
            state, message = FAILED, f"{type(exc).__name__}: {exc}"
        rec = read(cfg, kind) or {}
        rec.update({"state": state, "message": message, "finished_at": _now(),
                    "stopping": False})
        _write(cfg, kind, rec)
        _STOP.discard((cfg["name"], kind))

    threading.Thread(target=_run, name=f"job:{cfg['name']}:{kind}", daemon=True).start()
    return True
