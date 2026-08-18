"""SQLite-backed annotation store: per-cluster comments + starred markers/plots.

This is **user data**, so it lives OUTSIDE `.cache/` (which is regenerable) — by default at
`annotations/<name>.sqlite` (configurable via `annotations_dir`). Back it up like any other data.
A fresh connection is opened per call (SQLite is cheap and this keeps it thread-safe under Streamlit).
"""
from __future__ import annotations

import sqlite3
from contextlib import contextmanager
from datetime import datetime
from typing import Optional

import data as D

_SCHEMA = """
CREATE TABLE IF NOT EXISTS comments (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    cluster    INTEGER NOT NULL,
    text       TEXT    NOT NULL,
    created_at TEXT    NOT NULL
);
CREATE TABLE IF NOT EXISTS selected_markers (
    cluster    INTEGER NOT NULL,
    kind       TEXT    NOT NULL,   -- 'gene' | 'motif'
    feature    TEXT    NOT NULL,
    created_at TEXT    NOT NULL,
    PRIMARY KEY (cluster, kind, feature)
);
-- Legacy: the "mark plot as informative" feature was removed (it caused a DB write on every
-- feature-plot navigation). Table kept so existing annotation DBs stay valid and old data is
-- preserved; nothing reads or writes it anymore.
CREATE TABLE IF NOT EXISTS selected_plots (
    cluster    INTEGER NOT NULL,
    kind       TEXT    NOT NULL,   -- 'gene' | 'motif'
    feature    TEXT    NOT NULL,
    page       INTEGER,
    created_at TEXT    NOT NULL,
    PRIMARY KEY (cluster, kind, feature)
);
-- The annotation of record: one row per cluster, holding the decision the user is actually making.
-- Columns mirror the xlsx export template  so export
-- is a projection rather than a translation. Two template columns are deliberately absent:
-- key_marker_genes and select_marker_motifs are derived from `selected_markers` (the star
-- checkboxes), and `comments` is derived from the `comments` table — storing them twice would let
-- them drift.
CREATE TABLE IF NOT EXISTS annotations (
    cluster            INTEGER PRIMARY KEY,
    annot_origin       TEXT,       -- germ layer / lineage, e.g. "Ectoderm"
    annot_trajectory   TEXT,       -- free text
    annot_order        INTEGER,    -- renumbered cluster for figure order
    annot_abbrev       TEXT,       -- short label, e.g. "Fb"
    annot_type         TEXT,       -- full label, e.g. "Forebrain progenitors"
    reference_omg_type TEXT,       -- matching cell type in the reference atlas
    refs               TEXT,       -- e.g. "PMIDs:26371318,25820448"
    ncells             INTEGER,
    pct_cells          REAL,
    reviewed           INTEGER NOT NULL DEFAULT 0,
    updated_at         TEXT    NOT NULL
);
"""

# Editable annotation columns, in export order. Single source of truth for the form, the store and
# the exporter, so adding a column cannot leave one of the three behind.
ANNOTATION_FIELDS = ("annot_origin", "annot_trajectory", "annot_order", "annot_abbrev",
                     "annot_type", "reference_omg_type", "refs", "ncells", "pct_cells")


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


@contextmanager
def _conn(cfg: dict):
    con = sqlite3.connect(D.annotations_db(cfg))
    con.row_factory = sqlite3.Row
    try:
        con.executescript(_SCHEMA)
        yield con
        con.commit()
    finally:
        con.close()


# ---------------------------------------------------------------- comments
def add_comment(cfg: dict, cluster: int, text: str) -> None:
    text = (text or "").strip()
    if not text:
        return
    with _conn(cfg) as c:
        c.execute("INSERT INTO comments(cluster, text, created_at) VALUES (?, ?, ?)",
                  (int(cluster), text, _now()))


def list_comments(cfg: dict, cluster: int | None = None) -> list[dict]:
    """Comments in chronological order (oldest first). cluster=None -> all clusters."""
    with _conn(cfg) as c:
        if cluster is None:
            rows = c.execute("SELECT * FROM comments ORDER BY created_at, id").fetchall()
        else:
            rows = c.execute("SELECT * FROM comments WHERE cluster=? ORDER BY created_at, id",
                             (int(cluster),)).fetchall()
    return [dict(r) for r in rows]


def comment_count(cfg: dict, cluster: int) -> int:
    with _conn(cfg) as c:
        return c.execute("SELECT COUNT(*) FROM comments WHERE cluster=?", (int(cluster),)).fetchone()[0]


def delete_comment(cfg: dict, comment_id: int) -> None:
    with _conn(cfg) as c:
        c.execute("DELETE FROM comments WHERE id=?", (int(comment_id),))


# ---------------------------------------------------------------- starred markers
def selected_markers(cfg: dict, cluster: int, kind: str) -> set[str]:
    with _conn(cfg) as c:
        rows = c.execute("SELECT feature FROM selected_markers WHERE cluster=? AND kind=?",
                         (int(cluster), kind)).fetchall()
    return {r["feature"] for r in rows}


def set_markers(cfg: dict, cluster: int, kind: str, displayed: set[str], selected: set[str]) -> None:
    """Reconcile selection for the *displayed* features only: star `selected`, unstar the rest.

    Features not in `displayed` are left untouched, so filtering/Top-N never drops prior picks.
    Returns early when nothing changed, so merely re-rendering the table (a slider move, a tab
    switch) doesn't write to SQLite on every rerun.
    """
    current = selected_markers(cfg, cluster, kind)
    add = selected - current
    remove = (displayed - selected) & current
    if not add and not remove:
        return
    with _conn(cfg) as c:
        for f in add:
            c.execute("INSERT OR IGNORE INTO selected_markers(cluster, kind, feature, created_at) "
                      "VALUES (?, ?, ?, ?)", (int(cluster), kind, f, _now()))
        for f in remove:
            c.execute("DELETE FROM selected_markers WHERE cluster=? AND kind=? AND feature=?",
                      (int(cluster), kind, f))


# ---------------------------------------------------------------- annotation of record
def get_annotation(cfg: dict, cluster: int) -> dict:
    """The annotation row for a cluster, or {} if it has never been saved."""
    with _conn(cfg) as c:
        row = c.execute("SELECT * FROM annotations WHERE cluster=?", (int(cluster),)).fetchone()
    return dict(row) if row else {}


def all_annotations(cfg: dict) -> dict[int, dict]:
    with _conn(cfg) as c:
        rows = c.execute("SELECT * FROM annotations").fetchall()
    return {int(r["cluster"]): dict(r) for r in rows}


def set_annotation(cfg: dict, cluster: int, values: dict, reviewed: Optional[bool] = None) -> None:
    """Upsert the annotation for one cluster.

    Only keys in ANNOTATION_FIELDS are accepted — an unexpected key is a programming error, not
    something to silently write into a column that may not exist. `reviewed` is left unchanged
    when None so saving a field doesn't quietly un-review a cluster.
    """
    bad = set(values) - set(ANNOTATION_FIELDS)
    if bad:
        raise ValueError(f"unknown annotation field(s): {sorted(bad)}")
    cols = list(values)
    payload = [values[k] for k in cols]
    if reviewed is not None:
        cols.append("reviewed")
        payload.append(1 if reviewed else 0)
    assignments = ", ".join(f"{k}=excluded.{k}" for k in cols)
    with _conn(cfg) as c:
        c.execute(
            f"INSERT INTO annotations(cluster, {', '.join(cols)}, updated_at) "
            f"VALUES ({', '.join('?' * (len(cols) + 2))}) "
            f"ON CONFLICT(cluster) DO UPDATE SET {assignments}, updated_at=excluded.updated_at",
            [int(cluster), *payload, _now()])


def annotated_clusters(cfg: dict) -> set[int]:
    """Clusters with a usable annotation — a label is what makes a row worth exporting."""
    with _conn(cfg) as c:
        rows = c.execute("SELECT cluster FROM annotations WHERE "
                         "COALESCE(TRIM(annot_type), '') != '' "
                         "OR COALESCE(TRIM(annot_abbrev), '') != ''").fetchall()
    return {int(r["cluster"]) for r in rows}
