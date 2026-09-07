"""Persistencia ligera del historial de análisis (SQLite).

Sin servidor adicional: guarda un registro de cada análisis ejecutado
(sin datos sensibles) para que el usuario vea la actividad de la API.
"""

import os
import sqlite3
from datetime import datetime, timezone

DB_PATH = os.environ.get("BIOSEQ_DB_PATH", "history.db")


def _connect() -> sqlite3.Connection:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS analyses (
            id         INTEGER PRIMARY KEY AUTOINCREMENT,
            tool       TEXT NOT NULL,
            input_len  INTEGER NOT NULL,
            summary    TEXT,
            created_at TEXT NOT NULL
        )
        """
    )
    return conn


def record(tool: str, input_len: int, summary: str | None = None) -> None:
    conn = _connect()
    try:
        conn.execute(
            "INSERT INTO analyses (tool, input_len, summary, created_at) "
            "VALUES (?, ?, ?, ?)",
            (tool, input_len, summary, datetime.now(timezone.utc).isoformat()),
        )
        conn.commit()
    finally:
        conn.close()


def list_all(limit: int = 50) -> list[dict]:
    conn = _connect()
    try:
        rows = conn.execute(
            "SELECT * FROM analyses ORDER BY id DESC LIMIT ?", (limit,)
        ).fetchall()
        return [dict(row) for row in rows]
    finally:
        conn.close()


def count() -> int:
    conn = _connect()
    try:
        row = conn.execute("SELECT COUNT(*) AS c FROM analyses").fetchone()
        return int(row["c"])
    finally:
        conn.close()


def clear() -> int:
    conn = _connect()
    try:
        cur = conn.execute("DELETE FROM analyses")
        conn.commit()
        return cur.rowcount
    finally:
        conn.close()