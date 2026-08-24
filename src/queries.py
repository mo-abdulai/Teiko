"""Reusable SQL query definitions for the Teiko technical assessment."""

from __future__ import annotations

import os
import sqlite3


def _import_pandas_quietly():
    """Import pandas without leaking optional Arrow CPU probe diagnostics."""
    stderr_fd = os.dup(2)
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull_fd, 2)
        import pandas as pandas_module
    finally:
        os.dup2(stderr_fd, 2)
        os.close(stderr_fd)
        os.close(devnull_fd)
    return pandas_module


pd = _import_pandas_quietly()


def get_cell_counts_by_sample(connection: sqlite3.Connection) -> pd.DataFrame:
    """Return normalized cell-count records needed for Part 2 analysis."""
    return pd.read_sql_query(
        """
        SELECT
            s.sample_id AS sample,
            cp.name AS population,
            cc.count AS count
        FROM samples AS s
        INNER JOIN cell_counts AS cc
            ON cc.sample_id = s.sample_id
        INNER JOIN cell_populations AS cp
            ON cp.population_id = cc.population_id
        ORDER BY
            s.sample_id,
            cp.population_id
        """,
        connection,
    )


def get_sample_ids(connection: sqlite3.Connection) -> list[str]:
    """Return all sample IDs from the normalized database."""
    rows = connection.execute(
        """
        SELECT sample_id
        FROM samples
        ORDER BY sample_id
        """
    ).fetchall()
    return [row[0] for row in rows]


def get_population_names(connection: sqlite3.Connection) -> list[str]:
    """Return all immune-cell population names in database order."""
    rows = connection.execute(
        """
        SELECT name
        FROM cell_populations
        ORDER BY population_id
        """
    ).fetchall()
    return [row[0] for row in rows]
