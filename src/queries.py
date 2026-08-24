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


def get_baseline_melanoma_miraclib_pbmc_samples(
    connection: sqlite3.Connection,
) -> pd.DataFrame:
    """Return the Part 4 baseline melanoma/miraclib/PBMC cohort."""
    return pd.read_sql_query(
        """
        SELECT
            s.sample_id AS sample,
            subj.subject_id AS subject,
            p.project_name AS project,
            subj.condition AS condition,
            subj.sex AS sex,
            subj.treatment AS treatment,
            subj.response AS response,
            s.sample_type AS sample_type,
            s.time_from_treatment_start AS time_from_treatment_start
        FROM samples AS s
        INNER JOIN subjects AS subj
            ON subj.subject_id = s.subject_id
        INNER JOIN projects AS p
            ON p.project_id = subj.project_id
        WHERE subj.condition = ?
            AND subj.treatment = ?
            AND s.sample_type = ?
            AND s.time_from_treatment_start = ?
        ORDER BY
            p.project_name,
            subj.subject_id,
            s.sample_id
        """,
        connection,
        params=("melanoma", "miraclib", "PBMC", 0),
    )


def get_baseline_sample_counts_by_project(
    connection: sqlite3.Connection,
) -> pd.DataFrame:
    """Return sample counts by project for the Part 4 baseline cohort."""
    return pd.read_sql_query(
        """
        SELECT
            p.project_name AS project,
            COUNT(s.sample_id) AS sample_count
        FROM samples AS s
        INNER JOIN subjects AS subj
            ON subj.subject_id = s.subject_id
        INNER JOIN projects AS p
            ON p.project_id = subj.project_id
        WHERE subj.condition = ?
            AND subj.treatment = ?
            AND s.sample_type = ?
            AND s.time_from_treatment_start = ?
        GROUP BY p.project_name
        ORDER BY p.project_name
        """,
        connection,
        params=("melanoma", "miraclib", "PBMC", 0),
    )


def get_baseline_subject_counts_by_response(
    connection: sqlite3.Connection,
) -> pd.DataFrame:
    """Return distinct-subject response counts for the Part 4 baseline cohort."""
    return pd.read_sql_query(
        """
        SELECT
            COALESCE(subj.response, 'unknown') AS response,
            COUNT(DISTINCT subj.subject_id) AS subject_count
        FROM samples AS s
        INNER JOIN subjects AS subj
            ON subj.subject_id = s.subject_id
        WHERE subj.condition = ?
            AND subj.treatment = ?
            AND s.sample_type = ?
            AND s.time_from_treatment_start = ?
        GROUP BY subj.response
        ORDER BY
            CASE
                WHEN subj.response = 'yes' THEN 1
                WHEN subj.response = 'no' THEN 2
                WHEN subj.response IS NULL THEN 3
                ELSE 4
            END,
            subj.response
        """,
        connection,
        params=("melanoma", "miraclib", "PBMC", 0),
    )


def get_baseline_subject_counts_by_sex(connection: sqlite3.Connection) -> pd.DataFrame:
    """Return distinct-subject sex counts for the Part 4 baseline cohort."""
    return pd.read_sql_query(
        """
        SELECT
            subj.sex AS sex,
            COUNT(DISTINCT subj.subject_id) AS subject_count
        FROM samples AS s
        INNER JOIN subjects AS subj
            ON subj.subject_id = s.subject_id
        WHERE subj.condition = ?
            AND subj.treatment = ?
            AND s.sample_type = ?
            AND s.time_from_treatment_start = ?
        GROUP BY subj.sex
        ORDER BY subj.sex
        """,
        connection,
        params=("melanoma", "miraclib", "PBMC", 0),
    )


def get_baseline_melanoma_male_responder_b_cells(
    connection: sqlite3.Connection,
    *,
    male_value: str = "M",
) -> pd.DataFrame:
    """Return B-cell rows for melanoma male responders at baseline.

    This query intentionally includes all treatments and all sample types.
    """
    return pd.read_sql_query(
        """
        SELECT
            s.sample_id AS sample,
            subj.subject_id AS subject,
            p.project_name AS project,
            subj.condition AS condition,
            subj.sex AS sex,
            subj.treatment AS treatment,
            subj.response AS response,
            s.sample_type AS sample_type,
            s.time_from_treatment_start AS time_from_treatment_start,
            cp.name AS population,
            cc.count AS b_cell_count
        FROM samples AS s
        INNER JOIN subjects AS subj
            ON subj.subject_id = s.subject_id
        INNER JOIN projects AS p
            ON p.project_id = subj.project_id
        INNER JOIN cell_counts AS cc
            ON cc.sample_id = s.sample_id
        INNER JOIN cell_populations AS cp
            ON cp.population_id = cc.population_id
        WHERE subj.condition = ?
            AND subj.sex = ?
            AND subj.response = ?
            AND s.time_from_treatment_start = ?
            AND cp.name = ?
        ORDER BY
            p.project_name,
            subj.subject_id,
            s.sample_id
        """,
        connection,
        params=("melanoma", male_value, "yes", 0, "b_cell"),
    )


def get_baseline_melanoma_male_responder_b_cell_average(
    connection: sqlite3.Connection,
    *,
    male_value: str = "M",
) -> float | None:
    """Return the average B-cell count for melanoma male responders at baseline."""
    row = connection.execute(
        """
        SELECT AVG(cc.count) AS average_b_cell_count
        FROM samples AS s
        INNER JOIN subjects AS subj
            ON subj.subject_id = s.subject_id
        INNER JOIN cell_counts AS cc
            ON cc.sample_id = s.sample_id
        INNER JOIN cell_populations AS cp
            ON cp.population_id = cc.population_id
        WHERE subj.condition = ?
            AND subj.sex = ?
            AND subj.response = ?
            AND s.time_from_treatment_start = ?
            AND cp.name = ?
        """,
        ("melanoma", male_value, "yes", 0, "b_cell"),
    ).fetchone()
    return None if row is None or row[0] is None else float(row[0])
