"""Phase 1 tests for the SQLite data loader."""

from __future__ import annotations

import sqlite3

import pandas as pd
import pytest

from load_data import (
    CELL_POPULATIONS,
    DataValidationError,
    load_dataframe,
    validate_dataframe,
)
from src.database import connect_database, initialize_schema


def _fixture_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "project": "prj1",
                "subject": "sbj001",
                "condition": "melanoma",
                "age": 57,
                "sex": "M",
                "treatment": "miraclib",
                "response": None,
                "sample": "sample001",
                "sample_type": "PBMC",
                "time_from_treatment_start": 0,
                "b_cell": 100,
                "cd8_t_cell": 200,
                "cd4_t_cell": 300,
                "nk_cell": 150,
                "monocyte": 250,
            },
            {
                "project": "prj1",
                "subject": "sbj001",
                "condition": "melanoma",
                "age": 57,
                "sex": "M",
                "treatment": "miraclib",
                "response": None,
                "sample": "sample002",
                "sample_type": "PBMC",
                "time_from_treatment_start": 7,
                "b_cell": 110,
                "cd8_t_cell": 210,
                "cd4_t_cell": 310,
                "nk_cell": 160,
                "monocyte": 260,
            },
        ]
    )


def test_schema_tables_exist(tmp_path) -> None:
    db_path = tmp_path / "schema.db"
    connection = connect_database(db_path)
    try:
        initialize_schema(connection)
        tables = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            ).fetchall()
        }
    finally:
        connection.close()

    assert {
        "projects",
        "subjects",
        "samples",
        "cell_populations",
        "cell_counts",
    }.issubset(tables)


def test_population_count_and_sample_normalization(tmp_path) -> None:
    db_path = tmp_path / "loaded.db"
    load_dataframe(_fixture_dataframe(), db_path)

    connection = sqlite3.connect(db_path)
    try:
        population_count = connection.execute(
            "SELECT COUNT(*) FROM cell_populations"
        ).fetchone()[0]
        normalized_count = connection.execute(
            "SELECT COUNT(*) FROM cell_counts WHERE sample_id = 'sample001'"
        ).fetchone()[0]
        counts_by_name = dict(
            connection.execute(
                """
                SELECT cp.name, cc.count
                FROM cell_counts cc
                JOIN cell_populations cp ON cp.population_id = cc.population_id
                WHERE cc.sample_id = 'sample001'
                """
            ).fetchall()
        )
    finally:
        connection.close()

    assert population_count == 5
    assert normalized_count == 5
    assert counts_by_name == {
        "b_cell": 100,
        "cd8_t_cell": 200,
        "cd4_t_cell": 300,
        "nk_cell": 150,
        "monocyte": 250,
    }
    assert set(counts_by_name) == set(CELL_POPULATIONS)


def test_missing_response_is_stored_as_sql_null(tmp_path) -> None:
    db_path = tmp_path / "null_response.db"
    load_dataframe(_fixture_dataframe(), db_path)

    connection = sqlite3.connect(db_path)
    try:
        response = connection.execute(
            "SELECT response FROM subjects WHERE subject_id = 'sbj001'"
        ).fetchone()[0]
    finally:
        connection.close()

    assert response is None


def test_foreign_keys_are_enforced_and_check_passes(tmp_path) -> None:
    db_path = tmp_path / "foreign_keys.db"
    load_dataframe(_fixture_dataframe(), db_path)

    connection = connect_database(db_path)
    try:
        assert connection.execute("PRAGMA foreign_keys").fetchone()[0] == 1
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        with pytest.raises(sqlite3.IntegrityError):
            connection.execute(
                """
                INSERT INTO samples (
                    sample_id, subject_id, sample_type, time_from_treatment_start
                )
                VALUES ('orphan_sample', 'missing_subject', 'PBMC', 0)
                """
            )
    finally:
        connection.close()


def test_negative_cell_counts_are_rejected() -> None:
    df = _fixture_dataframe()
    df.loc[0, "b_cell"] = -1

    with pytest.raises(DataValidationError, match="non-negative"):
        validate_dataframe(df)


def test_inconsistent_subject_metadata_is_rejected() -> None:
    df = _fixture_dataframe()
    df.loc[1, "treatment"] = "phauximab"

    with pytest.raises(DataValidationError, match="treatment"):
        validate_dataframe(df)
