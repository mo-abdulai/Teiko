"""Lightweight tests for dashboard support helpers."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd

from dashboard import (
    count_metric_value,
    database_is_ready,
    metric_value,
    missing_required_files,
    significance_summary,
)
from src.database import connect_database, initialize_schema
from src.visualization import response_label


def test_response_label_mapping() -> None:
    assert response_label("yes") == "Responder"
    assert response_label("no") == "Non-responder"
    assert response_label("unknown") == "Unknown"


def test_significance_summary_reports_no_significant_populations() -> None:
    results = pd.DataFrame(
        {
            "population": ["b_cell", "cd8_t_cell"],
            "significant": ["False", False],
        }
    )

    assert "No immune-cell populations" in significance_summary(results)


def test_significance_summary_lists_significant_populations() -> None:
    results = pd.DataFrame(
        {
            "population": ["b_cell", "cd8_t_cell"],
            "significant": [True, False],
        }
    )

    assert "b_cell" in significance_summary(results)


def test_metric_value_extracts_requested_metric() -> None:
    metrics = pd.DataFrame(
        {
            "metric": [
                "matching_baseline_samples",
                "melanoma_male_responder_baseline_avg_b_cell",
            ],
            "value": ["485", "10206.15"],
        }
    )

    assert (
        metric_value(metrics, "melanoma_male_responder_baseline_avg_b_cell")
        == "10206.15"
    )
    assert metric_value(metrics, "missing") is None


def test_count_metric_value_formats_float_count_as_integer() -> None:
    metrics = pd.DataFrame(
        {
            "metric": ["matching_baseline_samples"],
            "value": [485.0],
        }
    )

    assert count_metric_value(metrics, "matching_baseline_samples") == "485"


def test_missing_required_files_reports_missing_inputs(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr("dashboard.get_database_path", lambda: tmp_path / "missing.db")

    missing = missing_required_files(Path(tmp_path))

    assert missing
    assert any(path.name == "missing.db" for path in missing)


def test_database_is_ready_rejects_empty_sqlite_file(tmp_path) -> None:
    db_path = tmp_path / "empty.db"
    sqlite3.connect(db_path).close()

    assert database_is_ready(db_path) is False


def test_database_is_ready_accepts_queryable_database(tmp_path) -> None:
    db_path = tmp_path / "ready.db"
    connection = connect_database(db_path)
    try:
        initialize_schema(connection)
        connection.execute(
            "INSERT INTO projects (project_id, project_name) VALUES (1, 'prj1')"
        )
        connection.execute(
            """
            INSERT INTO subjects (
                subject_id, project_id, condition, age, sex, treatment, response
            )
            VALUES ('subj1', 1, 'melanoma', 57, 'M', 'miraclib', 'yes')
            """
        )
        connection.execute(
            """
            INSERT INTO samples (
                sample_id, subject_id, sample_type, time_from_treatment_start
            )
            VALUES ('sample1', 'subj1', 'PBMC', 0)
            """
        )
        connection.executemany(
            "INSERT INTO cell_populations (population_id, name) VALUES (?, ?)",
            [
                (1, "b_cell"),
                (2, "cd8_t_cell"),
                (3, "cd4_t_cell"),
                (4, "nk_cell"),
                (5, "monocyte"),
            ],
        )
        connection.executemany(
            "INSERT INTO cell_counts (sample_id, population_id, count) VALUES (?, ?, ?)",
            [("sample1", population_id, 10) for population_id in range(1, 6)],
        )
        connection.commit()
    finally:
        connection.close()

    assert database_is_ready(db_path) is True
