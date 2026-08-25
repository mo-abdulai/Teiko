"""Lightweight tests for dashboard support helpers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from dashboard import (
    count_metric_value,
    metric_value,
    missing_required_files,
    significance_summary,
)
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
