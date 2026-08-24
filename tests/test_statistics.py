"""Phase 4 tests for response aggregation and statistical comparisons."""

from __future__ import annotations

import pandas as pd
import pytest

from src.analysis import (
    EXPECTED_POPULATIONS,
    SUBJECT_SUMMARY_COLUMNS,
    build_subject_level_response_summary,
)
from src.statistics import (
    STATISTICAL_RESULT_COLUMNS,
    StatisticalAnalysisError,
    apply_multiple_testing_correction,
    compare_response_groups,
)


def _frequency_rows(
    subject: str,
    response: str,
    *,
    percentages_by_timepoint: dict[int, dict[str, float]],
) -> list[dict[str, object]]:
    rows = []
    for timepoint, percentages in percentages_by_timepoint.items():
        for population in EXPECTED_POPULATIONS:
            rows.append(
                {
                    "subject": subject,
                    "sample": f"{subject}_{timepoint}",
                    "condition": "melanoma",
                    "treatment": "miraclib",
                    "sample_type": "PBMC",
                    "time_from_treatment_start": timepoint,
                    "response": response,
                    "population": population,
                    "count": 1,
                    "percentage": percentages[population],
                }
            )
    return rows


def _balanced_percentages(value: float) -> dict[str, float]:
    return {
        "b_cell": value,
        "cd8_t_cell": 5.0,
        "cd4_t_cell": 3.0,
        "nk_cell": 2.0,
        "monocyte": 100.0 - value - 10.0,
    }


def test_subject_level_aggregation_averages_longitudinal_timepoints() -> None:
    df = pd.DataFrame(
        _frequency_rows(
            "responder1",
            "yes",
            percentages_by_timepoint={
                0: _balanced_percentages(10.0),
                7: _balanced_percentages(20.0),
                14: _balanced_percentages(30.0),
            },
        )
        + _frequency_rows(
            "nonresponder1",
            "no",
            percentages_by_timepoint={
                0: _balanced_percentages(5.0),
                7: _balanced_percentages(5.0),
                14: _balanced_percentages(5.0),
            },
        )
    )

    summary = build_subject_level_response_summary(df)
    b_cell = summary[
        (summary["subject"] == "responder1") & (summary["population"] == "b_cell")
    ].iloc[0]

    assert list(summary.columns) == SUBJECT_SUMMARY_COLUMNS
    assert b_cell["mean_percentage"] == pytest.approx(20.0)
    assert b_cell["n_timepoints"] == 3


def test_primary_summary_has_one_row_per_subject_population() -> None:
    df = pd.DataFrame(
        _frequency_rows(
            "responder1",
            "yes",
            percentages_by_timepoint={
                0: _balanced_percentages(10.0),
                7: _balanced_percentages(20.0),
                14: _balanced_percentages(30.0),
            },
        )
        + _frequency_rows(
            "nonresponder1",
            "no",
            percentages_by_timepoint={
                0: _balanced_percentages(5.0),
                7: _balanced_percentages(10.0),
                14: _balanced_percentages(15.0),
            },
        )
    )

    summary = build_subject_level_response_summary(df)

    assert len(summary) == 2 * len(EXPECTED_POPULATIONS)
    assert not summary.duplicated(["subject", "population"]).any()


def test_baseline_summary_uses_only_time_zero() -> None:
    df = pd.DataFrame(
        _frequency_rows(
            "responder1",
            "yes",
            percentages_by_timepoint={
                0: _balanced_percentages(10.0),
                7: _balanced_percentages(90.0),
                14: _balanced_percentages(90.0),
            },
        )
        + _frequency_rows(
            "nonresponder1",
            "no",
            percentages_by_timepoint={
                0: _balanced_percentages(5.0),
                7: _balanced_percentages(90.0),
                14: _balanced_percentages(90.0),
            },
        )
    )

    summary = build_subject_level_response_summary(df, baseline_only=True)
    b_cell = summary[
        (summary["subject"] == "responder1") & (summary["population"] == "b_cell")
    ].iloc[0]

    assert b_cell["mean_percentage"] == pytest.approx(10.0)
    assert b_cell["n_timepoints"] == 1


def _statistical_summary() -> pd.DataFrame:
    rows = []
    for index in range(5):
        for population in EXPECTED_POPULATIONS:
            rows.append(
                {
                    "subject": f"yes_{index}_{population}",
                    "response": "yes",
                    "population": population,
                    "mean_percentage": 80.0 + index if population == "b_cell" else 50.0,
                    "n_timepoints": 3,
                }
            )
            rows.append(
                {
                    "subject": f"no_{index}_{population}",
                    "response": "no",
                    "population": population,
                    "mean_percentage": 10.0 + index if population == "b_cell" else 50.0,
                    "n_timepoints": 3,
                }
            )
    return pd.DataFrame(rows)


def test_mann_whitney_output_columns_counts_and_p_values() -> None:
    results = compare_response_groups(_statistical_summary())
    b_cell = results[results["population"] == "b_cell"].iloc[0]

    assert list(results.columns) == STATISTICAL_RESULT_COLUMNS
    assert b_cell["responder_n"] == 5
    assert b_cell["non_responder_n"] == 5
    assert b_cell["mann_whitney_u"] >= 0
    assert 0 <= b_cell["p_value"] <= 1


def test_multiple_testing_correction_adds_adjusted_p_values_and_significance() -> None:
    results = pd.DataFrame(
        {
            "population": EXPECTED_POPULATIONS,
            "p_value": [0.001, 0.02, 0.2, 0.8, 0.04],
        }
    )

    corrected = apply_multiple_testing_correction(results)

    assert "adjusted_p_value" in corrected.columns
    assert "significant" in corrected.columns
    assert corrected.loc[corrected["population"] == "b_cell", "significant"].iloc[0]


def test_non_significant_case_remains_non_significant() -> None:
    results = compare_response_groups(_statistical_summary())
    cd8 = results[results["population"] == "cd8_t_cell"].iloc[0]

    assert not bool(cd8["significant"])


def test_strong_difference_produces_small_p_value() -> None:
    results = compare_response_groups(_statistical_summary())
    b_cell = results[results["population"] == "b_cell"].iloc[0]

    assert b_cell["p_value"] < 0.05
    assert b_cell["median_difference"] > 0


def test_missing_group_raises_clear_error() -> None:
    summary = _statistical_summary()
    summary = summary[summary["response"] == "yes"]

    with pytest.raises(StatisticalAnalysisError, match="Both response groups"):
        compare_response_groups(summary)
