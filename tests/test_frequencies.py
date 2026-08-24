"""Phase 2 tests for cell-population relative-frequency calculations."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.analysis import (
    EXPECTED_POPULATIONS,
    FREQUENCY_COLUMNS,
    FrequencyValidationError,
    calculate_relative_frequencies,
)
from src.database import connect_database, initialize_schema
from src.queries import get_cell_counts_by_sample


def _sample_rows(sample: str, counts: dict[str, int]) -> list[dict[str, object]]:
    return [
        {"sample": sample, "population": population, "count": counts[population]}
        for population in EXPECTED_POPULATIONS
    ]


def _single_sample_df() -> pd.DataFrame:
    return pd.DataFrame(
        _sample_rows(
            "s1",
            {
                "b_cell": 100,
                "cd8_t_cell": 200,
                "cd4_t_cell": 300,
                "nk_cell": 150,
                "monocyte": 250,
            },
        )
    )


def test_calculates_correct_total_for_all_population_rows() -> None:
    result = calculate_relative_frequencies(_single_sample_df())

    assert set(result["total_count"]) == {1000}


def test_calculates_correct_percentages() -> None:
    result = calculate_relative_frequencies(_single_sample_df())
    percentages = dict(zip(result["population"], result["percentage"]))

    assert percentages["b_cell"] == pytest.approx(10.0)
    assert percentages["cd8_t_cell"] == pytest.approx(20.0)
    assert percentages["cd4_t_cell"] == pytest.approx(30.0)
    assert percentages["nk_cell"] == pytest.approx(15.0)
    assert percentages["monocyte"] == pytest.approx(25.0)


def test_percentages_sum_to_approximately_100_for_every_sample() -> None:
    result = calculate_relative_frequencies(
        pd.DataFrame(
            _sample_rows(
                "s1",
                {
                    "b_cell": 100,
                    "cd8_t_cell": 200,
                    "cd4_t_cell": 300,
                    "nk_cell": 150,
                    "monocyte": 250,
                },
            )
            + _sample_rows(
                "s2",
                {
                    "b_cell": 50,
                    "cd8_t_cell": 50,
                    "cd4_t_cell": 100,
                    "nk_cell": 100,
                    "monocyte": 200,
                },
            )
        )
    )

    sums = result.groupby("sample")["percentage"].sum()
    assert np.allclose(sums, 100.0)


def test_multiple_samples_are_calculated_independently() -> None:
    result = calculate_relative_frequencies(
        pd.DataFrame(
            _sample_rows(
                "s1",
                {
                    "b_cell": 100,
                    "cd8_t_cell": 200,
                    "cd4_t_cell": 300,
                    "nk_cell": 150,
                    "monocyte": 250,
                },
            )
            + _sample_rows(
                "s2",
                {
                    "b_cell": 1,
                    "cd8_t_cell": 1,
                    "cd4_t_cell": 1,
                    "nk_cell": 1,
                    "monocyte": 1,
                },
            )
        )
    )

    totals = dict(result.groupby("sample")["total_count"].first())
    assert totals == {"s1": 1000, "s2": 5}


def test_missing_population_fails_clearly() -> None:
    df = _single_sample_df()
    df = df[df["population"] != "monocyte"]

    with pytest.raises(FrequencyValidationError, match="expected populations"):
        calculate_relative_frequencies(df)


def test_duplicate_population_fails_clearly() -> None:
    df = pd.concat([_single_sample_df(), _single_sample_df().iloc[[0]]])

    with pytest.raises(FrequencyValidationError, match="duplicate"):
        calculate_relative_frequencies(df)


def test_zero_total_sample_fails_clearly() -> None:
    df = pd.DataFrame(
        _sample_rows(
            "s1",
            {
                "b_cell": 0,
                "cd8_t_cell": 0,
                "cd4_t_cell": 0,
                "nk_cell": 0,
                "monocyte": 0,
            },
        )
    )

    with pytest.raises(FrequencyValidationError, match="total_count > 0"):
        calculate_relative_frequencies(df)


def test_output_shape_and_column_order() -> None:
    result = calculate_relative_frequencies(_single_sample_df())

    assert list(result.columns) == FREQUENCY_COLUMNS
    assert len(result) == 5


def test_sqlite_query_and_frequency_calculation_integration(tmp_path) -> None:
    db_path = tmp_path / "phase2.db"
    connection = connect_database(db_path)
    try:
        initialize_schema(connection)
        connection.executescript(
            """
            INSERT INTO projects (project_id, project_name)
            VALUES (1, 'prj1');

            INSERT INTO subjects (
                subject_id, project_id, condition, age, sex, treatment, response
            )
            VALUES ('subject1', 1, 'melanoma', 50, 'female', 'miraclib', 'yes');

            INSERT INTO samples (
                sample_id, subject_id, sample_type, time_from_treatment_start
            )
            VALUES ('s1', 'subject1', 'PBMC', 0);

            INSERT INTO cell_populations (population_id, name)
            VALUES
                (1, 'b_cell'),
                (2, 'cd8_t_cell'),
                (3, 'cd4_t_cell'),
                (4, 'nk_cell'),
                (5, 'monocyte');

            INSERT INTO cell_counts (sample_id, population_id, count)
            VALUES
                ('s1', 1, 100),
                ('s1', 2, 200),
                ('s1', 3, 300),
                ('s1', 4, 150),
                ('s1', 5, 250);
            """
        )

        counts = get_cell_counts_by_sample(connection)
        result = calculate_relative_frequencies(counts)
    finally:
        connection.close()

    assert list(result.columns) == FREQUENCY_COLUMNS
    assert len(result) == 5
    assert set(result["total_count"]) == {1000}
    b_cell_percentage = result.loc[
        result["population"] == "b_cell", "percentage"
    ].iloc[0]
    assert b_cell_percentage == pytest.approx(10.0)
