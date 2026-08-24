"""Reusable analysis functions for immune-cell population data."""

from __future__ import annotations

import os

import numpy as np


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


EXPECTED_POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
FREQUENCY_COLUMNS = ["sample", "total_count", "population", "count", "percentage"]


class FrequencyValidationError(ValueError):
    """Raised when cell-frequency inputs or outputs fail validation."""


def calculate_relative_frequencies(cell_counts_df: pd.DataFrame) -> pd.DataFrame:
    """Calculate relative frequencies for each population in each sample.

    The input dataframe must contain one row per sample and population with
    columns: sample, population, and count.
    """
    _validate_cell_count_input(cell_counts_df)

    frequencies = cell_counts_df[["sample", "population", "count"]].copy()
    frequencies["count"] = pd.to_numeric(frequencies["count"], errors="raise")
    totals = frequencies.groupby("sample", sort=False)["count"].transform("sum")
    frequencies["total_count"] = totals

    zero_total_samples = frequencies.loc[
        frequencies["total_count"] <= 0, "sample"
    ].drop_duplicates()
    if not zero_total_samples.empty:
        raise FrequencyValidationError(
            "Samples must have total_count > 0. Affected samples: "
            + ", ".join(zero_total_samples.astype(str).head(10))
        )

    frequencies["percentage"] = (
        frequencies["count"] / frequencies["total_count"] * 100
    )
    frequencies["population"] = pd.Categorical(
        frequencies["population"],
        categories=EXPECTED_POPULATIONS,
        ordered=True,
    )
    frequencies = frequencies.sort_values(["sample", "population"]).reset_index(
        drop=True
    )
    frequencies["population"] = frequencies["population"].astype(str)
    frequencies = frequencies[FREQUENCY_COLUMNS]

    validate_relative_frequencies(frequencies)
    return frequencies


def validate_relative_frequencies(
    frequencies: pd.DataFrame,
    *,
    expected_samples: list[str] | None = None,
    expected_populations: list[str] | None = None,
    rtol: float = 1e-9,
    atol: float = 1e-9,
) -> None:
    """Validate the Part 2 frequency table before it is exported."""
    populations = expected_populations or EXPECTED_POPULATIONS

    if list(frequencies.columns) != FREQUENCY_COLUMNS:
        raise FrequencyValidationError(
            "Frequency output columns must be exactly: "
            + ", ".join(FREQUENCY_COLUMNS)
        )

    duplicate_rows = frequencies.duplicated(["sample", "population"])
    if duplicate_rows.any():
        examples = frequencies.loc[duplicate_rows, ["sample", "population"]].head(10)
        raise FrequencyValidationError(
            "Frequency output has duplicate sample/population rows: "
            + str(examples.to_dict(orient="records"))
        )

    if expected_samples is not None:
        actual_samples = set(frequencies["sample"])
        expected_sample_set = set(expected_samples)
        if actual_samples != expected_sample_set:
            missing = sorted(expected_sample_set - actual_samples)[:10]
            unexpected = sorted(actual_samples - expected_sample_set)[:10]
            raise FrequencyValidationError(
                "Frequency output samples do not match database samples. "
                f"Missing: {missing}; unexpected: {unexpected}"
            )

    expected_row_count = frequencies["sample"].nunique() * len(populations)
    if len(frequencies) != expected_row_count:
        raise FrequencyValidationError(
            f"Expected {expected_row_count:,} frequency rows, found {len(frequencies):,}."
        )

    _validate_population_structure(frequencies, populations)

    total_counts_per_sample = frequencies.groupby("sample")["total_count"].nunique()
    inconsistent_totals = total_counts_per_sample[total_counts_per_sample != 1]
    if not inconsistent_totals.empty:
        raise FrequencyValidationError(
            "All rows for each sample must share the same total_count. "
            "Affected samples: " + ", ".join(inconsistent_totals.index.astype(str)[:10])
        )

    count_sums = frequencies.groupby("sample")["count"].sum()
    sample_totals = frequencies.groupby("sample")["total_count"].first()
    if not count_sums.equals(sample_totals):
        mismatched = count_sums[count_sums != sample_totals].index.astype(str)[:10]
        raise FrequencyValidationError(
            "Sum of population counts must equal total_count. Affected samples: "
            + ", ".join(mismatched)
        )

    percentage_sums = frequencies.groupby("sample")["percentage"].sum()
    close_to_100 = np.isclose(percentage_sums.to_numpy(), 100.0, rtol=rtol, atol=atol)
    if not np.all(close_to_100):
        bad_samples = percentage_sums[~close_to_100].index.astype(str)[:10]
        raise FrequencyValidationError(
            "Percentages must sum to approximately 100 for every sample. "
            "Affected samples: " + ", ".join(bad_samples)
        )

    if (frequencies["count"] < 0).any():
        raise FrequencyValidationError("Cell counts must be non-negative.")
    if (frequencies["percentage"] < -atol).any() or (
        frequencies["percentage"] > 100 + atol
    ).any():
        raise FrequencyValidationError("Percentages must be between 0 and 100.")


def _validate_cell_count_input(cell_counts_df: pd.DataFrame) -> None:
    required_columns = {"sample", "population", "count"}
    missing_columns = required_columns - set(cell_counts_df.columns)
    if missing_columns:
        raise FrequencyValidationError(
            "Cell-count input is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    if cell_counts_df.empty:
        raise FrequencyValidationError("Cell-count input is empty.")

    duplicate_rows = cell_counts_df.duplicated(["sample", "population"])
    if duplicate_rows.any():
        examples = cell_counts_df.loc[duplicate_rows, ["sample", "population"]].head(10)
        raise FrequencyValidationError(
            "Input has duplicate sample/population rows: "
            + str(examples.to_dict(orient="records"))
        )

    counts = pd.to_numeric(cell_counts_df["count"], errors="raise")
    if (counts < 0).any():
        raise FrequencyValidationError("Cell counts must be non-negative.")

    _validate_population_structure(cell_counts_df, EXPECTED_POPULATIONS)


def _validate_population_structure(
    cell_counts_df: pd.DataFrame, expected_populations: list[str]
) -> None:
    expected_set = set(expected_populations)
    observed_populations = set(cell_counts_df["population"])
    unexpected = sorted(observed_populations - expected_set)
    if unexpected:
        raise FrequencyValidationError(
            "Unexpected immune-cell populations found: " + ", ".join(unexpected)
        )

    populations_by_sample = cell_counts_df.groupby("sample")["population"].agg(
        lambda values: set(values)
    )
    invalid_samples = populations_by_sample[populations_by_sample != expected_set]
    if not invalid_samples.empty:
        examples = []
        for sample, observed in invalid_samples.head(10).items():
            missing = sorted(expected_set - observed)
            extra = sorted(observed - expected_set)
            examples.append({"sample": sample, "missing": missing, "unexpected": extra})
        raise FrequencyValidationError(
            "Every sample must contain exactly the expected populations. "
            + str(examples)
        )
