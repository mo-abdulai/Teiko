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
RESPONSE_FREQUENCY_COLUMNS = [
    "subject",
    "sample",
    "condition",
    "treatment",
    "sample_type",
    "time_from_treatment_start",
    "response",
    "population",
    "count",
    "percentage",
]
SUBJECT_SUMMARY_COLUMNS = [
    "subject",
    "response",
    "population",
    "mean_percentage",
    "n_timepoints",
]


class FrequencyValidationError(ValueError):
    """Raised when cell-frequency inputs or outputs fail validation."""


class ResponseAnalysisValidationError(ValueError):
    """Raised when Part 3 response-analysis data fail validation."""


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


def build_response_frequency_data(cell_counts_df: pd.DataFrame) -> pd.DataFrame:
    """Attach reusable Part 2 percentages to the Part 3 response cohort."""
    required_columns = {
        "subject",
        "sample",
        "condition",
        "treatment",
        "sample_type",
        "time_from_treatment_start",
        "response",
        "population",
        "count",
    }
    missing_columns = required_columns - set(cell_counts_df.columns)
    if missing_columns:
        raise ResponseAnalysisValidationError(
            "Response cohort input is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )
    if cell_counts_df.empty:
        raise ResponseAnalysisValidationError("Response cohort input is empty.")

    frequencies = calculate_relative_frequencies(
        cell_counts_df[["sample", "population", "count"]]
    )
    metadata_columns = [
        "subject",
        "sample",
        "condition",
        "treatment",
        "sample_type",
        "time_from_treatment_start",
        "response",
    ]
    metadata = cell_counts_df[metadata_columns].drop_duplicates()
    if metadata["sample"].duplicated().any():
        examples = metadata.loc[metadata["sample"].duplicated(), "sample"].head(10)
        raise ResponseAnalysisValidationError(
            "Response cohort has conflicting sample metadata. Examples: "
            + ", ".join(examples.astype(str))
        )

    response_frequencies = metadata.merge(frequencies, on="sample", how="inner")
    response_frequencies = response_frequencies[
        RESPONSE_FREQUENCY_COLUMNS
    ].sort_values(
        ["subject", "time_from_treatment_start", "sample", "population"]
    )
    response_frequencies = response_frequencies.reset_index(drop=True)
    validate_response_frequency_data(response_frequencies)
    return response_frequencies


def validate_response_frequency_data(frequencies: pd.DataFrame) -> None:
    """Validate the Part 3 response cohort before statistical analysis."""
    required_columns = set(RESPONSE_FREQUENCY_COLUMNS)
    missing_columns = required_columns - set(frequencies.columns)
    if missing_columns:
        raise ResponseAnalysisValidationError(
            "Response frequency data are missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    expected_filters = {
        "condition": "melanoma",
        "treatment": "miraclib",
        "sample_type": "PBMC",
    }
    for column, expected in expected_filters.items():
        invalid = frequencies[column] != expected
        if invalid.any():
            bad_values = frequencies.loc[invalid, column].astype(str).unique()[:10]
            raise ResponseAnalysisValidationError(
                f"Response cohort has invalid {column} values: "
                + ", ".join(bad_values)
            )

    invalid_response = ~frequencies["response"].isin(["yes", "no"])
    if invalid_response.any():
        bad_values = frequencies.loc[invalid_response, "response"].astype(str).unique()
        raise ResponseAnalysisValidationError(
            "Response cohort must contain only yes/no responses. Found: "
            + ", ".join(bad_values[:10])
        )

    invalid_populations = set(frequencies["population"]) - set(EXPECTED_POPULATIONS)
    if invalid_populations:
        raise ResponseAnalysisValidationError(
            "Response cohort has unexpected populations: "
            + ", ".join(sorted(invalid_populations))
        )

    if (frequencies["percentage"] < -1e-9).any() or (
        frequencies["percentage"] > 100 + 1e-9
    ).any():
        raise ResponseAnalysisValidationError("Percentages must be between 0 and 100.")

    _validate_population_structure(frequencies, EXPECTED_POPULATIONS)


def build_subject_level_response_summary(
    frequencies: pd.DataFrame,
    *,
    baseline_only: bool = False,
) -> pd.DataFrame:
    """Aggregate response-cohort frequencies to one row per subject/population."""
    validate_response_frequency_data(frequencies)
    analysis_df = frequencies.copy()
    if baseline_only:
        analysis_df = analysis_df[analysis_df["time_from_treatment_start"] == 0]
        if analysis_df.empty:
            raise ResponseAnalysisValidationError(
                "Baseline response cohort contains no timepoint 0 rows."
            )

    response_counts = analysis_df.groupby("subject")["response"].nunique()
    conflicting_subjects = response_counts[response_counts > 1]
    if not conflicting_subjects.empty:
        raise ResponseAnalysisValidationError(
            "Subjects must not have conflicting response labels. Affected subjects: "
            + ", ".join(conflicting_subjects.index.astype(str)[:10])
        )

    summary = (
        analysis_df.groupby(["subject", "response", "population"], as_index=False)
        .agg(
            mean_percentage=("percentage", "mean"),
            n_timepoints=("time_from_treatment_start", "nunique"),
        )
        .sort_values(["subject", "population"])
        .reset_index(drop=True)
    )
    summary["population"] = pd.Categorical(
        summary["population"],
        categories=EXPECTED_POPULATIONS,
        ordered=True,
    )
    summary = summary.sort_values(["subject", "population"]).reset_index(drop=True)
    summary["population"] = summary["population"].astype(str)
    summary = summary[SUBJECT_SUMMARY_COLUMNS]
    validate_subject_level_response_summary(summary)
    return summary


def validate_subject_level_response_summary(summary: pd.DataFrame) -> None:
    """Validate one primary statistical observation per subject/population."""
    if list(summary.columns) != SUBJECT_SUMMARY_COLUMNS:
        raise ResponseAnalysisValidationError(
            "Subject summary columns must be exactly: "
            + ", ".join(SUBJECT_SUMMARY_COLUMNS)
        )
    duplicate_rows = summary.duplicated(["subject", "population"])
    if duplicate_rows.any():
        examples = summary.loc[duplicate_rows, ["subject", "population"]].head(10)
        raise ResponseAnalysisValidationError(
            "Subject summary must contain one row per subject/population. "
            + str(examples.to_dict(orient="records"))
        )
    if set(summary["population"]) != set(EXPECTED_POPULATIONS):
        raise ResponseAnalysisValidationError(
            "Subject summary must analyze exactly the expected populations."
        )
    if not {"yes", "no"}.issubset(set(summary["response"])):
        raise ResponseAnalysisValidationError(
            "Subject summary must contain both responder groups."
        )


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
