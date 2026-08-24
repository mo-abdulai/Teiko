"""Statistical analysis helpers for treatment response comparisons."""

from __future__ import annotations

import os

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

from src.analysis import EXPECTED_POPULATIONS


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


STATISTICAL_RESULT_COLUMNS = [
    "population",
    "responder_n",
    "non_responder_n",
    "responder_median_percentage",
    "non_responder_median_percentage",
    "median_difference",
    "mann_whitney_u",
    "p_value",
    "adjusted_p_value",
    "significant",
    "effect_size",
]


class StatisticalAnalysisError(ValueError):
    """Raised when statistical comparisons cannot be performed."""


def compare_response_groups(
    subject_summary: pd.DataFrame,
    *,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Compare responder and non-responder percentages by population."""
    _validate_statistical_input(subject_summary)

    rows = []
    for population in EXPECTED_POPULATIONS:
        population_df = subject_summary[subject_summary["population"] == population]
        responders = population_df.loc[
            population_df["response"] == "yes", "mean_percentage"
        ]
        non_responders = population_df.loc[
            population_df["response"] == "no", "mean_percentage"
        ]
        if responders.empty or non_responders.empty:
            raise StatisticalAnalysisError(
                f"Both response groups must contain data for {population}."
            )

        test_result = mannwhitneyu(
            responders,
            non_responders,
            alternative="two-sided",
        )
        responder_median = float(responders.median())
        non_responder_median = float(non_responders.median())
        u_statistic = float(test_result.statistic)
        rows.append(
            {
                "population": population,
                "responder_n": int(len(responders)),
                "non_responder_n": int(len(non_responders)),
                "responder_median_percentage": responder_median,
                "non_responder_median_percentage": non_responder_median,
                "median_difference": responder_median - non_responder_median,
                "mann_whitney_u": u_statistic,
                "p_value": float(test_result.pvalue),
                "effect_size": _rank_biserial_correlation(
                    u_statistic, len(responders), len(non_responders)
                ),
            }
        )

    results = pd.DataFrame(rows)
    results = apply_multiple_testing_correction(results, alpha=alpha)
    return results[STATISTICAL_RESULT_COLUMNS]


def apply_multiple_testing_correction(
    results: pd.DataFrame,
    *,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Apply Benjamini-Hochberg FDR correction to p-values."""
    if "p_value" not in results.columns:
        raise StatisticalAnalysisError("Results must contain a p_value column.")
    corrected = results.copy()
    reject, adjusted_p_values, _, _ = multipletests(
        corrected["p_value"].to_numpy(),
        alpha=alpha,
        method="fdr_bh",
    )
    corrected["adjusted_p_value"] = adjusted_p_values
    corrected["significant"] = reject
    corrected["population"] = pd.Categorical(
        corrected["population"],
        categories=EXPECTED_POPULATIONS,
        ordered=True,
    )
    corrected = corrected.sort_values("population").reset_index(drop=True)
    corrected["population"] = corrected["population"].astype(str)
    return corrected


def _rank_biserial_correlation(
    u_statistic: float, responder_n: int, non_responder_n: int
) -> float:
    """Return rank-biserial correlation oriented as responder minus non-responder."""
    return (2 * u_statistic / (responder_n * non_responder_n)) - 1


def _validate_statistical_input(subject_summary: pd.DataFrame) -> None:
    required_columns = {"subject", "response", "population", "mean_percentage"}
    missing_columns = required_columns - set(subject_summary.columns)
    if missing_columns:
        raise StatisticalAnalysisError(
            "Statistical input is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )
    if subject_summary.empty:
        raise StatisticalAnalysisError("Statistical input is empty.")
    duplicate_rows = subject_summary.duplicated(["subject", "population"])
    if duplicate_rows.any():
        examples = subject_summary.loc[
            duplicate_rows, ["subject", "population"]
        ].head(10)
        raise StatisticalAnalysisError(
            "Statistical input must contain one row per subject/population. "
            + str(examples.to_dict(orient="records"))
        )
    if set(subject_summary["population"]) != set(EXPECTED_POPULATIONS):
        raise StatisticalAnalysisError(
            "Statistical input must contain exactly the expected populations."
        )
    invalid_response = ~subject_summary["response"].isin(["yes", "no"])
    if invalid_response.any():
        raise StatisticalAnalysisError("Statistical input must only contain yes/no.")
