"""Run analysis phases for the Teiko technical assessment."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import plotly.express as px

from src.analysis import (
    build_response_frequency_data,
    build_subject_level_response_summary,
    calculate_relative_frequencies,
    validate_relative_frequencies,
)
from src.database import DATABASE_FILENAME, connect_database, get_database_path, get_repo_root
from src.queries import (
    get_baseline_melanoma_male_responder_b_cell_average,
    get_baseline_melanoma_male_responder_b_cells,
    get_baseline_melanoma_miraclib_pbmc_samples,
    get_baseline_sample_counts_by_project,
    get_baseline_subject_counts_by_response,
    get_baseline_subject_counts_by_sex,
    get_cell_counts_by_sample,
    get_melanoma_miraclib_pbmc_cell_counts,
    get_population_names,
    get_sample_ids,
)
from src.statistics import compare_response_groups


RELATIVE_FREQUENCIES_FILENAME = "relative_frequencies.csv"
B_CELL_AVERAGE_METRIC = "melanoma_male_responder_baseline_avg_b_cell"
RESPONSE_LABELS = {"yes": "Responder", "no": "Non-responder"}


def _format_int(value: int) -> str:
    return f"{value:,}"


def _format_float(value: float) -> str:
    return f"{value:.2f}"


def run_part2(
    connection: sqlite3.Connection, output_dir: Path, repo_root: Path
) -> None:
    """Generate the Part 2 relative-frequency output."""
    output_path = output_dir / RELATIVE_FREQUENCIES_FILENAME
    print("Running Part 2: cell population relative frequencies...")
    print("")

    sample_ids = get_sample_ids(connection)
    population_names = get_population_names(connection)
    cell_counts = get_cell_counts_by_sample(connection)

    frequencies = calculate_relative_frequencies(cell_counts)
    validate_relative_frequencies(
        frequencies,
        expected_samples=sample_ids,
        expected_populations=population_names,
    )

    frequencies.to_csv(output_path, index=False)

    print(f"Samples analyzed:      {_format_int(len(sample_ids)):>7}")
    print(f"Populations:          {_format_int(len(population_names)):>7}")
    print(f"Output rows:          {_format_int(len(frequencies)):>7}")
    print("")
    print("Validated percentage sums for all samples.")
    print(f"Wrote {output_path.relative_to(repo_root)}")
    print("")


def run_part3(
    connection: sqlite3.Connection, output_dir: Path, repo_root: Path
) -> dict[str, object]:
    """Generate Part 3 response statistical outputs and plots."""
    print("Running Part 3: treatment response analysis...")
    print("")

    raw_counts = get_melanoma_miraclib_pbmc_cell_counts(connection)
    response_frequencies = build_response_frequency_data(raw_counts)
    subject_summary = build_subject_level_response_summary(response_frequencies)
    baseline_summary = build_subject_level_response_summary(
        response_frequencies,
        baseline_only=True,
    )

    statistical_results = compare_response_groups(subject_summary)
    baseline_results = compare_response_groups(baseline_summary)

    statistical_results.to_csv(output_dir / "statistical_results.csv", index=False)
    baseline_results.to_csv(
        output_dir / "baseline_statistical_results.csv",
        index=False,
    )
    _write_response_boxplot(
        subject_summary,
        output_dir / "responder_boxplot.html",
        title="Subject-Level Immune-Cell Frequencies by Treatment Response",
    )
    _write_response_boxplot(
        baseline_summary,
        output_dir / "baseline_responder_boxplot.html",
        title="Baseline Immune-Cell Frequencies by Treatment Response",
    )

    subject_counts = (
        subject_summary[["subject", "response"]].drop_duplicates()["response"].value_counts()
    )
    sample_count = response_frequencies["sample"].nunique()
    timepoints = sorted(response_frequencies["time_from_treatment_start"].unique())
    significant = _significant_populations(statistical_results)
    baseline_significant = _significant_populations(baseline_results)

    print("Target cohort:")
    print(f"  Subjects: {_format_int(subject_summary['subject'].nunique())}")
    print(f"  Responders: {_format_int(int(subject_counts.get('yes', 0)))}")
    print(f"  Non-responders: {_format_int(int(subject_counts.get('no', 0)))}")
    print(f"  Samples: {_format_int(sample_count)}")
    print(f"  Timepoints: {', '.join(str(int(value)) for value in timepoints)}")
    print("")
    print("Primary subject-level analysis:")
    print(
        "  Significant populations after BH correction: "
        + (", ".join(significant) if significant else "none")
    )
    print("")
    print("Baseline analysis:")
    print(
        "  Significant populations after BH correction: "
        + (", ".join(baseline_significant) if baseline_significant else "none")
    )
    print("")
    print(f"Wrote {(output_dir / 'statistical_results.csv').relative_to(repo_root)}")
    print(
        f"Wrote {(output_dir / 'baseline_statistical_results.csv').relative_to(repo_root)}"
    )
    print(f"Wrote {(output_dir / 'responder_boxplot.html').relative_to(repo_root)}")
    print(
        f"Wrote {(output_dir / 'baseline_responder_boxplot.html').relative_to(repo_root)}"
    )
    print("")

    return {
        "response_frequencies": response_frequencies,
        "subject_summary": subject_summary,
        "baseline_summary": baseline_summary,
        "statistical_results": statistical_results,
        "baseline_results": baseline_results,
        "subject_counts": subject_counts.to_dict(),
        "sample_count": sample_count,
        "timepoints": timepoints,
    }


def _write_response_boxplot(summary, output_path: Path, *, title: str) -> None:
    plot_df = summary.copy()
    plot_df["response_label"] = plot_df["response"].map(RESPONSE_LABELS)
    fig = px.box(
        plot_df,
        x="population",
        y="mean_percentage",
        color="response_label",
        points="all",
        category_orders={
            "population": ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"],
            "response_label": ["Responder", "Non-responder"],
        },
        labels={
            "population": "Immune Cell Population",
            "mean_percentage": "Relative Frequency (%)",
            "response_label": "Treatment Response",
        },
        title=title,
    )
    fig.update_layout(boxmode="group")
    fig.write_html(output_path, include_plotlyjs="cdn")


def _significant_populations(results) -> list[str]:
    return results.loc[results["significant"], "population"].tolist()


def run_part4(
    connection: sqlite3.Connection, output_dir: Path, repo_root: Path
) -> dict[str, object]:
    """Generate the Part 4 baseline subset outputs."""
    print("Running Part 4: baseline subset analysis...")
    print("")

    baseline_samples = get_baseline_melanoma_miraclib_pbmc_samples(connection)
    project_counts = get_baseline_sample_counts_by_project(connection)
    response_counts = get_baseline_subject_counts_by_response(connection)
    sex_counts = get_baseline_subject_counts_by_sex(connection)
    b_cell_rows = get_baseline_melanoma_male_responder_b_cells(connection)
    b_cell_average = get_baseline_melanoma_male_responder_b_cell_average(connection)

    validate_part4_outputs(
        baseline_samples,
        project_counts,
        response_counts,
        sex_counts,
        b_cell_rows,
        b_cell_average,
    )

    formatted_average = _format_float(b_cell_average)
    average_output = _build_b_cell_average_output(
        matching_samples=len(b_cell_rows),
        formatted_average=formatted_average,
    )

    baseline_samples.to_csv(output_dir / "baseline_samples.csv", index=False)
    project_counts.to_csv(output_dir / "project_counts.csv", index=False)
    response_counts.to_csv(output_dir / "response_counts.csv", index=False)
    sex_counts.to_csv(output_dir / "sex_counts.csv", index=False)
    average_output.to_csv(output_dir / "baseline_b_cell_average.csv", index=False)

    print(
        "Baseline melanoma/miraclib/PBMC samples: "
        f"{_format_int(len(baseline_samples))}"
    )
    print("")
    print("Samples by project:")
    for row in project_counts.itertuples(index=False):
        print(f"  {row.project}: {_format_int(int(row.sample_count))}")
    print("")
    print("Subjects by response:")
    for row in response_counts.itertuples(index=False):
        print(f"  {row.response}: {_format_int(int(row.subject_count))}")
    print("")
    print("Subjects by sex:")
    for row in sex_counts.itertuples(index=False):
        print(f"  {row.sex}: {_format_int(int(row.subject_count))}")
    print("")
    print("Melanoma male responders at baseline:")
    print(f"  Matching samples: {_format_int(len(b_cell_rows))}")
    print(f"  Average B-cell count: {formatted_average}")
    print("")
    print("Part 4 outputs written successfully.")

    return {
        "baseline_samples": baseline_samples,
        "project_counts": project_counts,
        "response_counts": response_counts,
        "sex_counts": sex_counts,
        "b_cell_matching_samples": len(b_cell_rows),
        "b_cell_average": b_cell_average,
        "b_cell_average_formatted": formatted_average,
    }


def validate_part4_outputs(
    baseline_samples,
    project_counts,
    response_counts,
    sex_counts,
    b_cell_rows,
    b_cell_average: float | None,
) -> None:
    """Validate Part 4 query outputs before writing CSV files."""
    if baseline_samples["sample"].duplicated().any():
        examples = baseline_samples.loc[
            baseline_samples["sample"].duplicated(), "sample"
        ].head(10)
        raise RuntimeError(
            "Baseline cohort sample IDs must be unique. Examples: "
            + ", ".join(examples.astype(str))
        )

    expected_filters = {
        "condition": "melanoma",
        "treatment": "miraclib",
        "sample_type": "PBMC",
        "time_from_treatment_start": 0,
    }
    for column, expected in expected_filters.items():
        invalid = baseline_samples[column] != expected
        if invalid.any():
            bad_values = baseline_samples.loc[invalid, column].astype(str).unique()[:10]
            raise RuntimeError(
                f"Baseline cohort has invalid {column} values: "
                + ", ".join(bad_values)
            )

    baseline_subjects = baseline_samples["subject"].nunique()
    if int(project_counts["sample_count"].sum()) != len(baseline_samples):
        raise RuntimeError("Project sample counts do not reconcile to baseline samples.")
    if int(response_counts["subject_count"].sum()) != baseline_subjects:
        raise RuntimeError(
            "Response subject counts do not reconcile to distinct baseline subjects."
        )
    if int(sex_counts["subject_count"].sum()) != baseline_subjects:
        raise RuntimeError(
            "Sex subject counts do not reconcile to distinct baseline subjects."
        )

    _validate_b_cell_rows(b_cell_rows)
    if b_cell_rows.empty:
        raise RuntimeError("Final B-cell average query returned no matching samples.")
    calculated_average = float(b_cell_rows["b_cell_count"].mean())
    if b_cell_average is None or abs(calculated_average - b_cell_average) > 1e-9:
        raise RuntimeError("Final B-cell average does not match underlying rows.")


def _validate_b_cell_rows(b_cell_rows) -> None:
    checks = {
        "condition": "melanoma",
        "sex": "M",
        "response": "yes",
        "time_from_treatment_start": 0,
        "population": "b_cell",
    }
    for column, expected in checks.items():
        invalid = b_cell_rows[column] != expected
        if invalid.any():
            bad_values = b_cell_rows.loc[invalid, column].astype(str).unique()[:10]
            raise RuntimeError(
                f"Final B-cell query has invalid {column} values: "
                + ", ".join(bad_values)
            )


def _build_b_cell_average_output(
    *, matching_samples: int, formatted_average: str
):
    import pandas as pd

    return pd.DataFrame(
        [
            {
                "metric": "matching_baseline_samples",
                "value": str(matching_samples),
            },
            {
                "metric": B_CELL_AVERAGE_METRIC,
                "value": formatted_average,
            },
        ]
    )


def main() -> None:
    """Run the current analysis workflow."""
    repo_root = get_repo_root()
    db_path = get_database_path()
    output_dir = repo_root / "outputs"

    if not db_path.exists():
        raise FileNotFoundError(
            f"Could not find {DATABASE_FILENAME} at {db_path}. "
            "Run `python load_data.py` before running analysis."
        )

    output_dir.mkdir(exist_ok=True)

    connection = connect_database(db_path)
    try:
        run_part2(connection, output_dir, repo_root)
        run_part3(connection, output_dir, repo_root)
        run_part4(connection, output_dir, repo_root)
    finally:
        connection.close()


if __name__ == "__main__":
    main()
