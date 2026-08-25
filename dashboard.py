"""Streamlit dashboard for the Teiko technical assessment."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from src.analysis import (
    build_response_frequency_data,
    build_subject_level_response_summary,
)
from src.database import get_database_path, get_repo_root
from src.queries import (
    get_melanoma_miraclib_pbmc_cell_counts,
    get_overview_counts,
    get_sample_counts_by_condition,
    get_sample_counts_by_sample_type,
    get_sample_counts_by_timepoint,
    get_sample_counts_by_treatment,
)
from src.visualization import POPULATION_ORDER, create_response_boxplot, response_label


REQUIRED_OUTPUTS = [
    "relative_frequencies.csv",
    "statistical_results.csv",
    "baseline_statistical_results.csv",
    "baseline_samples.csv",
    "project_counts.csv",
    "response_counts.csv",
    "sex_counts.csv",
    "baseline_b_cell_average.csv",
]


def expected_output_paths(repo_root: Path) -> dict[str, Path]:
    """Return expected pipeline output paths keyed by filename."""
    output_dir = repo_root / "outputs"
    return {filename: output_dir / filename for filename in REQUIRED_OUTPUTS}


def missing_required_files(repo_root: Path) -> list[Path]:
    """Return required dashboard inputs that are currently missing."""
    paths = [get_database_path(), *expected_output_paths(repo_root).values()]
    return [path for path in paths if not path.exists()]


def significance_summary(results: pd.DataFrame) -> str:
    """Return a concise significance interpretation for a results table."""
    if results.empty or "significant" not in results.columns:
        return "No statistical result rows are available."
    significant_mask = results["significant"].map(_is_significant)
    significant = results[significant_mask]
    if significant.empty:
        return (
            "No immune-cell populations showed a statistically significant "
            "responder/non-responder difference after Benjamini-Hochberg "
            "correction at alpha = 0.05."
        )
    populations = ", ".join(significant["population"].astype(str).tolist())
    return (
        "Significant populations after Benjamini-Hochberg correction at "
        f"alpha = 0.05: {populations}."
    )


def _is_significant(value: object) -> bool:
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return bool(value)


def metric_value(metrics: pd.DataFrame, metric: str) -> str | None:
    """Return a metric value from a two-column metric/value dataframe."""
    matches = metrics.loc[metrics["metric"] == metric, "value"]
    if matches.empty:
        return None
    return str(matches.iloc[0])


def count_metric_value(metrics: pd.DataFrame, metric: str) -> str | None:
    """Return a metric count formatted without a decimal point."""
    value = metric_value(metrics, metric)
    if value is None:
        return None
    return f"{int(float(value)):,}"


@st.cache_data
def load_csv(path: str) -> pd.DataFrame:
    """Load a CSV file with Streamlit caching."""
    return pd.read_csv(path)


@st.cache_data
def load_overview_data(db_path: str) -> dict[str, object]:
    """Load dashboard overview metrics and distributions from SQLite."""
    connection = sqlite3.connect(db_path)
    try:
        return {
            "counts": get_overview_counts(connection),
            "conditions": get_sample_counts_by_condition(connection),
            "treatments": get_sample_counts_by_treatment(connection),
            "sample_types": get_sample_counts_by_sample_type(connection),
            "timepoints": get_sample_counts_by_timepoint(connection),
        }
    finally:
        connection.close()


@st.cache_data
def load_response_summaries(db_path: str) -> dict[str, pd.DataFrame]:
    """Load Part 3 primary and baseline subject-level summaries."""
    connection = sqlite3.connect(db_path)
    try:
        raw_counts = get_melanoma_miraclib_pbmc_cell_counts(connection)
    finally:
        connection.close()

    response_frequencies = build_response_frequency_data(raw_counts)
    primary = build_subject_level_response_summary(response_frequencies)
    baseline = build_subject_level_response_summary(
        response_frequencies,
        baseline_only=True,
    )
    return {
        "response_frequencies": response_frequencies,
        "primary": primary,
        "baseline": baseline,
    }


def main() -> None:
    """Render the Streamlit dashboard."""
    st.set_page_config(page_title="Teiko Immune Cell Analytics", layout="wide")

    repo_root = get_repo_root()
    db_path = get_database_path()
    missing_files = missing_required_files(repo_root)

    st.title("Teiko Immune Cell Analytics")
    st.caption(
        "Interactive exploration of immune-cell population data from the "
        "clinical-trial dataset."
    )

    with st.sidebar:
        st.header("Pipeline Status")
        if missing_files:
            st.error("Pipeline outputs missing")
        else:
            st.success("Pipeline data available")
        st.caption("Run `make pipeline` to rebuild the database and outputs.")

    if missing_files:
        st.error("Run `make pipeline` before launching the dashboard.")
        st.write("Missing required files:")
        for path in missing_files:
            st.code(str(path.relative_to(repo_root)))
        st.stop()

    paths = expected_output_paths(repo_root)
    data = {
        "relative_frequencies": load_csv(str(paths["relative_frequencies.csv"])),
        "statistical_results": load_csv(str(paths["statistical_results.csv"])),
        "baseline_statistical_results": load_csv(
            str(paths["baseline_statistical_results.csv"])
        ),
        "baseline_samples": load_csv(str(paths["baseline_samples.csv"])),
        "project_counts": load_csv(str(paths["project_counts.csv"])),
        "response_counts": load_csv(str(paths["response_counts.csv"])),
        "sex_counts": load_csv(str(paths["sex_counts.csv"])),
        "b_cell_average": load_csv(str(paths["baseline_b_cell_average.csv"])),
    }

    overview_data = load_overview_data(str(db_path))
    response_data = load_response_summaries(str(db_path))

    tabs = st.tabs(
        ["Overview", "Cell Frequencies", "Treatment Response", "Baseline Cohort"]
    )
    with tabs[0]:
        render_overview_tab(overview_data)
    with tabs[1]:
        render_cell_frequencies_tab(data["relative_frequencies"])
    with tabs[2]:
        render_treatment_response_tab(
            response_data,
            data["statistical_results"],
            data["baseline_statistical_results"],
        )
    with tabs[3]:
        render_baseline_cohort_tab(
            data["baseline_samples"],
            data["project_counts"],
            data["response_counts"],
            data["sex_counts"],
            data["b_cell_average"],
        )


def render_overview_tab(overview_data: dict[str, object]) -> None:
    """Render dataset overview metrics and distribution charts."""
    counts = overview_data["counts"]
    st.subheader("Dataset Overview")
    cols = st.columns(4)
    cols[0].metric("Samples", f"{counts['samples']:,}")
    cols[1].metric("Subjects", f"{counts['subjects']:,}")
    cols[2].metric("Projects", f"{counts['projects']:,}")
    cols[3].metric("Cell Populations", f"{counts['cell_populations']:,}")

    left, right = st.columns(2)
    with left:
        st.plotly_chart(
            _bar_chart(
                overview_data["conditions"],
                x="condition",
                y="sample_count",
                title="Samples by Condition",
                x_label="Condition",
                y_label="Sample Count",
            ),
            use_container_width=True,
        )
        st.plotly_chart(
            _bar_chart(
                overview_data["sample_types"],
                x="sample_type",
                y="sample_count",
                title="Samples by Sample Type",
                x_label="Sample Type",
                y_label="Sample Count",
            ),
            use_container_width=True,
        )
    with right:
        st.plotly_chart(
            _bar_chart(
                overview_data["treatments"],
                x="treatment",
                y="sample_count",
                title="Samples by Treatment",
                x_label="Treatment",
                y_label="Sample Count",
            ),
            use_container_width=True,
        )
        st.plotly_chart(
            _bar_chart(
                overview_data["timepoints"],
                x="time_from_treatment_start",
                y="sample_count",
                title="Samples by Timepoint",
                x_label="Days from Treatment Start",
                y_label="Sample Count",
            ),
            use_container_width=True,
        )


def render_cell_frequencies_tab(frequencies: pd.DataFrame) -> None:
    """Render Part 2 relative-frequency exploration."""
    st.subheader("Cell Population Relative Frequencies")
    selected_populations = st.multiselect(
        "Population",
        options=POPULATION_ORDER,
        default=POPULATION_ORDER,
    )
    sample_search = st.text_input("Sample search", placeholder="e.g. sample00000")
    filtered = frequencies[frequencies["population"].isin(selected_populations)]
    if sample_search:
        filtered = filtered[
            filtered["sample"].astype(str).str.contains(sample_search, case=False)
        ]

    if filtered.empty:
        st.warning("No frequency rows match the selected filters.")
        return

    overall, selected = st.tabs(["Overall Distribution", "Selected Sample"])
    with overall:
        fig = px.box(
            filtered,
            x="population",
            y="percentage",
            points=False,
            category_orders={"population": POPULATION_ORDER},
            labels={
                "population": "Immune Cell Population",
                "percentage": "Relative Frequency (%)",
            },
            title="Relative Frequency Distribution Across Samples",
        )
        st.plotly_chart(fig, use_container_width=True)

    with selected:
        sample_options = filtered["sample"].drop_duplicates().head(200).tolist()
        selected_sample = st.selectbox("Sample", options=sample_options)
        sample_df = filtered[filtered["sample"] == selected_sample]
        fig = px.bar(
            sample_df,
            x="population",
            y="percentage",
            category_orders={"population": POPULATION_ORDER},
            labels={
                "population": "Immune Cell Population",
                "percentage": "Relative Frequency (%)",
            },
            title=f"Selected Sample Composition: {selected_sample}",
        )
        st.plotly_chart(fig, use_container_width=True)

    st.dataframe(
        filtered,
        width="stretch",
        column_config={
            "percentage": st.column_config.NumberColumn(
                "percentage",
                format="%.2f",
            )
        },
    )


def render_treatment_response_tab(
    response_data: dict[str, pd.DataFrame],
    primary_results: pd.DataFrame,
    baseline_results: pd.DataFrame,
) -> None:
    """Render Part 3 responder/non-responder analysis."""
    st.subheader("Treatment Response Analysis")
    st.markdown(
        "`condition = melanoma`, `treatment = miraclib`, `sample_type = PBMC`, "
        "`response = yes/no`"
    )

    response_frequencies = response_data["response_frequencies"]
    subject_rows = response_data["primary"][["subject", "response"]].drop_duplicates()
    subject_counts = subject_rows["response"].value_counts()
    cols = st.columns(4)
    cols[0].metric("Subjects", f"{subject_rows['subject'].nunique():,}")
    cols[1].metric("Responders", f"{int(subject_counts.get('yes', 0)):,}")
    cols[2].metric("Non-responders", f"{int(subject_counts.get('no', 0)):,}")
    cols[3].metric("Samples", f"{response_frequencies['sample'].nunique():,}")

    analysis_choice = st.radio(
        "Analysis view",
        options=[
            "Primary longitudinal subject-level analysis",
            "Baseline-only analysis",
        ],
        horizontal=True,
    )
    is_baseline = analysis_choice.startswith("Baseline")
    summary = response_data["baseline"] if is_baseline else response_data["primary"]
    results = baseline_results if is_baseline else primary_results
    title = (
        "Baseline Immune-Cell Frequencies by Treatment Response"
        if is_baseline
        else "Subject-Level Immune-Cell Frequencies by Treatment Response"
    )

    if is_baseline:
        st.info(
            "The baseline-only view is the more appropriate comparison when "
            "considering potential pre-treatment markers of response, because "
            "post-treatment measurements are not available before therapy begins."
        )
    else:
        st.info(
            "The primary view averages each subject's available longitudinal "
            "PBMC/miraclib measurements before comparing response groups."
        )

    st.plotly_chart(
        create_response_boxplot(summary, title=title),
        use_container_width=True,
    )
    st.write(significance_summary(results))
    st.dataframe(
        results,
        width="stretch",
        column_config={
            "responder_median_percentage": st.column_config.NumberColumn(format="%.3f"),
            "non_responder_median_percentage": st.column_config.NumberColumn(
                format="%.3f"
            ),
            "median_difference": st.column_config.NumberColumn(format="%.3f"),
            "mann_whitney_u": st.column_config.NumberColumn(format="%.1f"),
            "p_value": st.column_config.NumberColumn(format="%.4g"),
            "adjusted_p_value": st.column_config.NumberColumn(format="%.4g"),
            "effect_size": st.column_config.NumberColumn(format="%.3f"),
        },
    )


def render_baseline_cohort_tab(
    baseline_samples: pd.DataFrame,
    project_counts: pd.DataFrame,
    response_counts: pd.DataFrame,
    sex_counts: pd.DataFrame,
    b_cell_average: pd.DataFrame,
) -> None:
    """Render Part 4 baseline subset analysis."""
    st.subheader("Baseline Cohort")
    st.markdown(
        "`condition = melanoma`, `treatment = miraclib`, `sample_type = PBMC`, "
        "`time_from_treatment_start = 0`"
    )

    cols = st.columns(3)
    cols[0].metric("Baseline Samples", f"{len(baseline_samples):,}")
    cols[1].metric("Unique Subjects", f"{baseline_samples['subject'].nunique():,}")
    cols[2].metric("Projects Represented", f"{baseline_samples['project'].nunique():,}")

    b_cell_value = metric_value(
        b_cell_average,
        "melanoma_male_responder_baseline_avg_b_cell",
    )
    matching_samples = count_metric_value(b_cell_average, "matching_baseline_samples")
    st.metric("Average Baseline B-cell Count", b_cell_value or "Unavailable")
    st.caption(
        "Final B-cell filter: melanoma males, response = yes, time = 0, "
        "all sample types, all treatment types."
    )
    if matching_samples:
        st.caption(f"Matching baseline samples used: {matching_samples}")

    left, middle, right = st.columns(3)
    with left:
        st.plotly_chart(
            _bar_chart(
                project_counts,
                x="project",
                y="sample_count",
                title="Baseline Samples by Project",
                x_label="Project",
                y_label="Sample Count",
            ),
            use_container_width=True,
        )
    with middle:
        response_plot = response_counts.copy()
        response_plot["response_label"] = response_plot["response"].map(response_label)
        st.plotly_chart(
            _bar_chart(
                response_plot,
                x="response_label",
                y="subject_count",
                title="Baseline Subjects by Response",
                x_label="Treatment Response",
                y_label="Distinct Subject Count",
            ),
            use_container_width=True,
        )
    with right:
        st.plotly_chart(
            _bar_chart(
                sex_counts,
                x="sex",
                y="subject_count",
                title="Baseline Subjects by Sex",
                x_label="Sex",
                y_label="Distinct Subject Count",
            ),
            use_container_width=True,
        )

    st.markdown("#### Baseline Sample Table")
    filtered = baseline_samples.copy()
    projects = st.multiselect(
        "Project",
        options=sorted(filtered["project"].dropna().unique()),
        default=sorted(filtered["project"].dropna().unique()),
    )
    responses = st.multiselect(
        "Response",
        options=sorted(filtered["response"].dropna().unique()),
        default=sorted(filtered["response"].dropna().unique()),
    )
    sexes = st.multiselect(
        "Sex",
        options=sorted(filtered["sex"].dropna().unique()),
        default=sorted(filtered["sex"].dropna().unique()),
    )
    filtered = filtered[
        filtered["project"].isin(projects)
        & filtered["response"].isin(responses)
        & filtered["sex"].isin(sexes)
    ]
    st.dataframe(filtered, width="stretch")


def _bar_chart(df: pd.DataFrame, *, x: str, y: str, title: str, x_label: str, y_label: str):
    fig = px.bar(
        df,
        x=x,
        y=y,
        labels={x: x_label, y: y_label},
        title=title,
    )
    fig.update_layout(showlegend=False)
    return fig


if __name__ == "__main__":
    main()
