"""Reusable Plotly visualizations for the Teiko assessment."""

from __future__ import annotations

import plotly.express as px


POPULATION_ORDER = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
RESPONSE_LABELS = {"yes": "Responder", "no": "Non-responder"}


def response_label(value: str) -> str:
    """Return a human-readable treatment-response label."""
    return RESPONSE_LABELS.get(value, "Unknown")


def create_response_boxplot(summary_df, *, title: str):
    """Create a responder/non-responder boxplot from subject-level data."""
    plot_df = summary_df.copy()
    plot_df["response_label"] = plot_df["response"].map(response_label)
    fig = px.box(
        plot_df,
        x="population",
        y="mean_percentage",
        color="response_label",
        points="all",
        category_orders={
            "population": POPULATION_ORDER,
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
    return fig
