"""Phase 3 tests for Part 4 database subset queries."""

from __future__ import annotations

import pytest

from src.database import connect_database, initialize_schema
from src.queries import (
    get_baseline_melanoma_male_responder_b_cell_average,
    get_baseline_melanoma_male_responder_b_cells,
    get_baseline_melanoma_miraclib_pbmc_samples,
    get_baseline_sample_counts_by_project,
    get_baseline_subject_counts_by_response,
    get_baseline_subject_counts_by_sex,
    get_melanoma_miraclib_pbmc_cell_counts,
)


POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def _insert_subject(
    connection,
    *,
    subject: str,
    project_id: int,
    condition: str,
    sex: str,
    treatment: str,
    response: str | None,
) -> None:
    connection.execute(
        """
        INSERT INTO subjects (
            subject_id, project_id, condition, age, sex, treatment, response
        )
        VALUES (?, ?, ?, 50, ?, ?, ?)
        """,
        (subject, project_id, condition, sex, treatment, response),
    )


def _insert_sample(
    connection,
    *,
    sample: str,
    subject: str,
    sample_type: str,
    timepoint: int,
    b_cell_count: int,
    other_count: int = 1,
) -> None:
    connection.execute(
        """
        INSERT INTO samples (
            sample_id, subject_id, sample_type, time_from_treatment_start
        )
        VALUES (?, ?, ?, ?)
        """,
        (sample, subject, sample_type, timepoint),
    )
    population_counts = {
        "b_cell": b_cell_count,
        "cd8_t_cell": other_count,
        "cd4_t_cell": other_count,
        "nk_cell": other_count,
        "monocyte": other_count,
    }
    connection.executemany(
        """
        INSERT INTO cell_counts (sample_id, population_id, count)
        VALUES (?, ?, ?)
        """,
        [
            (sample, index + 1, population_counts[population])
            for index, population in enumerate(POPULATIONS)
        ],
    )


@pytest.fixture()
def part4_connection(tmp_path):
    db_path = tmp_path / "part4.db"
    connection = connect_database(db_path)
    initialize_schema(connection)
    connection.executemany(
        "INSERT INTO projects (project_id, project_name) VALUES (?, ?)",
        [(1, "prj1"), (2, "prj2"), (3, "prj3")],
    )
    connection.executemany(
        "INSERT INTO cell_populations (population_id, name) VALUES (?, ?)",
        [(index + 1, population) for index, population in enumerate(POPULATIONS)],
    )

    subjects = [
        ("subj_yes", 1, "melanoma", "M", "miraclib", "yes"),
        ("subj_yes_repeat", 1, "melanoma", "F", "miraclib", "yes"),
        ("subj_no", 2, "melanoma", "F", "miraclib", "no"),
        ("subj_null", 2, "melanoma", "M", "miraclib", None),
        ("subj_wrong_condition", 1, "carcinoma", "M", "miraclib", "yes"),
        ("subj_wrong_treatment", 1, "melanoma", "M", "phauximab", "yes"),
        ("subj_all_treat_none", 2, "melanoma", "M", "none", "yes"),
        ("subj_all_treat_phaux", 3, "melanoma", "M", "phauximab", "yes"),
        ("subj_female_responder", 3, "melanoma", "F", "none", "yes"),
        ("subj_male_nonresponder", 3, "melanoma", "M", "none", "no"),
    ]
    for subject, project_id, condition, sex, treatment, response in subjects:
        _insert_subject(
            connection,
            subject=subject,
            project_id=project_id,
            condition=condition,
            sex=sex,
            treatment=treatment,
            response=response,
        )

    samples = [
        ("sample_yes", "subj_yes", "PBMC", 0, 100),
        ("sample_yes_repeat_a", "subj_yes_repeat", "PBMC", 0, 110),
        ("sample_yes_repeat_b", "subj_yes_repeat", "PBMC", 0, 115),
        ("sample_no", "subj_no", "PBMC", 0, 200),
        ("sample_null", "subj_null", "PBMC", 0, 300),
        ("sample_wrong_condition", "subj_wrong_condition", "PBMC", 0, 400),
        ("sample_wrong_treatment", "subj_wrong_treatment", "PBMC", 0, 500),
        ("sample_wrong_type", "subj_yes", "WB", 0, 600),
        ("sample_wrong_time", "subj_yes", "PBMC", 7, 700),
        ("sample_none_pbmc", "subj_all_treat_none", "PBMC", 0, 800),
        ("sample_phaux_wb", "subj_all_treat_phaux", "WB", 0, 1000),
        ("sample_female_responder", "subj_female_responder", "PBMC", 0, 1200),
        ("sample_male_nonresponder", "subj_male_nonresponder", "PBMC", 0, 1400),
    ]
    for sample, subject, sample_type, timepoint, b_cell_count in samples:
        _insert_sample(
            connection,
            sample=sample,
            subject=subject,
            sample_type=sample_type,
            timepoint=timepoint,
            b_cell_count=b_cell_count,
            other_count=9999,
        )

    try:
        yield connection
    finally:
        connection.close()


def test_baseline_filtering_selects_only_target_cohort(part4_connection) -> None:
    cohort = get_baseline_melanoma_miraclib_pbmc_samples(part4_connection)

    assert set(cohort["sample"]) == {
        "sample_yes",
        "sample_yes_repeat_a",
        "sample_yes_repeat_b",
        "sample_no",
        "sample_null",
    }
    assert set(cohort["condition"]) == {"melanoma"}
    assert set(cohort["treatment"]) == {"miraclib"}
    assert set(cohort["sample_type"]) == {"PBMC"}
    assert set(cohort["time_from_treatment_start"]) == {0}


def test_project_counts_are_sample_counts(part4_connection) -> None:
    counts = get_baseline_sample_counts_by_project(part4_connection)

    assert dict(zip(counts["project"], counts["sample_count"])) == {
        "prj1": 3,
        "prj2": 2,
    }


def test_response_counts_use_distinct_subjects_not_samples(part4_connection) -> None:
    counts = get_baseline_subject_counts_by_response(part4_connection)

    assert dict(zip(counts["response"], counts["subject_count"])) == {
        "yes": 2,
        "no": 1,
        "unknown": 1,
    }


def test_response_counts_keep_null_distinct_from_no(part4_connection) -> None:
    counts = get_baseline_subject_counts_by_response(part4_connection)
    result = dict(zip(counts["response"], counts["subject_count"]))

    assert result["no"] == 1
    assert result["unknown"] == 1


def test_sex_counts_group_distinct_subjects(part4_connection) -> None:
    counts = get_baseline_subject_counts_by_sex(part4_connection)

    assert dict(zip(counts["sex"], counts["subject_count"])) == {"F": 2, "M": 2}


def test_final_b_cell_filter_includes_only_melanoma_male_responders_at_time0(
    part4_connection,
) -> None:
    rows = get_baseline_melanoma_male_responder_b_cells(part4_connection)

    assert set(rows["subject"]) == {
        "subj_yes",
        "subj_wrong_treatment",
        "subj_all_treat_none",
        "subj_all_treat_phaux",
    }
    assert set(rows["condition"]) == {"melanoma"}
    assert set(rows["sex"]) == {"M"}
    assert set(rows["response"]) == {"yes"}
    assert set(rows["time_from_treatment_start"]) == {0}
    assert set(rows["population"]) == {"b_cell"}


def test_final_b_cell_average_uses_all_treatments_and_sample_types(
    part4_connection,
) -> None:
    rows = get_baseline_melanoma_male_responder_b_cells(part4_connection)
    average = get_baseline_melanoma_male_responder_b_cell_average(part4_connection)

    assert set(rows["treatment"]) == {"miraclib", "none", "phauximab"}
    assert set(rows["sample_type"]) == {"PBMC", "WB"}
    assert set(rows["b_cell_count"]) == {100, 500, 600, 800, 1000}
    assert average == pytest.approx((100 + 500 + 600 + 800 + 1000) / 5)


def test_final_b_cell_average_ignores_other_population_values(
    part4_connection,
) -> None:
    average = get_baseline_melanoma_male_responder_b_cell_average(part4_connection)

    assert average == pytest.approx(600.0)


def test_two_decimal_formatting_does_not_change_underlying_average(
    part4_connection,
) -> None:
    average = get_baseline_melanoma_male_responder_b_cell_average(part4_connection)

    assert average == pytest.approx(600.0)
    assert f"{average:.2f}" == "600.00"


def test_part3_response_cohort_excludes_missing_response_and_wrong_filters(
    part4_connection,
) -> None:
    rows = get_melanoma_miraclib_pbmc_cell_counts(part4_connection)

    assert set(rows["subject"]) == {"subj_yes", "subj_yes_repeat", "subj_no"}
    assert set(rows["condition"]) == {"melanoma"}
    assert set(rows["treatment"]) == {"miraclib"}
    assert set(rows["sample_type"]) == {"PBMC"}
    assert set(rows["response"]) == {"yes", "no"}
    assert rows["sample"].nunique() == 5
    assert set(rows["time_from_treatment_start"]) == {0, 7}
    assert len(rows) == 25
