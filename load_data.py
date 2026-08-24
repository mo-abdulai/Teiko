"""Load cell-count.csv into a normalized SQLite database.

Running this script rebuilds the known project database, ``cell_counts.db``,
from scratch. The rebuild behavior is intentional for a reproducible Phase 1
pipeline and avoids duplicate data on repeated runs.
"""

from __future__ import annotations

import os
import sqlite3
from pathlib import Path
from typing import Any


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

from src.database import (
    DATABASE_FILENAME,
    connect_database,
    get_database_path,
    get_repo_root,
    initialize_schema,
    recreate_database,
)


CSV_FILENAME = "cell-count.csv"
CELL_POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
SUBJECT_LEVEL_COLUMNS = [
    "project",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
]
REQUIRED_COLUMNS = [
    *SUBJECT_LEVEL_COLUMNS,
    "subject",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    *CELL_POPULATIONS,
]
EXPECTED_COUNTS = {
    "projects": 3,
    "subjects": 3500,
    "samples": 10500,
    "cell_populations": 5,
    "cell_counts": 52500,
}


class DataValidationError(ValueError):
    """Raised when the CSV does not satisfy Phase 1 loading assumptions."""


def _format_int(value: int) -> str:
    return f"{value:,}"


def _sample_values(values: pd.Series, limit: int = 5) -> str:
    sample = values.astype(str).head(limit).tolist()
    return ", ".join(sample)


def _is_integer_compatible(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    return numeric.notna() & (numeric % 1 == 0)


def read_dataset(csv_path: Path | str) -> pd.DataFrame:
    """Read the assessment CSV file."""
    return pd.read_csv(csv_path)


def validate_dataframe(df: pd.DataFrame) -> None:
    """Validate the source dataframe before any database inserts occur."""
    missing_columns = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing_columns:
        raise DataValidationError(
            "CSV is missing required columns: " + ", ".join(missing_columns)
        )

    if df["subject"].isna().any() or (df["subject"].astype(str).str.strip() == "").any():
        raise DataValidationError("Subject IDs must not be missing.")
    if df["sample"].isna().any() or (df["sample"].astype(str).str.strip() == "").any():
        raise DataValidationError("Sample IDs must not be missing.")

    columns_expected_non_null = [column for column in REQUIRED_COLUMNS if column != "response"]
    null_columns = [
        column for column in columns_expected_non_null if df[column].isna().any()
    ]
    if null_columns:
        raise DataValidationError(
            "Only response may contain missing values; found missing values in: "
            + ", ".join(null_columns)
        )

    duplicate_samples = df.loc[df["sample"].duplicated(), "sample"]
    if not duplicate_samples.empty:
        raise DataValidationError(
            "Sample IDs must be unique. Examples: " + _sample_values(duplicate_samples)
        )

    combo_columns = ["subject", "time_from_treatment_start", "sample_type"]
    duplicate_combos = df.loc[df.duplicated(combo_columns), combo_columns]
    if not duplicate_combos.empty:
        examples = duplicate_combos.head(5).to_dict(orient="records")
        raise DataValidationError(
            "Duplicate (subject, time_from_treatment_start, sample_type) combinations: "
            + str(examples)
        )

    for column in ["age", "time_from_treatment_start"]:
        if not _is_integer_compatible(df[column]).all():
            bad_values = df.loc[~_is_integer_compatible(df[column]), column]
            raise DataValidationError(
                f"{column} must be numeric and integer-compatible. Examples: "
                + _sample_values(bad_values)
            )

    count_values = df[CELL_POPULATIONS].apply(pd.to_numeric, errors="coerce")
    if count_values.isna().any().any():
        bad_columns = count_values.columns[count_values.isna().any()].tolist()
        raise DataValidationError(
            "Cell counts must be numeric. Invalid columns: " + ", ".join(bad_columns)
        )
    if (count_values < 0).any().any():
        bad_columns = count_values.columns[(count_values < 0).any()].tolist()
        raise DataValidationError(
            "Cell counts must be non-negative. Invalid columns: " + ", ".join(bad_columns)
        )
    integer_counts = count_values.apply(lambda column: column % 1 == 0)
    if not integer_counts.all().all():
        bad_columns = integer_counts.columns[~integer_counts.all()].tolist()
        raise DataValidationError(
            "Cell counts must be integer-compatible. Invalid columns: "
            + ", ".join(bad_columns)
        )

    for column in SUBJECT_LEVEL_COLUMNS:
        distinct_per_subject = df.groupby("subject", dropna=False)[column].nunique(
            dropna=False
        )
        conflicting = distinct_per_subject[distinct_per_subject > 1]
        if not conflicting.empty:
            raise DataValidationError(
                f"Subject-level field '{column}' has conflicting values. "
                f"Affected subjects: {_sample_values(conflicting.index.to_series())}"
            )


def prepare_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with types normalized for SQLite insertion."""
    prepared = df.copy()
    prepared["subject"] = prepared["subject"].astype(str)
    prepared["sample"] = prepared["sample"].astype(str)
    prepared["age"] = pd.to_numeric(prepared["age"], errors="raise").astype(int)
    prepared["time_from_treatment_start"] = pd.to_numeric(
        prepared["time_from_treatment_start"], errors="raise"
    ).astype(int)
    for column in CELL_POPULATIONS:
        prepared[column] = pd.to_numeric(prepared[column], errors="raise").astype(int)
    prepared["response"] = prepared["response"].where(prepared["response"].notna(), None)
    return prepared


def load_dataframe(df: pd.DataFrame, db_path: Path | str | None = None) -> dict[str, int]:
    """Validate and load a dataframe into a fresh SQLite database."""
    validate_dataframe(df)
    prepared = prepare_dataframe(df)

    target_path = recreate_database(db_path)
    connection = connect_database(target_path)
    try:
        with connection:
            initialize_schema(connection)
            counts = _insert_normalized_records(connection, prepared)
        validate_loaded_data(connection, prepared, counts)
    finally:
        connection.close()

    return counts


def _insert_normalized_records(
    connection: sqlite3.Connection, df: pd.DataFrame
) -> dict[str, int]:
    projects = sorted(df["project"].unique().tolist())
    connection.executemany(
        "INSERT INTO projects (project_name) VALUES (?)",
        [(project,) for project in projects],
    )
    project_lookup = dict(
        connection.execute("SELECT project_name, project_id FROM projects").fetchall()
    )

    subject_df = (
        df[["subject", *SUBJECT_LEVEL_COLUMNS]]
        .drop_duplicates(subset=["subject"])
        .sort_values("subject")
    )
    subject_records = [
        (
            row.subject,
            project_lookup[row.project],
            row.condition,
            int(row.age),
            row.sex,
            row.treatment,
            _nullable(row.response),
        )
        for row in subject_df.itertuples(index=False)
    ]
    connection.executemany(
        """
        INSERT INTO subjects (
            subject_id, project_id, condition, age, sex, treatment, response
        )
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        subject_records,
    )

    sample_df = df[["sample", "subject", "sample_type", "time_from_treatment_start"]]
    sample_records = [
        (
            row.sample,
            row.subject,
            row.sample_type,
            int(row.time_from_treatment_start),
        )
        for row in sample_df.itertuples(index=False)
    ]
    connection.executemany(
        """
        INSERT INTO samples (
            sample_id, subject_id, sample_type, time_from_treatment_start
        )
        VALUES (?, ?, ?, ?)
        """,
        sample_records,
    )

    population_records = [(index + 1, name) for index, name in enumerate(CELL_POPULATIONS)]
    connection.executemany(
        "INSERT INTO cell_populations (population_id, name) VALUES (?, ?)",
        population_records,
    )
    population_lookup = dict(
        connection.execute("SELECT name, population_id FROM cell_populations").fetchall()
    )

    long_counts = df.melt(
        id_vars=["sample"],
        value_vars=CELL_POPULATIONS,
        var_name="population",
        value_name="count",
    )
    count_records = [
        (row.sample, population_lookup[row.population], int(row.count))
        for row in long_counts.itertuples(index=False)
    ]
    connection.executemany(
        """
        INSERT INTO cell_counts (sample_id, population_id, count)
        VALUES (?, ?, ?)
        """,
        count_records,
    )

    return {
        "projects": len(projects),
        "subjects": len(subject_records),
        "samples": len(sample_records),
        "cell_populations": len(population_records),
        "cell_counts": len(count_records),
    }


def _nullable(value: Any) -> Any:
    return None if pd.isna(value) else value


def validate_loaded_data(
    connection: sqlite3.Connection, df: pd.DataFrame, expected_counts: dict[str, int]
) -> None:
    """Validate database row counts and relational integrity after loading."""
    table_counts = {
        "projects": _query_count(connection, "projects"),
        "subjects": _query_count(connection, "subjects"),
        "samples": _query_count(connection, "samples"),
        "cell_populations": _query_count(connection, "cell_populations"),
        "cell_counts": _query_count(connection, "cell_counts"),
    }
    if table_counts != expected_counts:
        raise RuntimeError(
            f"Loaded row counts do not match inserted counts: {table_counts} != {expected_counts}"
        )

    sample_count_problems = connection.execute(
        """
        SELECT sample_id, COUNT(*) AS records
        FROM cell_counts
        GROUP BY sample_id
        HAVING records != ?
        LIMIT 5
        """,
        (len(CELL_POPULATIONS),),
    ).fetchall()
    if sample_count_problems:
        raise RuntimeError(
            "Every sample must have exactly five cell-count records. Examples: "
            + str(sample_count_problems)
        )

    csv_sample_count = df["sample"].nunique()
    if table_counts["samples"] != csv_sample_count:
        raise RuntimeError(
            f"Expected {csv_sample_count} CSV samples in SQLite, found {table_counts['samples']}."
        )

    missing_samples = connection.execute(
        """
        SELECT COUNT(*)
        FROM samples s
        LEFT JOIN cell_counts cc ON cc.sample_id = s.sample_id
        WHERE cc.sample_id IS NULL
        """
    ).fetchone()[0]
    if missing_samples:
        raise RuntimeError(f"{missing_samples} samples have no cell-count records.")

    foreign_key_violations = connection.execute("PRAGMA foreign_key_check;").fetchall()
    if foreign_key_violations:
        raise RuntimeError(f"Foreign key violations found: {foreign_key_violations}")


def _query_count(connection: sqlite3.Connection, table_name: str) -> int:
    return int(connection.execute(f"SELECT COUNT(*) FROM {table_name}").fetchone()[0])


def main() -> None:
    repo_root = get_repo_root()
    csv_path = repo_root / CSV_FILENAME
    db_path = get_database_path()

    if not csv_path.exists():
        raise FileNotFoundError(f"Could not find {CSV_FILENAME} at {csv_path}")

    print(f"Loading {CSV_FILENAME}...")
    df = read_dataset(csv_path)
    validate_dataframe(df)
    print(f"Validated {_format_int(len(df))} input rows.")

    counts = load_dataframe(df, db_path)
    print(f"Created {DATABASE_FILENAME}.")
    print("")
    print("Loaded:")
    print(f"  Projects:       {_format_int(counts['projects']):>7}")
    print(f"  Subjects:       {_format_int(counts['subjects']):>7}")
    print(f"  Samples:        {_format_int(counts['samples']):>7}")
    print(f"  Populations:    {_format_int(counts['cell_populations']):>7}")
    print(f"  Cell counts:    {_format_int(counts['cell_counts']):>7}")
    print("")
    _print_expected_count_note(counts)
    print("Database integrity checks passed.")


def _print_expected_count_note(counts: dict[str, int]) -> None:
    mismatches = {
        key: (counts[key], expected)
        for key, expected in EXPECTED_COUNTS.items()
        if counts.get(key) != expected
    }
    if not mismatches:
        print("Expected Phase 1 row counts produced.")
        return

    print("Loaded counts differ from the known dataset contract:")
    for key, (actual, expected) in mismatches.items():
        print(f"  {key}: actual {_format_int(actual)}, expected {_format_int(expected)}")


if __name__ == "__main__":
    main()
