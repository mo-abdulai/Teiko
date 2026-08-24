"""SQLite helpers for the Teiko technical assessment."""

from __future__ import annotations

import sqlite3
from pathlib import Path


DATABASE_FILENAME = "cell_counts.db"
REPO_ROOT = Path(__file__).resolve().parents[1]
DATABASE_PATH = REPO_ROOT / DATABASE_FILENAME


def get_repo_root() -> Path:
    """Return the repository root inferred from this module location."""
    return REPO_ROOT


def get_database_path() -> Path:
    """Return the canonical project database path."""
    return DATABASE_PATH


def connect_database(path: Path | str | None = None) -> sqlite3.Connection:
    """Connect to SQLite with foreign key enforcement enabled."""
    db_path = Path(path) if path is not None else DATABASE_PATH
    connection = sqlite3.connect(db_path)
    connection.execute("PRAGMA foreign_keys = ON;")
    return connection


def recreate_database(path: Path | str | None = None) -> Path:
    """Remove and recreate the assessment database.

    The default behavior is intentionally idempotent for Phase 1: only the
    known project database path is removed before a fresh load.
    """
    db_path = Path(path) if path is not None else DATABASE_PATH
    if db_path.exists():
        db_path.unlink()
    return db_path


def initialize_schema(connection: sqlite3.Connection) -> None:
    """Create the normalized cell-count schema and indexes."""
    connection.executescript(
        """
        CREATE TABLE IF NOT EXISTS projects (
            project_id INTEGER PRIMARY KEY,
            project_name TEXT NOT NULL UNIQUE
        );

        CREATE TABLE IF NOT EXISTS subjects (
            subject_id TEXT PRIMARY KEY,
            project_id INTEGER NOT NULL,
            condition TEXT NOT NULL,
            age INTEGER NOT NULL,
            sex TEXT NOT NULL,
            treatment TEXT NOT NULL,
            response TEXT NULL,
            FOREIGN KEY(project_id) REFERENCES projects(project_id)
        );

        CREATE TABLE IF NOT EXISTS samples (
            sample_id TEXT PRIMARY KEY,
            subject_id TEXT NOT NULL,
            sample_type TEXT NOT NULL,
            time_from_treatment_start INTEGER NOT NULL
                CHECK(time_from_treatment_start >= 0),
            FOREIGN KEY(subject_id) REFERENCES subjects(subject_id)
        );

        CREATE TABLE IF NOT EXISTS cell_populations (
            population_id INTEGER PRIMARY KEY,
            name TEXT NOT NULL UNIQUE
        );

        CREATE TABLE IF NOT EXISTS cell_counts (
            sample_id TEXT NOT NULL,
            population_id INTEGER NOT NULL,
            count INTEGER NOT NULL CHECK(count >= 0),
            PRIMARY KEY(sample_id, population_id),
            FOREIGN KEY(sample_id) REFERENCES samples(sample_id),
            FOREIGN KEY(population_id) REFERENCES cell_populations(population_id)
        );

        CREATE INDEX IF NOT EXISTS idx_subjects_project_id
            ON subjects(project_id);
        CREATE INDEX IF NOT EXISTS idx_subjects_condition
            ON subjects(condition);
        CREATE INDEX IF NOT EXISTS idx_subjects_treatment
            ON subjects(treatment);
        CREATE INDEX IF NOT EXISTS idx_subjects_response
            ON subjects(response);
        CREATE INDEX IF NOT EXISTS idx_subjects_sex
            ON subjects(sex);
        CREATE INDEX IF NOT EXISTS idx_subjects_cohort
            ON subjects(condition, treatment, response);

        CREATE INDEX IF NOT EXISTS idx_samples_subject_id
            ON samples(subject_id);
        CREATE INDEX IF NOT EXISTS idx_samples_sample_type
            ON samples(sample_type);
        CREATE INDEX IF NOT EXISTS idx_samples_time
            ON samples(time_from_treatment_start);
        CREATE INDEX IF NOT EXISTS idx_samples_type_time
            ON samples(sample_type, time_from_treatment_start);

        CREATE INDEX IF NOT EXISTS idx_cell_counts_population_id
            ON cell_counts(population_id);
        """
    )
