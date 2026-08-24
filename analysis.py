"""Run analysis phases for the Teiko technical assessment.

Phase 2 implements Part 2 only: cell-population relative frequencies.
"""

from __future__ import annotations

from src.analysis import calculate_relative_frequencies, validate_relative_frequencies
from src.database import DATABASE_FILENAME, connect_database, get_database_path, get_repo_root
from src.queries import get_cell_counts_by_sample, get_population_names, get_sample_ids


OUTPUT_FILENAME = "relative_frequencies.csv"


def _format_int(value: int) -> str:
    return f"{value:,}"


def main() -> None:
    """Run the current analysis workflow."""
    repo_root = get_repo_root()
    db_path = get_database_path()
    output_dir = repo_root / "outputs"
    output_path = output_dir / OUTPUT_FILENAME

    if not db_path.exists():
        raise FileNotFoundError(
            f"Could not find {DATABASE_FILENAME} at {db_path}. "
            "Run `python load_data.py` before running analysis."
        )

    print("Running Part 2: cell population relative frequencies...")
    print("")

    connection = connect_database(db_path)
    try:
        sample_ids = get_sample_ids(connection)
        population_names = get_population_names(connection)
        cell_counts = get_cell_counts_by_sample(connection)
    finally:
        connection.close()

    frequencies = calculate_relative_frequencies(cell_counts)
    validate_relative_frequencies(
        frequencies,
        expected_samples=sample_ids,
        expected_populations=population_names,
    )

    output_dir.mkdir(exist_ok=True)
    frequencies.to_csv(output_path, index=False)

    print(f"Samples analyzed:      {_format_int(len(sample_ids)):>7}")
    print(f"Populations:          {_format_int(len(population_names)):>7}")
    print(f"Output rows:          {_format_int(len(frequencies)):>7}")
    print("")
    print("Validated percentage sums for all samples.")
    print(f"Wrote {output_path.relative_to(repo_root)}")


if __name__ == "__main__":
    main()
