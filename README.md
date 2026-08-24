# Teiko Technical Assessment

I am analyzing immune-cell population data from clinical-trial samples using a reproducible SQLite-backed Python pipeline and an interactive Streamlit dashboard. The current pipeline validates and normalizes `cell-count.csv` into SQLite, calculates immune-cell relative frequencies, and runs the Part 4 baseline subset analysis from the relational database. The responder statistical analysis and full dashboard remain planned for later phases.

## Assessment Goals

### Part 1 - Data Management

The dataset is modeled in a relational SQLite database using a normalized schema. The root-level `load_data.py` script initializes the database, creates the schema, and loads `cell-count.csv` in a reproducible way.

Run Phase 1 with:

```bash
python load_data.py
```

This creates `cell_counts.db` in the repository root. The loader is idempotent for this assessment phase: if the known project database already exists, it is removed and rebuilt from the source CSV.

### Part 2 - Cell Population Frequencies

For every sample, I calculate:

```text
total_count = sum of the five immune-cell counts
percentage = population_count / total_count * 100
```

The expected output columns are:

```text
sample
total_count
population
count
percentage
```

The dataset has 10,500 samples and 5 immune-cell populations, so this step produces 52,500 long-format frequency rows in `outputs/relative_frequencies.csv`.

### Part 3 - Treatment Response Analysis

My target cohort is:

```text
condition = melanoma
treatment = miraclib
sample_type = PBMC
response in {yes, no}
```

I will compare responder and non-responder relative frequencies for each immune-cell population. The dataset is longitudinal:

```text
3 samples per subject
days 0, 7, 14
```

Because the same subjects are measured repeatedly, I will avoid automatically treating every sample as an independent observation.

My planned primary statistical approach is to create subject-level immune-cell frequency summaries where appropriate, compare responder and non-responder distributions with the Mann-Whitney U test, correct the five population tests using Benjamini-Hochberg multiple-testing correction, use adjusted `p < 0.05` as the statistical significance threshold, and visualize the distributions with Plotly boxplots.

I also plan to include a baseline-only view because measurements taken after treatment should not be interpreted as pre-treatment predictors when discussing potential response prediction.

### Part 4 - Baseline Subset Analysis

This part is implemented. My baseline cohort definition is:

```text
condition = melanoma
treatment = miraclib
sample_type = PBMC
time_from_treatment_start = 0
```

For this cohort, I report:

```text
samples by project -> sample count
response breakdown -> distinct subject count
sex breakdown -> distinct subject count
```

The generated outputs are:

```text
outputs/baseline_samples.csv
outputs/project_counts.csv
outputs/response_counts.csv
outputs/sex_counts.csv
```

The separate final calculation is:

```text
condition = melanoma
sex = M
response = yes
time_from_treatment_start = 0
population = b_cell
```

The wording requests all sample types and all treatment types, so I do not apply the PBMC or miraclib filters to this final calculation. The result is written to:

```text
outputs/baseline_b_cell_average.csv
```

## Dataset

```text
Rows: 10,500
Subjects: 3,500
Samples per subject: 3
Timepoints: 0, 7, 14
Projects: 3
Immune populations: 5
```

The exact CSV columns are:

```text
project
subject
condition
age
sex
treatment
response
sample
sample_type
time_from_treatment_start
b_cell
cd8_t_cell
cd4_t_cell
nk_cell
monocyte
```

Known dataset characteristics:

- no exact duplicate rows
- no duplicate sample IDs
- every subject has repeated longitudinal samples
- each subject currently has measurements at days 0, 7, and 14
- `response` is the only known column containing missing values
- 1,422 response values are missing
- I will preserve missing responses as null rather than imputing them as `"no"`

A missing response does not necessarily mean a patient was a non-responder, so preserving null values avoids introducing an unsupported assumption.

## Implemented Database Schema

```text
projects
--------
project_id INTEGER PRIMARY KEY
project_name TEXT NOT NULL UNIQUE

subjects
--------
subject_id TEXT PRIMARY KEY
project_id INTEGER NOT NULL FK
condition TEXT NOT NULL
age INTEGER NOT NULL
sex TEXT NOT NULL
treatment TEXT NOT NULL
response TEXT NULLABLE

samples
-------
sample_id TEXT PRIMARY KEY
subject_id TEXT NOT NULL FK
sample_type TEXT NOT NULL
time_from_treatment_start INTEGER NOT NULL

cell_populations
----------------
population_id INTEGER PRIMARY KEY
name TEXT NOT NULL UNIQUE

cell_counts
-----------
sample_id TEXT NOT NULL FK
population_id INTEGER NOT NULL FK
count INTEGER NOT NULL CHECK(count >= 0)
PRIMARY KEY (sample_id, population_id)
```

The five CSV cell-count columns are normalized into long-form `cell_counts` rows linked to `cell_populations` definitions:

```text
b_cell
cd8_t_cell
cd4_t_cell
nk_cell
monocyte
```

Missing `response` values are preserved as SQL `NULL`; they are not converted to `"no"` or a string sentinel.

I am using this design because it reduces unnecessary duplication of subject-level metadata across repeated samples, models longitudinal measurements explicitly, represents immune-cell populations as data rows rather than permanent schema columns, and allows a future cell population to be added without changing the table structure. Normalized tables also support useful indexes and efficient cohort queries, and the model can scale to hundreds of projects, many thousands or millions of samples, additional cell populations, and future analytics.

The loader validates required columns, unique sample IDs, unique `(subject, time_from_treatment_start, sample_type)` combinations, integer-compatible timepoints, numeric non-negative cell counts, non-missing subject/sample IDs, and null-aware consistency for subject-level metadata (`project`, `condition`, `age`, `sex`, `treatment`, `response`).

## Project Structure

```text
teiko-technical/
|-- README.md
|-- Makefile
|-- requirements.txt
|-- .gitignore
|-- load_data.py
|-- analysis.py
|-- dashboard.py
|-- src/
|   |-- __init__.py
|   |-- database.py
|   |-- queries.py
|   |-- analysis.py
|   `-- statistics.py
|-- outputs/
|   `-- .gitkeep
`-- tests/
    |-- __init__.py
    |-- test_loader.py
    |-- test_frequencies.py
    `-- test_queries.py
```

Root-level orchestration scripts provide the required command-line entry points for the assessment. The `src/` package contains reusable application modules. The `outputs/` directory will hold generated analysis artifacts. The `tests/` directory will contain automated validation.

`load_data.py` must remain at the repository root because this is required by the assessment grader.

## Reproducibility

My intended execution workflow is:

```bash
make setup
make pipeline
make dashboard
```

### `make setup`

Installs all project dependencies.

### `make pipeline`

Currently executes the implemented reproducible pipeline:

```bash
python load_data.py
python analysis.py
```

The pipeline currently:

1. initialize/rebuild the SQLite database
2. load `cell-count.csv`
3. generate the Part 2 relative-frequency table
4. run Part 4 database subset queries
5. generate the implemented required tables

Part 3 responder/non-responder statistical analysis and responder boxplots are not implemented yet.


### `make dashboard`

Starts my Streamlit dashboard. The current dashboard is a minimal placeholder while implementation is in progress.

## Outputs

Generated artifacts currently include:

```text
outputs/relative_frequencies.csv
outputs/baseline_samples.csv
outputs/project_counts.csv
outputs/response_counts.csv
outputs/sex_counts.csv
outputs/baseline_b_cell_average.csv
```

Planned later artifacts for Part 3 may include statistical result tables and Plotly responder visualizations.

## Dashboard

My intended Streamlit dashboard sections are:

```text
Overview
Cell Frequencies
Treatment Response
Baseline Cohort
```

I plan to use interactive filters and Plotly visualizations.

Dashboard: _To be added after deployment._

## Testing Strategy

The test suite currently covers database row-count validation, primary-key and uniqueness constraints, correct population-frequency calculations, percentages summing to approximately 100% for every sample, SQL cohort filtering, correct handling of null response values, Part 4 sample-vs-subject counting semantics, and the final B-cell average filter.

## Design Principles

The implementation will be guided by reproducibility, separation of concerns, normalized relational modeling, SQL-first filtering where the assessment specifically asks for database queries, statistically defensible analysis, explicit generated outputs, and simple execution in GitHub Codespaces.
