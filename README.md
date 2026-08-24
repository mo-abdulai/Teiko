# Teiko Technical

I am analyzing immune-cell population data from clinical-trial samples using a reproducible SQLite-backed Python pipeline and an interactive Streamlit dashboard. The project is currently in the initial scaffold stage: the repository structure, execution targets, and implementation plan are in place, while the database loader, analysis pipeline, statistics, and full dashboard will be implemented next.

## Assessment Goals

### Part 1 - Data Management

I will model the dataset in a relational SQLite database using a normalized schema. The root-level `load_data.py` script will initialize the database, create the schema, and load `cell-count.csv` in a reproducible way.

### Part 2 - Cell Population Frequencies

For every sample, I will calculate:

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

The dataset has 10,500 samples and 5 immune-cell populations, so this step will produce 52,500 long-format frequency rows.

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

My target cohort is:

```text
condition = melanoma
sample_type = PBMC
treatment = miraclib
time_from_treatment_start = 0
```

I will report sample count by project, responder/non-responder subject counts, and male/female subject counts.

The separate final calculation is:

```text
condition = melanoma
sex = male
response = yes
time_from_treatment_start = 0
```

The wording requests all sample types and all treatment types, so I will not automatically apply the PBMC or miraclib filters to this final calculation.

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

## Planned Database Schema

```text
projects
--------
project_id PK
project_name UNIQUE

subjects
--------
subject_id PK
project_id FK
condition
age
sex
treatment
response NULLABLE

samples
-------
sample_id PK
subject_id FK
sample_name UNIQUE
sample_type
time_from_treatment_start

cell_populations
----------------
population_id PK
name UNIQUE

cell_counts
-----------
sample_id FK
population_id FK
count
PRIMARY KEY (sample_id, population_id)
```

I am using this design because it reduces unnecessary duplication of subject-level metadata across repeated samples, models longitudinal measurements explicitly, represents immune-cell populations as data rows rather than permanent schema columns, and allows a future cell population to be added without changing the table structure. Normalized tables also support useful indexes and efficient cohort queries, and the model can scale to hundreds of projects, many thousands or millions of samples, additional cell populations, and future analytics.

I may adjust the schema slightly during implementation if dataset validation shows that a field I assumed to be subject-level actually varies across samples.

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

Will eventually execute my complete reproducible pipeline:

1. initialize/rebuild the SQLite database
2. load `cell-count.csv`
3. generate the Part 2 relative-frequency table
4. run Part 3 responder/non-responder statistical analysis
5. run Part 4 database subset queries
6. generate required tables and plots

This repository is currently in the initial scaffold stage, so the full pipeline is not implemented yet.

### `make dashboard`

Starts my Streamlit dashboard. The current dashboard is a minimal placeholder while implementation is in progress.

## Planned Outputs

Expected generated artifacts include:

```text
outputs/relative_frequencies.csv
outputs/statistical_results.csv
outputs/baseline_samples.csv
outputs/project_counts.csv
outputs/response_counts.csv
outputs/sex_counts.csv
outputs/responder_boxplot.html
```

Exact filenames may be refined as I implement the pipeline.

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

I plan to include tests for database row-count validation, primary-key and uniqueness constraints, correct population-frequency calculations, percentages summing to approximately 100% for every sample, correct SQL cohort filtering, correct handling of null response values, and correct Part 4 aggregations.

## Design Principles

The implementation will be guided by reproducibility, separation of concerns, normalized relational modeling, SQL-first filtering where the assessment specifically asks for database queries, statistically defensible analysis, explicit generated outputs, and simple execution in GitHub Codespaces.
