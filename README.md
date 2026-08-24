# Teiko Technical Assessment

I am analyzing immune-cell population data from clinical-trial samples using a reproducible SQLite-backed Python pipeline and an interactive Streamlit dashboard. The current pipeline validates and normalizes `cell-count.csv` into SQLite, calculates immune-cell relative frequencies, runs responder/non-responder statistical comparisons, creates Plotly boxplots, runs the Part 4 baseline subset analysis from the relational database, and presents the results in a local Streamlit dashboard.

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

This part is implemented. My target cohort is:

```text
condition = melanoma
treatment = miraclib
sample_type = PBMC
response = yes/no
```

I compare responder and non-responder relative frequencies for each immune-cell population. The dataset is longitudinal:

```text
3 samples per subject
days 0, 7, 14
```

Because the same subjects are measured repeatedly, the primary test uses subject-level mean relative frequencies. For each `subject x population`, I average the available PBMC/miraclib timepoints before comparing response groups, so each subject contributes at most one primary statistical observation per population.

The statistical method is:

```text
Mann-Whitney U test
two-sided alternative
Benjamini-Hochberg FDR correction
adjusted p < 0.05 significance threshold
```

I also run a baseline-only secondary analysis at `time_from_treatment_start = 0` because post-treatment measurements should not be interpreted as pre-treatment predictors when discussing potential response prediction.

Generated outputs:

```text
outputs/statistical_results.csv
outputs/baseline_statistical_results.csv
outputs/responder_boxplot.html
outputs/baseline_responder_boxplot.html
```

I include rank-biserial correlation as a simple Mann-Whitney effect size. Positive values mean responder percentages tend to be higher than non-responder percentages; negative values mean they tend to be lower.

#### Part 3 Results

The target cohort contains 656 subjects, including 331 responders and 325 non-responders, with 1,968 samples across timepoints 0, 7, and 14.

Primary subject-level longitudinal mean comparison:

| population | responder median % | non-responder median % | median difference | p-value | adjusted p-value | significant |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| b_cell | 9.6714 | 9.8446 | -0.1732 | 0.3458 | 0.4322 | no |
| cd8_t_cell | 24.8969 | 25.0097 | -0.1128 | 0.6221 | 0.6221 | no |
| cd4_t_cell | 30.2098 | 29.8225 | 0.3873 | 0.0124 | 0.0621 | no |
| nk_cell | 14.7397 | 14.9598 | -0.2201 | 0.1267 | 0.3169 | no |
| monocyte | 19.7945 | 20.2768 | -0.4823 | 0.2645 | 0.4322 | no |

Baseline-only comparison:

| population | responder median % | non-responder median % | median difference | p-value | adjusted p-value | significant |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| b_cell | 9.7850 | 9.7582 | 0.0269 | 0.5485 | 0.8853 | no |
| cd8_t_cell | 24.3957 | 24.6010 | -0.2053 | 0.5140 | 0.8853 | no |
| cd4_t_cell | 29.6338 | 29.5305 | 0.1033 | 0.7964 | 0.8853 | no |
| nk_cell | 14.9973 | 14.8905 | 0.1069 | 0.8853 | 0.8853 | no |
| monocyte | 19.6056 | 20.2860 | -0.6804 | 0.2114 | 0.8853 | no |

No immune-cell population is significant after Benjamini-Hochberg correction in either analysis. The unadjusted `cd4_t_cell` primary comparison is nominally small, but it does not meet the adjusted significance threshold, so I do not treat it as a statistically significant responder-associated population. These results describe group associations only; they do not prove predictive performance.

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
|   |-- visualization.py
|   `-- statistics.py
|-- outputs/
|   `-- .gitkeep
`-- tests/
    |-- __init__.py
    |-- test_loader.py
    |-- test_dashboard.py
    |-- test_frequencies.py
    |-- test_queries.py
    `-- test_statistics.py
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
4. run Part 3 responder/non-responder statistical analysis
5. generate Part 3 Plotly boxplots
6. run Part 4 database subset queries
7. generate the implemented required tables


### `make dashboard`

Starts the implemented Streamlit dashboard:

```bash
streamlit run dashboard.py
```

The dashboard expects the pipeline outputs to exist. If they are missing, it shows a clear prompt to run `make pipeline`.

## Outputs

Generated artifacts currently include:

```text
outputs/relative_frequencies.csv
outputs/statistical_results.csv
outputs/baseline_statistical_results.csv
outputs/responder_boxplot.html
outputs/baseline_responder_boxplot.html
outputs/baseline_samples.csv
outputs/project_counts.csv
outputs/response_counts.csv
outputs/sex_counts.csv
outputs/baseline_b_cell_average.csv
```

Static PNG plot export is not generated because Plotly HTML export satisfies the current requirement without adding extra image-export dependencies.

## Dashboard

The Streamlit dashboard has four sections:

```text
Overview
Cell Frequencies
Treatment Response
Baseline Cohort
```

The Overview section derives top-level KPIs and sample distributions from SQLite. Cell Frequencies loads `outputs/relative_frequencies.csv` and provides population/sample filtering, an overall distribution boxplot, a selected-sample composition chart, and an interactive table.

The Treatment Response section presents the Part 3 target cohort, primary subject-level comparison, baseline-only comparison, Mann-Whitney U results, Benjamini-Hochberg adjusted p-values, significance interpretation, and responder/non-responder boxplots. In the current dataset, no immune-cell populations were statistically significant after correction.

The Baseline Cohort section presents the Part 4 baseline filter, baseline sample KPIs, project sample counts, response subject counts, sex subject counts, a filterable baseline sample table, and the final B-cell answer. The final B-cell metric is explicitly shown as all sample types and all treatment types, not PBMC/miraclib-only.

Dashboard: Run locally or in GitHub Codespaces with `make dashboard`.

## Testing Strategy

The test suite currently covers database row-count validation, primary-key and uniqueness constraints, correct population-frequency calculations, percentages summing to approximately 100% for every sample, SQL cohort filtering, correct handling of null response values, Part 4 sample-vs-subject counting semantics, the final B-cell average filter, subject-level response aggregation, Mann-Whitney output, Benjamini-Hochberg correction, non-significant and strong-difference statistical cases, baseline-only filtering, and lightweight dashboard helper behavior.

## Design Principles

The implementation will be guided by reproducibility, separation of concerns, normalized relational modeling, SQL-first filtering where the assessment specifically asks for database queries, statistically defensible analysis, explicit generated outputs, and simple execution in GitHub Codespaces.
