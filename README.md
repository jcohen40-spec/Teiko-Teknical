# Loblaw Bio — Immune Cell Population Analysis

Analysis pipeline and interactive dashboard for exploring how immune cell
populations relate to treatment response in Bob Loblaw's clinical trial
data (`cell-count.csv`).

## Quickstart (GitHub Codespaces or local)

```bash
make setup      # install dependencies
make pipeline   # build the database and run all analyses (Parts 1-4)
make dashboard  # launch the interactive dashboard
```

`make dashboard` starts a Streamlit server on `localhost:8501`. In
GitHub Codespaces, a popup/notification will offer to forward that port
and open it in your browser automatically; if not, open the **Ports**
tab and click the forwarded `8501` link.

**Dashboard link:** _(add your deployed/forwarded URL here once running)_

## Project layout

```
.
├── load_data.py                 # Part 1: builds cell_counts.db from cell-count.csv
├── cell-count.csv                # source data
├── cell_counts.db                 # generated SQLite database (committed for convenience;
│                                   #   regenerated deterministically by load_data.py)
├── analysis/
│   ├── frequencies.py            # Part 2: relative frequency per sample/population
│   ├── response_comparison.py    # Part 3: responder vs non-responder stats + plots
│   ├── baseline_summary.py       # Part 4: baseline cohort breakdowns
│   ├── melanoma_male_responder_bcell.py   # ad-hoc question script
│   └── output/                   # generated CSVs and plots (committed, regenerable)
├── dashboard/
│   └── app.py                    # Streamlit dashboard (Parts 2-4, interactive)
├── requirements.txt
├── Makefile
└── README.md
```

## Database schema

`load_data.py` builds a normalized SQLite database with four tables:

```
projects(project_id PK)
  └─< subjects(subject_id PK, project_id FK, condition, age, sex, treatment, response)
        └─< samples(sample_id PK, subject_id FK, sample_type, time_from_treatment_start)
              └─1 cell_counts(sample_id PK/FK, b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte)
```

**Design rationale:**

- **One row per real-world entity, not per CSV row.** The source CSV repeats
  each subject's demographic/clinical fields (`condition`, `age`, `sex`,
  `treatment`, `response`) on every sample row. Those fields don't change
  per sample — they describe the *subject* — so they're stored once per
  subject rather than duplicated. This removes redundancy and prevents the
  data from ever going inconsistent (e.g. one sample row saying `age=57`
  and another for the same subject saying `age=58` by data-entry error).
- **`cell_counts` is one-to-one with `samples`**, keyed on `sample_id`
  directly (used as both primary key and foreign key), because every
  sample has exactly one set of cell counts. Cell population counts are
  stored **wide** (one column per population — `b_cell`, `cd8_t_cell`,
  etc.) rather than as separate long-format rows, mirroring the shape of
  the source data directly and keeping the common case (retrieve all 5
  counts for a sample) a single-row lookup.
- **Foreign keys enforce referential integrity**: a sample can't
  reference a subject that doesn't exist, and cell counts can't reference
  a sample that doesn't exist. This is checked automatically by SQLite
  (`PRAGMA foreign_keys = ON`) rather than trusted to application code.
- **Indexes** on the foreign key columns (`subjects.project_id`,
  `samples.subject_id`) keep the joins used throughout `analysis/` fast
  even as the dataset grows.

### How this would scale

At the current size (3 projects, 3,500 subjects, 10,500 samples) a single
SQLite file with in-Python joins is fast and simple. If this grew to
**hundreds of projects, thousands of samples, and a wider variety of
analyses**, a few things would change:

- **Database engine**: SQLite is single-writer and file-based, which is
  fine for a local analysis tool but becomes a bottleneck with concurrent
  writers (e.g. multiple pipelines loading new trial data at once) or
  when the dashboard needs to be served to many simultaneous users.
  I'd move to a client-server database (Postgres) at that point — the
  schema itself wouldn't need to change, since it's already in standard
  normalized relational form.
- **Cell population columns → a long-format table.** The wide
  `cell_counts` schema (one column per population) works well for a
  fixed, known set of 5 populations. If Loblaw Bio started measuring
  new populations over time, or different panels per project, a wide
  schema would require schema migrations (`ALTER TABLE ADD COLUMN`)
  every time a new population is introduced. A long-format
  `cell_counts(sample_id FK, population, count)` table (which I used in
  an earlier draft of this schema) avoids that entirely — new
  populations are just new rows, and it also makes it trivial to write
  population-agnostic aggregate queries (`GROUP BY population`) that
  don't need to be rewritten as the panel grows.
- **Partitioning / indexing strategy.** With thousands of samples per
  project, queries that filter heavily by `project_id`, `condition`, or
  `time_from_treatment_start` (as Part 3 and Part 4 already do) would
  benefit from composite indexes on those common filter combinations,
  and potentially partitioning large tables by project or by trial
  phase.
- **Separating raw data from derived/analytic tables.** Right now,
  Parts 2-4 recompute their summaries from `cell_counts` on every run.
  At scale, with expensive or frequently-reused analyses, I'd
  materialize commonly-needed aggregates (e.g. per-sample relative
  frequencies) into their own table or a scheduled batch job, so the
  dashboard doesn't recompute the same joins/aggregations from scratch
  on every page load.
- **A metadata/analysis-registry table.** With "various types of
  analytics," it becomes worth tracking *what* analyses have been run,
  on what data version, with what parameters — essentially light
  experiment tracking — so results are reproducible and auditable as
  the number of analyses grows beyond what any one person can hold in
  their head.

## Code structure

- **`load_data.py`** (Part 1) is deliberately dependency-free (standard
  library only: `csv`, `sqlite3`, `pathlib`) and self-contained, since
  it's the very first thing that has to run — no environment setup
  required before the database can be built.
- **`analysis/`** holds one script per analytical question (Parts 2-4,
  plus the ad-hoc question script), each following the same pattern:
  - a small, pure, importable function that does the actual query/
    computation and returns a DataFrame (e.g. `compute_frequencies()`,
    `run_statistical_tests()`, `get_baseline_cohort()`)
  - a `main()` that calls those functions, prints a human-readable
    summary, and saves results to `analysis/output/`
  - an `if __name__ == "__main__":` guard so each script can be run
    standalone (`python analysis/frequencies.py`) **or** imported
    without side effects.

  This split matters because the **dashboard imports these same
  functions directly** rather than re-implementing the logic or
  shelling out to the scripts — so the numbers shown in the dashboard
  and the numbers in the CSV/plot outputs can never drift out of sync
  with each other; there is exactly one implementation of each
  calculation.
- **`dashboard/app.py`** is a thin presentation layer on top of
  `analysis/` — it contains no analytical logic of its own, only
  Streamlit UI code (tables, charts, filters) wired to the functions
  described above.
- **Statistical choices (Part 3)**: an independent two-sample t-test
  (Welch's, unequal variance) is used to compare responders vs.
  non-responders per population — with exactly two groups this is
  equivalent to one-way ANOVA, so the simpler/more familiar test is used
  directly. Bonferroni correction is applied across the 5 simultaneous
  per-population tests to control the family-wise false positive rate.
  Effect size (Cohen's d) and the raw mean difference are reported
  alongside the p-value, since with ~1,000 samples per group even a
  small, practically negligible difference can be statistically
  significant — reporting effect size avoids over-interpreting a small
  p-value as a large biological effect.

## Reproducing the outputs

`make pipeline` runs, in order:

1. `python load_data.py` — builds `cell_counts.db` from `cell-count.csv`
2. `python analysis/frequencies.py` — writes `analysis/output/frequencies.csv`
3. `python analysis/response_comparison.py` — writes
   `analysis/output/response_comparison_stats.csv`,
   `response_comparison_boxplot.png`, and
   `response_comparison_mean_errorbar.png`
4. `python analysis/baseline_summary.py` — writes
   `analysis/output/baseline_samples_per_project.csv`,
   `baseline_subjects_by_response.csv`, and `baseline_subjects_by_sex.csv`

All four steps are idempotent and can be re-run any time the source CSV
changes.
