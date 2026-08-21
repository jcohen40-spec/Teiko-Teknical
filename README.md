# Jonathan Cohen - Teiko Teknical

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

Running this line in the Teiko Teknical workspace should also start the dashboard if it isn't up and running.

```bash
python -m streamlit run dashboard/app.py
```

**Dashboard link:** https://psychic-broccoli-jjxjv6r9464vc5gv4-8501.app.github.dev/

Note - the above link is only functional if the dashboard is active.

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

-I grouped in the way I found most intuitive. Project ID is the umbrella variable beneath which all other variables
would fall. Based on my rationale, I found that project ID should operate as a primary key. Then, the subjects variable
would be nested within with a foreign key mapping back to project ID and subject ID as a primary key. I also grouped
variables that fell beneath subject here, which were condition, age, sex, treatment, and response. Next, I nested sample
beneath subject, as each subject produced several samples, but the specifics of each sample are nested within subject. Thus, sample ID
operates as the primary key, and subject ID is the foreign key, structurally. I also placed the sample type and time from treatment
variables here, as they were sample-specific. Lastly, I nested cell counts within sample, using sample ID as a foreign key. In terms 
of my preferences, I thought it best not to create a primary key for cell counts, as it didn't seem like a useful variable in a structural sense.
Of course, for later parts of this analysis, it was helpful to pivot the cell counts to be in a long format as opposed to a wide one because
it was easier to calculate proportions accordingly. However, in terms of the SQLite database, it felt more natural to have separate variables
for each cell type, so one could easily, by sample, have all cell counts present as variables attached. This was my general rationale, 
and it mostly was produced based on how I understood each variable fell beneath, above, or amongst others in terms of hierarchy.

### How this would scale

At the current size (3 projects, 3,500 subjects, 10,500 samples), a single
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
  fixed, known set of 5 populations. If we started measuring
  new populations over time, or different panels per project, a wide
  schema would require schema migrations (`ALTER TABLE ADD COLUMN`)
  every time a new population is introduced. A long-format
  `cell_counts(sample_id FK, population, count)` table avoids that entirely — new
  populations are just new rows, and it also makes it trivial to write
  population-agnostic aggregate queries (`GROUP BY population`) that
  don't need to be rewritten as the panel grows. Moving to a long format
  might be less structurally intuitive, but it would be more analytically efficient
  with more cell types or a larger dataset. 
- **Partitioning/indexing strategy.** If we find that we are frequently filtering
  by the same grouping methods, we might consider adding unique indexes that reflect
  rows that are often grouped.
- **Separating raw data from derived/analytic tables.** Right now,
  Parts 2-4 recompute their summaries from `cell_counts` on every run.
  At scale, saving these values to avoid duplicate recalculation would be
  beneficial and save on computational energy and time.
- **A metadata/analysis-registry table.** It might be beneficial, at larger scale,
  to develop a storage method to understand how the data has been grouped, queried, analyzed, etc
  such that it is easier to go back and review previous processes that were run.
## Code structure

- **`load_data.py`** (Part 1) is deliberately dependency-free (standard
  library only: `csv`, `sqlite3`, `pathlib`) and self-contained, since
  it's the very first thing that has to run — no environment setup
  required before the database can be built.
- **`analysis/`** holds one script per analytical question (Parts 2-4, each following the same pattern:
  - a small, pure, importable function that does the actual query/
    computation and returns a DataFrame (e.g. `compute_frequencies()`,
    `run_statistical_tests()`, `get_baseline_cohort()`)
  - a `main()` that calls those functions, prints a 
    summary, and saves results to `analysis/output/`
  - an `if __name__ == "__main__":` guard so each script can be run
    standalone (`python analysis/frequencies.py`) **or** imported
    without side effects.

  This split matters because the **dashboard imports these same
  functions directly** rather than re-implementing the logic or
  shelling out to the scripts. Thus, the numbers shown in the dashboard
  and the numbers in the CSV/plot outputs can never drift out of sync
  with each other.
- **`dashboard/app.py`** doesn't contain any analytical code, but rather
  is used to produce an interactive visual of the analysis.
- **Statistical choices (Part 3)**: an independent two-sample t-test
  (Welch's, unequal variance) is used to compare responders vs.
  non-responders per population, as we are looking for a difference in proportions across two response groups.
- I applied a Bonferroni correction to account for five statistical tests being conducted back to back
- I also reported effect size (Cohen's d) and the raw mean difference
  with the p-value as a sanity check. It seemed, in the boxplots produced, that a major difference
  was not present among the response groups for the cd4_t cells, despite a significant p-value.
  I wanted to ensure that the significance was not a false positive, so I found these two values
  helpful in identifying that significance as true. I produced a mean with error bars plot as well to
  help look more closely at the differences across the response groups for each cell type. This
  was more helpful for confirming the significance, as there was no overlap in the error bars for the cd4_t
  cells while overlap was visible for all other cell types for the two response groups. Reporting
  effect size helped understand whether the significance was substantive or not.

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
