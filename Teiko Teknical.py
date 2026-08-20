"""

Part 2: Initial Analysis - Data Overview

"What is the frequency of each cell type in each sample?"

For every sample, computes the total cell count (sum across all five
populations) and the relative frequency (%) of each population within
that sample and also reports cell count by type in each sample.

Usage
~~~~~

Produces a long-format table with one row per (sample, population):

columns:    sample | total_count | population | count | percentage

and writes it to analysis/output/frequencies.csv, printing a preview
to stdout.

Can also be imported and used programmatically, e.g. by the dashboard:

    from analysis.frequencies import compute_frequencies
    df = compute_frequencies("cell_counts.db")
"""

import sqlite3
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent
DB_PATH = REPO_ROOT / "cell_counts.db"
OUTPUT_PATH = Path(__file__).resolve().parent / "output" / "frequencies.csv"

# The five cell population columns stored
# on cell_counts (wide format).
POPULATION_COLUMNS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]


def compute_frequencies(db_path: Path = DB_PATH) -> pd.DataFrame:
    """
    Compute per-sample, per-population relative frequencies.

    Returns a long-format DataFrame with columns:
        sample, total_count, population, count, percentage
    """
    conn = sqlite3.connect(db_path)
    try:
        # cell_counts is stored wide (one row per sample, one column per
        # population). Pull it into a DataFrame and reshape to long format
        # so we get one row per (sample, population), as required.
        wide = pd.read_sql_query(
            f"SELECT sample_id AS sample, {', '.join(POPULATION_COLUMNS)} "
            "FROM cell_counts",
            conn,
        )
    finally:
        conn.close()

    wide["total_count"] = wide[POPULATION_COLUMNS].sum(axis=1)

    long_df = wide.melt(
        id_vars=["sample", "total_count"],
        value_vars=POPULATION_COLUMNS,
        var_name="population",
        value_name="count",
    )

    long_df["percentage"] = (long_df["count"] / long_df["total_count"] * 100).round(2)

    # Order rows by sample, then keep population order consistent with the
    # source columns for readability.
    long_df["population"] = pd.Categorical(
        long_df["population"], categories=POPULATION_COLUMNS, ordered=True
    )
    long_df = long_df.sort_values(["sample", "population"]).reset_index(drop=True)

    return long_df[["sample", "total_count", "population", "count", "percentage"]]


# I did not print the entire table due to its size (52500 rows), but the option to print the entire table is available
# by removing .head(10) from print(df.head(10).to_string(index=False))

def main() -> None:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Could not find {DB_PATH}. Run `python load_data.py` from the "
            "repo root first."
        )

    df = compute_frequencies(DB_PATH)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_PATH, index=False)

    print(f"Computed frequencies for {df['sample'].nunique()} samples "
          f"({len(df)} sample-population rows).")
    print(f"Saved to: {OUTPUT_PATH}")
    print()
    print(df.head(10).to_string(index=False))

    num_rows = df.shape[0]
    print("Number of rows:", num_rows)

if __name__ == "__main__":
    main()

"""
Part 3: 

Filters applied (per spec):
    - condition == 'melanoma'
    - treatment == 'miraclib'
    - sample_type == 'PBMC'

For each of the 5 cell types, compares the relative frequency (%)
between responders (response == 'yes') and non-responders (response == 'no')

    - A boxplot (one panel per cell population) showing the two groups side by
      side, saved as a single figure.
    - An independent two-sample t-test (Welch's t-test, which does not
      assume the two groups have equal variance) comparing the two
      groups per cell population. 
    -Bonferronia correction used to account for multiple statistical tests 
     at once which can inflate odds of false positive. Multiply p-values by 5
     since that is the number of tests. Simple correction, easily explained.

Usage
-----

Outputs:
    analysis/output/response_comparison_boxplot.png
    analysis/output/response_comparison_stats.csv

Can also be imported and used programmatically, e.g. by the dashboard:

    from analysis.response_comparison import (
        get_filtered_frequencies, run_statistical_tests
    )
"""

import sqlite3
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parent
DB_PATH = REPO_ROOT / "cell_counts.db"
OUTPUT_DIR = Path(__file__).resolve().parent / "output"

CONDITION = "melanoma"
TREATMENT = "miraclib"
SAMPLE_TYPE = "PBMC"


def get_filtered_frequencies(db_path: Path = DB_PATH) -> pd.DataFrame:
    """
    Return the Part 2 long-format frequency table, restricted to
    melanoma + miraclib + PBMC samples, with each row's response
    ('yes'/'no') attached.

    Columns: sample, total_count, population, count, percentage, response
    """
    freq = compute_frequencies(db_path)

    conn = sqlite3.connect(db_path)
    try:
        meta = pd.read_sql_query(
            """
            SELECT s.sample_id AS sample, s.sample_type,
                   sub.condition, sub.treatment, sub.response
            FROM samples s
            JOIN subjects sub ON s.subject_id = sub.subject_id
            WHERE sub.condition = ?
              AND sub.treatment = ?
              AND s.sample_type = ?
            """,
            conn,
            params=(CONDITION, TREATMENT, SAMPLE_TYPE),
        )
    finally:
        conn.close()

    merged = freq.merge(meta[["sample", "response"]], on="sample", how="inner")
    return merged


def run_statistical_tests(filtered: pd.DataFrame) -> pd.DataFrame:
    """
    For each population, run an independent two-sample t-test (Welch's,
    unequal variance) comparing relative frequency (%) between responders
    and non-responders. Applies Bonferroni correction across the 5
    populations tested.

    Returns a DataFrame with columns:
          population, n_responder, n_non_responder,
          mean_responder_pct, mean_non_responder_pct, mean_diff, cohens_d,
          t_statistic, p_value, p_value_adj, significant
    """
    rows = []
    for pop in POPULATION_COLUMNS:
        sub = filtered[filtered["population"] == pop]
        resp = sub.loc[sub["response"] == "yes", "percentage"]
        non_resp = sub.loc[sub["response"] == "no", "percentage"]

        t_stat, p_val = stats.ttest_ind(resp, non_resp, equal_var=False)

        mean_diff = resp.mean() - non_resp.mean()
        pooled_std = ((resp.std() ** 2 + non_resp.std() ** 2) / 2) ** 0.5
        cohens_d = mean_diff / pooled_std if pooled_std > 0 else float("nan")

        rows.append(
            {
                "population": pop,
                "n_responder": len(resp),
                "n_non_responder": len(non_resp),
                "mean_responder_pct": round(resp.mean(), 2),
                "mean_non_responder_pct": round(non_resp.mean(), 2),
                "mean_diff": round(mean_diff, 3),
                "cohens_d": round(cohens_d, 3),
                "t_statistic": round(t_stat, 2),
                "p_value": p_val,
            }
        )

    results = pd.DataFrame(rows)

    # Bonferroni correction: multiply each raw p-value by the number of
    # tests (5), capped at 1.0. Simple and conservative -- easy to explain.
    n_tests = len(results)
    results["p_value_adj"] = (results["p_value"] * n_tests).clip(upper=1.0)
    results["significant"] = results["p_value_adj"] < 0.05

    return results.sort_values("p_value_adj").reset_index(drop=True)
    

def make_boxplot(filtered: pd.DataFrame, output_path: Path) -> None:
    """Boxplot of relative frequency (%) by response, one panel per population."""
    fig, axes = plt.subplots(1, len(POPULATION_COLUMNS), figsize=(4 * len(POPULATION_COLUMNS), 5), sharey=False)

    for ax, pop in zip(axes, POPULATION_COLUMNS):
        sub = filtered[filtered["population"] == pop]
        groups = [
            sub.loc[sub["response"] == "no", "percentage"],
            sub.loc[sub["response"] == "yes", "percentage"],
        ]
        ax.boxplot(groups, tick_labels=["Non-responder", "Responder"], showmeans=True)
        ax.set_title(pop)
        ax.set_ylabel("Relative frequency (%)")
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle(
        f"Cell population frequencies: responders vs. non-responders\n"
        f"({CONDITION}, {TREATMENT}, {SAMPLE_TYPE} samples)"
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def make_mean_errorbar_plot(filtered: pd.DataFrame, output_path: Path) -> None:
    """
    Plot mean relative frequency (%) +/- standard error of the mean (SEM)
    for responders vs. non-responders, one pair of points per population.

    Unlike the boxplot (which shows the full spread of individual samples),
    this zooms in on just the means and their uncertainty -- much easier
    to see small-but-real shifts like cd4_t_cell's, which get visually
    swamped by outliers/spread on a boxplot.
    """
    x = range(len(POPULATION_COLUMNS))
    offset = 0.08

    non_resp_means, non_resp_errs = [], []
    resp_means, resp_errs = [], []

    for pop in POPULATION_COLUMNS:
        sub = filtered[filtered["population"] == pop]
        resp = sub.loc[sub["response"] == "yes", "percentage"]
        non_resp = sub.loc[sub["response"] == "no", "percentage"]

        resp_means.append(resp.mean())
        resp_errs.append(resp.sem())  # SEM = std / sqrt(n)
        non_resp_means.append(non_resp.mean())
        non_resp_errs.append(non_resp.sem())

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.errorbar(
        [i - offset for i in x], non_resp_means, yerr=non_resp_errs,
        fmt="o", capsize=4, label="Non-responder", color="tab:orange",
    )
    ax.errorbar(
        [i + offset for i in x], resp_means, yerr=resp_errs,
        fmt="o", capsize=4, label="Responder", color="tab:blue",
    )

    ax.set_xticks(list(x))
    ax.set_xticklabels(POPULATION_COLUMNS)
    ax.set_ylabel("Mean relative frequency (%) \u00b1 SEM")
    ax.set_title(
        f"Mean cell population frequency \u00b1 SEM: responders vs. non-responders\n"
        f"({CONDITION}, {TREATMENT}, {SAMPLE_TYPE} samples)"
    )
    ax.legend()
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main() -> None:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Could not find {DB_PATH}. Run `python load_data.py` from the "
            "repo root first."
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    filtered = get_filtered_frequencies(DB_PATH)
    print(
        f"Filtered to {filtered['sample'].nunique()} samples "
        f"({CONDITION}, {TREATMENT}, {SAMPLE_TYPE})."
    )
    print(filtered.groupby("response")["sample"].nunique().rename("n_samples"))
    print()

    results = run_statistical_tests(filtered)
    stats_path = OUTPUT_DIR / "response_comparison_stats.csv"
    results.to_csv(stats_path, index=False)

    print("Independent t-test results (Bonferroni-corrected across 5 populations):")
    print(results.to_string(index=False))
    print(f"\nSaved statistics to: {stats_path}")

    boxplot_path = OUTPUT_DIR / "response_comparison_boxplot.png"
    make_boxplot(filtered, boxplot_path)
    print(f"Saved boxplot to: {boxplot_path}")

    errorbar_path = OUTPUT_DIR / "response_comparison_mean_errorbar.png"
    make_mean_errorbar_plot(filtered, errorbar_path)
    print(f"Saved mean +/- SEM plot to: {errorbar_path}")

    sig = results[results["significant"]]["population"].tolist()
    if sig:
        print(f"\nSignificant populations (Bonferroni-adjusted p < 0.05): {sig}")
    else:
        print("\nNo populations reached significance after Bonferroni correction.")


if __name__ == "__main__":
    main()


"""
Part 4:

Identifies all melanoma, PBMC samples at baseline
(time_from_treatment_start == 0) from patients treated with miraclib,
then breaks that cohort down by:

    - number of samples per project
    - number of subjects who are responders vs. non-responders
    - number of subjects who are male vs. female

Usage
-----

Outputs (printed to stdout, and saved to analysis/output/):
    baseline_samples_per_project.csv
    baseline_subjects_by_response.csv
    baseline_subjects_by_sex.csv

Can also be imported and used programmatically, e.g. by the dashboard:

    from analysis.baseline_summary import (
        get_baseline_cohort, samples_per_project,
        subjects_by_response, subjects_by_sex,
    )
"""

import sqlite3
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent
DB_PATH = REPO_ROOT / "cell_counts.db"
OUTPUT_DIR = Path(__file__).resolve().parent / "output"

CONDITION = "melanoma"
TREATMENT = "miraclib"
SAMPLE_TYPE = "PBMC"
BASELINE_TIME = 0


def get_baseline_cohort(db_path: Path = DB_PATH) -> pd.DataFrame:
    """
    Query the database for all melanoma, PBMC, baseline (time == 0)
    samples from subjects treated with miraclib.

    Returns one row per sample, with subject-level fields (project,
    response, sex) attached via join.
    """
    conn = sqlite3.connect(db_path)
    try:
        df = pd.read_sql_query(
            """
            SELECT
                s.sample_id,
                s.sample_type,
                s.time_from_treatment_start,
                sub.subject_id,
                sub.project_id,
                sub.condition,
                sub.treatment,
                sub.response,
                sub.sex
            FROM samples s
            JOIN subjects sub ON s.subject_id = sub.subject_id
            WHERE sub.condition = ? 
              AND sub.treatment = ? 
              AND s.sample_type = ? 
              AND s.time_from_treatment_start = ?
            """,
            conn,
            params=(CONDITION, TREATMENT, SAMPLE_TYPE, BASELINE_TIME),
        )
    finally:
        conn.close()
    return df


def samples_per_project(cohort: pd.DataFrame) -> pd.DataFrame:
    # Sample counts by project
    return (
        cohort.groupby("project_id")["sample_id"]
        .nunique()
        .rename("n_samples")
        .reset_index()
        .sort_values("project_id")
        .reset_index(drop=True)
    )


def subjects_by_response(cohort: pd.DataFrame) -> pd.DataFrame:
    # Subject counts by response (Y/N)
    subj = cohort.drop_duplicates("subject_id")
    return (
        subj.groupby("response")["subject_id"]
        .nunique()
        .rename("n_subjects")
        .reset_index()
        .sort_values("response")
        .reset_index(drop=True)
    )


def subjects_by_sex(cohort: pd.DataFrame) -> pd.DataFrame:
    # Subject counts by sex (M/F)
    subj = cohort.drop_duplicates("subject_id")
    return (
        subj.groupby("sex")["subject_id"]
        .nunique()
        .rename("n_subjects")
        .reset_index()
        .sort_values("sex")
        .reset_index(drop=True)
    )


def main() -> None:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Could not find {DB_PATH}. Run `python load_data.py` from the "
            "repo root first."
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    cohort = get_baseline_cohort(DB_PATH)
    n_samples = cohort["sample_id"].nunique()
    n_subjects = cohort["subject_id"].nunique()
    print(
        f"Baseline cohort ({CONDITION}, {TREATMENT}, {SAMPLE_TYPE}, "
        f"time_from_treatment_start={BASELINE_TIME}):"
    )
    print(f"  {n_samples} samples from {n_subjects} subjects\n")

    by_project = samples_per_project(cohort)
    by_project.to_csv(OUTPUT_DIR / "baseline_samples_per_project.csv", index=False)
    print("Samples per project:")
    print(by_project.to_string(index=False))
    print()

    by_response = subjects_by_response(cohort)
    by_response.to_csv(OUTPUT_DIR / "baseline_subjects_by_response.csv", index=False)
    print("Subjects by response:")
    print(by_response.to_string(index=False))
    print()

    by_sex = subjects_by_sex(cohort)
    by_sex.to_csv(OUTPUT_DIR / "baseline_subjects_by_sex.csv", index=False)
    print("Subjects by sex:")
    print(by_sex.to_string(index=False))

    print(f"\nSaved breakdowns to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()


# Additional Question

"""
"Considering melanoma males of all sample and treatment types, what is
the average number of B cells for responders at time=0?"

Filters applied:
    - condition == 'melanoma'
    - sex == 'M'
    - response == 'yes'
    - time_from_treatment_start == 0

Usage
-----

Can also be imported and used programmatically, e.g. by the dashboard:

    from analysis.melanoma_male_responder_bcell import average_b_cell_count
    avg = average_b_cell_count("cell_counts.db")
"""

REPO_ROOT = Path(__file__).resolve().parent
DB_PATH = REPO_ROOT / "cell_counts.db"

CONDITION = "melanoma"
SEX = "M"
RESPONSE = "yes"
BASELINE_TIME = 0


def get_matching_samples(db_path: Path = DB_PATH) -> pd.DataFrame:
    """
    Return every sample (any sample_type, any treatment) matching:
    melanoma, male, responder, baseline (time_from_treatment_start == 0).

    One row per sample, including its b_cell count.
    """
    conn = sqlite3.connect(db_path)
    try:
        df = pd.read_sql_query(
            """
            SELECT
                s.sample_id,
                s.sample_type,
                s.time_from_treatment_start,
                sub.subject_id,
                sub.condition,
                sub.treatment,
                sub.response,
                sub.sex,
                cc.b_cell
            FROM samples s
            JOIN subjects sub ON s.subject_id = sub.subject_id
            JOIN cell_counts cc ON cc.sample_id = s.sample_id
            WHERE sub.condition = ?
              AND sub.sex = ?
              AND sub.response = ?
              AND s.time_from_treatment_start = ?
            """,
            conn,
            params=(CONDITION, SEX, RESPONSE, BASELINE_TIME),
        )
    finally:
        conn.close()
    return df


def average_b_cell_count(db_path: Path = DB_PATH) -> float:
    #Average b cell count given specified parameters
    df = get_matching_samples(db_path)
    return round(df["b_cell"].mean(), 2)


def main() -> None:
    if not DB_PATH.exists():
        raise FileNotFoundError(
            f"Could not find {DB_PATH}. Run `python load_data.py` from the "
            "repo root first."
        )

    df = get_matching_samples(DB_PATH)
    avg = round(df["b_cell"].mean(), 2)

    print(
        f"Melanoma, male, responders, time_from_treatment_start=0 "
        f"(all sample types, all treatments):"
    )
    print(f"  {len(df)} matching samples")
    print(f"  Average B cell count: {avg:.2f}")


if __name__ == "__main__":
    main()