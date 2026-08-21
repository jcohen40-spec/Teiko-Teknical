#!/usr/bin/env python3
"""
analysis/response_comparison.py

Part 3: Responder vs. Non-Responder Comparison

Question: Which immune cell populations differ between melanoma patients
on miraclib who respond to treatment vs. those who don't?

Filters applied (per spec):
    - condition == 'melanoma'
    - treatment == 'miraclib'
    - sample_type == 'PBMC'

For each of the 5 cell populations, compares the relative frequency (%)
between responders (response == 'yes') and non-responders (response == 'no')
using:
    - A boxplot (one panel per population) showing the two groups side by
      side, saved as a single figure.
    - An independent two-sample t-test (Welch's t-test, which does not
      assume the two groups have equal variance) comparing the two
      groups per population. With exactly two groups, this is
      equivalent to one-way ANOVA -- ANOVA's F-statistic reduces to the
      t-statistic squared (F = t^2) and gives the identical p-value, so
      the t-test is used directly since it's the more standard choice
      for a two-group comparison.
    - Bonferroni correction across the 5 simultaneous tests (one per
      population), since testing multiple populations at once inflates
      the chance of a false positive by chance alone. Bonferroni simply
      multiplies each raw p-value by the number of tests (5) -- the
      simplest and most conservative standard correction, and
      straightforward to explain to a non-statistician.

Usage
-----
    python analysis/response_comparison.py

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

matplotlib.use("Agg")  # headless-safe backend for script/dashboard use
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

try:
    from analysis.frequencies import compute_frequencies, POPULATION_COLUMNS
except ImportError:
    # Allows running this file directly as a script (python analysis/response_comparison.py)
    from frequencies import compute_frequencies, POPULATION_COLUMNS

REPO_ROOT = Path(__file__).resolve().parent.parent
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

    mean_diff is (mean_responder - mean_non_responder), in percentage
    points. cohens_d is a standardized effect size (mean_diff divided by
    the pooled standard deviation) -- useful alongside the p-value since
    with large sample sizes even a very small, practically negligible
    difference can be statistically significant. Rough interpretation:
    ~0.2 = small, ~0.5 = medium, ~0.8 = large effect.
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
