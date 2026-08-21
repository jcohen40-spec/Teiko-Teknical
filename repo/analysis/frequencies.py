#!/usr/bin/env python3
"""
analysis/frequencies.py

Part 2: Initial Analysis - Data Overview

Answers: "What is the frequency of each cell type in each sample?"

For every sample, computes the total cell count (sum across all five
populations) and the relative frequency (%) of each population within
that sample.

Usage
-----
    python analysis/frequencies.py

Produces a long-format table with one row per (sample, population):

    sample | total_count | population | count | percentage

and writes it to analysis/output/frequencies.csv, printing a preview
to stdout.

Can also be imported and used programmatically, e.g. by the dashboard:

    from analysis.frequencies import compute_frequencies
    df = compute_frequencies("cell_counts.db")
"""

import sqlite3
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
DB_PATH = REPO_ROOT / "cell_counts.db"
OUTPUT_PATH = Path(__file__).resolve().parent / "output" / "frequencies.csv"

# The five cell population columns stored on cell_counts (wide format).
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


if __name__ == "__main__":
    main()
