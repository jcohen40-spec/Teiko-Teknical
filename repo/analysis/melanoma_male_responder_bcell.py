#!/usr/bin/env python3
"""
analysis/melanoma_male_responder_bcell.py

Ad-hoc question:
"Considering melanoma males of all sample and treatment types, what is
the average number of B cells for responders at time=0?"

Filters applied:
    - condition == 'melanoma'
    - sex == 'M'
    - response == 'yes'
    - time_from_treatment_start == 0
    (sample_type and treatment are NOT filtered -- "all types" per the
    question)

Usage
-----
    python analysis/melanoma_male_responder_bcell.py

Can also be imported and used programmatically, e.g. by the dashboard:

    from analysis.melanoma_male_responder_bcell import average_b_cell_count
    avg = average_b_cell_count("cell_counts.db")
"""

import sqlite3
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
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
    """Average b_cell count across the matching samples, rounded to 2 dp."""
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
