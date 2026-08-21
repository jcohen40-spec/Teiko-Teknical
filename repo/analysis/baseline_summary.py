#!/usr/bin/env python3
"""
analysis/baseline_summary.py

Part 4: Baseline Melanoma / Miraclib / PBMC Cohort Summary

Identifies all melanoma, PBMC samples at baseline
(time_from_treatment_start == 0) from patients treated with miraclib,
then breaks that cohort down by:

    - number of samples per project
    - number of subjects who are responders vs. non-responders
    - number of subjects who are male vs. female

Usage
-----
    python analysis/baseline_summary.py

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

REPO_ROOT = Path(__file__).resolve().parent.parent
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
    """Count of baseline samples, broken down by project."""
    return (
        cohort.groupby("project_id")["sample_id"]
        .nunique()
        .rename("n_samples")
        .reset_index()
        .sort_values("project_id")
        .reset_index(drop=True)
    )


def subjects_by_response(cohort: pd.DataFrame) -> pd.DataFrame:
    """
    Count of unique subjects broken down by responder / non-responder.
    Counts subjects (not samples) since response is a subject-level
    attribute -- one subject shouldn't be double counted even if this
    cohort happened to include more than one sample per subject.
    """
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
    """Count of unique subjects broken down by sex (M / F)."""
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
