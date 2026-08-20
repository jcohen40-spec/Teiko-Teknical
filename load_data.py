#!/usr/bin/env python3
"""
load_data.py

Create an SQLite database which reads in data from cell-counts.csv

Usage:
    python load_data.py

This will create (or overwrite) `cell_counts.db` in the same directory
as this script.

Schema
~~~~~~
projects
    project_id   TEXT (PRIMARY KEY)

subjects (organizes subjects by project and relevant vars)
    subject_id   TEXT (PRIMARY KEY)
    project_id   TEXT (FOREIGN KEY THAT REFERENCES projects(project_id))
    condition    TEXT
    age          INTEGER
    sex          TEXT
    treatment    TEXT
    response     TEXT        can be "none"

samples (organizes samples by subject of origin and a couple other vars)
    sample_id                  TEXT (PRIMARY KEY)
    subject_id                 TEXT (FOREIGN KEY THAT REFERENCES subjects(subject_id))
    sample_type                TEXT
    time_from_treatment_start  INTEGER

cell_counts (organizes cells by type and sample they belong to)
    sample_id     TEXT (PRIMARY KEY and FOREIGN KEY referencing samples(sample_id))
    b_cell        INTEGER
    cd8_t_cell    INTEGER
    cd4_t_cell    INTEGER
    nk_cell       INTEGER
    monocyte      INTEGER

Design notes
------------
- Subject-level clinical/demographic fields (condition, age, sex,
  treatment, response) are constant across all samples belonging to a
  subject in the source CSV, so they are stored once per subject rather
  than repeated on every sample/cell-count row. This removes redundancy
  and guards against inconsistent updates.
- Cell population counts are stored in "wide" format: one row per
  sample, with one column per cell population. This mirrors the
  structure of the source CSV directly and keeps a 1:1 relationship
  between samples and cell_counts rows (sample_id is both the primary
  key in cell_counts and the foreign key back to samples).
"""

import csv
import sqlite3
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
CSV_PATH = SCRIPT_DIR / "cell-count.csv"
DB_PATH = SCRIPT_DIR / "cell_counts.db"

SCHEMA_SQL = """
PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id   TEXT PRIMARY KEY,
    name         TEXT
);

CREATE TABLE subjects (
    subject_id   TEXT PRIMARY KEY,
    project_id   TEXT NOT NULL REFERENCES projects(project_id),
    condition    TEXT,
    age          INTEGER,
    sex          TEXT,
    treatment    TEXT,
    response     TEXT
);

CREATE TABLE samples (
    sample_id                  TEXT PRIMARY KEY,
    subject_id                 TEXT NOT NULL REFERENCES subjects(subject_id),
    sample_type                TEXT,
    time_from_treatment_start  INTEGER
);

CREATE TABLE cell_counts (
    sample_id     TEXT PRIMARY KEY REFERENCES samples(sample_id),
    b_cell        INTEGER NOT NULL,
    cd8_t_cell    INTEGER NOT NULL,
    cd4_t_cell    INTEGER NOT NULL,
    nk_cell       INTEGER NOT NULL,
    monocyte      INTEGER NOT NULL
);

CREATE INDEX idx_subjects_project ON subjects(project_id);
CREATE INDEX idx_samples_subject ON samples(subject_id);
"""


def init_db(conn: sqlite3.Connection) -> None:
    """Create the database schema (dropping any existing tables first)."""
    conn.executescript(SCHEMA_SQL)
    conn.commit()


def load_csv(conn: sqlite3.Connection, csv_path: Path) -> None:
    """Read cell-count.csv and populate the normalized tables."""
    cur = conn.cursor()

    seen_projects = set()
    seen_subjects = set()

    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            project_id = row["project"]
            subject_id = row["subject"]
            sample_id = row["sample"]

            if project_id not in seen_projects:
                cur.execute(
                    "INSERT OR IGNORE INTO projects (project_id, name) VALUES (?, ?)",
                    (project_id, project_id),
                )
                seen_projects.add(project_id)

            if subject_id not in seen_subjects:
                cur.execute(
                    """
                    INSERT OR IGNORE INTO subjects
                        (subject_id, project_id, condition, age, sex, treatment, response)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        subject_id,
                        project_id,
                        row["condition"],
                        int(row["age"]),
                        row["sex"],
                        row["treatment"],
                        row["response"] if row["response"] not in (None, "") else None,
                    ),
                )
                seen_subjects.add(subject_id)

            cur.execute(
                """
                INSERT INTO samples
                    (sample_id, subject_id, sample_type, time_from_treatment_start)
                VALUES (?, ?, ?, ?)
                """,
                (
                    sample_id,
                    subject_id,
                    row["sample_type"],
                    int(row["time_from_treatment_start"]),
                ),
            )

            cur.execute(
                """
                INSERT INTO cell_counts
                    (sample_id, b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    sample_id,
                    int(row["b_cell"]),
                    int(row["cd8_t_cell"]),
                    int(row["cd4_t_cell"]),
                    int(row["nk_cell"]),
                    int(row["monocyte"]),
                ),
            )

    conn.commit()


def main() -> None:
    if not CSV_PATH.exists():
        raise FileNotFoundError(
            f"Could not find {CSV_PATH}. Make sure cell-count.csv is in the "
            "same directory as load_data.py."
        )

    # Start fresh each run.
    if DB_PATH.exists():
        DB_PATH.unlink()

    conn = sqlite3.connect(DB_PATH)
    try:
        init_db(conn)
        load_csv(conn, CSV_PATH)

        # Quick sanity check / summary printed to stdout.
        cur = conn.cursor()
        counts = {}
        for table in ("projects", "subjects", "samples", "cell_counts"):
            cur.execute(f"SELECT COUNT(*) FROM {table}")
            counts[table] = cur.fetchone()[0]

        print(f"Database created at: {DB_PATH}")
        for table, n in counts.items():
            print(f"  {table}: {n} rows")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
