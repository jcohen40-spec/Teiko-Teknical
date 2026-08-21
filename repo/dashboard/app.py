#!/usr/bin/env python3
"""
dashboard/app.py

Interactive Streamlit dashboard for Bob's immune cell population
analysis (Parts 2-4). Reads directly from cell_counts.db and reuses
the same analysis functions as the standalone scripts in analysis/,
so the dashboard and the CLI scripts can never disagree with each
other -- there is exactly one implementation of each calculation.

Usage
-----
    streamlit run dashboard/app.py

(or, from the repo root: `make dashboard`)
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

# Make the repo root importable so `import analysis...` works regardless
# of the working directory Streamlit is launched from.
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from analysis.frequencies import compute_frequencies, POPULATION_COLUMNS
from analysis.response_comparison import (
    get_filtered_frequencies,
    run_statistical_tests,
    CONDITION as RESP_CONDITION,
    TREATMENT as RESP_TREATMENT,
    SAMPLE_TYPE as RESP_SAMPLE_TYPE,
)
from analysis.baseline_summary import (
    get_baseline_cohort,
    samples_per_project,
    subjects_by_response,
    subjects_by_sex,
)

from analysis.melanoma_male_responder_bcell import average_b_cell_count

DB_PATH = REPO_ROOT / "cell_counts.db"

st.set_page_config(page_title="Loblaw Bio | Immune Cell Dashboard", layout="wide")


@st.cache_data
def load_frequencies():
    return compute_frequencies(DB_PATH)


@st.cache_data
def load_response_comparison():
    filtered = get_filtered_frequencies(DB_PATH)
    stats_df = run_statistical_tests(filtered)
    return filtered, stats_df


@st.cache_data
def load_baseline_cohort():
    return get_baseline_cohort(DB_PATH)


def main():
    st.title("Immune Cell Population Dashboard")
    st.caption("Loblaw Bio — miraclib / phauximab clinical trial cell count analysis")

    if not DB_PATH.exists():
        st.error(
            f"Could not find `{DB_PATH.name}`. Run `python load_data.py` "
            "(or `make pipeline`) from the repo root first, then reload."
        )
        st.stop()

    tab1, tab2, tab3, tab4 = st.tabs(
        [
            "Part 2 — Cell Frequencies",
            "Part 3 — Responder vs. Non-Responder",
            "Part 4 — Baseline Cohort",
            "Part 5 — Additional Question",
        ]
    )

    # ---------------------------------------------------------------
    # Part 2: relative frequency table
    # ---------------------------------------------------------------
    with tab1:
        st.header("Relative frequency of each cell population, per sample")
        st.write(
            "For every sample, the total cell count (sum across all five "
            "populations) and each population's share of that total."
        )

        freq = load_frequencies()

        samples = sorted(freq["sample"].unique())
        selected_samples = st.multiselect(
            "Filter by sample (leave empty to show all)",
            options=samples,
            default=[],
        )
        table = freq if not selected_samples else freq[freq["sample"].isin(selected_samples)]

        st.dataframe(table, use_container_width=True, height=400)
        st.caption(f"{table['sample'].nunique()} samples, {len(table)} rows shown.")

        st.download_button(
            "Download full table as CSV",
            data=freq.to_csv(index=False),
            file_name="frequencies.csv",
            mime="text/csv",
        )

        st.subheader("Average relative frequency by population")
        avg_by_pop = (
            freq.groupby("population", observed=True)["percentage"]
            .mean()
            .reindex(POPULATION_COLUMNS)
        )
        st.bar_chart(avg_by_pop)

    # ---------------------------------------------------------------
    # Part 3: responder vs non-responder
    # ---------------------------------------------------------------
    with tab2:
        st.header("Responders vs. non-responders")
        st.write(
            f"Filtered to **{RESP_CONDITION}**, **{RESP_TREATMENT}**, "
            f"**{RESP_SAMPLE_TYPE}** samples only."
        )

        filtered, stats_df = load_response_comparison()

        n_resp = filtered.loc[filtered["response"] == "yes", "sample"].nunique()
        n_non_resp = filtered.loc[filtered["response"] == "no", "sample"].nunique()
        col1, col2 = st.columns(2)
        col1.metric("Responders (n)", n_resp)
        col2.metric("Non-responders (n)", n_non_resp)

        st.subheader("Statistical comparison (t-test, Bonferroni-corrected)")
        st.dataframe(
            stats_df.style.format(
                {
                    "mean_responder_pct": "{:.2f}",
                    "mean_non_responder_pct": "{:.2f}",
                    "mean_diff": "{:.3f}",
                    "cohens_d": "{:.3f}",
                    "t_statistic": "{:.2f}",
                    "p_value": "{:.4f}",
                    "p_value_adj": "{:.4f}",
                }
            ),
            use_container_width=True,
        )
        sig = stats_df.loc[stats_df["significant"], "population"].tolist()
        if sig:
            st.success(f"Significant after correction (p_adj < 0.05): {', '.join(sig)}")
        else:
            st.info("No populations reached significance after correction.")

        st.subheader("Boxplot: relative frequency by response")
        selected_pop = st.selectbox(
            "Cell population", options=POPULATION_COLUMNS, index=2  # default cd4_t_cell
        )
        sub = filtered[filtered["population"] == selected_pop]
        fig, ax = plt.subplots(figsize=(5, 4))
        groups = [
            sub.loc[sub["response"] == "no", "percentage"],
            sub.loc[sub["response"] == "yes", "percentage"],
        ]
        ax.boxplot(groups, tick_labels=["Non-responder", "Responder"], showmeans=True)
        ax.set_ylabel("Relative frequency (%)")
        ax.set_title(selected_pop)
        ax.grid(axis="y", alpha=0.3)
        st.pyplot(fig)
        plt.close(fig)

        st.subheader("Mean \u00b1 SEM, all populations")
        means = (
            filtered.groupby(["population", "response"], observed=True)["percentage"]
            .agg(["mean", "sem"])
            .reset_index()
        )
        fig2, ax2 = plt.subplots(figsize=(9, 4))
        offset = 0.08
        for i, pop in enumerate(POPULATION_COLUMNS):
            for resp_val, color, off in [("no", "tab:orange", -offset), ("yes", "tab:blue", offset)]:
                row = means[(means["population"] == pop) & (means["response"] == resp_val)]
                if not row.empty:
                    ax2.errorbar(
                        i + off, row["mean"].values[0], yerr=row["sem"].values[0],
                        fmt="o", capsize=4, color=color,
                        label=("Non-responder" if resp_val == "no" else "Responder") if i == 0 else None,
                    )
        ax2.set_xticks(range(len(POPULATION_COLUMNS)))
        ax2.set_xticklabels(POPULATION_COLUMNS)
        ax2.set_ylabel("Mean relative frequency (%) \u00b1 SEM")
        ax2.legend()
        ax2.grid(axis="y", alpha=0.3)
        st.pyplot(fig2)
        plt.close(fig2)

    # ---------------------------------------------------------------
    # Part 4: baseline cohort
    # ---------------------------------------------------------------
    with tab3:
        st.header("Baseline melanoma / miraclib / PBMC cohort")
        st.write("Samples at `time_from_treatment_start == 0` only.")

        cohort = load_baseline_cohort()
        n_samples = cohort["sample_id"].nunique()
        n_subjects = cohort["subject_id"].nunique()
        st.metric("Baseline samples", n_samples)
        st.caption(f"from {n_subjects} subjects")

        col1, col2, col3 = st.columns(3)

        with col1:
            st.subheader("Samples per project")
            by_project = samples_per_project(cohort)
            st.dataframe(by_project, use_container_width=True, hide_index=True)
            st.bar_chart(by_project.set_index("project_id")["n_samples"])

        with col2:
            st.subheader("Subjects by response")
            by_response = subjects_by_response(cohort)
            st.dataframe(by_response, use_container_width=True, hide_index=True)
            st.bar_chart(by_response.set_index("response")["n_subjects"])

        with col3:
            st.subheader("Subjects by sex")
            by_sex = subjects_by_sex(cohort)
            st.dataframe(by_sex, use_container_width=True, hide_index=True)
            st.bar_chart(by_sex.set_index("sex")["n_subjects"])

    # ---------------------------------------------------------------
    # Part 5: melanoma male responder B-cell question
    # ---------------------------------------------------------------
    with tab4:
        st.header("Additional Question")
        st.write(
            "Considering melanoma males of all sample and treatment types, "
            "what is the average number of B cells for responders at time=0?"
        )

        avg_b_cells = average_b_cell_count(DB_PATH)

        st.metric(
            "Average B cell count",
            f"{avg_b_cells:.2f}",
        )

        st.caption(
            "Filters: melanoma | male | responder | "
            "time_from_treatment_start = 0 | all sample types | all treatments"
        )

if __name__ == "__main__":
    main()


