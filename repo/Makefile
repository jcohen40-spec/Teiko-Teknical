.PHONY: setup pipeline dashboard

# Installs all dependencies needed to run the pipeline and dashboard.
setup:
	pip install -r requirements.txt || pip install --break-system-packages -r requirements.txt

# Runs the full data pipeline end to end, with no manual steps:
#   1. Initializes the SQLite database and loads cell-count.csv (Part 1)
#   2. Computes per-sample relative cell population frequencies (Part 2)
#   3. Compares responders vs. non-responders on miraclib, with stats
#      and plots (Part 3)
#   4. Summarizes the baseline melanoma/miraclib/PBMC cohort (Part 4)
pipeline:
	python load_data.py
	python analysis/frequencies.py
	python analysis/response_comparison.py
	python analysis/baseline_summary.py

# Starts the interactive dashboard (Streamlit) on localhost.
dashboard:
	streamlit run dashboard/app.py
