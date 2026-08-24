.PHONY: setup pipeline dashboard

setup:
	python -m pip install -r requirements.txt

pipeline:
	# TODO: Replace this placeholder with the full reproducible pipeline.
	@echo "Pipeline placeholder: this will eventually execute:"
	@echo "  python load_data.py"
	@echo "  python analysis.py"

dashboard:
	streamlit run dashboard.py
