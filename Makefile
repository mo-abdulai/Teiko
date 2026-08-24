.PHONY: setup pipeline dashboard

setup:
	python -m pip install -r requirements.txt

pipeline:
	python load_data.py
	python analysis.py
	@echo "Pipeline complete. Run 'make dashboard' to open the Streamlit dashboard."

dashboard:
	streamlit run dashboard.py
