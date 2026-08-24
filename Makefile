.PHONY: setup pipeline dashboard

setup:
	python -m pip install -r requirements.txt

pipeline:
	python load_data.py
	python analysis.py
	@echo "TODO: Add Part 3 responder statistics in a later phase."

dashboard:
	streamlit run dashboard.py
