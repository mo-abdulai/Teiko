.PHONY: setup pipeline dashboard

setup:
	python -m pip install -r requirements.txt

pipeline:
	python load_data.py
	python analysis.py
	@echo "TODO: Add Part 3 statistics and Part 4 subset analysis in later phases."

dashboard:
	streamlit run dashboard.py
