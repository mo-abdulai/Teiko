.PHONY: setup pipeline dashboard

setup:
	python -m pip install -r requirements.txt

pipeline:
	python load_data.py
	@echo "Analysis stages will be implemented in subsequent phases."

dashboard:
	streamlit run dashboard.py
