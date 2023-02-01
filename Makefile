
all : clean integration_tests

.PHONY: all clean integration_tests

clean:
	rm -rf output

integration_tests:
	bash integration_tests/00_help/run.sh
	bash integration_tests/01_qc_parser/run.sh
	bash integration_tests/02_read_filter/run.sh
	bash integration_tests/03_fusion_parser/run.sh
	bash integration_tests/06_requantify/run.sh
