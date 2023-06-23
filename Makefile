all : clean test

.PHONY: all test clean

clean:
	rm -rf test/output
	# uncomment this line to avoid resume
	#rm -rf work
	rm -f report.html*
	rm -f timeline.html*
	rm -f trace.txt*
	rm -f dag.dot*
	rm -f .nextflow.log*
	rm -rf .nextflow*


test:
	bash test/scripts/test_00.sh
	bash test/scripts/test_01.sh
	bash test/scripts/test_02.sh

