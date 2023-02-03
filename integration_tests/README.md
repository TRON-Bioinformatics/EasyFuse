# Integration tests

This folder contains scripts that run integration tests over easy fuse.

These tests, when possible, are integrated in the continuous integration environment.

## Requirements

Requirements to run the tests:
- To run these tests easy fuse Python package must be installed and available from the command line
- Shell scripts (ie: `run.sh`) must be called from the root of the repository as they contain relative paths

## Test cases

The test `00_test_case` runs the full pipeline of easy fuse. It requires a `config.ini` file and a functional set of
references. This test is not fully automated. This test does not run in the continuous integration environment.
The expected results of this run are provided under the folder `00_test_case/expected`.

The rest of the test cases test smaller parts of the pipeline that do not have any major dependencies of software or 
references.

The test `06_summarize_data` requires R 4.2.0 and the R libraries `optparse`, `dplyr`, `tidyr`, `tidyselect`, `readr`, 
`stringr`, `XML` and `randomForest` to be installed.

The subcommands `skewer-wrapper`, `soapfuse-wrapper`, `annotation` and `star-index` do not have tests due to the 
dependencies.

To run all automated test cases run: `make integration_tests`
