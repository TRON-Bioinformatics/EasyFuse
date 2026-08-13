# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.2.1] - 2026-08-14

### Feature Additions
- EasyFuse now supports mouse index

### Code Refactoring
- Validation subworkflow have been updated to add an extra argument `--species` which is set to `human` by default
- Updated readme and know-how for using easyfuse with mouse index
- Updated pipeline config
- Updated pipeline schema

## [2.2.0] - 2026-08-10

#### Fixed

- Added `fusion_protein_sequence` and `fusion_protein_sequence_bp` to the merged final output table.

#### Code Refactoring
- **Modules**: the earlier modules have been restructed into per tool module file along with the environment.yml file, this reduces maintainence overhead
    - Each tool is now encapsulated in its own module directory under `modules/`, containing:
        - `main.nf` — Nextflow process definition
        - `environment.yml` — per-tool Conda environment specification
        - `tests/` — nf-test test cases
- **Subworkflows**: closely related modules are clubbed together into subworkflows to reduce the final pipeline code
- **main.nf**: the main pipeline script now only has the named subworkflows that are part of the pipeline flow diagram in the main README file.

#### Technical Improvements
- **Input validation**: easyfuse now uses the input validation subworkflow to validate the input samplesheet and parameters
- **Configurable fusion tools**: Users can now select which fusion detection tools to run (Arriba, STARFusion, FusionCatcher); all three are run by default
- **Improved resource management**: New config files with retry logic to automatically increase allocated resources on failure
- **Multi-container support**: Added Singularity and Docker profiles alongside the existing Conda profile
- Three execution profiles are now fully supported:
    - **`conda`** — builds per-process Conda environments from `environment.yml` files
    - **`singularity`** — pulls Seqera Wave containers; supports local cache via `NXF_SINGULARITY_CACHEDIR`
    - **`docker`** — pulls Docker images from Seqera Wave registry
    - **`slurm`** — can be combined with any container profile for HPC cluster execution

#### Infrastructure & Dependencies
 - Process resource labels (`process_single`, `process_low`, `process_medium`, `process_high`) now scale linearly with `task.attempt`, so resources are automatically increased on retry.
 - Exit codes `130–145`, `104`, and `175` trigger a retry; all other failures terminate the task.

#### Development & Testing
- **CI pipeline**: Added GitHub Actions workflows running nf-tests and Python unittests, testing individual modules with data from `TRON-Bioinformatics/test-datasets` using both Conda and Singularity profiles
- **nf-tests**: Added nf-test to test every module that runs in EasyFuse (testdata stored in the TRON-Bioinformatics/testdatasets repo)
- **Python unit tests**: CI pipeline also has Python unit tests that were previously run locally
- - Triggers on pull requests and releases
    - Supports self-hosted runners for `TRON-Private` organisation branches; uses GitHub-hosted `ubuntu-24.04` runners for `main`/`dev`
    - Tests with Nextflow `25.04.0` (required) and `latest-everything` (non-blocking)
    - Tests with the `singularity` profile against the test config
    - A `confirm-pass` job enforces that all required checks pass before merging


## [2.1.0] - 2026-02-13

### Functional changes

- No changes to pipeline output or results; this release focuses on internal refactoring and maintainability improvements.
- In rare cases, `no_frame` or `neo_frame` fusion annotations may change due to improved handling of reference annotations.

### Internal & Technical Changes

#### Code Refactoring

- Refactored fusion annotation to modular fusionannotator.py and supporting components
- Refactored fusion parsing; fusiontoolparser.py now consumes standardized per-tool CSVs via repeated --tool and writes Detected_Fusions.csv with consistent headers
- Nextflow workflow reorganized: parsing moved to modules/04_fusionparsing.nf, annotation to modules/05_fusionannotation.nf; downstream stages renumbered
- Retrained model due to slight changes in annotations

#### Technical Improvements

- Robust header-based range extraction in read_selection.py for wildtype ranges
- Stricter chromosome filtering to primary contigs and strand handling in tool parsers
- Pipeline outputs aligned for ARRIBA (only fusions.tsv) and downstream consumers


#### Infrastructure & Dependencies

- Updated environments: conda channels switched to nodefaults
- Removed logzero dependency and related logging calls; reduced log output across scripts
- Removed legacy monolithic fusion annotation script (replaced by fusionannotator.py)
- Removed ARRIBA discarded output from pipeline (structural cleanup; discarded fusions were not used downstream before, only high-confidence calls proceed)

#### Development & Testing

- New utilities for annotation: gff3_to_db.py (build gffutils DB) and gtf2tsl.py (extract TSL)
- New modular fusion parsing with per-tool parsers (Arriba, STAR-Fusion, FusionCatcher, InFusion, MapSplice, SOAPfuse) producing standardized CSV via parse_tool.py
- New Nextflow parsing processes (PARSE_ARRIBA, PARSE_STAR_FUSION, PARSE_FUSION_CATCHER) and dedicated conda env (environments/fusionparsing.yml)
- Unit tests and test runner for fusion annotation module

## [2.0.4] - 2024-12-09

### Added

- Added full length protein sequence to the final output
- Specify computational requirements via predefined labels: single, low and medium

### Changed

- Updated NextflowVersion to 24.10.1
- Updated resource management
- Fixed exon count in final output
- Fixed tool_frac column in final output
- Updated prediction model based on new results

## [2.0.3] - 2024-04-04

### Added

- Arriba v2.4.0 high confidence calls as fusion candidates
- [easyquant] (https://github.com/TRON-Bioinformatics/Easyquant) v0.5.2 for read support requantification
- Unit/integration tests using pytest

### Changed

- Fixed issue with gene names in fusion annotation script
- Updated prediction model based on new results
- Moved conversion, parsing and annotation code from the easyfuse-src package
- Removed unnecessary columns from final output

## [2.0.2] - 2023-11-24

### Changed

- Upgraded pipeline to Ensembl v110
- Updated to FusionCatcher v1.33


## [2.0.1] - 2023-08-11

### Changed

- Simplified installation and dependency management through migration of EasyFuse package to Bioconda

### Fixed

- Fixed bug in QC workflow


## [2.0.0] - 2023-07-07

### Added

- EasyFuse as NextFlow pipeline for increased usability, stability, and scalability
- Python code as python package outsourced to separate repository
- Internal detection tools were reduced to StarFusion and FusionCatcher
- Prediction model has been changed to EF_requant_type to not rely on specific tool features
- Overall reduced detection performance in sensitivity and precision compared to EasyFuse 1.3.7


## [1.3.7] - 2022-12-15

### Changed

- Updated models and provide additional models for feature subsets
- Cleaned code and made it more robust
- Updated error handling
- Cleaned up Dockerfile and made versioning more strict

### Fixed

- Fixed bugs related to Python compatibility
- Fixed read counts from tools in final results table


## [1.3.6] - 2022-07-21

### Added

- Add support for Singularity
- Update example output files
- Make Dockerfile more flexible
- Update README


## [1.3.5] - 2022-06-20

### Added

- used a breakpoint-specific identifier (BPID) for joined annotation and in output
- new output file names
- separate output files for predicted fusions .pred.csv and all candidates .all.csv
- new output format including column BPID
- retrained model on new output column format
- cleaned up R code and updated R dependencies
- added Docker example scripts with test data and run_test.sh script
- added support for INI and JSON config files and make them more user-friendly

### Fixed

- fixed content of columns <tool>_detected, tool_count, and tool_frac
- fixed several bugs in input file/folder parsing


## [1.3.4] - 2021-05-29

### Changed

- updated prediction models EF_full (default) and EF_full_rep (with replicates).

### Fixed

- Fixed bug in fastqc parsing from different sources


## [1.3.3] - 2021-05-29

### Added

- Added option to select Fusioncatcher index

### Fixed

- Fixed some bugs for non-queueing environments


## [1.3.2] - 2021-05-29

### Added

- include updated prediction model
- remove the creation of unused files


## [1.3.1] - 2021-05-29

### Added

- Restructured code
- Switched to Python config file as it provides more convenience and better usability

### Fixed 

- Fixed some bugs in summarization scripts
- Fixed some bugs when no fusions are found


## [1.3.0] - 2021-05-21

### Added

- Initial release on GitHub


[Unreleased]: https://github.com/TRON-Bioinformatics/EasyFuse/v2.0.3...dev
[2.0.4]: https://github.com/TRON-Bioinformatics/EasyFuse/v2.0.3...v2.0.4
[2.0.3]: https://github.com/TRON-Bioinformatics/EasyFuse/v2.0.2...v2.0.3
[2.0.2]: https://github.com/TRON-Bioinformatics/EasyFuse/v2.0.1...v2.0.2
[2.0.1]: https://github.com/TRON-Bioinformatics/EasyFuse/v2.0.0...v2.0.1
[2.0.0]: https://github.com/TRON-Bioinformatics/EasyFuse/v1.3.7...v2.0.0
[1.3.7]: https://github.com/TRON-Bioinformatics/EasyFuse/v1.3.6...v1.3.7
[1.3.6]: https://github.com/TRON-Bioinformatics/EasyFuse/v1.3.5...v1.3.6
[1.3.5]: https://github.com/TRON-Bioinformatics/EasyFuse/v1.3.4...v1.3.5
[1.3.4]: https://github.com/TRON-Bioinformatics/EasyFuse/v1.3.3...v1.3.4
[1.3.3]: https://github.com/TRON-Bioinformatics/EasyFuse/v1.3.2...v1.3.3
[1.3.2]: https://github.com/TRON-Bioinformatics/EasyFuse/v1.3.1...v1.3.2
[1.3.1]: https://github.com/TRON-Bioinformatics/EasyFuse/v1.3.0...v1.3.1
[1.3.0]: https://github.com/TRON-Bioinformatics/EasyFuse/releases/tag/v1.3.0
