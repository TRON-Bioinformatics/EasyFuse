# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


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
