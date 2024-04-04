# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.3] - 2024-04-04

### Added

- Arriba v2.4.0 for fusion prediction
- High confidence filter for Arriba
- Conversion, parsing and annotation scripts
- Easyquant v0.5.2 for requantification
- Unit/integration tests using `pytest`

### Changed

- rm `tool_frac` column in R script
- Fixed issue in fusion annotation script
- Updated prediction model based on new results

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


[2.1.0]: https://github.com/TRON-Bioinformatics/EasyFuse/v2.0.2...v2.1.0
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
