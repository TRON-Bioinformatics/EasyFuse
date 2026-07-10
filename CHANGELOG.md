# TRON-Bioinformatics/easyfuse: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.2.0 - [date]

Release of TRON-Bioinformatics/easyfuse.

### `Added`
- **Multi-container support**: Added Singularity and Docker profiles alongside the existing Conda profile
- **Configurable fusion tools**: Users can now select which fusion detection tools to run (Arriba, STARFusion, FusionCatcher); all three are run by default
- **Improved resource management**: New config files with retry logic to automatically increase allocated resources on failure
- **CI pipeline**: Added GitHub Actions workflows running nf-tests and Python unittests, testing individual modules with data from `TRON-Bioinformatics/test-datasets` using both Conda and Singularity profiles
- **nf-tests**: Added nf-test to test every module that runs in EasyFuse (testdata stored in the TRON-Bioinformatics/testdatasets repo)
- **Python unit tests**: CI pipeline also has Python unit tests that were previously run locally

### `Fixed`
 - Minor improvements in nextflow code to follow best practices

### `Dependencies`

### `Deprecated`
