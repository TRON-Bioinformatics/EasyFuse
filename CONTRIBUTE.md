# Developer guide


## Run the tests

After the EasyFuse Python package is installed as described in the README, run the integration tests with `make`.

You may need to edit the test scripts to use an appropriate reference folder in `--reference`.

## Prepare a new release

Increase the version in the property `VERSION` in the file `nextflow.config`.

## Test the integration of easyfuse-src

The Python and R code are distributed as a Python package and this is imported in 
EasyFuse pipeline via conda. But, for development purposes the import of this
package can be skipped while all other dependencies are still managed via conda.

Use the flag `--disable_pyeasyfuse_conda`, this will rely on whatever version of
pyeasyfuse has been installed on your environment.
