# Developer guide

## Set up the environment

EasyFuse requires Python 3.7.
The EasyFuse package is described in the file `pyproject.toml` as described in [PEP621](https://packaging.python.org/en/latest/specifications/declaring-project-metadata/#declaring-project-metadata). 

We use Poetry (`pip install poetry`) to manage the development environment, you only need to run:
```
poetry install
```

This will install the right Python interpreter and all of its dependencies. 
To avoid automatic updates of the dependencies we use the poetry lock which automatically records the entire dependency tree.
If you need to add a new dependency, add it to the `pyproject.toml` file and run `poetry lock` to update the lock file.

After installation EasyFuse will be available in the command line as `poetry run easy-fuse --help`.

## Run unit tests

Run:
```
poetry run python -m unittest discover tests
```

## Run the integration tests

See the [integration tests README](integration_tests/README.md).

## Build, distribute and install

To build the package into a wheel file run: `poetry build`

This will create a wheel file under the `dist` folder such as `easy_fuse-x.y.z-py3-none-any.whl`

This wheel file can be distributed and installed with `pip install easy_fuse-x.y.z-py3-none-any.whl`

After installation, EasyFuse should be available in the command line: `easy-fuse --help`

## Cleaning the code

PyCharm provides a useful tool to reformat code under `code -> reformat`.
This allows to organize imports and enforce basic style guidelines.

Black provides a more advanced tool to reformat code. It can be installed with `pip install black`. 
Then just run `black easy_fuse` to reformat the code.

Vulture enables inspecting the code base for dead code. Install vulture with `pip install vulture`. 
Then run as follows: `vulture easy_fuse`. Do not trust a vulture blindly.

Autoflake is a tool to remove unused imports. Install it with `pip install autoflake`. 
Then run as follows: `autoflake --in-place --remove-all-unused-imports --remove-unused-variables -r easy_fuse`

