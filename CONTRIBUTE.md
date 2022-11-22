# Developer guide

## Set up the environment

EasyFuse requires Python 3.7.

Install all Python dependencies as follows: `pip install -r requirements.txt`

## Run unit tests

Run:
```
python -m unittest discover tests
```

There is no need to install EasyFuse to run the tests.

## Build and install

To build the package into a wheel file first install if not already done: `pip install wheel`

Then build: `python setup.py bdist_wheel`

This will create a wheel file under the dist folder such as `easy_fuse-x.y.z-py3-none-any.whl`

Install the wheel file: `pip install dist/easy_fuse-x.y.z-py3-none-any.whl`

Now, EasyFuse should be available in the command line: `easy-fuse --help`

## Cleaning the code

PyCharm provides a useful tool to reformat code under `code -> reformat`.
This allows to organize imports and enforce basic style guidelines.

Vulture enables inspecting the code base for dead code. Install vulture with `pip install vulture`. 
Then run as follows: `vulture easy_fuse`. Do not trust a vulture blindly.

