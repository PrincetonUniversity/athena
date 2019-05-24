# Overview

Regression tests are comparisons of code output to accepted solutions. They are used to ensure that parts of the code that work are not broken in the course of development. There should ultimately be a suite of stringent tests utilizing different parts of the code. Modifications to the code should not be pushed to the repository if they break regression tests.

This guide will describe how to use regression tests, as well as how to write new ones. See also the [condensed presentation][1].

All regression testing is confined to `athena/tst/regression/`; running these tests should have no effect outside this directory.

The regression tests for Athena++ are written in Python. They should be compatible with any Python 2.7 or 3 distribution. If they break in *later* versions of Python, please file a bug report and we can try to fix any incompatibilities.

# Usage

From the top-level `athena/` directory go to the regression test directory:
```ShellSession
cd tst/regression
```
To run all test scripts that have been written, simply call
```ShellSession
python run_tests.py
```

Instead of running all tests, selected tests can be run with
```ShellSession
python run_tests.py <name_1> [<name_2> ...]
```
where each name is either of a suite of tests or an individual test within a suite. For example, the `scripts/tests/gr/` directory contains the general relativity suite, which can be run with
```ShellSession
python run_tests.py gr
```
One specific test in the suite is defined by the file `scripts/tests/gr/mhd_shocks_hlld.py`, and this can be run with
```ShellSession
python run_tests.py gr/mhd_shocks_hlld
```

# Writing new tests

This section is written from the perspective of the `athena/tst/regression/` directory.

## Test locations

If a new suite of tests is desired (for example, there should be separate suites for relativity, radiation, AMR, etc.), create a new directory under `scripts/tests/`. In order for Python to recognize the directory as containing scripts, it must have a file named `__init__.py`. This file can and probably should be left empty.

In order to create a new regression test, start with a new file `<test_name>.py` in `scripts/tests/<suite_name>/`. At this time, the regression test framework does not support any further nested directories. No other file need be modified; in particular, the test does not need to be enrolled anywhere. In fact, the top-level script that is called, `run_tests.py`, should *not* be modified.

A new test may very well need a new saved dataset to compare against. All such data files are kept in `data/`. New files, usually Athena++ outputs, can be added here if need be. Data files in this directory are considered part of the repository: new files should be committed, and other developers' files should not be tampered with. Please refrain from including excessively large and unwieldy files.

**Note**: Storing reference solutions and comparing them to results of regression tests should only be performed for a limited set of tests. Exact comparison is impossible since that may depend on compiler/processor/library version, and inexact comparison is subjective. Analytic solutions or error convergence are preferred metrics for analyzing the regression test results. 

One may want to write common functions for use in other tests. These can be placed in `scripts/utils/` as detailed below.

## Definition of a test

The file `<test_name>.py` is a valid test if it defines three valid Python functions, `prepare()`, `run()`, and `analyze()`, that take no inputs, where the last returns either `True` or `False` depending on whether the test passes or not.

There are no other requirements for a test to work, but of course tests should behave themselves and should not modify/delete anything that could affect the rest of the code.

See a preexisting file, especially [`example.py`][2], for a working example.

## Side effects

Running regression tests *should* leave the files outside the regression directory unchanged. This is tricky, since the configure script overwrites `athena/Makefile` and `athena/src/defs.hpp`. To work around this, the main regression script `run_tests.py` (which should not be modified) saves copies of these files, should they exist, before running any tests. Upon completing tests or encountering an exception, the regression script writes the original files back (or deletes whatever is there if those files did not exist before).

Because calling `make` on the `Makefile` output by the configure script would overwrite the objects and executable in the main directory, one must be careful in calling `make`. Instead of calling it directly, the `make()` function in `scripts/utils/athena.py` overrides the definitions in the `Makefile` to place the executable and all objects inside the regression directory. The executable will then place its outputs in the same location.

## Utilities

Because most regression tests will perform similar actions, a set of common utility functions is provided in `scripts/utils/`. This includes for example functionality for taking L1 error norms.

To use these functions in a script, one must import the appropriate file. For example, from the test implementation in `scripts/tests/<test_name>.py`, one can call the `make()` function in `scripts/utils/athena.py` by writing
```Python
import scripts.utils.athena as athena
athena.make()
```

### Athena++ interfacing

The file `scripts/utils/athena.py` defines Python routines for configuring, compiling, and running Athena++. The publicly useful functions are detailed here.

##### `configure(*args, **kwargs)`

- Inputs:
  - Positional: any number of strings, to be prefaced with `-` when calling `athena/configure.py`; must go before keywords
  - Keywords: any number of `key='val'` assignments, to be used as `--key=val` when calling `athena/configure.py`; must go after positional arguments
- Outputs: (none)
- Example: `configure('omp', prob='shock_tube')` will result in the code being configured as though one called `python configure.py -omp --prob=shock_tube` from the home directory.

##### `make()`

- Inputs: (none)
- Outputs: (none)
- Example: Simply call `make()` and the code will be compiled into `obj/` and linked into `bin/`. Note these are the subdirectories of `athena/tst/regression/`, not those of the same names in `athena/`.

##### `run(input_filename, arguments)`

- Inputs:
  - `input_filename`: string containing the path to the desired input file, relative to `athena/inputs/`
  - `arguments`: list of strings, one for each command-line argument to be fed to Athena++
- Outputs: (none)
- Example: The following code snippet runs Athena++ with a relativistic shock tube, ensuring the CFL number, simulation time, and number of grid points are correct:
```Python
import scripts.utils.athena as athena
arguments = ['time/cfl_number=0.4', 'time/tlim=0.4', 'mesh/nx1=400']
athena.run('hydro_sr/athinput.mb_1', arguments)
```

### Comparisons between datasets

Some utilities for comparing one dataset to another are provided in `scripts/utils/comparison.py`. The available functions are documented here.

##### `l1_norm(faces, vals)`

- Inputs:
  - `faces`: 1D array of N+1 interface locations
  - `vals`: 1D array of N volume-averaged values of a function `f`
- Outputs: Returns the integral of the absolute value of `f` over the domain.

##### `l1_diff(faces_1, vals_1, faces_2, vals_2)`

- Inputs:
  - `faces_1`: 1D array of N1+1 interface locations
  - `vals_1`: 1D array of N1 volume-averaged values of a function `f`
  - `faces_2`: 1D array of N2+1 interface locations
  - `vals_2`: 1D array of N2 volume-averaged values of a function `g`
- Outputs: Returns the integral of the absolute value of `f-g` over the domain. This is calculated by taking the union of all N1+N2+2 interface locations to define a refined grid, pushing each of the N1 cells' values for `f` into the corresponding new cells, and doing likewise for `g`. The procedure of `l1_norm()` is applied to the refined grid and the unambiguously calculated `f-g` values on that grid. Note no interpolation is involved; the answer is mathematically exact to the extent `f` and `g` are represented by step functions.

### Reading data into Python

See the [[documentation|Reading Data into Python]] on the general-purpose Python scripts for reading Athena++ data.

  [1]: presentations/regression_testing.pdf
  [2]: https://github.com/PrincetonUniversity/athena-public-version/blob/master/tst/regression/scripts/tests/example.py
