### Python utilities
Although Athena++ is primarily a C++ code, Python files serve several important auxiliary functions in the project:
* [[Configuring]] the solver with `./configure.py` in the root directory
* [[Regression Testing]] with scripts in the `tst/regression/` folder
* [[Reading Data Into Python]] and visualizing with scripts in the `vis/` folder

### PEP 8 style guidelines
The ~7,000 lines of Python code in the project (*as of 6/11/18*) therefore merit a set of policies to ensure their maintainability. When modifying the `.py` files in the above locations or adding new Python to the repository, follow [PEP 8 (Style Guide for Python Code)](https://www.python.org/dev/peps/pep-0008/) and [PEP 257 (Docstring Conventions)](https://www.python.org/dev/peps/pep-0257/) guidelines. Some of the most important rules, taken directly from PEP and [pycodestyle Error Codes](http://pycodestyle.pycqa.org/en/latest/intro.html#error-codes), to follow are: (commentary/reasoning in italics and Athena++ modifications to the rules in bold)
* Use 4 spaces per indentation level. Use spaces, not tabs. `[E101, E111]`
  * *The 4 space indent accentuates nested blocks, which can prevent nasty bugs in a whitespace-sensitive language like Python.*
* The Python standard library is conservative and requires limiting lines to 79 characters (and docstrings/comments to 72). [...] it is okay to increase the nominal line length from 80 to 100 characters `[E501]`
  * **Athena++ change**: Relaxed to 90 characters to be consistent with above C++ style
* Function names should be lowercase, with words separated by underscores as necessary to improve readability.
* Whitespace:
  * Always surround these binary operators with a single space on either side: assignment (=), augmented assignment (+=, -= etc.), comparisons (==, <, >, !=, <>, <=, >=, in, not in, is, is not), Booleans (and, or, not). The exception to that is when `=` is used to set named parameters. `[E225]`
  * **Athena++ change**: "Use whitespace around arithmetic operator" `[E226]` is disabled so that array indexing with expressions, such as `rho[nx3+2, nx2-1, nx1+i]`, can save whitespace.
  * Use whitespace after ‘,’ , ‘;’, or ‘:’ `[E231]`
* Blank lines:
  * Surround top-level function and class definitions with two blank lines. `[E302]`
  * Method definitions inside a class are surrounded by a single blank line. `[E301]`
* The preferred way of wrapping long lines is by using Python's implied line continuation inside parentheses, brackets and braces. Long lines can be broken over multiple lines by wrapping expressions in parentheses. These should be used in preference to using a backslash for line continuation. `[E122, E124, E127, E128]`
  * *There is a potential danger when using `\` for line continuation when there are stray trailing whitespace characters. Using `(...)` to wrap an expression over multiple lines is more expressive.*
* Continuation lines should align wrapped elements either vertically using Python's implicit line joining inside parentheses, brackets and braces, or using a hanging indent. When using a hanging indent the following should be considered; there should be no arguments on the first line and further indentation should be used to clearly distinguish itself as a continuation line. `[E125, E131, E133]`
  * *The use of hanging indents is encouraged. `[E121, E126]` error codes are suppressed, so any amount of additional indentation (not just 4 spaces) is allowed as long as it is distinguishable from required indentation (such as from an `if` statement). For example:*
```python
# Set cell center functions for preset coordinates
if center_func_1 is None:
    if (coord == 'cartesian' or coord == 'minkowski' or coord == 'tilted'
           or coord == 'sinusoidal' or coord == 'kerr-schild'):
        def center_func_1(xm, xp): return 0.5 * (xm + xp)
```

<!-- The PEP 8 guide provides a mix of clear-cut rules and recommendations; the ambiguity of choice may be confusing for some users. This also explains why certain pycodestyle error codes are DEFAULT_IGNORE -->

For the most up to date Python style rules for Athena++, refer to the `setup.cfg` file in the root project directory.

### Checking Python style
<!-- mention/consider alternatives such as PyLint? -->
The [`flake8`](http://flake8.pycqa.org/en/latest/) Python module is the tool used for checking Python style in Athena++. This popular command line utility wraps the PyFlakes, `pycodestyle` (formerly `pep8`), and [McCabe](https://github.com/PyCQA/mccabe) (for checking code complexity) tools. It can be installed from `pip` and run directly from the root Athena++ directory without any options.  If `flake8` was installed to a directory in the shell's search `PATH`, it can be executed via:
```
flake8
```
Or, the module can be invoked as a script from the specific version of Python it was installed for:
```
python2.7 -m flake8 ...
```
The tool will check every `*.py` file in the subdirectories and output the line numbers, error codes, and description for every violation of a style rule. Analogous to `CPPLINT.cfg` for the [[C++ Style]] checker `cpplint.py`, the `flake8` linter refers to a [project configuration file](http://flake8.pycqa.org/en/latest/user/configuration.html), `setup.cfg`, in the root Athena++ directory. This version-controlled file stores the suppressed error codes and other linter options under the `[flake8]` section.
<!--  Note, Google's `cpplint.py` stored in `tst/style/` should be excluded from the linting process. -->

In a similar fashion as the C++ style linter, the [[Continuous Integration (CI)]] servers will execute `flake8` before running any [[Regression Testing]] and will report the build as failed whenever Python style conventions are broken. Athena++ developers should use this tool locally when modifying or adding Python source code; if the `flake8` output is not empty, correct any style errors before pushing your local changes.

The following comment ("no quality assurance") can be added to the end of any line to disable [a subset of `pycodestyle` error codes](http://pycodestyle.pycqa.org/en/latest/intro.html#error-codes) for that line (such as maximum line length):
```python
# noqa
```
Its use should be reserved for special cases, such as long strings whose formatting must not be modified (like the XDMF headers) or where dividing the strings might greatly diminish readability.

[`autopep8`](https://github.com/hhatto/autopep8) can be used on the command line to automatically fix most style inconsistencies. The tool can be installed from `pip`, and it uses the local project configuration defined in `setup.cfg`'s `[flake8]` section to decide which style conventions to apply, by default. To modify all the `*.py` files in-place from the root project directory, use:
```
autopep8 --in-place --recursive --aggressive --aggressive --global-config=./setup.cfg .
```
(The use of `--global-config` is redundant here). The converted files may need to be manually edited to fix any style violations that `autopep8` was unable to automatically correct.
