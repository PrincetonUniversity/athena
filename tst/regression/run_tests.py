#!/usr/bin/env python

"""
Regression test script.

Usage: From this directory, call this script with python:
      python run_tests.py

Notes:
  - Requires Python 2.7+. (compliant with Python 3)
  - This file should not be modified when adding new scripts.
  - To add a new script, create a new .py file in scripts/tests/ subdirectory.
  - See scripts/tests/example.py for an example.
    - Example can be forced to run, but does not run by default in full test.
  - For more information, check online regression test documentation.
"""
from __future__ import print_function

# Python modules
import argparse
import os
from pkgutil import iter_modules
import traceback

# Prevent generation of .pyc files
# This should be set before importing any user modules
import sys
sys.dont_write_bytecode = True

# Athena++ modules
import scripts.utils.athena as athena # noqa


# Main function
def main(**kwargs):

    # Save existing files
    athena.save_files()

    # Make list of tests to run
    tests = kwargs.pop('tests')
    test_names = []
    # Get MPI run syntax and other flags
    mpirun_cmd = kwargs.pop('mpirun')
    mpirun_opts = kwargs.pop('mpirun_opts')
    silent_opt = kwargs.pop('global_silent')
    # Get args to pass to scripts.utils.athena as list of strings
    athena_config_args = kwargs.pop('config')
    athena_run_args = kwargs.pop('run')

    if len(tests) == 0:  # run all tests
        for _, directory, ispkg in iter_modules(path=['scripts/tests']):
            if ispkg:
                dir_test_names = [name for _, name, _ in
                                  iter_modules(path=['scripts/tests/' +
                                                     directory],
                                               prefix=directory + '.')]
                test_names.extend(dir_test_names)
    else:  # run selected tests
        for test in tests:
            if test == 'example':
                test_names.append(test)
                continue
            if test[-1] == '/':
                test = test[:-1]  # remove trailing slash
            if '/' in test:  # specific test specified
                test_names.append(test.replace('/', '.'))
            else:  # test suite specified
                dir_test_names = [name for _, name, _ in
                                  iter_modules(path=['scripts/tests/' +
                                                     test],
                                               prefix=test + '.')]
                test_names.extend(dir_test_names)
    test_names = list(set(test_names))

    # Run tests
    current_dir = os.getcwd()
    test_results = []
    test_errors = []
    try:
        for name in test_names:
            try:
                name_full = 'scripts.tests.' + name
                module = __import__(name_full, globals(), locals(),
                                    fromlist=['prepare', 'run', 'analyze'])
                os.system('rm -rf {0}/bin'.format(current_dir))

                # insert arguments to athena.run and athena.configure
                # by changing global values through module
                module.athena.global_config_args = athena_config_args
                module.athena.global_run_args = athena_run_args
                module.athena.global_silent = silent_opt

                try:
                    module.prepare(**kwargs)
                except Exception:
                    traceback.print_exc()
                    test_errors.append('prepare()')
                    raise TestError(name_full.replace('.', '/') + '.py')
                try:
                    module.run(mpirun_cmd=mpirun_cmd, mpirun_opts=mpirun_opts)
                except Exception:
                    traceback.print_exc()
                    test_errors.append('run()')
                    raise TestError(name_full.replace('.', '/') + '.py')
                try:
                    result = module.analyze()
                except Exception:
                    traceback.print_exc()
                    test_errors.append('analyze()')
                    raise TestError(name_full.replace('.', '/') + '.py')
            except TestError as err:
                test_results.append(False)
                print('---> Error in ' + str(err))
            else:
                test_results.append(result)
                test_errors.append(None)
            finally:
                os.system('rm -rf {0}/bin'.format(current_dir))
            # For CI, print after every individual test has finished
            print('{} test: prepare(), run(), analyze() finished'.format(name))

    # Restore any previously-existing files once runs are complete
    finally:
        athena.restore_files()

    # Report test results
    print('\nResults:')
    for name, result, error in zip(test_names, test_results, test_errors):
        result_string = 'passed' if result else 'failed'
        error_string = ' -- unexpected failure in {0} stage'.format(error) \
                       if error is not None else ''
        print('    {0}: {1}{2}'.format(name, result_string, error_string))
    print('')
    num_tests = len(test_results)
    num_passed = test_results.count(True)
    test_string = 'test' if num_tests == 1 else 'tests'
    print('Summary: {0} out of {1} {2} passed\n'.format(num_passed, num_tests,
                                                        test_string))
    # For CI calling scripts, explicitly raise error if not all tests passed
    if num_passed == num_tests:
        return 0
    else:
        raise TestError()


# Exception for unexpected behavior by individual tests
class TestError(RuntimeError):
    pass


# Execute main function
if __name__ == '__main__':
    help_msg = ('names of tests to run, relative to scripts/tests/,'
                'excluding .py')
    parser = argparse.ArgumentParser()
    parser.add_argument('tests',
                        type=str,
                        default=None,
                        nargs='*',
                        help=help_msg)

    parser.add_argument('--mpirun',
                        default='mpirun',
                        choices=['mpirun', 'srun', 'mpiexec'],
                        help='select MPI run command')

    parser.add_argument('--mpirun_opts',
                        default=[],
                        action='append',
                        help='select MPI run command')

    parser.add_argument('--silent', '-s',
                        dest='global_silent',
                        default=False,
                        action='store_true',
                        help='redirect stdout of make to devnull (for CI logs)')

    parser.add_argument("--config", "-c",
                        default=[],
                        action='append',
                        help=('arguments to pass to athena.configure'))

    parser.add_argument("--run", "-r",
                        default=[],
                        action='append',
                        help=('arguments to pass to athena.run'))

    args = parser.parse_args()
    main(**vars(args))
