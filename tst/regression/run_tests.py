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

# Python modules
from __future__ import print_function
import argparse
import os
from collections import OrderedDict
from importlib import reload
import logging
import logging.config
from pkgutil import iter_modules
from timeit import default_timer as timer

# Prevent generation of .pyc files
# This should be set before importing any user modules
import sys
sys.dont_write_bytecode = True

# Athena++ modules
import scripts.utils.athena as athena  # noqa

logger = logging.getLogger('athena')


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
    coverage_cmd = kwargs.pop('coverage')
    # Get args to pass to scripts.utils.athena as list of strings
    athena_config_args = kwargs.pop('config')
    athena_run_args = kwargs.pop('run')

    if len(tests) == 0:  # run all tests
        for _, directory, ispkg in iter_modules(path=['scripts/tests']):
            if ispkg:
                dir_test_names = [name for _, name, _ in
                                  iter_modules(path=['scripts/tests/'
                                                     + directory],
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
                                  iter_modules(path=['scripts/tests/'
                                                     + test],
                                               prefix=test + '.')]
                test_names.extend(dir_test_names)
    # Remove duplicate test entries while preserving the original order
    test_names = list(OrderedDict.fromkeys(test_names))

    # Run tests
    current_dir = os.getcwd()
    test_times = []
    test_results = []
    test_errors = []
    try:
        # Check that required modules are installed for all test dependencies
        deps_installed = True
        for name in test_names:
            try:
                name_full = 'scripts.tests.' + name
                module = __import__(name_full, globals(), locals(),
                                    fromlist=['prepare', 'run', 'analyze'])
            except ImportError as e:
                if sys.version_info >= (3, 6, 0):  # ModuleNotFoundError subclass
                    missing_module = e.name
                else:
                    missing_module = e.message.split(' ')[-1]
                logger.warning('Unable to import "{:}".'.format(missing_module))
                deps_installed = False
        if not deps_installed:
            logger.warning('##########################################################')
            logger.warning('# WARNING! Not all required Python mdules are available. #')
            logger.warning('##########################################################')
        # Run each test
        for name in test_names:
            t0 = timer()
            try:
                name_full = 'scripts.tests.' + name
                module = __import__(name_full, globals(), locals(),
                                    fromlist=['prepare', 'run', 'analyze'])
                reload(module)
                os.system('rm -rf {0}/bin'.format(current_dir))
                os.system('rm -rf {0}/obj'.format(current_dir))

                # Change formatting of full single test name, e.g. gr.compile_kerr-schild:
                # (Lcov test names may only contain letters, numbers, and '_')
                lcov_test_name = name.replace('.', '_').replace('-', '_')
                # (optional) build test-dependent code coverage command
                if coverage_cmd is not None:
                    module.athena.global_test_name = lcov_test_name
                    module.athena.global_coverage_cmd = coverage_cmd

                # insert arguments for athena.run() and athena.configure()
                # by changing global values through module
                module.athena.global_config_args = athena_config_args
                module.athena.global_run_args = athena_run_args
                module.athena.global_silent = silent_opt

                try:
                    module.prepare(**kwargs)
                except Exception:
                    logger.error("Exception occurred", exc_info=True)
                    test_errors.append('prepare()')
                    raise TestError(name_full.replace('.', '/') + '.py')
                try:
                    run_ret = module.run(mpirun_cmd=mpirun_cmd, mpirun_opts=mpirun_opts)
                except Exception:
                    logger.error("Exception occurred", exc_info=True)
                    test_errors.append('run()')
                    raise TestError(name_full.replace('.', '/') + '.py')
                else:
                    # (optional) if test_name.run() completes w/o error, perform code
                    # coverage analysis after the final athena.run() call, producing
                    # test_name.info (no suffix). Function only actually executes on cmd
                    # line if --coverage=CMD is passed to this run_tests.py script:
                    if run_ret is None:
                        athena.analyze_code_coverage(lcov_test_name, '')
                try:
                    result = module.analyze()
                except Exception:
                    logger.error("Exception occurred", exc_info=True)
                    test_errors.append('analyze()')
                    raise TestError(name_full.replace('.', '/') + '.py')
            except TestError as err:
                test_results.append(False)
                logger.error('---> Error in ' + str(err))
                # do not measure runtime for failed/incomplete tests
                test_times.append(None)
            else:
                test_times.append(timer() - t0)
                msg = 'Test {0} took {1:.3g} seconds to complete.'
                msg = msg.format(name, test_times[-1])
                logging.getLogger('athena.tests.' + name).debug(msg)
                test_results.append(result)
                test_errors.append(None)
            finally:
                os.system('rm -rf {0}/bin'.format(current_dir))
                os.system('rm -rf {0}/obj'.format(current_dir))
            # For CI, print after every individual test has finished
            logger.info('{} test: prepare(), run(), analyze() finished'.format(name))

    # Restore any previously-existing files once ALL runs are complete
    finally:
        athena.restore_files()

    # Report test results
    logger.info('\nResults:')
    for name, result, error, time in zip(test_names, test_results, test_errors,
                                         test_times):
        result_string = 'passed' if result else 'failed'
        error_string = ' -- unexpected failure in {0} stage'.format(error) \
                       if error is not None else '; time elapsed: {0:.3g} s'.format(time)
        logger.info('    {0}: {1}{2}'.format(name, result_string, error_string))
    logger.info('')
    num_tests = len(test_results)
    num_passed = test_results.count(True)
    test_string = 'test' if num_tests == 1 else 'tests'
    logger.info('Summary: {0} out of {1} {2} passed\n'.format(num_passed, num_tests,
                                                              test_string))
    # For CI calling scripts, explicitly raise error if not all tests passed
    if num_passed == num_tests:
        return 0
    else:
        raise TestError()


# Exception for unexpected behavior by individual tests
class TestError(RuntimeError):
    pass


class CriticalExceptionFilter(logging.Filter):
    """Filter out critical exceptions"""
    def filter(self, record):
        return not record.exc_info or record.levelno != logging.CRITICAL


class MakeFilter(logging.Filter):
    """Filter out make output"""
    def filter(self, record):
        return 'athena.make' != record.name[:len('athena.make')]


def log_init(args):
    """Initialize log"""
    kwargs = vars(args)
    logging.basicConfig(level=0)  # setting this to zero gives output control to handler
    logger.propagate = False  # don't use default handler
    c_handler = logging.StreamHandler()  # console/terminal handler
    c_handler.setLevel(kwargs.pop('loglevel'))
    c_handler.addFilter(CriticalExceptionFilter())  # let stderr print errors to screen
    c_handler.setFormatter(logging.Formatter('%(message)s'))  # only show the message
    if kwargs.pop('verbose'):
        c_handler.setLevel(0)
        c_handler.setFormatter(logging.Formatter('%(levelname)s:%(name)s: %(message)s'))
    if kwargs.pop('hide_make'):
        c_handler.addFilter(MakeFilter())
    logger.addHandler(c_handler)
    # setup logfile
    log_fn = kwargs.pop('logfile')
    if log_fn:
        f_handler = logging.FileHandler(log_fn)
        f_handler.setLevel(0)  # log everything
        f_format = \
            logging.Formatter('%(asctime)s|%(levelname)s:%(name)s: %(message)s')
        f_handler.setFormatter(f_format)
        logger.addHandler(f_handler)
    # setup runtime diagnostics file
    if kwargs.pop('diagnostics_file'):
        rtd_handler = logging.FileHandler('diagnostics.log')
        rtd_handler.setLevel(0)  # log everything
        rtd_fmt = logging.Formatter('%(asctime)s|%(levelname)s:%(name)s: %(message)s')
        rtd_handler.setFormatter(rtd_fmt)
        record = logger.makeRecord('athena.tests', logging.INFO, None, None,
                                   'Starting new tests', None, None)
        rtd_handler.emit(record)
        for log in ['athena.run', 'athena.tests']:
            logging.getLogger(log).addHandler(rtd_handler)

    logger.debug('Starting Athena++ regression tests')


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
                        # 2x MPI, Slurm, PBS/Torque, LSF, Cray ALPS
                        choices=['mpirun', 'mpiexec', 'srun', 'qsub', 'lsrun', 'aprun'],
                        help='change MPI run wrapper command (e.g. for job schedulers)')

    parser.add_argument('--mpirun_opts',
                        default=[],
                        action='append',
                        help='add options to mpirun command')

    parser.add_argument('--silent', '-s',
                        dest='global_silent',
                        default=False,
                        action='store_true',
                        help='redirect stdout of make to devnull (for CI logs)')

    parser.add_argument("--config", "-c",
                        default=[],
                        action='append',
                        help='arguments to pass to athena.configure')

    parser.add_argument("--run", "-r",
                        default=[],
                        action='append',
                        help='arguments to pass to athena.run')

    parser.add_argument("--coverage", "-cov",
                        type=str,
                        default=None,
                        help=('code coverage command to run after a successful test;'
                              ' automatically passes -coverage to configure.py.'
                              ' Currently, assumes that Lcov is being used and appends '
                              ' -t and -o options w/ reformatted test name to COVERAGE.'))
    parser.add_argument('-d', '--debug',
                        help="print debugging information",
                        action="store_const",
                        dest="loglevel",
                        const=logging.DEBUG,
                        default=logging.INFO)
    parser.add_argument('-v', '--verbose',
                        default=False,
                        action='store_true',
                        help="print all output, timestamps and logging information")
    parser.add_argument('--logfile',
                        type=str,
                        default=None,
                        help='set filename of logfile')
    parser.add_argument('--diagnostics_file',
                        default=False,
                        action='store_true',
                        help='Write runtime diagnostics to runtime.log')
    parser.add_argument('--hide_make',
                        default=False,
                        action='store_true',
                        help='Hide output from make command')

    args = parser.parse_args()
    log_init(args)

    try:
        logger.debug('args: ' + str(vars(args)))
        main(**vars(args))
    except Exception:
        logger.critical('', exc_info=True)
        raise
