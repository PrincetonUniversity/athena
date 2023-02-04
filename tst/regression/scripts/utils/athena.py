# Functions for interfacing with Athena++ during testing

# Modules
import logging
import os
import sys
import subprocess
from timeit import default_timer as timer
from .log_pipe import LogPipe

# Global variables
athena_rel_path = '../../'
saved_filenames = ['src/defs.hpp', 'Makefile']
saved_files = []
global_config_args = []
global_run_args = []
global_test_name = None
global_coverage_cmd = None
global_silent = False


# Function for configuring Athena++
def configure(*args, **kwargs):
    current_dir = os.getcwd()
    os.chdir(athena_rel_path)
    try:
        pybin = sys.executable if sys.executable is not None else 'python'
        configure_command = [pybin, 'configure.py']
        for arg in args:
            configure_command.append('-{0}'.format(arg))
        for key, val in kwargs.items():
            if val:
                configure_command.append('--{0}={1}'.format(key, val))
        if global_coverage_cmd is not None:
            configure_command.append('-coverage')
        out_log = LogPipe('athena.configure', logging.INFO)
        err_log = LogPipe('athena.configure', logging.ERROR)
        cmd = ' '.join(['Executing: '] + configure_command + global_config_args)
        logging.getLogger('athena.configure').debug(cmd)
        try:
            subprocess.check_call(configure_command + global_config_args,
                                  stdout=out_log, stderr=err_log)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
        finally:
            out_log.close()
            err_log.close()
    finally:
        os.chdir(current_dir)


# Function for compiling Athena++
def make(clean_first=True, obj_only=False):
    current_dir = os.getcwd()
    os.chdir(athena_rel_path)
    logger = logging.getLogger('athena.make')
    out_log = open(os.devnull, 'w') if global_silent else LogPipe('athena.make',
                                                                  logging.INFO)
    try:
        exe_dir = 'EXE_DIR:={0}/bin/'.format(current_dir)
        obj_dir = 'OBJ_DIR:={0}/obj/'.format(current_dir)
        clean_command = ['make', 'clean', exe_dir, obj_dir]
        if obj_only:
            # used in pgen_compile.py to save expensive linking time for Intel Compiler:
            make_command = ['make', '-j8', 'objs']
        else:
            # disable parallel GNU Make execution for Lcov (issues w/ Jenkins filesystem)
            if (global_coverage_cmd is not None):
                make_command = ['make']
            else:
                make_command = ['make', '-j8']
        make_command += [exe_dir, obj_dir]
        try:
            if clean_first:
                logger.debug('Executing: ' + ' '.join(clean_command))
                subprocess.check_call(clean_command, stdout=out_log)
            logger.debug('Executing: ' + ' '.join(make_command))
            t0 = timer()
            subprocess.check_call(make_command, stdout=out_log)
            logger.debug('Compilation took {0:.3g} seconds.'.format(timer() - t0))
        except subprocess.CalledProcessError as err:
            logger.error("Something bad happened", exc_info=True)
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
    finally:
        out_log.close()
        os.chdir(current_dir)


# Functions for running Athena++
def run(input_filename, arguments, lcov_test_suffix=None):
    current_dir = os.getcwd()
    os.chdir('bin')
    out_log = LogPipe('athena.run', logging.INFO)
    try:
        input_filename_full = '../' + athena_rel_path + 'inputs/' + \
                              input_filename
        run_command = ['./athena', '-i', input_filename_full]
        try:
            cmd = run_command + arguments + global_run_args
            logging.getLogger('athena.run').debug('Executing: ' + ' '.join(cmd))
            subprocess.check_call(cmd, stdout=out_log)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
        else:
            os.chdir(current_dir)
            # (optional) if execution completes without error, and a lcov_test_suffix is
            # explicitly passed, process Lcov tracefile immediately after run_command
            analyze_code_coverage(global_test_name, lcov_test_suffix)
    finally:
        out_log.close()
        os.chdir(current_dir)


def restart(input_filename, arguments):
    current_dir = os.getcwd()
    os.chdir('bin')
    out_log = LogPipe('athena.make', logging.INFO)
    try:
        run_command = ['./athena', '-r', input_filename]
        try:
            cmd = run_command + arguments
            logging.getLogger('athena.run').debug('Executing (restart): ' + ' '.join(cmd))
            subprocess.check_call(cmd, stdout=out_log)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
    finally:
        out_log.close()
        os.chdir(current_dir)


def mpirun(mpirun_cmd, mpirun_opts, nproc, input_filename, arguments,
           lcov_test_suffix=None):
    current_dir = os.getcwd()
    os.chdir('bin')
    out_log = LogPipe('athena.run', logging.INFO)
    try:
        input_filename_full = '../' + athena_rel_path + 'inputs/' + \
                              input_filename
        run_command = [mpirun_cmd] + mpirun_opts + ['-n', str(nproc), './athena', '-i',
                                                    input_filename_full]
        run_command = list(filter(None, run_command))  # remove any empty strings
        try:
            cmd = run_command + arguments + global_run_args
            logging.getLogger('athena.run').debug('Executing (mpirun): ' + ' '.join(cmd))
            subprocess.check_call(cmd, stdout=out_log)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
        else:
            os.chdir(current_dir)
            # (optional) if execution completes without error, and a lcov_test_suffix is
            # explicitly passed, process Lcov tracefile immediately after run_command
            analyze_code_coverage(global_test_name, lcov_test_suffix)
    finally:
        out_log.close()
        os.chdir(current_dir)


# Function for saving configure-generated files that may already exist
def save_files():
    global saved_files
    for filename in saved_filenames:
        rel_path_to_file = athena_rel_path + filename
        if os.path.isfile(rel_path_to_file):
            with open(rel_path_to_file, 'r') as current_file:
                saved_files.append(current_file.read())
        else:
            saved_files.append(None)


# Function for restoring configure-generated files that previously existed
def restore_files():
    for filename, saved_file in zip(saved_filenames, saved_files):
        rel_path_to_file = athena_rel_path + filename
        if saved_file is None:
            os.system('rm -f ' + rel_path_to_file)
        else:
            with open(rel_path_to_file, 'w') as current_file:
                current_file.write(saved_file)


# Function for analyzing code coverage after athena.run() or athena.mpirun()
def analyze_code_coverage(test_name, lcov_test_suffix=None):
    # Only run Lcov if a string is passed to optional lcov_test_suffix argument (most
    # regression tests only need to run Lcov ONCE after test.analyze() is complete).

    # Regression tests with multiple athena.make() calls and binaries require that Lcov is
    # executed in test.run() after intermediate athena.mpi/run() calls, before the test is
    # complete. Test author must ensure that the obj/ directory contains the
    # correct/matching files when athena.run() is called with lcov_test_suffix=string
    if (lcov_test_suffix is not None and global_coverage_cmd is not None):
        # Empty string --> use base test name for output file
        if lcov_test_suffix == '':
            lcov_test_name = global_test_name
        # Append nonempty string suffix to base test name with an underscore
        else:
            lcov_test_name = '_'.join([global_test_name, lcov_test_suffix])
        # For now, assumes Lcov flags for adding test-dependent info (name, output file):
        test_lcov_cmd = (
            global_coverage_cmd
            + ' --test-name {0} -output-file {0}.info'.format(lcov_test_name)
        )
        # is this usage of os.system() safe?
        os.system(test_lcov_cmd)


# General exception class for these functions
class AthenaError(RuntimeError):
    pass
