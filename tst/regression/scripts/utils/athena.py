# Functions for interfacing with Athena++ during testing

# Modules
import logging
import os
import subprocess
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
        configure_command = ['python', 'configure.py']
        for arg in args:
            configure_command.append('-{0}'.format(arg))
        for key, val in kwargs.items():
            if val:
                configure_command.append('--{0}={1}'.format(key, val))
        if global_coverage_cmd is not None:
            configure_command.append('-coverage')
        stdout_f = LogPipe('athena.configure', logging.INFO)
        stderr_f = LogPipe('athena.configure', logging.ERROR)
        try:
            subprocess.check_call(configure_command + global_config_args,
                                  stdout=stdout_f, stderr=stderr_f)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
        finally:
            stdout_f.close()
            stderr_f.close()
    finally:
        os.chdir(current_dir)


# Function for compiling Athena++
def make(clean_first=True, obj_only=False):
    current_dir = os.getcwd()
    os.chdir(athena_rel_path)
    stdout_f = open(os.devnull, 'w') if global_silent else LogPipe('athena.make',
                                                                   logging.INFO)
    try:
        exe_dir = 'EXE_DIR:={0}/bin/'.format(current_dir)
        obj_dir = 'OBJ_DIR:={0}/obj/'.format(current_dir)
        clean_command = ['make', 'clean', exe_dir, obj_dir]
        if obj_only:
            # used in pgen_compile.py to save expensive linking time for Intel Compiler:
            make_command = ['make', '-j8', 'objs', exe_dir, obj_dir]
        else:
            # KGF: temporarily disable parallel compilation to debug Jenkins+Gcov issues
            # (could disable only if --coverage was used)
            make_command = ['make',  # '-j8',
                            exe_dir, obj_dir]
        try:
            if clean_first:
                subprocess.check_call(clean_command, stdout=stdout_f)
            # KGF: temporarily ignore "--silent" option for devnull redirection
            # (what about stderr?) stdout, stderr default behavior:
            # "with the default settings of None, no redirection will occur; the child's
            # file handles will be inherited from the parent."
            subprocess.check_call(make_command, stdout=stdout_f)
        except subprocess.CalledProcessError as err:
            logging.getLogger().error("Something bad happened", exc_info=True)
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
    finally:
        stdout_f.close()
        os.chdir(current_dir)


# Functions for running Athena++
def run(input_filename, arguments, lcov_test_suffix=None):
    current_dir = os.getcwd()
    os.chdir('bin')
    stdout_f = LogPipe('athena.run', logging.INFO)
    try:
        input_filename_full = '../' + athena_rel_path + 'inputs/' + \
                              input_filename
        run_command = ['./athena', '-i', input_filename_full]
        try:
            subprocess.check_call(run_command + arguments + global_run_args,
                                  stdout=stdout_f)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
        else:
            os.chdir(current_dir)
            # (optional) if execution completes without error, and a lcov_test_suffix is
            # explicitly passed, process Lcov tracefile immediately after run_command
            analyze_code_coverage(global_test_name, lcov_test_suffix)
    finally:
        stdout_f.close()
        os.chdir(current_dir)


def restart(input_filename, arguments):
    current_dir = os.getcwd()
    os.chdir('bin')
    stdout_f = LogPipe('athena.make', logging.INFO)
    try:
        run_command = ['./athena', '-r', input_filename]
        try:
            subprocess.check_call(run_command + arguments)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
    finally:
        stdout_f.close()
        os.chdir(current_dir)


def mpirun(mpirun_cmd, mpirun_opts, nproc, input_filename, arguments,
           lcov_test_suffix=None):
    current_dir = os.getcwd()
    os.chdir('bin')
    stdout_f = LogPipe('athena.run', logging.INFO)
    try:
        input_filename_full = '../' + athena_rel_path + 'inputs/' + \
                              input_filename
        run_command = [mpirun_cmd] + mpirun_opts + ['-n', str(nproc), './athena', '-i',
                                                    input_filename_full]
        run_command = list(filter(None, run_command))  # remove any empty strings
        try:
            subprocess.check_call(run_command + arguments, stdout=stdout_f)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
        else:
            os.chdir(current_dir)
            # (optional) if execution completes without error, and a lcov_test_suffix is
            # explicitly passed, process Lcov tracefile immediately after run_command
            analyze_code_coverage(global_test_name, lcov_test_suffix)
    finally:
        stdout_f.close()
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
