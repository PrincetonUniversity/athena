# Functions for interfacing with Athena++ during testing

# Modules
import os
import sys
import subprocess

# Global variables
athena_rel_path = '../../'
saved_filenames = ['src/defs.hpp', 'Makefile']
saved_files = []
global_config_args = []
global_run_args = []
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
        try:
            subprocess.check_call(configure_command + global_config_args)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
    finally:
        os.chdir(current_dir)


# Function for compiling Athena++
def make(clean_first=True, obj_only=False):
    current_dir = os.getcwd()
    os.chdir(athena_rel_path)
    try:
        exe_dir = 'EXE_DIR:={0}/bin/'.format(current_dir)
        obj_dir = 'OBJ_DIR:={0}/obj/'.format(current_dir)
        clean_command = ['make', 'clean', exe_dir, obj_dir]
        if obj_only:
            make_command = ['make', '-j', 'objs', exe_dir, obj_dir]
        else:
            make_command = ['make', '-j', exe_dir, obj_dir]
        try:
            stdout_f = open(os.devnull, 'w') if global_silent else sys.stdout
            if clean_first:
                subprocess.check_call(clean_command, stdout=stdout_f)
            subprocess.check_call(make_command, stdout=stdout_f)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
        finally:
            if stdout_f is not sys.stdout:
                stdout_f.close()
    finally:
        os.chdir(current_dir)


# Functions for running Athena++
def run(input_filename, arguments):
    current_dir = os.getcwd()
    os.chdir('bin')
    try:
        input_filename_full = '../' + athena_rel_path + 'inputs/' + \
                              input_filename
        run_command = ['./athena', '-i', input_filename_full]
        try:
            subprocess.check_call(run_command + arguments + global_run_args)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
    finally:
        os.chdir(current_dir)


def restart(input_filename, arguments):
    current_dir = os.getcwd()
    os.chdir('bin')
    try:
        run_command = ['./athena', '-r', input_filename]
        try:
            subprocess.check_call(run_command + arguments)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
    finally:
        os.chdir(current_dir)


def mpirun(mpirun_cmd, mpirun_opts, nproc, input_filename, arguments):
    current_dir = os.getcwd()
    os.chdir('bin')
    try:
        input_filename_full = '../' + athena_rel_path + 'inputs/' + \
                              input_filename
        run_command = [mpirun_cmd] + mpirun_opts + ['-n', str(nproc), './athena', '-i',
                                                    input_filename_full]
        run_command = list(filter(None, run_command))  # remove any empty strings
        try:
            subprocess.check_call(run_command + arguments)
        except subprocess.CalledProcessError as err:
            raise AthenaError('Return code {0} from command \'{1}\''
                              .format(err.returncode, ' '.join(err.cmd)))
    finally:
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


# General exception class for these functions
class AthenaError(RuntimeError):
    pass
