# Functions for interfacing with Athena++ during testing

# Modules
import numpy as np
import os

# Global variables
athena_rel_path = '../../'
saved_filenames = ['src/defs.hpp', 'Makefile']
saved_files = []

# Function for configuring Athena++
def configure(*args, **kwargs):
  current_dir = os.getcwd()
  os.chdir(athena_rel_path)
  try:
    configure_string = 'python configure.py'
    for arg in args:
      configure_string += ' -{0}'.format(arg)
    for key,val in kwargs.iteritems():
      configure_string += ' --{0}={1}'.format(key,val)
    os.system(configure_string)
  finally:
    os.chdir(current_dir)

# Function for compiling Athena++
def make():
  current_dir = os.getcwd()
  os.chdir(athena_rel_path)
  try:
    exe_dir = current_dir + '/bin/'
    obj_dir = current_dir + '/obj/'
    os.system('make clean EXE_DIR:={0} OBJ_DIR:={1}'.format(exe_dir,obj_dir))
    os.system('make EXE_DIR:={0} OBJ_DIR:={1}'.format(exe_dir,obj_dir))
  finally:
    os.chdir(current_dir)

# Function for running Athena++
def run(input_filename, *args):
  current_dir = os.getcwd()
  os.chdir('bin')
  try:
    input_filename_full = '../' + athena_rel_path + 'inputs/' + input_filename
    run_string = './athena -i {0}'.format(input_filename_full)
    for arg in args:
      run_string += ' {0}'.format(arg)
    os.system(run_string)
  finally:
    os.chdir(current_dir)

# Function for reading Athena++ data
def read(filename, headings=None):
  with open(filename, 'r') as data_file:
    raw_data = data_file.readlines()
  data_array = []
  for line in raw_data:
    if line.split()[0][0] == '#':
      continue
    row = []
    for val in line.split()[1:]:
      row.append(float(val))
    data_array.append(row)
  data_array = np.array(data_array)
  if headings is not None:
    data_dict = {}
    for i in range(len(headings)):
      data_dict[headings[i]] = data_array[:,i]
    return data_dict
  else:
    return data_array

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
  for filename,saved_file in zip(saved_filenames,saved_files):
    rel_path_to_file = athena_rel_path + filename
    if saved_file is None:
      os.system('rm -f ' + rel_path_to_file)
    else:
      with open(rel_path_to_file, 'w') as current_file:
        current_file.write(saved_file)
