# Functions for interfacing with Athena++ during testing

# Modules
import numpy as np
import os
import struct
import subprocess

# Global variables
athena_rel_path = '../../'
saved_filenames = ['src/defs.hpp', 'Makefile']
saved_files = []

# Function for configuring Athena++
def configure(*args, **kwargs):
  current_dir = os.getcwd()
  os.chdir(athena_rel_path)
  try:
    configure_command = ['python', 'configure.py']
    for arg in args:
      configure_command.append('-{0}'.format(arg))
    for key,val in kwargs.iteritems():
      configure_command.append('--{0}={1}'.format(key,val))
    try:
      subprocess.check_call(configure_command)
    except subprocess.CalledProcessError as err:
      raise AthenaError('Return code {0} from command \'{1}\''\
          .format(err.returncode,' '.join(err.cmd)))
  finally:
    os.chdir(current_dir)

# Function for compiling Athena++
def make():
  current_dir = os.getcwd()
  os.chdir(athena_rel_path)
  try:
    exe_dir = 'EXE_DIR:={0}/bin/'.format(current_dir)
    obj_dir = 'OBJ_DIR:={0}/obj/'.format(current_dir)
    clean_command = ['make', 'clean', exe_dir, obj_dir]
    make_command = ['make', exe_dir, obj_dir]
    try:
      subprocess.check_call(clean_command)
      subprocess.check_call(make_command)
    except subprocess.CalledProcessError as err:
      raise AthenaError('Return code {0} from command \'{1}\''\
          .format(err.returncode,' '.join(err.cmd)))
  finally:
    os.chdir(current_dir)

# Function for running Athena++
def run(input_filename, arguments):
  current_dir = os.getcwd()
  os.chdir('bin')
  try:
    input_filename_full = '../' + athena_rel_path + 'inputs/' + input_filename
    run_command = ['./athena', '-i', input_filename_full]
    try:
      subprocess.check_call(run_command+arguments)
    except subprocess.CalledProcessError as err:
      raise AthenaError('Return code {0} from command \'{1}\''\
          .format(err.returncode,' '.join(err.cmd)))
  finally:
    os.chdir(current_dir)

# Function for reading Athena++ tabular data
def read_tab(filename, headings=None):
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

# Function for reading Athena++ VTK data
def read_vtk(filename):

  # Read raw data
  with open(filename, 'r') as data_file:
    raw_data = data_file.read()

  # Skip header
  current_index = 0
  current_char = raw_data[current_index]
  while current_char == '#':
    while current_char != '\n':
      current_index += 1
      current_char = raw_data[current_index]
    current_index += 1
    current_char = raw_data[current_index]

  # Function for skipping though the file
  def skip_string(expected_string):
    expected_string_len = len(expected_string)
    if raw_data[current_index:current_index+expected_string_len] != expected_string:
      raise AthenaError('File not formatted as expected')
    return current_index+expected_string_len

  # Read metadata
  current_index = skip_string('BINARY\nDATASET RECTILINEAR_GRID\nDIMENSIONS ')
  end_of_line_index = current_index + 1
  while raw_data[end_of_line_index] != '\n':
    end_of_line_index += 1
  face_dimensions = map(int, raw_data[current_index:end_of_line_index].split(' '))
  current_index = end_of_line_index + 1

  # Function for reading interface locations
  def read_faces(letter, num_faces):
    identifier_string = '{0}_COORDINATES {1} float\n'.format(letter,num_faces)
    begin_index = skip_string(identifier_string)
    format_string = '>' + 'f'*num_faces
    end_index = begin_index + 4*num_faces
    vals = struct.unpack(format_string, raw_data[begin_index:end_index])
    return vals,end_index+1

  # Read interface locations
  x_faces,current_index = read_faces('X', face_dimensions[0])
  y_faces,current_index = read_faces('Y', face_dimensions[1])
  z_faces,current_index = read_faces('Z', face_dimensions[2])

  # Prepare to read quantities defined on grid
  cell_dimensions = np.array([max(dim-1,1) for dim in face_dimensions])
  num_cells = cell_dimensions.prod()
  current_index = skip_string('CELL_DATA {0} \n'.format(num_cells))  # note trailing space
  data = {}

  # Function for reading scalar data
  def read_cell_scalars():
    begin_index = skip_string('SCALARS ')
    end_of_word_index = begin_index + 1
    while raw_data[end_of_word_index] != ' ':
      end_of_word_index += 1
    array_name = raw_data[begin_index:end_of_word_index]
    string_to_skip = 'SCALARS {0} float\nLOOKUP_TABLE default\n'.format(array_name)
    begin_index = skip_string(string_to_skip)
    format_string = '>' + 'f'*num_cells
    end_index = begin_index + 4*num_cells
    data[array_name] = struct.unpack(format_string, raw_data[begin_index:end_index])
    dimensions = tuple(cell_dimensions[::-1])
    data[array_name] = np.array(data[array_name]).reshape(dimensions)
    return end_index+1

  # Function for reading vector data
  def read_cell_vectors():
    begin_index = skip_string('VECTORS ')
    end_of_word_index = begin_index + 1
    while raw_data[end_of_word_index] != '\n':
      end_of_word_index += 1
    array_name = raw_data[begin_index:end_of_word_index]
    string_to_skip = 'VECTORS {0}\n'.format(array_name)
    array_name = array_name[:-6]  # remove ' float'
    begin_index = skip_string(string_to_skip)
    format_string = '>' + 'f'*num_cells*3
    end_index = begin_index + 4*num_cells*3
    data[array_name] = struct.unpack(format_string, raw_data[begin_index:end_index])
    dimensions = tuple(np.append(cell_dimensions[::-1],3))
    data[array_name] = np.array(data[array_name]).reshape(dimensions)
    return end_index+1

  # Read quantities defined on grid
  while current_index < len(raw_data):
    expected_string = 'SCALARS'
    expected_string_len = len(expected_string)
    if raw_data[current_index:current_index+expected_string_len] == expected_string:
      current_index = read_cell_scalars()
      continue
    expected_string = 'VECTORS'
    expected_string_len = len(expected_string)
    if raw_data[current_index:current_index+expected_string_len] == expected_string:
      current_index = read_cell_vectors()
      continue
    raise AthenaError('File not formatted as expected')
  return x_faces,y_faces,z_faces,data

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

# General exception class for these functions
class AthenaError(RuntimeError):
  pass
