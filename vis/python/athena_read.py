"""
Read Athena++ output data files.
"""

# Python Modules
import numpy as np
import h5py
import struct

#=======================================================================================

def tab(filename, headings=None, dimensions=1):
  """Read .tab files and return dict or array."""

  # Check for valid number of dimensions
  if dimensions != 1 and dimensions !=2 and dimensions != 3:
    raise AthenaError('Improper number of dimensions')

  # Read raw data
  with open(filename, 'r') as data_file:
    raw_data = data_file.readlines()

  # Organize data into array of numbers
  data_array = []
  first_line = True
  last_line_number = len(raw_data)
  line_number = 0
  for line in raw_data:
    line_number += 1
    if line.split()[0][0] == '#':  # comment line
      continue
    row = []
    col = 0
    for val in line.split():
      col += 1
      if col == 1:
        if first_line:
          i_min = int(val)
        if line_number == last_line_number:
          i_max = int(val)
      elif col == 3 and dimensions >= 2:
        if first_line:
          j_min = int(val)
        if line_number == last_line_number:
          j_max = int(val)
      elif col == 5 and dimensions == 3:
        if first_line:
          k_min = int(val)
        if line_number == last_line_number:
          j_max = int(val)
      else:
        row.append(float(val))
    first_line = False
    data_array.append(row)

  # Reshape array based on number of dimensions
  if dimensions == 1:
    j_min = j_max = 0
  if dimensions <= 2:
    k_min = k_max = 0
  array_shape = (k_max-k_min+1,j_max-j_min+1,i_max-i_min+1,len(row))
  data_array = np.reshape(data_array, array_shape)

  # Store separate variables as dictionary entries if headings given
  if headings is not None:
    data_dict = {}
    for n in range(len(headings)):
      data_dict[headings[n]] = data_array[:,:,:,n]
    return data_dict
  else:
    return data_array

#=======================================================================================

def vtk(filename):
  """Read .vtk files and return dict of arrays of data."""

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
    vals = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
    return vals,end_index+1

  # Read interface locations
  x_faces,current_index = read_faces('X', face_dimensions[0])
  y_faces,current_index = read_faces('Y', face_dimensions[1])
  z_faces,current_index = read_faces('Z', face_dimensions[2])

  # Prepare to read quantities defined on grid
  cell_dimensions = np.array([max(dim-1,1)
      for dim in face_dimensions])
  num_cells = cell_dimensions.prod()
  current_index = skip_string('CELL_DATA {0}\n'.format(num_cells))
  if raw_data[current_index:current_index+1] == '\n':
    current_index = skip_string('\n')  # extra newline inserted by join script
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

#=======================================================================================

def athdf(filename, data=None, quantities=None):
  """Read .athdf files and populate dict of arrays of data."""

  # Open file
  with h5py.File(filename, 'r') as f:

    # Create list of all quantities if none given
    if data is not None:
      quantities = data.values()
    elif quantities is None:
      quantities = f[u'MeshBlock0'].keys()
      quantities = [q for q in quantities \
          if q != u'x1f' and q != u'x2f' and q != u'x3f']

    # Get block count, dimensions, and sizes
    num_blocks = len(f.keys())
    dims = 0
    block_size = []
    coords = [u'x1f',u'x2f',u'x3f']
    for key in coords:
      if key in f[u'MeshBlock0'].keys():
        dims += 1
        block_size.append(len(f[u'MeshBlock0'][key][:]) - 1)
    coords = coords[:dims]

    # Order blocks
    edges = np.empty((num_blocks,dims))
    for block_num,block_name in zip(range(num_blocks),f.keys()):
      for dim,coord in zip(range(dims),coords):
        edges[block_num,dim] = f[block_name][coord][0]
    edges_unique = []
    for dim in range(dims):
      edges_unique.append(set(edges[:,dim]))
    indices = np.empty((num_blocks,3,2), dtype=int)
    for block_num in range(num_blocks):
      for dim in range(dims):
        num_prior = sum(edge < edges[block_num,dim] for edge in edges_unique[dim])
        indices[block_num,dim,0] = num_prior * block_size[dim]
        indices[block_num,dim,1] = (num_prior+1) * block_size[dim]
      for dim in range(dims,3):
        indices[block_num,dim,0] = 0
        indices[block_num,dim,1] = 1

    # Prepare arrays if needed
    if data is None:
      nx1 = block_size[0] * len(edges_unique[0])
      nx2 = block_size[1] * len(edges_unique[1]) if dims >= 2 else 1
      nx3 = block_size[2] * len(edges_unique[2]) if dims >= 3 else 1
      data = {}
      data[u'x1f'] = np.empty(nx1+1)
      if dims >= 2:
        data[u'x2f'] = np.empty(nx2+1)
      if dims >= 3:
        data[u'x3f'] = np.empty(nx3+1)
      for q in quantities:
        data[q] = np.empty((nx3,nx2,nx1))

    # Read interface data
    for n,block_name in zip(range(num_blocks),f.keys()):
      for dim,coord in zip(range(dims),coords):
        need_interfaces = True
        for dim_other in range(dims):
          if dim_other == dim:
            continue
          if indices[n,dim_other,0] != 0:
            need_interfaces = False
        if not need_interfaces:
          continue
        data[coord][indices[n,dim,0]:indices[n,dim,1]] = f[block_name][coord][:-1]
        if indices[n,dim,1] == block_size[dim] * len(edges_unique[dim]):
          data[coord][indices[n,dim,1]] = f[block_name][coord][-1]

    # Read value data
    for n,block_name in zip(range(num_blocks),f.keys()):
      for q in quantities:
        data[q][indices[n,2,0]:indices[n,2,1],indices[n,1,0]:indices[n,1,1],\
            indices[n,0,0]:indices[n,0,1]] = f[block_name][q][:]
  return data

#=======================================================================================

class AthenaError(RuntimeError):
  """General exception class for Athena++ read functions."""
  pass
