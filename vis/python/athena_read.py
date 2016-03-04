"""
Read Athena++ output data files.
"""

# Python modules
import numpy as np

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

  # Python module
  import struct

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

def athdf(filename, data=None, quantities=None, level=0, vol_func=None):
  """Read .athdf files and populate dict of arrays of data."""

  # Python module for reading hdf5 files
  import h5py

  # Set naive volume function assuming Cartesian coordinates
  if vol_func is None:
    vol_func = lambda x1m,x1p,x2m,x2p,x3m,x3p: (x1p-x1m) * (x2p-x2m) * (x3p-x3m)

  # Open file
  with h5py.File(filename, 'r') as f:

    # Create list of all quantities if none given
    if data is not None:
      quantities = data.values()
    elif quantities is None:
      quantities = f[u'MeshBlock0'].keys()
    quantities = [str(q) for q in quantities \
        if q != 'x1f' and q != 'x2f' and q != 'x3f']

    # Extract size information
    block_size = f.attrs['MeshBlockSize']
    root_grid_size = f.attrs['RootGridSize']
    nx1 = root_grid_size[0] * 2**level
    nx2 = root_grid_size[1] * 2**level
    nx3 = root_grid_size[2] * 2**level
    lx1 = nx1 / block_size[0]
    lx2 = nx2 / block_size[1]
    lx3 = nx3 / block_size[2]
    if nx3 > 1:
      dim = 3
    elif nx2 > 1:
      dim = 2
    else:
      dim = 1

    # Prepare arrays
    if data is not None:
      for q in quantities:
        data[q].fill(0.0)
    else:
      data = {}
      for q in quantities:
        data[q] = np.zeros((nx3,nx2,nx1))
    data['x1f'] = np.empty(nx1+1)
    data['x2f'] = np.empty(nx2+1)
    data['x3f'] = np.empty(nx3+1)
    restricted_data = np.zeros((lx3,lx2,lx1), dtype=bool)

    # Go through blocks in data file
    for block in f.itervalues():

      # Extract location information
      block_level = block.attrs['Level'][0]
      block_location = block.attrs['LogicalLocation']

      # Prolongate coarse data and copy same-level data
      if block_level <= level:
        s = 2**(level-block_level)
        il = block_location[0] * s * block_size[0]
        jl = block_location[1] * s * block_size[1]
        kl = block_location[2] * s * block_size[2]
        iu = il + s * block_size[0]
        ju = jl + s * block_size[1]
        ku = kl + s * block_size[2]
        for q in quantities:
          for ko in range(s):
            for jo in range(s):
              for io in range(s):
                data[q][kl+ko:ku+ko:s,jl+jo:ju+jo:s,il+io:iu+io:s] = block[q][:]

      # Restrict fine data
      else:
        s = 2**(block_level-level)
        ir_vals = np.arange(block_size[0])
        jr_vals = np.arange(block_size[1])
        kr_vals = np.arange(block_size[2])
        i_vals = (ir_vals + block_location[0] * block_size[0]) / s
        j_vals = (jr_vals + block_location[1] * block_size[1]) / s
        k_vals = (kr_vals + block_location[2] * block_size[2]) / s
        for k,kr in zip(k_vals,kr_vals):
          x3m = block['x3f'][kr]
          x3p = block['x3f'][kr+1]
          for j,jr in zip(j_vals,jr_vals):
            x2m = block['x2f'][jr]
            x2p = block['x2f'][jr+1]
            for i,ir in zip(i_vals,ir_vals):
              x1m = block['x1f'][ir]
              x1p = block['x1f'][ir+1]
              vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
              for q in quantities:
                data[q][k,j,i] += block[q][:][kr,jr,ir] * vol
        loc1 = block_location[0] / s
        loc2 = block_location[1] / s
        loc3 = block_location[2] / s
        restricted_data[loc3,loc2,loc1] = True

    # Record interface locations
    for d,nx in zip(np.arange(dim)+1,(nx1,nx2,nx3)[:dim]):
      xf_string = 'x' + str(d) + 'f'
      x_0 = f['MeshBlock0'][xf_string][:][0]
      x_1 = f['MeshBlock0'][xf_string][:][1]
      x_2 = f['MeshBlock0'][xf_string][:][2]
      ratio_block0 = (x_2-x_1) / (x_1-x_0)
      level_block0 = f['MeshBlock0'].attrs['Level'][0]
      ratio = ratio_block0 ** (1.0 / 2**(level-level_block0))
      ratio_powers = ratio ** np.arange(nx)
      ratio_powers_sum = np.cumsum(ratio_powers)
      data[xf_string][0] = x_0
      data[xf_string][1:] = x_0 + (x_1-x_0) * ratio_powers_sum
      data[xf_string][-1] \
          = f['MeshBlock'+str(f.attrs['TotalMeshBlock'][0]-1)][xf_string][:][-1]

  # Remove volume factors from restricted data
  for loc3 in range(lx3):
    for loc2 in range(lx2):
      for loc1 in range(lx1):
        if restricted_data[loc3,loc2,loc1]:
          il = loc1 * block_size[0]
          jl = loc2 * block_size[1]
          kl = loc3 * block_size[2]
          iu = il + block_size[0]
          ju = jl + block_size[1]
          ku = kl + block_size[2]
          for k in range(kl,ku):
            x3m = data['x3f'][k]
            x3p = data['x3f'][k+1]
            for j in range(jl,ju):
              x2m = data['x2f'][j]
              x2p = data['x2f'][j+1]
              for i in range(il,iu):
                x1m = data['x1f'][i]
                x1p = data['x1f'][i+1]
                vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
                for q in quantities:
                  data[q][k,j,i] /= vol
  return data

#=======================================================================================

class AthenaError(RuntimeError):
  """General exception class for Athena++ read functions."""
  pass
