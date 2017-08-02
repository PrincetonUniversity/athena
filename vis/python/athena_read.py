"""
Read Athena++ output data files.
"""

# Python modules
import numpy as np
import re
import struct
import sys
import warnings

#=========================================================================================

def hst(filename, raw=False):
  """Read .hst files and return dict of 1D arrays."""

  # Read data
  with open(filename, 'r') as data_file:

    # Find header
    header_found = False
    multiple_headers = False
    header_location = None
    line = data_file.readline()
    while len(line) > 0:
      if line == '# Athena++ history data\n':
        if header_found:
          multiple_headers = True
        else:
          header_found = True
        header_location = data_file.tell()
      line = data_file.readline()
    if multiple_headers:
      warnings.warn('Multiple headers found; using most recent data', AthenaWarning)
    if header_location is None:
      raise AthenaError('No header found')

    # Parse header
    data_file.seek(header_location)
    header = data_file.readline()
    data_names = re.findall(r'\[\d+\]=(\S+)', header)
    if len(data_names) == 0:
      raise AthenaError('Header could not be parsed')

    # Prepare dictionary of results
    data = {}
    for name in data_names:
      data[name] = []

    # Read data
    for line in data_file:
      for name,val in zip(data_names, line.split()):
        data[name].append(float(val))

  # Finalize data
  for key,val in data.iteritems():
    data[key] = np.array(val)
  if not raw:
    if data_names[0] != 'time':
      raise AthenaError('Cannot remove spurious data because time column could not be' \
          + ' identified')
    branches_removed = False
    while not branches_removed:
      branches_removed = True
      for n in range(1, len(data['time'])):
        if data['time'][n] <= data['time'][n-1]:
          branch_index = np.where(data['time'][:n] >= data['time'][n])[0][0]
          for key,val in data.iteritems():
            data[key] = np.concatenate((val[:branch_index], val[n:]))
          branches_removed = False
          break
  return data

#=========================================================================================

def tab(filename, raw=False, dimensions=None):
  """Read .tab files and return dict or array."""

  # Check for valid number of dimensions
  if raw and not (dimensions == 1 or dimensions == 2 or dimensions == 3):
    raise AthenaError('Improper number of dimensions')
  if not raw and dimensions is not None:
    warnings.warn('Ignoring specified number of dimensions', AthenaWarning)

  # Parse header
  if not raw:
    data_dict = {}
    with open(filename, 'r') as data_file:
      line = data_file.readline()
      attributes = re.search(r'time=(\S+)\s+cycle=(\S+)\s+variables=(\S+)', line)
      line = data_file.readline()
      headings = line.split()[1:]
    data_dict['time'] = float(attributes.group(1))
    data_dict['cycle'] = int(attributes.group(2))
    data_dict['variables'] = attributes.group(3)
    if headings[0] == 'i' and headings[2] == 'j' and headings[4] == 'k':
      headings = headings[1:2] + headings[3:4] + headings[5:]
      dimensions = 3
    elif headings[0] == 'i' and headings[2] == 'j':
      headings = headings[1:2] + headings[3:]
      dimensions = 2
    elif headings[0] == 'i':
      headings = headings[1:]
      dimensions = 1
    else:
      raise AthenaError('Could not parse header')

  # Go through lines
  data_array = []
  with open(filename, 'r') as data_file:
    first_line = True
    for line in data_file:

      # Skip comments
      if line.split()[0][0] == '#':
        continue

      # Extract cell indices
      vals = line.split()
      if first_line:
        i_min = i_max = int(vals[0])
        if dimensions == 2 or dimensions == 3:
          j_min = j_max = int(vals[2])
        if dimensions == 3:
          k_min = k_max = int(vals[4])
        num_entries = len(vals) - dimensions
        first_line = False
      else:
        i_max = max(i_max, int(vals[0]))
        if dimensions == 2 or dimensions == 3:
          j_max = max(j_max, int(vals[2]))
        if dimensions == 3:
          k_max = max(k_max, int(vals[4]))

      # Extract cell values
      if dimensions == 1:
        vals = vals[1:]
      if dimensions == 2:
        vals = vals[1:2] + vals[3:]
      if dimensions == 3:
        vals = vals[1:2] + vals[3:4] + vals[5:]
      data_array.append([float(val) for val in vals])

  # Reshape array based on number of dimensions
  if dimensions == 1:
    array_shape = (i_max-i_min+1, num_entries)
    array_transpose = (1, 0)
  if dimensions == 2:
    array_shape = (j_max-j_min+1, i_max-i_min+1, num_entries)
    array_transpose = (2, 0, 1)
  if dimensions == 3:
    array_shape = (k_max-k_min+1, j_max-j_min+1, i_max-i_min+1, num_entries)
    array_transpose = (3, 0, 1, 2)
  data_array = np.transpose(np.reshape(data_array, array_shape), array_transpose)

  # Finalize data
  if raw:
    return data_array
  else:
    for n,heading in enumerate(headings):
      data_dict[heading] = data_array[n,...]
    return data_dict

#=========================================================================================

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
    identifier_string = '{0}_COORDINATES {1} float\n'.format(letter, num_faces)
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

#=========================================================================================

def athdf(filename, data=None, quantities=None, dtype=np.float32, level=None,
    subsample=False, fast_restrict=False, x1_min=None, x1_max=None, x2_min=None,
    x2_max=None, x3_min=None, x3_max=None, vol_func=None, vol_params=None,
    face_func_1=None, face_func_2=None, face_func_3=None, center_func_1=None,
    center_func_2=None, center_func_3=None):
  """Read .athdf files and populate dict of arrays of data."""

  # Python module
  import h5py

  # Prepare dictionary for results
  if data is None:
    data = {}
    new_data = True
  else:
    new_data = False

  # Open file
  with h5py.File(filename, 'r') as f:

    # Extract size information
    max_level = f.attrs['MaxLevel']
    if level is None:
      level = max_level
    block_size = f.attrs['MeshBlockSize']
    root_grid_size = f.attrs['RootGridSize']
    levels = f['Levels'][:]
    logical_locations = f['LogicalLocations'][:]
    nx_vals = []
    for d in range(3):
      if block_size[d] == 1 and root_grid_size[d] > 1:  # sum or slice
        other_locations = [location for location in \
            zip(levels, logical_locations[:,(d+1)%3], logical_locations[:,(d+2)%3])]
        if len(set(other_locations)) == len(other_locations):  # effective slice
          nx_vals.append(1)
        else:  # nontrivial sum
          nx_vals.append(2**level)
      elif block_size[d] == 1:  # singleton dimension
        nx_vals.append(1)
      else:  # normal case
        nx_vals.append(root_grid_size[d] * 2**level)
    nx1 = nx_vals[0]
    nx2 = nx_vals[1]
    nx3 = nx_vals[2]
    lx1 = nx1 / block_size[0]
    lx2 = nx2 / block_size[1]
    lx3 = nx3 / block_size[2]
    num_extended_dims = 0
    for nx in nx_vals:
      if nx > 1:
        num_extended_dims += 1

    # Set volume function for preset coordinates if needed
    coord = f.attrs['Coordinates']
    if level < max_level and not subsample and not fast_restrict and vol_func is None:
      x1_rat = f.attrs['RootGridX1'][2]
      x2_rat = f.attrs['RootGridX2'][2]
      x3_rat = f.attrs['RootGridX3'][2]
      if coord == 'cartesian' or coord == 'minkowski' or coord == 'tilted' \
          or coord == 'sinusoidal':
        if (nx1 == 1 or x1_rat == 1.0) and (nx2 == 1 or x2_rat == 1.0) and \
            (nx3 == 1 or x3_rat == 1.0):
          fast_restrict = True
        else:
          vol_func = lambda xm,xp,ym,yp,zm,zp: (xp-xm) * (yp-ym) * (zp-zm)
      elif coord == 'cylindrical':
        if nx1 == 1 and (nx2 == 1 or x2_rat == 1.0) and (nx3 == 1 or x3_rat == 1.0):
          fast_restrict = True
        else:
          vol_func = lambda rm,rp,phim,phip,zm,zp: (rp**2-rm**2) * (phip-phim) * (zp-zm)
      elif coord == 'spherical_polar' or coord == 'schwarzschild':
        if nx1 == 1 and nx2 == 1 and (nx3 == 1 or x3_rat == 1.0):
          fast_restrict = True
        else:
          vol_func = lambda rm,rp,thetam,thetap,phim,phip: \
              (rp**3-rm**3) * abs(np.cos(thetam)-np.cos(thetap)) * (phip-phim)
      elif coord == 'kerr-schild':
        if nx1 == 1 and nx2 == 1 and (nx3 == 3 or x3_rat == 1.0):
          fast_restrict = True
        else:
          a = vol_params[0]
          def vol_func(rm, rp, thetam, thetap, phim, phip):
            cosm = np.cos(thetam)
            cosp = np.cos(thetap)
            return \
                ((rp**3-rm**3) * abs(cosm-cosp) + a**2 * (rp-rm) * abs(cosm**3-cosp**3)) \
                * (phip-phim)
      else:
        raise AthenaError('Coordinates not recognized')

    # Set cell center functions for preset coordinates
    if center_func_1 is None:
      if coord == 'cartesian' or coord == 'minkowski' or coord == 'tilted' \
          or coord == 'sinusoidal' or coord == 'kerr-schild':
        center_func_1 = lambda xm,xp : 0.5*(xm+xp)
      elif coord == 'cylindrical':
        center_func_1 = lambda xm,xp : 2.0/3.0 * (xp**3-xm**3) / (xp**2-xm**2)
      elif coord == 'spherical_polar':
        center_func_1 = lambda xm,xp : 3.0/4.0 * (xp**4-xm**4) / (xp**3-xm**3)
      elif coord == 'schwarzschild':
        center_func_1 = lambda xm,xp : (0.5*(xm**3+xp**3)) ** 1.0/3.0
      else:
        raise AthenaError('Coordinates not recognized')
    if center_func_2 is None:
      if coord == 'cartesian' or coord == 'cylindrical' or coord == 'minkowski' \
          or coord == 'tilted' or coord == 'sinusoidal' or coord == 'kerr-schild':
        center_func_2 = lambda xm,xp : 0.5*(xm+xp)
      elif coord == 'spherical_polar':
        def center_func_2(xm, xp):
          sm = np.sin(xm)
          cm = np.cos(xm)
          sp = np.sin(xp)
          cp = np.cos(xp)
          return (sp-xp*cp - sm+xm*cm) / (cm-cp)
      elif coord == 'schwarzschild':
        center_func_2 = lambda xm,xp : np.arccos(0.5*(np.cos(xm)+np.cos(xp)))
      else:
        raise AthenaError('Coordinates not recognized')
    if center_func_3 is None:
      if coord == 'cartesian' or coord == 'cylindrical' or coord == 'spherical_polar' \
          or coord == 'minkowski' or coord == 'tilted' or coord == 'sinusoidal' \
          or coord == 'schwarzschild' or coord == 'kerr-schild':
        center_func_3 = lambda xm,xp : 0.5*(xm+xp)
      else:
        raise AthenaError('Coordinates not recognized')

    # Check output level compared to max level in file
    if level < max_level and not subsample and not fast_restrict:
      warnings.warn('Exact restriction being used: performance severely affected; see' \
          + ' documentation', AthenaWarning)
      sys.stderr.flush()
    if level > max_level:
      warnings.warn('Requested refinement level higher than maximum level in file: all' \
          + ' cells will be prolongated', AthenaWarning)
      sys.stderr.flush()

    # Check that subsampling and/or fast restriction will work if needed
    if level < max_level and (subsample or fast_restrict):
      max_restrict_factor = 2**(max_level-level)
      for current_block_size in block_size:
        if current_block_size != 1 and current_block_size%max_restrict_factor != 0:
          raise AthenaError('Block boundaries at finest level must be cell boundaries' \
              + ' at desired level for subsampling or fast restriction to work')

    # Create list of all quantities if none given
    var_quantities = f.attrs['VariableNames'][:]
    coord_quantities = ('x1f', 'x2f', 'x3f', 'x1v', 'x2v', 'x3v')
    attr_quantities = [key for key in f.attrs]
    if not new_data:
      quantities = data.values()
    elif quantities is None:
      quantities = var_quantities
    else:
      for q in quantities:
        if q not in var_quantities and q not in coord_quantities:
          possibilities = '", "'.join(var_quantities)
          possibilities = '"' + possibilities + '"'
          error_string = 'Quantity not recognized: file does not include "{0}" but does' \
              + ' include {1}'
          raise AthenaError(error_string.format(q, possibilities))
    quantities = [str(q) for q in quantities \
        if q not in coord_quantities and q not in attr_quantities]

    # Store file attribute metadata
    for key in attr_quantities:
      data[str(key)] = f.attrs[key]

    # Get metadata describing file layout
    num_blocks = f.attrs['NumMeshBlocks']
    dataset_names = f.attrs['DatasetNames'][:]
    dataset_sizes = f.attrs['NumVariables'][:]
    dataset_sizes_cumulative = np.cumsum(dataset_sizes)
    variable_names = f.attrs['VariableNames'][:]
    quantity_datasets = []
    quantity_indices = []
    for q in quantities:
      var_num = np.where(variable_names == q)[0][0]
      dataset_num = np.where(dataset_sizes_cumulative > var_num)[0][0]
      if dataset_num == 0:
        dataset_index = var_num
      else:
        dataset_index = var_num - dataset_sizes_cumulative[dataset_num-1]
      quantity_datasets.append(dataset_names[dataset_num])
      quantity_indices.append(dataset_index)

    # Locate fine block for coordinates in case of slice
    fine_block = np.where(levels == max_level)[0][0]
    x1m = f['x1f'][fine_block,0]
    x1p = f['x1f'][fine_block,1]
    x2m = f['x2f'][fine_block,0]
    x2p = f['x2f'][fine_block,1]
    x3m = f['x3f'][fine_block,0]
    x3p = f['x3f'][fine_block,1]

    # Populate coordinate arrays
    face_funcs = (face_func_1, face_func_2, face_func_3)
    center_funcs = (center_func_1, center_func_2, center_func_3)
    for d,nx,face_func,center_func in zip(range(1, 4), nx_vals, face_funcs, center_funcs):
      if nx == 1:
        xm = (x1m, x2m, x3m)[d-1]
        xp = (x1p, x2p, x3p)[d-1]
        data['x'+repr(d)+'f'] = np.array([xm, xp])
      else:
        xmin = f.attrs['RootGridX'+repr(d)][0]
        xmax = f.attrs['RootGridX'+repr(d)][1]
        xrat_root = f.attrs['RootGridX'+repr(d)][2]
        if (face_func is not None):
          data['x'+repr(d)+'f'] = face_func(xmin, xmax, xrat_root, nx+1)
        elif (xrat_root == 1.0):
          data['x'+repr(d)+'f'] = np.linspace(xmin, xmax, nx+1)
        else:
          xrat = xrat_root ** (1.0 / 2**level)
          data['x'+repr(d)+'f'] = \
              xmin + (1.0-xrat**np.arange(nx+1)) / (1.0-xrat**nx) * (xmax-xmin)
      data['x'+repr(d)+'v'] = np.empty(nx)
      for i in range(nx):
        data['x'+repr(d)+'v'][i] = \
            center_func(data['x'+repr(d)+'f'][i], data['x'+repr(d)+'f'][i+1])

    # Account for selection
    x1_select = False
    x2_select = False
    x3_select = False
    i_min = j_min = k_min = 0
    i_max = nx1
    j_max = nx2
    k_max = nx3
    error_string = '{0} must be {1} than {2} in order to intersect data range'
    if x1_min is not None and x1_min >= data['x1f'][1]:
      if x1_min >= data['x1f'][-1]:
        raise AthenaError(error_string.format('x1_min', 'less', data['x1f'][-1]))
      x1_select = True
      i_min = np.where(data['x1f'] <= x1_min)[0][-1]
    if x1_max is not None and x1_max <= data['x1f'][-2]:
      if x1_max <= data['x1f'][0]:
        raise AthenaError(error_string.format('x1_max', 'greater', data['x1f'][0]))
      x1_select = True
      i_max = np.where(data['x1f'] >= x1_max)[0][0]
    if x2_min is not None and x2_min >= data['x2f'][1]:
      if x2_min >= data['x2f'][-1]:
        raise AthenaError(error_string.format('x2_min', 'less', data['x2f'][-1]))
      x2_select = True
      j_min = np.where(data['x2f'] <= x2_min)[0][-1]
    if x2_max is not None and x2_max <= data['x2f'][-2]:
      if x2_max <= data['x2f'][0]:
        raise AthenaError(error_string.format('x2_max', 'greater', data['x2f'][0]))
      x2_select = True
      j_max = np.where(data['x2f'] >= x2_max)[0][0]
    if x3_min is not None and x3_min >= data['x3f'][1]:
      if x3_min >= data['x3f'][-1]:
        raise AthenaError(error_string.format('x3_min', 'less', data['x3f'][-1]))
      x3_select = True
      k_min = np.where(data['x3f'] <= x3_min)[0][-1]
    if x3_max is not None and x3_max <= data['x3f'][-2]:
      if x3_max <= data['x3f'][0]:
        raise AthenaError(error_string.format('x3_max', 'greater', data['x3f'][0]))
      x3_select = True
      k_max = np.where(data['x3f'] >= x3_max)[0][0]

    # Adjust coordinates if selection made
    if x1_select:
      data['x1f'] = data['x1f'][i_min:i_max+1]
      data['x1v'] = data['x1v'][i_min:i_max]
    if x2_select:
      data['x2f'] = data['x2f'][j_min:j_max+1]
      data['x2v'] = data['x2v'][j_min:j_max]
    if x3_select:
      data['x3f'] = data['x3f'][k_min:k_max+1]
      data['x3v'] = data['x3v'][k_min:k_max]

    # Prepare arrays for data and bookkeeping
    if new_data:
      for q in quantities:
        data[q] = np.zeros((k_max-k_min, j_max-j_min, i_max-i_min), dtype=dtype)
    else:
      for q in quantities:
        data[q].fill(0.0)
    if not subsample and not fast_restrict and max_level > level:
      restricted_data = np.zeros((lx3, lx2, lx1), dtype=bool)

    # Go through blocks in data file
    for block_num in range(num_blocks):

      # Extract location information
      block_level = levels[block_num]
      block_location = logical_locations[block_num,:]

      # Prolongate coarse data and copy same-level data
      if block_level <= level:

        # Calculate scale (number of copies per dimension)
        s = 2 ** (level - block_level)

        # Calculate destination indices, without selection
        il_d = block_location[0] * block_size[0] * s if nx1 > 1 else 0
        jl_d = block_location[1] * block_size[1] * s if nx2 > 1 else 0
        kl_d = block_location[2] * block_size[2] * s if nx3 > 1 else 0
        iu_d = il_d + block_size[0] * s if nx1 > 1 else 1
        ju_d = jl_d + block_size[1] * s if nx2 > 1 else 1
        ku_d = kl_d + block_size[2] * s if nx3 > 1 else 1

        # Calculate (prolongated) source indices, with selection
        il_s = max(il_d, i_min) - il_d
        jl_s = max(jl_d, j_min) - jl_d
        kl_s = max(kl_d, k_min) - kl_d
        iu_s = min(iu_d, i_max) - il_d
        ju_s = min(ju_d, j_max) - jl_d
        ku_s = min(ku_d, k_max) - kl_d
        if il_s >= iu_s or jl_s >= ju_s or kl_s >= ku_s:
          continue

        # Account for selection in destination indices
        il_d = max(il_d, i_min) - i_min
        jl_d = max(jl_d, j_min) - j_min
        kl_d = max(kl_d, k_min) - k_min
        iu_d = min(iu_d, i_max) - i_min
        ju_d = min(ju_d, j_max) - j_min
        ku_d = min(ku_d, k_max) - k_min

        # Assign values
        for q,dataset,index in zip(quantities, quantity_datasets, quantity_indices):
          block_data = f[dataset][index,block_num,:]
          if s > 1:
            if nx1 > 1:
              block_data = np.repeat(block_data, s, axis=2)
            if nx2 > 1:
              block_data = np.repeat(block_data, s, axis=1)
            if nx3 > 1:
              block_data = np.repeat(block_data, s, axis=0)
          data[q][kl_d:ku_d,jl_d:ju_d,il_d:iu_d] = \
              block_data[kl_s:ku_s,jl_s:ju_s,il_s:iu_s]

      # Restrict fine data
      else:

        # Calculate scale
        s = 2 ** (block_level - level)

        # Calculate destination indices, without selection
        il_d = block_location[0] * block_size[0] / s if nx1 > 1 else 0
        jl_d = block_location[1] * block_size[1] / s if nx2 > 1 else 0
        kl_d = block_location[2] * block_size[2] / s if nx3 > 1 else 0
        iu_d = il_d + block_size[0] / s if nx1 > 1 else 1
        ju_d = jl_d + block_size[1] / s if nx2 > 1 else 1
        ku_d = kl_d + block_size[2] / s if nx3 > 1 else 1

        # Calculate (restricted) source indices, with selection
        il_s = max(il_d, i_min) - il_d
        jl_s = max(jl_d, j_min) - jl_d
        kl_s = max(kl_d, k_min) - kl_d
        iu_s = min(iu_d, i_max) - il_d
        ju_s = min(ju_d, j_max) - jl_d
        ku_s = min(ku_d, k_max) - kl_d
        if il_s >= iu_s or jl_s >= ju_s or kl_s >= ku_s:
          continue

        # Account for selection in destination indices
        il_d = max(il_d, i_min) - i_min
        jl_d = max(jl_d, j_min) - j_min
        kl_d = max(kl_d, k_min) - k_min
        iu_d = min(iu_d, i_max) - i_min
        ju_d = min(ju_d, j_max) - j_min
        ku_d = min(ku_d, k_max) - k_min

        # Account for restriction in source indices
        if nx1 > 1:
          il_s *= s
          iu_s *= s
        if nx2 > 1:
          jl_s *= s
          ju_s *= s
        if nx3 > 1:
          kl_s *= s
          ku_s *= s

        # Apply subsampling
        if subsample:

          # Calculate fine-level offsets (nearest cell at or below center)
          o1 = s/2 - 1 if nx1 > 1 else 0
          o2 = s/2 - 1 if nx2 > 1 else 0
          o3 = s/2 - 1 if nx3 > 1 else 0

          # Assign values
          for q,dataset,index in zip(quantities, quantity_datasets, quantity_indices):
            data[q][kl_d:ku_d,jl_d:ju_d,il_d:iu_d] = \
                f[dataset][index,block_num,kl_s+o3:ku_s:s,jl_s+o2:ju_s:s,il_s+o1:iu_s:s]

        # Apply fast (uniform Cartesian) restriction
        elif fast_restrict:

          # Calculate fine-level offsets
          io_vals = range(s) if nx1 > 1 else (0,)
          jo_vals = range(s) if nx2 > 1 else (0,)
          ko_vals = range(s) if nx3 > 1 else (0,)

          # Assign values
          for q,dataset,index in zip(quantities, quantity_datasets, quantity_indices):
            for ko in ko_vals:
              for jo in jo_vals:
                for io in io_vals:
                  data[q][kl_d:ku_d,jl_d:ju_d,il_d:iu_d] += \
                      f[dataset]\
                      [index,block_num,kl_s+ko:ku_s:s,jl_s+jo:ju_s:s,il_s+io:iu_s:s]
            data[q][kl_d:ku_d,jl_d:ju_d,il_d:iu_d] /= s ** num_extended_dims

        # Apply exact (volume-weighted) restriction
        else:

          # Calculate sets of indices
          i_s_vals = range(il_s, iu_s)
          j_s_vals = range(jl_s, ju_s)
          k_s_vals = range(kl_s, ku_s)
          i_d_vals = range(il_d, iu_d)
          j_d_vals = range(jl_d, ju_d)
          k_d_vals = range(kl_d, ku_d)
          if nx1 > 1:
            i_d_vals = np.repeat(i_d_vals, s)
          if nx2 > 1:
            j_d_vals = np.repeat(j_d_vals, s)
          if nx3 > 1:
            k_d_vals = np.repeat(k_d_vals, s)

          # Accumulate values
          for k_s,k_d in zip(k_s_vals, k_d_vals):
            if nx3 > 1:
              x3m = f['x3f'][block_num,k_s]
              x3p = f['x3f'][block_num,k_s+1]
            for j_s,j_d in zip(j_s_vals, j_d_vals):
              if nx2 > 1:
                x2m = f['x2f'][block_num,j_s]
                x2p = f['x2f'][block_num,j_s+1]
              for i_s,i_d in zip(i_s_vals, i_d_vals):
                if nx1 > 1:
                  x1m = f['x1f'][block_num,i_s]
                  x1p = f['x1f'][block_num,i_s+1]
                vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
                for q,dataset,index in \
                    zip(quantities, quantity_datasets, quantity_indices):
                  data[q][k_d,j_d,i_d] += f[dataset][index,block_num,k_s,j_s,i_s] * vol
          loc1 = (nx1 > 1 ) * block_location[0] / s
          loc2 = (nx2 > 1 ) * block_location[1] / s
          loc3 = (nx3 > 1 ) * block_location[2] / s
          restricted_data[loc3,loc2,loc1] = True

  # Remove volume factors from restricted data
  if level < max_level and not subsample and not fast_restrict:
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
            il = max(il, i_min) - i_min
            jl = max(jl, j_min) - j_min
            kl = max(kl, k_min) - k_min
            iu = min(iu, i_max) - i_min
            ju = min(ju, j_max) - j_min
            ku = min(ku, k_max) - k_min
            for k in range(kl, ku):
              if nx3 > 1:
                x3m = data['x3f'][k]
                x3p = data['x3f'][k+1]
              for j in range(jl, ju):
                if nx2 > 1:
                  x2m = data['x2f'][j]
                  x2p = data['x2f'][j+1]
                for i in range(il, iu):
                  if nx1 > 1:
                    x1m = data['x1f'][i]
                    x1p = data['x1f'][i+1]
                  vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
                  for q in quantities:
                    data[q][k,j,i] /= vol

  # Return dictionary containing requested data arrays
  return data

#=========================================================================================

class AthenaError(RuntimeError):
  """General exception class for Athena++ read functions."""
  pass

class AthenaWarning(RuntimeWarning):
  """General warning class for Athena++ read functions."""
  pass
