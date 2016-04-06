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

def athdf(filename, data=None, quantities=None, level=0, subsample=False,
    fast_restrict=False, vol_func=None, coord='cartesian'):
  """Read .athdf files and populate dict of arrays of data."""

  # Python module for reading hdf5 files
  import h5py

  # Set volume function for preset coordinates if needed
  if not subsample and not fast_restrict and vol_func is None:
    if coord == 'cartesian':
      vol_func = lambda xm,xp,ym,yp,zm,zp: (xp-xm) * (yp-ym) * (zp-zm)
    elif coord == 'cylindrical':
      vol_func = lambda rm,rp,phim,phip,zm,zp: 0.5*(rp**2-rm**2) * (phip-phim) * (zp-zm)
    elif coord == 'spherical':
      vol_func = lambda rm,rp,thetam,thetap,phim,phip: \
          1.0/3.0*(rp**3-rm**3) * abs(np.cos(thetam)-np.cos(thetap)) * (phip-phim)
    else:
      raise AthenaError('Coordinates not recognized')

  # Open file
  with h5py.File(filename, 'r') as f:

    # Create list of all quantities if none given
    if data is not None:
      quantities = data.values()
    elif quantities is None:
      quantities = f['MeshBlock0'].keys()
    quantities = [str(q) for q in quantities \
        if q != 'x1f' and q != 'x2f' and q != 'x3f']

    # Extract size information
    block_size = f.attrs['MeshBlockSize']
    root_grid_size = f.attrs['RootGridSize']
    nx1 = root_grid_size[0] * 2**level
    nx2 = root_grid_size[1] * 2**level if root_grid_size[1] > 1 else 1
    nx3 = root_grid_size[2] * 2**level if root_grid_size[2] > 1 else 1
    lx1 = nx1 / block_size[0]
    lx2 = nx2 / block_size[1]
    lx3 = nx3 / block_size[2]
    if nx3 > 1:
      dim = 3
    elif nx2 > 1:
      dim = 2
    else:
      dim = 1

    # Check that subsampling and/or fast restriction will work if needed
    max_level = f.attrs['MaxLevel'][0]
    if subsample or fast_restrict:
      max_restrict_factor = 2**(max_level-level)
      for current_block_size in block_size[:dim]:
        if current_block_size % max_restrict_factor != 0:
          raise AthenaError('Block boundaries at finest level must be cell ' \
              + 'boundaries at desired level for\nsubsampling or fast restriction to ' \
              + 'work')

    # Prepare arrays to be returned
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

    # Prepare bookkeeping arrays
    x1f_level = np.full_like(data['x1f'], level+1, dtype=int)
    x2f_level = np.full_like(data['x2f'], level+1, dtype=int)
    x3f_level = np.full_like(data['x3f'], level+1, dtype=int)
    if not subsample and not fast_restrict and max_level > level:
      restricted_data = np.zeros((lx3,lx2,lx1), dtype=bool)

    # Account for singleton dimensions in arrays of face locations
    if nx2 == 1:
      data['x2f'] = f['MeshBlock0']['x2f'][:]
    if nx3 == 1:
      data['x3f'] = f['MeshBlock0']['x3f'][:]

    # Go through blocks in data file
    for block in f.itervalues():

      # Extract location information
      block_level = block.attrs['Level'][0]
      block_location = block.attrs['LogicalLocation']

      # Populate interface arrays if appropriate
      for d in np.arange(dim)+1:

        # Extract basics about this block and the direction
        loc = block_location[d-1]
        size = block_size[d-1]
        nx = (nx1,nx2,nx3)[d-1]
        xf_block = block['x'+str(d)+'f'][:]
        xf = data['x'+str(d)+'f']
        xf_level = (x1f_level,x2f_level,x3f_level)[d-1]

        # Refine coarse positions if no better positions are yet known to exist
        if block_level < level:
          level_diff = level - block_level
          s = 2**level_diff
          index_low = loc * size * s
          index_high = index_low + size * s
          if np.any(xf_level[index_low+1:index_high] <= level_diff):
            continue
          xf[index_low:index_high+1:s] = xf_block
          ratio_block = ((xf_block[-1]-xf_block[-2]) / (xf_block[1]-xf_block[0])) \
              ** (1.0/(size-1))
          for l in range(level_diff):
            ss = 2**(level_diff-l)
            ratio = ratio_block ** (1.0/2**(l+1))
            xf_low = xf[index_low:index_high:ss]
            xf_high = xf[index_low+ss:index_high+1:ss]
            xf[index_low+ss/2:index_high:ss] = xf_low + (xf_high-xf_low) / (1.0+ratio)
          xf_level[index_low:index_high+1] = \
              np.minimum(xf_level[index_low:index_high+1], level_diff)

        # Copy exact values from sufficiently refined block if values have not been set
        else:
          level_diff = block_level - level
          s = 2**level_diff
          index_low = loc * size
          index_high = index_low + size
          if index_low%s == 0:
            index_first_aligned = index_low
          else:
            index_first_aligned = index_low + s - index_low%s
          if np.all(xf_level[index_first_aligned/s:index_high/s+1] == 0):
            continue
          xf[index_first_aligned/s:index_high/s+1] = \
              xf_block[index_first_aligned-index_low:index_high+1-index_low:s]
          xf_level[index_first_aligned/s:index_high/s+1] = 0

      # Prolongate coarse data and copy same-level data
      if block_level <= level:

        # Calculate scale (number of copies per dimension)
        s = 2**(level-block_level)

        # Calculate fine-level begin indices
        il = block_location[0] * block_size[0] * s
        jl = block_location[1] * block_size[1] * s if dim >= 2 else 0
        kl = block_location[2] * block_size[2] * s if dim == 3 else 0

        # Calculate fine-level end indices
        iu = il + block_size[0] * s
        ju = jl + block_size[1] * s if dim >= 2 else 1
        ku = kl + block_size[2] * s if dim == 3 else 1

        # Calculate fine-level offsets
        io_vals = range(s)
        jo_vals = range(s) if dim >= 2 else (0,)
        ko_vals = range(s) if dim == 3 else (0,)

        # Assign values
        for q in quantities:
          for ko in ko_vals:
            for jo in jo_vals:
              for io in io_vals:
                data[q][kl+ko:ku+ko:s,jl+jo:ju+jo:s,il+io:iu+io:s] = block[q][:]

      # Restrict fine data
      else:

        # Apply subsampling
        if subsample:

          # Calculate scale (fine-level stride)
          s = 2**(block_level-level)

          # Calculate coarse-level begin indices
          il = block_location[0] * block_size[0] / s
          jl = block_location[1] * block_size[1] / s if dim >= 2 else 0
          kl = block_location[2] * block_size[2] / s if dim == 3 else 0

          # Calculate coarse-level end indices
          iu = il + block_size[0] / s
          ju = jl + block_size[1] / s if dim >= 2 else 1
          ku = kl + block_size[2] / s if dim == 3 else 1

          # Calculate fine-level offset (nearest cell at or below center)
          o = s/2 - 1

          # Assign values
          for q in quantities:
            data[q][kl:ku,jl:ju,il:iu] = block[q][o::s,o::s,o::s]

        # Apply fast (uniform Cartesian) restriction
        elif fast_restrict:

          # Calculate scale (fine-level stride)
          s = 2**(block_level-level)

          # Calculate coarse-level begin indices
          il = block_location[0] * block_size[0] / s
          jl = block_location[1] * block_size[1] / s if dim >= 2 else 0
          kl = block_location[2] * block_size[2] / s if dim == 3 else 0

          # Calculate coarse-level end indices
          iu = il + block_size[0] / s
          ju = jl + block_size[1] / s if dim >= 2 else 1
          ku = kl + block_size[2] / s if dim == 3 else 1

          # Calculate fine-level offsets
          io_vals = range(s)
          jo_vals = range(s) if dim >= 2 else (0,)
          ko_vals = range(s) if dim == 3 else (0,)

          # Assign values
          for q in quantities:
            for ko in ko_vals:
              for jo in jo_vals:
                for io in io_vals:
                  data[q][kl:ku,jl:ju,il:iu] += block[q][ko::s,jo::s,io::s]
            data[q][kl:ku,jl:ju,il:iu] /= s**dim

        # Apply exact (volume-weighted) restriction
        else:

          # Calculate scale
          s = 2**(block_level-level)

          # Calculate fine-level indices
          ir_vals = np.arange(block_size[0])
          jr_vals = np.arange(block_size[1])
          kr_vals = np.arange(block_size[2])

          # Calculate coarse-level indices
          i_vals = (ir_vals + block_location[0] * block_size[0]) / s
          j_vals = (jr_vals + block_location[1] * block_size[1]) / s
          k_vals = (kr_vals + block_location[2] * block_size[2]) / s

          # Accumulate values
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
                  data[q][k,j,i] += block[q][kr,jr,ir] * vol
          loc1 = block_location[0] / s
          loc2 = block_location[1] / s
          loc3 = block_location[2] / s
          restricted_data[loc3,loc2,loc1] = True

    # Remove volume factors from restricted data
    if not subsample and not fast_restrict and max_level > level:
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
