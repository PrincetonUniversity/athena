"""
Read Athena++ output data files.
"""

# Python modules
import re
import struct
import sys
import warnings
from io import open  # Consistent binary I/O from Python 2 and 3

# Other Python modules
import numpy as np

# Load HDF5 reader
try:
    import h5py
except ImportError:
    pass


# ========================================================================================


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
            for name, val in zip(data_names, line.split()):
                data[name].append(float(val))

    # Finalize data
    for key, val in data.iteritems():
        data[key] = np.array(val)
    if not raw:
        if data_names[0] != 'time':
            raise AthenaError(
                    'Cannot remove spurious data because time column could not be' +
                    ' identified')
        branches_removed = False
        while not branches_removed:
            branches_removed = True
            for n in range(1, len(data['time'])):
                if data['time'][n] <= data['time'][n-1]:
                    branch_index = np.where(data['time'][:n] >= data['time'][n])[0][0]
                    for key, val in data.iteritems():
                        data[key] = np.concatenate((val[:branch_index], val[n:]))
                    branches_removed = False
                    break
    return data


# ========================================================================================


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
        elif ((headings[0] == 'i' and headings[2] == 'j') or
                (headings[0] == 'i' and headings[2] == 'k') or
                (headings[0] == 'j' and headings[2] == 'k')):
            headings = headings[1:2] + headings[3:]
            dimensions = 2
        elif headings[0] == 'i' or headings[0] == 'j' or headings[0] == 'k':
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
        for n, heading in enumerate(headings):
            data_dict[heading] = data_array[n, ...]
        return data_dict


# ========================================================================================


def vtk(filename):
    """Read .vtk files and return dict of arrays of data."""

    # Read raw data
    with open(filename, 'rb') as data_file:
        raw_data = data_file.read()
    raw_data_ascii = raw_data.decode('ascii', 'replace')

    # Skip header
    current_index = 0
    current_char = raw_data_ascii[current_index]
    while current_char == '#':
        while current_char != '\n':
            current_index += 1
            current_char = raw_data_ascii[current_index]
        current_index += 1
        current_char = raw_data_ascii[current_index]

    # Function for skipping though the file
    def skip_string(expected_string):
        expected_string_len = len(expected_string)
        if raw_data_ascii[current_index:current_index +
                          expected_string_len] != expected_string:
            raise AthenaError('File not formatted as expected')
        return current_index+expected_string_len

    # Read metadata
    current_index = skip_string('BINARY\nDATASET RECTILINEAR_GRID\nDIMENSIONS ')
    end_of_line_index = current_index + 1
    while raw_data_ascii[end_of_line_index] != '\n':
        end_of_line_index += 1
    face_dimensions = list(map(
            int, raw_data_ascii[current_index:end_of_line_index].split(' ')))
    current_index = end_of_line_index + 1

    # Function for reading interface locations
    def read_faces(letter, num_faces):
        identifier_string = '{0}_COORDINATES {1} float\n'.format(letter, num_faces)
        begin_index = skip_string(identifier_string)
        format_string = '>' + 'f'*num_faces
        end_index = begin_index + 4*num_faces
        vals = np.array(struct.unpack(format_string, raw_data[begin_index:end_index]))
        return vals, end_index+1

    # Read interface locations
    x_faces, current_index = read_faces('X', face_dimensions[0])
    y_faces, current_index = read_faces('Y', face_dimensions[1])
    z_faces, current_index = read_faces('Z', face_dimensions[2])

    # Prepare to read quantities defined on grid
    cell_dimensions = np.array([max(dim-1, 1)
                                for dim in face_dimensions])
    num_cells = cell_dimensions.prod()
    current_index = skip_string('CELL_DATA {0}\n'.format(num_cells))
    if raw_data_ascii[current_index:current_index+1] == '\n':
        current_index = skip_string('\n')  # extra newline inserted by join script
    data = {}

    # Function for reading scalar data
    def read_cell_scalars():
        begin_index = skip_string('SCALARS ')
        end_of_word_index = begin_index + 1
        while raw_data_ascii[end_of_word_index] != ' ':
            end_of_word_index += 1
        array_name = raw_data_ascii[begin_index:end_of_word_index]
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
        while raw_data_ascii[end_of_word_index] != '\n':
            end_of_word_index += 1
        array_name = raw_data_ascii[begin_index:end_of_word_index]
        string_to_skip = 'VECTORS {0}\n'.format(array_name)
        array_name = array_name[:-6]  # remove ' float'
        begin_index = skip_string(string_to_skip)
        format_string = '>' + 'f'*num_cells*3
        end_index = begin_index + 4*num_cells*3
        data[array_name] = struct.unpack(format_string, raw_data[begin_index:end_index])
        dimensions = tuple(np.append(cell_dimensions[::-1], 3))
        data[array_name] = np.array(data[array_name]).reshape(dimensions)
        return end_index+1

    # Read quantities defined on grid
    while current_index < len(raw_data):
        expected_string = 'SCALARS'
        expected_string_len = len(expected_string)
        if raw_data_ascii[current_index:current_index +
                          expected_string_len] == expected_string:
            current_index = read_cell_scalars()
            continue
        expected_string = 'VECTORS'
        expected_string_len = len(expected_string)
        if raw_data_ascii[current_index:current_index +
                          expected_string_len] == expected_string:
            current_index = read_cell_vectors()
            continue
        raise AthenaError('File not formatted as expected')
    return x_faces, y_faces, z_faces, data


# ========================================================================================


class athdf(dict):
    """Read .athdf files and populate dict of arrays of data."""

    # Initialization
    def __init__(
            self,
            filename,
            data=None,
            quantities=None,
            dtype=np.float32,
            level=None,
            return_levels=False,
            subsample=False,
            fast_restrict=False,
            x1_min=None,
            x1_max=None,
            x2_min=None,
            x2_max=None,
            x3_min=None,
            x3_max=None,
            vol_func=None,
            vol_params=None,
            face_func_1=None,
            face_func_2=None,
            face_func_3=None,
            center_func_1=None,
            center_func_2=None,
            center_func_3=None):

        # Import necessary module for reading HDF5 files
        try:
            h5py
        except NameError:
            raise ImportError(
                    'athdf could not be executed because h5py could not be imported.')

        # Prepare dictionary for results
        if data is None:
            data = {}
            self.new_data = True
            self._existing_keys = []
        else:
            self.new_data = False
            self._existing_keys = list(data.keys())
        super(athdf, self).__init__(data)
        self.filename = filename
        self.quantities = quantities
        self.dtype = dtype
        self.level = level
        self.return_levels = return_levels
        self.subsample = subsample
        self.fast_restrict = fast_restrict
        self.x1_min = x1_min
        self.x1_max = x1_max
        self.x2_min = x2_min
        self.x2_max = x2_max
        self.x3_min = x3_min
        self.x3_max = x3_max
        self.vol_func = vol_func
        self.vol_params = vol_params
        self.face_func_1 = face_func_1
        self.face_func_2 = face_func_2
        self.face_func_3 = face_func_3
        self.center_func_1 = center_func_1
        self.center_func_2 = center_func_2
        self.center_func_3 = center_func_3

        # Open file
        with h5py.File(filename, 'r') as f:

            # Extract size information
            self.max_level = f.attrs['MaxLevel']
            if self.level is None:
                self.level = self.max_level
            self.block_size = f.attrs['MeshBlockSize']
            root_grid_size = f.attrs['RootGridSize']
            self.levels = f['Levels'][:]
            self.logical_locations = f['LogicalLocations'][:]

            # Determine sums and slices
            nx_vals = []
            for d in range(3):

                # Sum or slice
                if self.block_size[d] == 1 and root_grid_size[d] > 1:
                    other_locations = [
                            location for location in zip(
                                self.levels,
                                self.logical_locations[:, (d+1) % 3],
                                self.logical_locations[:, (d+2) % 3])]

                    # Effective slice
                    if len(set(other_locations)) == len(other_locations):
                        nx_vals.append(1)

                    # Nontrivial sum
                    else:
                        num_blocks_this_dim = 0
                        for (level_this_dim, loc_this_dim) in zip(
                                self.levels, self.logical_locations[:, d]):
                            if level_this_dim <= self.level:
                                num_blocks_this_dim = max(
                                        num_blocks_this_dim, (loc_this_dim + 1)
                                        * 2 ** (self.level-level_this_dim))
                            else:
                                num_blocks_this_dim = max(
                                        num_blocks_this_dim, (loc_this_dim + 1)
                                        / 2 ** (level_this_dim-self.level))
                        nx_vals.append(num_blocks_this_dim)

                # Singleton dimension
                elif self.block_size[d] == 1:
                    nx_vals.append(1)

                # Normal case
                else:
                    nx_vals.append(root_grid_size[d] * 2 ** self.level)

            # Set dimensions
            self.nx1 = nx_vals[0]
            self.nx2 = nx_vals[1]
            self.nx3 = nx_vals[2]
            self.lx1 = self.nx1 / self.block_size[0]
            self.lx2 = self.nx2 / self.block_size[1]
            self.lx3 = self.nx3 / self.block_size[2]
            self.num_extended_dims = 0
            for nx in nx_vals:
                if nx > 1:
                    self.num_extended_dims += 1

            # Set volume function for preset coordinates if needed
            coord = f.attrs['Coordinates'].decode('ascii', 'replace')
            if (self.level < self.max_level and not self.subsample
                    and not self.fast_restrict and self.vol_func is None):
                x1_rat = f.attrs['RootGridX1'][2].decode('ascii', 'replace')
                x2_rat = f.attrs['RootGridX2'][2].decode('ascii', 'replace')
                x3_rat = f.attrs['RootGridX3'][2].decode('ascii', 'replace')
                if (coord == 'cartesian' or coord == 'minkowski' or coord == 'tilted'
                        or coord == 'sinusoidal'):
                    if (
                            (self.nx1 == 1 or x1_rat == 1.0)
                            and (self.nx2 == 1 or x2_rat == 1.0)
                            and (self.nx3 == 1 or x3_rat == 1.0)):
                        self.fast_restrict = True
                    else:
                        self.vol_func = lambda xm, xp, ym, yp, zm, zp: (
                                xp - xm) * (yp - ym) * (zp - zm)
                elif coord == 'cylindrical':
                    if (self.nx1 == 1 and (self.nx2 == 1 or x2_rat == 1.0)
                            and (self.nx3 == 1 or x3_rat == 1.0)):
                        self.fast_restrict = True
                    else:
                        self.vol_func = lambda rm, rp, phim, phip, zm, zp: (
                                rp**2 - rm**2) * (phip - phim) * (zp - zm)
                elif coord == 'spherical_polar' or coord == 'schwarzschild':
                    if self.nx1 == 1 and self.nx2 == 1 and (
                            self.nx3 == 1 or x3_rat == 1.0):
                        self.fast_restrict = True
                    else:
                        self.vol_func = lambda rm, rp, thetam, thetap, phim, phip: ((
                                rp**3 - rm**3) * abs(np.cos(thetam) - np.cos(thetap))
                                * (phip - phim))
                elif coord == 'kerr-schild':
                    if self.nx1 == 1 and self.nx2 == 1 and (
                            self.nx3 == 1 or x3_rat == 1.0):
                        self.fast_restrict = True
                    else:
                        a = vol_params[0]

                        def vol_func(rm, rp, thetam, thetap, phim, phip):
                            cosm = np.cos(thetam)
                            cosp = np.cos(thetap)
                            return ((rp**3 - rm**3) * abs(cosm - cosp) + a**2 *
                                    (rp - rm) * abs(cosm**3 - cosp**3)) * (phip - phim)
                else:
                    raise AthenaError('Coordinates not recognized')

            # Set cell center functions for preset coordinates
            if center_func_1 is None:
                if (coord == 'cartesian' or coord == 'minkowski' or coord == 'tilted'
                        or coord == 'sinusoidal' or coord == 'kerr-schild'):
                    def center_func_1(xm, xp): return 0.5 * (xm + xp)
                elif coord == 'cylindrical':
                    def center_func_1(xm, xp): return (
                            2.0/3.0 * (xp**3 - xm**3) / (xp**2 - xm**2))
                elif coord == 'spherical_polar':
                    def center_func_1(xm, xp): return (
                            3.0/4.0 * (xp**4 - xm**4) / (xp**3 - xm**3))
                elif coord == 'schwarzschild':
                    def center_func_1(xm, xp): return (0.5 * (xm**3 + xp**3)) ** (1.0/3.0)
                else:
                    raise AthenaError('Coordinates not recognized')
            if center_func_2 is None:
                if (coord == 'cartesian' or coord == 'cylindrical' or coord == 'minkowski'
                        or coord == 'tilted' or coord == 'sinusoidal'
                        or coord == 'kerr-schild'):
                    def center_func_2(xm, xp): return 0.5 * (xm + xp)
                elif coord == 'spherical_polar':
                    def center_func_2(xm, xp):
                        sm = np.sin(xm)
                        cm = np.cos(xm)
                        sp = np.sin(xp)
                        cp = np.cos(xp)
                        return (sp - xp*cp - sm + xm*cm) / (cm-cp)
                elif coord == 'schwarzschild':
                    def center_func_2(xm, xp): return np.arccos(
                            0.5 * (np.cos(xm) + np.cos(xp)))
                else:
                    raise AthenaError('Coordinates not recognized')
            if center_func_3 is None:
                if (coord == 'cartesian' or coord == 'cylindrical'
                        or coord == 'spherical_polar' or coord == 'minkowski'
                        or coord == 'tilted' or coord == 'sinusoidal'
                        or coord == 'schwarzschild' or coord == 'kerr-schild'):
                    def center_func_3(xm, xp): return 0.5 * (xm + xp)
                else:
                    raise AthenaError('Coordinates not recognized')

            # Check output level compared to max level in file
            if (self.level < self.max_level and not self.subsample
                    and not self.fast_restrict):
                warnings.warn(
                        'Exact restriction being used: performance severely affected;' +
                        ' see documentation',
                        AthenaWarning)
                sys.stderr.flush()
            if self.level > self.max_level:
                warnings.warn(
                        'Requested refinement level higher than maximum level in file:' +
                        ' all cells will be prolongated',
                        AthenaWarning)
                sys.stderr.flush()

            # Check that subsampling and/or fast restriction will work if needed
            if self.level < self.max_level and (self.subsample or self.fast_restrict):
                max_restrict_factor = 2 ** (self.max_level - self.level)
                for current_block_size in self.block_size:
                    if (current_block_size != 1 and current_block_size %
                            max_restrict_factor != 0):
                        raise AthenaError(
                                'Block boundaries at finest level must be cell' +
                                ' boundaries at desired level for subsampling or fast' +
                                ' restriction to work')

            # Create list of all quantities if none given
            var_quantities = np.array([x.decode('ascii', 'replace')
                                       for x in f.attrs['VariableNames'][:]])
            coord_quantities = ('x1f', 'x2f', 'x3f', 'x1v', 'x2v', 'x3v')
            attr_quantities = [key for key in f.attrs]
            other_quantities = ('Levels',)
            if not self.new_data:
                self.quantities = self.values()
            elif self.quantities is None:
                self.quantities = var_quantities
            else:
                for q in self.quantities:
                    if q not in var_quantities and q not in coord_quantities:
                        possibilities = '", "'.join(var_quantities)
                        possibilities = '"' + possibilities + '"'
                        error_string = (
                                'Quantity not recognized: file does not include "{0}" but'
                                ' does include {1}')
                        raise AthenaError(error_string.format(q, possibilities))
            self.quantities = [str(q) for q in self.quantities
                               if q not in coord_quantities
                               and q not in attr_quantities
                               and q not in other_quantities]

            # Store file attribute metadata
            for key in attr_quantities:
                self[str(key)] = f.attrs[key]

            # Get metadata describing file layout
            self.num_blocks = f.attrs['NumMeshBlocks']
            dataset_names = np.array([x.decode('ascii', 'replace')
                                      for x in f.attrs['DatasetNames'][:]])
            dataset_sizes = f.attrs['NumVariables'][:]
            dataset_sizes_cumulative = np.cumsum(dataset_sizes)
            variable_names = np.array([x.decode('ascii', 'replace')
                                       for x in f.attrs['VariableNames'][:]])
            self.quantity_datasets = {}
            self.quantity_indices = {}
            for q in self.quantities:
                var_num = np.where(variable_names == q)[0][0]
                dataset_num = np.where(dataset_sizes_cumulative > var_num)[0][0]
                if dataset_num == 0:
                    dataset_index = var_num
                else:
                    dataset_index = var_num - dataset_sizes_cumulative[dataset_num-1]
                self.quantity_datasets[q] = dataset_names[dataset_num]
                self.quantity_indices[q] = dataset_index

            # Locate fine block for coordinates in case of slice
            fine_block = np.where(self.levels == self.max_level)[0][0]
            self.x1m = f['x1f'][fine_block, 0]
            self.x1p = f['x1f'][fine_block, 1]
            self.x2m = f['x2f'][fine_block, 0]
            self.x2p = f['x2f'][fine_block, 1]
            self.x3m = f['x3f'][fine_block, 0]
            self.x3p = f['x3f'][fine_block, 1]

            # Populate coordinate arrays
            face_funcs = (face_func_1, face_func_2, face_func_3)
            center_funcs = (center_func_1, center_func_2, center_func_3)
            for d, nx, face_func, center_func in zip(
                    range(1, 4), nx_vals, face_funcs, center_funcs):
                xf = 'x' + repr(d) + 'f'
                xv = 'x' + repr(d) + 'v'
                if nx == 1:
                    xm = (self.x1m, self.x2m, self.x3m)[d-1]
                    xp = (self.x1p, self.x2p, self.x3p)[d-1]
                    self[xf] = np.array([xm, xp])
                else:
                    xmin = f.attrs['RootGridX'+repr(d)][0]
                    xmax = f.attrs['RootGridX'+repr(d)][1]
                    xrat_root = f.attrs['RootGridX'+repr(d)][2]
                    if xrat_root == -1.0 and face_func is None:
                        raise AthenaError(
                                'Must specify user-defined face_func_{0}'.format(d))
                    elif face_func is not None:
                        self[xf] = face_func(xmin, xmax, xrat_root, nx + 1)
                    elif xrat_root == 1.0:
                        if np.all(self.levels == self.level):
                            self[xf] = np.empty(nx + 1)
                            for n_block in range(int(nx / self.block_size[d-1])):
                                sample_location = [0, 0, 0]
                                sample_location[d-1] = n_block
                                sample_block = np.where(np.all(
                                    self.logical_locations == sample_location,
                                    axis=1))[0][0]
                                index_low = n_block * self.block_size[d-1]
                                index_high = index_low + self.block_size[d-1] + 1
                                self[xf][index_low:index_high] = f[xf][sample_block, :]
                        else:
                            self[xf] = np.linspace(xmin, xmax, nx + 1)
                    else:
                        xrat = xrat_root ** (1.0 / 2 ** self.level)
                        self[xf] = (
                                xmin + (1.0 - xrat**np.arange(nx+1)) / (1.0 - xrat**nx)
                                * (xmax - xmin))
                self[xv] = np.empty(nx)
                for i in range(nx):
                    self[xv][i] = center_func(self[xf][i], self[xf][i+1])

            # Account for selection
            x1_select = False
            x2_select = False
            x3_select = False
            self.i_min = self.j_min = self.k_min = 0
            self.i_max = self.nx1
            self.j_max = self.nx2
            self.k_max = self.nx3
            error_string = '{0} must be {1} than {2} in order to intersect data range'
            if x1_min is not None and x1_min >= self['x1f'][1]:
                if x1_min >= self['x1f'][-1]:
                    raise AthenaError(error_string.format(
                            'x1_min', 'less', self['x1f'][-1]))
                x1_select = True
                self.i_min = np.where(self['x1f'] <= x1_min)[0][-1]
            if x1_max is not None and x1_max <= self['x1f'][-2]:
                if x1_max <= self['x1f'][0]:
                    raise AthenaError(
                            error_string.format('x1_max', 'greater', self['x1f'][0]))
                x1_select = True
                self.i_max = np.where(self['x1f'] >= x1_max)[0][0]
            if x2_min is not None and x2_min >= self['x2f'][1]:
                if x2_min >= self['x2f'][-1]:
                    raise AthenaError(error_string.format(
                            'x2_min', 'less', self['x2f'][-1]))
                x2_select = True
                self.j_min = np.where(self['x2f'] <= x2_min)[0][-1]
            if x2_max is not None and x2_max <= self['x2f'][-2]:
                if x2_max <= self['x2f'][0]:
                    raise AthenaError(
                            error_string.format('x2_max', 'greater', self['x2f'][0]))
                x2_select = True
                self.j_max = np.where(self['x2f'] >= x2_max)[0][0]
            if x3_min is not None and x3_min >= self['x3f'][1]:
                if x3_min >= self['x3f'][-1]:
                    raise AthenaError(error_string.format(
                            'x3_min', 'less', self['x3f'][-1]))
                x3_select = True
                self.k_min = np.where(self['x3f'] <= x3_min)[0][-1]
            if x3_max is not None and x3_max <= self['x3f'][-2]:
                if x3_max <= self['x3f'][0]:
                    raise AthenaError(
                            error_string.format('x3_max', 'greater', self['x3f'][0]))
                x3_select = True
                self.k_max = np.where(self['x3f'] >= x3_max)[0][0]

            # Adjust coordinates if selection made
            if x1_select:
                self['x1f'] = self['x1f'][self.i_min:self.i_max+1]
                self['x1v'] = self['x1v'][self.i_min:self.i_max]
            if x2_select:
                self['x2f'] = self['x2f'][self.j_min:self.j_max+1]
                self['x2v'] = self['x2v'][self.j_min:self.j_max]
            if x3_select:
                self['x3f'] = self['x3f'][self.k_min:self.k_max+1]
                self['x3v'] = self['x3v'][self.k_min:self.k_max]

            # Set to None any quantities not set by initialization
            for i in self.quantities:
                if i not in self:
                    self[i] = None

    # Function for contingently accessing data
    def __getitem__(self, item):
        if self._need_to_read(item):
            self._grab_quantities([item])
        return super(athdf, self).__getitem__(item)

    # Function for contingently setting data
    def __setitem__(self, key, value):
        if not self.new_data:
            try:
                self._existing_keys.remove(key)
            except ValueError:
                pass
        return super(athdf, self).__setitem__(key, value)

    # Function for returning constant shape of all 3D arrays
    def _shape(self):
        return (self.k_max - self.k_min, self.j_max - self.j_min, self.i_max - self.i_min)

    # Function for determining if a quantity must be set
    def _need_to_read(self, quantity):
        if not self.new_data:
            if quantity in self._existing_keys:
                return True
        try:
            if super(athdf, self).__getitem__(quantity) is None:
                return True
        except KeyError:
            if quantity in self.quantities:
                return True
        return False

    # Function for setting all needed quantities
    def _grab_quantities(self, quantities):

        # Create list of quantities to be set
        quantities = [q for q in quantities if self._need_to_read(q)]

        # Open file
        with h5py.File(self.filename, 'r') as f:

            # Prepare arrays for data and bookkeeping
            if self.new_data:
                for q in quantities:
                    self[q] = np.zeros((self._shape()), dtype=self.dtype)
                if self.return_levels:
                    self['Levels'] = np.empty((self._shape()), dtype=np.int32)
            else:
                for q in quantities:
                    self[q].fill(0.0)
            if (not self.subsample and not self.fast_restrict
                    and self.max_level > self.level):
                restricted_data = np.zeros((self.lx3, self.lx2, self.lx1), dtype=bool)

            # Go through blocks in data file
            for block_num in range(self.num_blocks):

                # Extract location information
                block_level = self.levels[block_num]
                block_location = self.logical_locations[block_num, :]

                # Prolongate coarse data and copy same-level data
                if block_level <= self.level:

                    # Calculate scale (number of copies per dimension)
                    s = 2 ** (self.level - block_level)

                    # Calculate destination indices, without selection
                    il_d = (block_location[0] *
                            self.block_size[0] * s) if self.nx1 > 1 else 0
                    jl_d = (block_location[1] *
                            self.block_size[1] * s) if self.nx2 > 1 else 0
                    kl_d = (block_location[2] *
                            self.block_size[2] * s) if self.nx3 > 1 else 0
                    iu_d = il_d + self.block_size[0] * s if self.nx1 > 1 else 1
                    ju_d = jl_d + self.block_size[1] * s if self.nx2 > 1 else 1
                    ku_d = kl_d + self.block_size[2] * s if self.nx3 > 1 else 1

                    # Calculate (prolongated) source indices, with selection
                    il_s = max(il_d, self.i_min) - il_d
                    jl_s = max(jl_d, self.j_min) - jl_d
                    kl_s = max(kl_d, self.k_min) - kl_d
                    iu_s = min(iu_d, self.i_max) - il_d
                    ju_s = min(ju_d, self.j_max) - jl_d
                    ku_s = min(ku_d, self.k_max) - kl_d
                    if il_s >= iu_s or jl_s >= ju_s or kl_s >= ku_s:
                        continue

                    # Account for selection in destination indices
                    il_d = max(il_d, self.i_min) - self.i_min
                    jl_d = max(jl_d, self.j_min) - self.j_min
                    kl_d = max(kl_d, self.k_min) - self.k_min
                    iu_d = min(iu_d, self.i_max) - self.i_min
                    ju_d = min(ju_d, self.j_max) - self.j_min
                    ku_d = min(ku_d, self.k_max) - self.k_min

                    # Assign values
                    for q in quantities:
                        dataset = self.quantity_datasets[q]
                        index = self.quantity_indices[q]
                        block_data = f[dataset][index, block_num, :]
                        if s > 1:
                            if self.nx1 > 1:
                                block_data = np.repeat(block_data, s, axis=2)
                            if self.nx2 > 1:
                                block_data = np.repeat(block_data, s, axis=1)
                            if self.nx3 > 1:
                                block_data = np.repeat(block_data, s, axis=0)
                        self[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = (
                                block_data[kl_s:ku_s, jl_s:ju_s, il_s:iu_s])

                # Restrict fine data
                else:

                    # Calculate scale
                    s = 2 ** (block_level - self.level)

                    # Calculate destination indices, without selection
                    il_d = (block_location[0] *
                            self.block_size[0] / s) if self.nx1 > 1 else 0
                    jl_d = (block_location[1] *
                            self.block_size[1] / s) if self.nx2 > 1 else 0
                    kl_d = (block_location[2] *
                            self.block_size[2] / s) if self.nx3 > 1 else 0
                    iu_d = il_d + self.block_size[0] / s if self.nx1 > 1 else 1
                    ju_d = jl_d + self.block_size[1] / s if self.nx2 > 1 else 1
                    ku_d = kl_d + self.block_size[2] / s if self.nx3 > 1 else 1

                    # Calculate (restricted) source indices, with selection
                    il_s = max(il_d, self.i_min) - il_d
                    jl_s = max(jl_d, self.j_min) - jl_d
                    kl_s = max(kl_d, self.k_min) - kl_d
                    iu_s = min(iu_d, self.i_max) - il_d
                    ju_s = min(ju_d, self.j_max) - jl_d
                    ku_s = min(ku_d, self.k_max) - kl_d
                    if il_s >= iu_s or jl_s >= ju_s or kl_s >= ku_s:
                        continue

                    # Account for selection in destination indices
                    il_d = max(il_d, self.i_min) - self.i_min
                    jl_d = max(jl_d, self.j_min) - self.j_min
                    kl_d = max(kl_d, self.k_min) - self.k_min
                    iu_d = min(iu_d, self.i_max) - self.i_min
                    ju_d = min(ju_d, self.j_max) - self.j_min
                    ku_d = min(ku_d, self.k_max) - self.k_min

                    # Account for restriction in source indices
                    if self.nx1 > 1:
                        il_s *= s
                        iu_s *= s
                    if self.nx2 > 1:
                        jl_s *= s
                        ju_s *= s
                    if self.nx3 > 1:
                        kl_s *= s
                        ku_s *= s

                    # Apply subsampling
                    if self.subsample:

                        # Calculate fine-level offsets (nearest cell at or below center)
                        o1 = s/2 - 1 if self.nx1 > 1 else 0
                        o2 = s/2 - 1 if self.nx2 > 1 else 0
                        o3 = s/2 - 1 if self.nx3 > 1 else 0

                        # Assign values
                        for q in quantities:
                            dataset = self.quantity_datasets[q]
                            index = self.quantity_indices[q]
                            self[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = (f[dataset][
                                index, block_num, kl_s+o3:ku_s:s,
                                jl_s+o2:ju_s:s, il_s+o1:iu_s:s]
                            )

                    # Apply fast (uniform Cartesian) restriction
                    elif self.fast_restrict:

                        # Calculate fine-level offsets
                        io_vals = range(s) if self.nx1 > 1 else (0,)
                        jo_vals = range(s) if self.nx2 > 1 else (0,)
                        ko_vals = range(s) if self.nx3 > 1 else (0,)

                        # Assign values
                        for q in quantities:
                            dataset = self.quantity_datasets[q]
                            index = self.quantity_indices[q]
                            for ko in ko_vals:
                                for jo in jo_vals:
                                    for io in io_vals:
                                        self[q][kl_d:ku_d,
                                                jl_d:ju_d,
                                                il_d:iu_d] += f[dataset][
                                                    index, block_num,
                                                    kl_s+ko:ku_s:s,
                                                    jl_s+jo:ju_s:s,
                                                    il_s+io:iu_s:s]
                            self[q][kl_d:ku_d, jl_d:ju_d,
                                    il_d:iu_d] /= s ** self.num_extended_dims

                    # Apply exact (volume-weighted) restriction
                    else:

                        # Calculate sets of indices
                        i_s_vals = range(il_s, iu_s)
                        j_s_vals = range(jl_s, ju_s)
                        k_s_vals = range(kl_s, ku_s)
                        i_d_vals = range(il_d, iu_d)
                        j_d_vals = range(jl_d, ju_d)
                        k_d_vals = range(kl_d, ku_d)
                        if self.nx1 > 1:
                            i_d_vals = np.repeat(i_d_vals, s)
                        if self.nx2 > 1:
                            j_d_vals = np.repeat(j_d_vals, s)
                        if self.nx3 > 1:
                            k_d_vals = np.repeat(k_d_vals, s)

                        # Accumulate values
                        for k_s, k_d in zip(k_s_vals, k_d_vals):
                            if self.nx3 > 1:
                                self.x3m = f['x3f'][block_num, k_s]
                                self.x3p = f['x3f'][block_num, k_s+1]
                            for j_s, j_d in zip(j_s_vals, j_d_vals):
                                if self.nx2 > 1:
                                    self.x2m = f['x2f'][block_num, j_s]
                                    self.x2p = f['x2f'][block_num, j_s+1]
                                for i_s, i_d in zip(i_s_vals, i_d_vals):
                                    if self.nx1 > 1:
                                        self.x1m = f['x1f'][block_num, i_s]
                                        self.x1p = f['x1f'][block_num, i_s+1]
                                    vol = self.vol_func(
                                            self.x1m, self.x1p, self.x2m, self.x2p,
                                            self.x3m, self.x3p)
                                    for q in quantities:
                                        dataset = self.quantity_datasets[q]
                                        index = self.quantity_indices[q]
                                        self[q][k_d, j_d, i_d] += (
                                            vol
                                            * f[dataset][index, block_num, k_s, j_s, i_s])
                        loc1 = (self.nx1 > 1) * block_location[0] / s
                        loc2 = (self.nx2 > 1) * block_location[1] / s
                        loc3 = (self.nx3 > 1) * block_location[2] / s
                        restricted_data[loc3, loc2, loc1] = True

                # Set level information for cells in this block
                if self.return_levels:
                    self['Levels'][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_level

        # Remove volume factors from restricted data
        if self.level < self.max_level and not self.subsample and not self.fast_restrict:
            for loc3 in range(self.lx3):
                for loc2 in range(self.lx2):
                    for loc1 in range(self.lx1):
                        if restricted_data[loc3, loc2, loc1]:
                            il = loc1 * self.block_size[0]
                            jl = loc2 * self.block_size[1]
                            kl = loc3 * self.block_size[2]
                            iu = il + self.block_size[0]
                            ju = jl + self.block_size[1]
                            ku = kl + self.block_size[2]
                            il = max(il, self.i_min) - self.i_min
                            jl = max(jl, self.j_min) - self.j_min
                            kl = max(kl, self.k_min) - self.k_min
                            iu = min(iu, self.i_max) - self.i_min
                            ju = min(ju, self.j_max) - self.j_min
                            ku = min(ku, self.k_max) - self.k_min
                            for k in range(kl, ku):
                                if self.nx3 > 1:
                                    self.x3m = self['x3f'][k]
                                    self.x3p = self['x3f'][k+1]
                                for j in range(jl, ju):
                                    if self.nx2 > 1:
                                        self.x2m = self['x2f'][j]
                                        self.x2p = self['x2f'][j+1]
                                    for i in range(il, iu):
                                        if self.nx1 > 1:
                                            self.x1m = self['x1f'][i]
                                            self.x1p = self['x1f'][i+1]
                                        vol = self.vol_func(
                                                self.x1m, self.x1p, self.x2m, self.x2p,
                                                self.x3m, self.x3p)
                                        for q in quantities:
                                            self[q][k, j, i] /= vol


# ========================================================================================


def restrict_like(vals, levels, vols=None):
    """Average cell values according to given mesh refinement scheme."""

    # Determine maximum amount of restriction
    nx3, nx2, nx1 = vals.shape
    max_level = np.max(levels)
    if nx3 > 1 and nx3 % 2**max_level != 0:
        raise AthenaError('x3-dimension wrong size to be restricted')
    if nx2 > 1 and nx2 % 2**max_level != 0:
        raise AthenaError('x2-dimension wrong size to be restricted')
    if nx1 % 2**max_level != 0:
        raise AthenaError('x1-dimension wrong size to be restricted')

    # Construct volume weighting
    if vols is None:
        vols = np.ones_like(vals)
    else:
        if vols.shape != vals.shape:
            raise AthenaError('Array of volumes must match cell values in size')

    # Restrict data
    vals_restricted = np.copy(vals)
    for level in range(max_level):
        level_difference = max_level - level
        stride = 2 ** level_difference
        if nx3 > 1:
            vals_level = np.reshape(vals * vols, (nx3/stride, stride, nx2/stride, stride,
                                                  nx1/stride, stride))
            vols_level = np.reshape(
                    vols, (nx3/stride, stride, nx2/stride, stride, nx1/stride, stride))
            vals_sum = np.sum(np.sum(np.sum(vals_level, axis=5), axis=3), axis=1)
            vols_sum = np.sum(np.sum(np.sum(vols_level, axis=5), axis=3), axis=1)
            vals_level = np.repeat(
                    np.repeat(np.repeat(vals_sum / vols_sum, stride, axis=0), stride,
                              axis=1), stride, axis=2)
        elif nx2 > 1:
            vals_level = np.reshape(vals * vols, (nx2/stride, stride, nx1/stride, stride))
            vols_level = np.reshape(vols, (nx2/stride, stride, nx1/stride, stride))
            vals_sum = np.sum(np.sum(vals_level, axis=3), axis=1)
            vols_sum = np.sum(np.sum(vols_level, axis=3), axis=1)
            vals_level = np.repeat(
                    np.repeat(vals_sum / vols_sum, stride, axis=0), stride, axis=1)
        else:
            vals_level = np.reshape(vals * vols, (nx1/stride, stride))
            vols_level = np.reshape(vols, (nx1/stride, stride))
            vals_sum = np.sum(vals_level, axis=1)
            vols_sum = np.sum(vols_level, axis=1)
            vals_level = np.repeat(vals_sum / vols_sum, stride, axis=0)
        vals_restricted = np.where(levels == level, vals_level, vals_restricted)
    return vals_restricted


# ========================================================================================


def athinput(filename):
    """Read athinput file and returns a dictionary of dictionaries."""

    # Read data, removing comments, extra whitespace, and empty lines
    with open(filename, 'r') as athinput:
        lines = filter(None, [i.split('#')[0].strip() for i in athinput.readlines()])

    # Split into blocks, first element will be empty
    blocks = ('\n'.join(lines)).split('<')[1:]

    # Function for interpreting strings numerically
    def typecast(x):
        try:
            return int(x)
        except ValueError:
            pass
        try:
            return float(x)
        except ValueError:
            pass
        try:
            return complex(x)
        except ValueError:
            pass
        return x

    # Function for parsing assignment based on first '='
    def parse_line(line):
        out = [i.strip() for i in line.split('=')]
        out[1] = '='.join(out[1:])
        out[1] = typecast(out[1])
        return out[:2]

    # Assign values into dictionaries
    data = {}
    for block in blocks:
        info = list(filter(None, block.split('\n')))
        key = info.pop(0)[:-1]  # last character is '>'
        data[key] = dict(map(parse_line, info))
    return data


# ========================================================================================


class AthenaError(RuntimeError):
    """General exception class for Athena++ read functions."""
    pass


class AthenaWarning(RuntimeWarning):
    """General warning class for Athena++ read functions."""
    pass
