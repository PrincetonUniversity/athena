### Overview

The file [`athena/vis/python/athena_read.py`][1] contains functions for reading Athena++ output data into a Python script. This guide documents these functions and how to call them.

### Importing the module

A Python script must import a file as a module before using it. This can be a bit unintuitive for files located in arbitrary locations. In your script, use the standard `sys` module to add the `athena/vis/python` directory to Python's path, then load the desired module:
```Python
import sys
sys.path.insert(0, '<path/to/>athena/vis/python')
import athena_read
```

### Available functions

Currently, there are readers for HDF5, VTK, tabular, and history output, all tailored for Athena++. These read a single file at a time, so if for instance you have many VTK files you can either read them one at a time or join them into a single file using the [join script][2]. The HDF5 reader supports files produced with mesh refinement. The function specifications are below.

##### `athdf(filename, ...)`

- Inputs:
  - `filename`: string containing the path to the desired Athena++ HDF5 output file to be read.
  - `raw` (optional, default `False`): flag indicating data should be returned in its raw state (i.e. with blocks not joined into a contiguous array). Overrides all other options. Output dictionary will collate all cell-centered variables into 4D arrays (with the first index corresponding to the MeshBlock).
  - `data` (optional, default `None`): dictionary similar in structure to output, containing arrays of the correct size to be overwritten by the read-in data. Values will only be written corresponding to existing keys, overriding `quantities` argument. If `None`, a new dictionary is created.
  - `quantities` (optional, default `None`): list of strings naming desired quantities to be read. The interface locations `x1f`, `x2f`, and `x3f` will be read automatically and need not be given explicitly. Similarly the cell centers `x1v`, `x2v`, and `x3v` will be calculated automatically. Overridden by keys present in `data` argument if given.
  - `dtype` (optional, default `np.float32`): data type for cell data (excludes coordinate arrays). The default is designed to match the type written to the HDF5 file by Athena++. Overridden if `data` argument given.
  - `level` (optional, default `None`): nonnegative integer specifying refinement level relative to root grid (defined to be 0) to which data should be cast. All data must be cast to the same level. See notes for how this is accomplished. If `None`, the maximum refinement level present in the file is used.
  - `return_levels` (optional, default `False`): Boolean indicating the returned dictionary should contain the array `'Levels'`, which has the same shape as the other cell data arrays and contains the refinement level of each cell.
  - `subsample` (optional, default `False`): Boolean indicating restriction should be accomplished via subsampling (defined in notes). If `True`, `fast_restrict` and `vol_func` are ignored.
  - `fast_restrict` (optional, default `False`): Boolean indicating restriction should be accomplished via an approximate fast method. Ignored if `subsample` is `True`. If `True`, `vol_func` is ignored.
  - `x1_min`/`x1_max`/`x2_min`/`x2_max`/`x3_min`/`x3_max` (optional, default `None`): Numeric values indicating what selection of the full domain to return. The returned data will be truncated at these values (actually at the smallest cell boundaries at the desired refinement level containing these values) if given. For example, a constant x2-slice can be extracted from a 3D dataset by specifying `x2_min` and `x2_max` to be the same value.
  - `vol_func` (optional, default `None`): function taking 6 numeric inputs and returning a numeric output. The inputs should be the lower and upper x1 coordinates, lower and upper x2 coordinates, and lower and upper x3 coordinates of a constant-coordinate-surface region, and the output should be the volume of the region. Used only for exact restriction; ignored if `subsample` or `fast_restrict` is `True`. If `subsample` and `fast_restrict` are `False` and `vol_func` is `None`, the appropriate volume function will be used based on the coordinates specified in the file.
  - `vol_params` (optional, default `None`): arguments to be passed into predefined volume functions, in cases where data not stored in the file is needed. Presently this only applies to Kerr-Schild coordinates, where `vol_params[0]` is expected to be the spin of the black hole (with the same units as mass). Only used if `subsample` and `fast_restrict` are `False`, `vol_func` is `None`, and the file specifies `'kerr-schild'` for its `Coordinates` attribute.
  - `face_func_1`/`face_func_2`/`face_func_3` (optional, default `None`): functions taking three floats and an integer and returning an array of floats. The inputs should be the minimum value of the coordinate, the maximum value, the cell width ratio (at the root grid level), and number of faces (at the output level), and the output should be an array containing the face locations. If `None`, the appropriate uniform or geometric ratio spacing will be used.
  - `center_func_1`/`center_func_2`/`center_func_3` (optional, default `None`): functions taking 2 numeric inputs and returning a numeric output. The inputs should be the lower and upper coordinates of a cell in the appropriate direction, and the output should be the appropriate coordinate of the cell center. If `None`, the appropriate functions will be used based on the coordinates specified in the file.
  - `num_ghost` (optional, default `0`): Number of ghost cells in HDF5 file data (i.e., the configure parameter `--nghost`, but only nonzero of the output was produced with the input parameter `ghost_zones` set to `true`). Must be given if file has ghost zones, and must be `0` if file does not have ghost zones. Array will only have one copy (not necessarily from the active zone) of any cells belonging to internal ghost zones.
- Outputs: Dictionary with the following entries:
  - `'x1f'`/`'x2f'`/`'x3f'`: 1D arrays of `nx1+1`/`nx2+1`/`nx3+1` x1/x2/x3-interfaces.
  - `'x1v'`/`'x2v'`/`'x3v'`: 1D arrays of `nx1`/`nx2`/`nx3` x1/x2/x3-coordinates of cell centers.
  - 3D array for each desired quantity, with indices ordered k, j, i (x3, x2, x1). Note all data is scalar data. Dictionary keys correspond to array names as stored in the HDF5 file.
  - Each file-level attribute. For example, the key `'Time'` will correspond to a float containing the simulation time, and the key `'VariableNames'` will correspond to an array of strings.
- Notes:
  - This function requires the `h5py` module to be installed.
  - Where coarse data needs to be prolongated to a finer output refinement level, values in coarse cells are copied to every fine cell overlapping in volume.
  - Finer data will be restricted according to the specified method:
    - Subsampling: Each output coarse cell will copy values from a single fine cell from the data file, the chosen fine cell being closest (on the low side) to the logical center of the coarse cell. For example, if a block in the data file is 3 levels (a factor of 8) too refined, its cells with (0-indexed) indices 3, 11, 17, ... will be copied. This is generally the fastest restriction method. If the maximum amount of coarsening needed is k refinement levels, this method requires all nonsingleton MeshBlock sizes be divisible by 2<sup>k</sup>; equivalently the most refined MeshBlock boundaries in the data file must lie on cell boundaries in the output grid.
    - Fast restriction: Each output coarse cell will be an unweighted average of all fine cells overlapping it. This is equivalent to exact, volume-weighted restriction in the case of uniform grid spacing and Cartesian coordinates. This is generally slightly slower than subsampling, but much faster than exact, volume-weighted averaging. If the maximum amount of coarsening needed is k refinement levels, this method requires all nonsingleton MeshBlock sizes be divisible by 2<sup>k</sup>; equivalently the most refined MeshBlock boundaries in the data file must lie on cell boundaries in the output grid.
    - Exact, volume-weighted restriction: Each output coarse cell will be a volume-weighted average of the fine cells overlapping it. **Caution:** this method can be extremely slow.
- Example: Suppose one wants to read the conserved quantities from `my_output.athdf`. The following code snippet does this, describing the total energy and x-momentum in the cell with coordinates `i=10`, `j=20`, and `k=30`.
```Python
import sys
sys.path.insert(0, '<path/to/>athena/vis/python')
import athena_read
data = athena_read.athdf('my_output.athdf')
print('x-bounds: {0} to {1}'.format(data['x1f'][10], data['x1f'][11]))
print('y-bounds: {0} to {1}'.format(data['x2f'][20], data['x2f'][21]))
print('z-bounds: {0} to {1}'.format(data['x3f'][30], data['x3f'][31]))
print('total energy: {0}'.format(data['Etot'][30,20,10]))
print('x-momentum: {0}'.format(data['mom1'][30,20,10]))
```
Note that Athena++'s HDF5 format stores even vector quantities as separate scalars.

##### `vtk(filename)`

- Inputs:
  - `filename`: string containing the path to the desired Athena++ VTK output file to be read.
- Outputs:
  1. 1D array of `nx1+1` x-interfaces.
  2. 1D array of `nx2+1` y-interfaces.
  3. 1D array of `nx3+1` z-interfaces.
  4. Dictionary of quantities defined on grid, given as either 3D arrays (for scalars; indices are ordered k, j, i, that is x3, x2, x1) or 4D arrays (for vectors; last index indexes the component of the vector). Dictionary keys correspond to array names as stored in the VTK file.
- Example: Suppose one wants to read the conserved quantities from `my_output.vtk`. The following code snippet does this, describing the total energy and x-momentum in the cell with coordinates `i=10`, `j=20`, and `k=30`:
```Python
import sys
sys.path.insert(0, '<path/to/>athena/vis/python')
import athena_read
x,y,z,data = athena_read.vtk('my_output.vtk')
print('x-bounds: {0} to {1}'.format(x[10], x[11]))
print('y-bounds: {0} to {1}'.format(y[20], y[21]))
print('z-bounds: {0} to {1}'.format(z[30], z[31]))
print('total energy: {0}'.format(data['Etot'][30,20,10]))
print('x-momentum: {0}'.format(data['mom'][30,20,10,0]))
```
Note that Athena++'s VTK format stores vector quantities in 4D arrays, with the last index corresponding to the vector component.

##### `tab(filename, raw=False, dimensions=None)`

- Inputs:
  - `filename`: string containing the path to the desired Athena++ tabular output file to be read.
  - `raw` (optional, default `False`): flag indicating any metadata and header data should be ignored, returning an appropriately shaped array containing all floating point cell data. Must specify `dimensions` in this case.
  - `dimensions` (optional, default `None`): integer specifying if the file has 1D, 2D, or 3D data. Only used if `raw` flag is `True`.
- Outputs: If `raw` is `True`, a 2D/3D/4D numpy array containing the data from the file, excluding the cell index columns, indexed with column number, k (if data is 3D), j (if data is 2D or 3D), and i. Otherwise a dictionary whose keywords are the headers found in the file and whose values are the corresponding 1D/2D/3D arrays read from the file, again in k, j, i (x3, x2, x1) order. In this latter case, the dictionary will also include entries for the file `time`, `cycle`, and `variables` metadata.
- Example: Suppose one wants to read the 1D primitive hydrodynamics data from `my_output.tab`. The following code snippet does this, printing the values from the third column in the file:
```Python
import sys
sys.path.insert(0, '<path/to/>athena/vis/python')
import athena_read
data = athena_read.tab('my_output.tab')
data_raw = athena_read.tab('my_output.tab', raw=True, dimensions=1)
print(data['rho'])
print(data_raw[1,:])
```

##### `hst(filename, raw=False)`

- Inputs:
  - `filename`: string containing the path to the desired Athena++ history output file to be read.
  - `raw` (optional, default `False`): flag indicating all data should be returned, even if the times are not monotonically increasing (as for instance might happen after a restart). If `False`, any time a simulation time is found that is less than or equal to the previous time, all previous data with times less than or equal to this time is ignored.
- Outputs: A dictionary whose keys are the column headings in the file and whose values are the corresponding 1D arrays of numbers.
- Example: Suppose one wants to report the total conserved mass four history time intervals after the start of the simulation. This might be done as:
```Python
import sys
sys.path.insert(0, '<path/to/>athena/vis/python')
import athena_read
data = athena_read.hst('my_output.hst')
print('mass at time {0}: {1}'.format(data['time'][4], data['mass'][4]))
```

##### `restrict_like(vals, levels, vols=None)`

- Inputs:
  - `vals`: array of cell values, assumed to be cast to maximum refinement level in file.
  - `levels`: array of refinement levels, matching `vals` in size.
  - `vols` (optional, default `None`): if given, an array of cell volumes matching `vals` in size. If omitted, all cell volumes are assumed identical (as in uniform Cartesian grids).
- Outputs: Array of the same size as input in which values have been averaged over cells in the refinement scheme specified by `levels`.
- Example: Suppose a simulation is run both with refinement and at a uniform resolution matching the highest refinement level. In order to properly evaluate differences between these results, either the refinement results should be prolongated in a nontrivial (and not built-in) way to the uniform resolution, or else the uniform data must be restricted in a volume-weighted way to the grid used with refinement. This latter option is done below:
```Python
import sys
sys.path.insert(0, '<path/to/>athena/vis/python')
import athena_read
data_refinement = athena_read.athdf('refinement.athdf', return_levels=True)
data_uniform = athena_read.athdf('uniform.athdf')
rho_uniform_restricted = athena_read.restrict_like(data_uniform['rho'], data_refinement['Levels'])
rho_error = np.mean(abs(data_refinement['rho'] - rho_uniform_restricted))
```

  [1]: https://github.com/PrincetonUniversity/athena-public-version/blob/master/vis/python/athena_read.py
  [2]: https://github.com/PrincetonUniversity/athena-public-version/blob/master/vis/vtk/join_vtk%2B%2B.c