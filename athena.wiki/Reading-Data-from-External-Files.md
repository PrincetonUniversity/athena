Some applications require having on-hand arrays of precomputed values, such as equation-of-state tables for non-ideal gases or opacity tables for computing radiation transfer. The `src/inputs/` directory contains functions designed to facilitate reading such data into Athena. For example, these functions can be called during problem set up to load data into user-allocated arrays.

Currently, only HDF5 data files are supported.

### HDF5

`hdf5_reader.hpp` declares the following function:

#### `void HDF5ReadRealArray(...)`

- Inputs:
  - `const char *filename`: name of data file, including any necessary path.
  - `const char *dataset_name`: name of dataset to load from, including any necessary path (i.e. HDF5 groups).
  - `int rank_file`: number of dimensions in file dataset.
  - `const int *start_file`: array of length `rank_file` indexing starting position (starting from 0's) from which to read data.
  - `const int *count_file`: array of length `rank_file` giving number of consecutive elements to read in each dimension. The product of these numbers will be the total number of elements read.
  - `int rank_mem`: number of dimensions in `array`.
  - `const int *start_mem`: array of length `rank_mem` indexing starting position (starting from 0's) where read data should be placed.
  - `const int *count_mem`: array of length `rank_mem` giving number of consecutive elements to read in each dimension. The product of these numbers must match the product of the `rank_file` elements of `count_file`.
  - `AthenaArray<Real> &array`: AthenaArray, already allocated, into which results should be placed.
  - `bool collective=false`: When set to `true` and the code is compiled with MPI, the file reads will be done collectively. By default, reads will be serial and thus safe (though possibly slow when running on many MPI ranks). If using collective reads, each MPI rank (i.e. Mesh) must make the same number of calls to this function.
  - `bool noop=false`: When set to `true`, no data will be read. This is intended to be used with collective calls in cases where different ranks would otherwise make different numbers of calls to this function.
- Outputs: 
  - `array` will contain the values from the file.
- Notes:
  - This function is designed to read floating-point data.
  - The data in the file will be cast to the length of `Real` when placed in the output array.
  - The function will automatically determine the extent of each dimension for both the input dataset and the output array.
- Example: Suppose `example.hdf5` contains a 12&times;15 array named `vals` in group `data_group` from which we want to extract just the last 2 rows into a 2&times;15 array. A problem generator could do this based on the following snippet:
```C++
#include "../inputs/hdf5_reader.hpp"
int start_file[2] = {10, 0};
int count_file[2] = {2, 15};
int start_mem[2] = {0, 0};
int count_mem[2] = {2, 15};
AthenaArray<Real> array;
array.NewAthenaArray(2, 15);
HDF5ReadRealArray("example.hdf5", "/data_group/vals", 2, start_file, count_file,
    2, start_mem, count_mem, array);
```