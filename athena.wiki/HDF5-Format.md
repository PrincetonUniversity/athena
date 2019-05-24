Athena++ uses its own HDF5 format when writing `.athdf` files. Below are the full specifications of the file format within the HDF5 framework. Both 32-bit and 64-bit numeric datatypes (ints and floats) are used, and they are always stored in big-endian order. "Strings" here refer to variable length char arrays with at most 20 characters (excluding terminating null).

### File-level attributes

- `NumCycles`
  - scalar 32-bit int
  - cycle number at which file is written
- `Time`
  - scalar 64-bit float
  - simulation time at which file is written
- `Coordinates`
  - scalar string
  - name of coordinate system (matches configure option and `.cpp` file)
- `RootGridX1`/`RootGridX2`/`RootGridX3`
  - triples of 64-bit floats
  - minimum, maximum, and geometric ratio in x1/x2/x3-direction
- `RootGridSize`
  - triple of 32-bit ints
  - numbers of cells at root level in x1-, x2-, and x3-directions
- `NumMeshBlocks`
  - scalar 32-bit int
  - total number of MeshBlocks in the simulation
- `MeshBlockSize`
  - triple of 32-bit ints
  - numbers of cells in each MeshBlock in x1-, x2-, and x3-directions
- `MaxLevel`
  - scalar 32-bit int
  - highest level of mesh refinement present, with root being 0
- `NumVariables`
  - array of 32-bit ints
  - length is number of cell-centered datasets in file
  - each entry is number of scalar variables in corresponding dataset
- `DatasetNames`
  - array of strings
  - length matches that of `NumVariables`
  - each entry is name of cell-centered dataset
  - order matches that of `NumVariables`
- `VariableNames`
  - array of strings
  - length equals sum of entries in `NumVariables`
  - each entry is name of cell-centered variable
  - order matches that of `DatasetNames`, with variables ordered by their index within the dataset

### Datasets
Let "NBlocks" be the value of the `NumMeshBlocks` attribute. Let "nx1," "nx2," and "nx3" be the values in the `MeshBlockSize` attribute.

- `Levels`
  - (NBlocks) array of 32-bit ints
  - Refinement levels of MeshBlocks, with root being 0
- `LogicalLocations`
  - (NBlocks)&times;(3) array of 64-bit ints
  - For each MeshBlock, the offsets in the x1-, x2-, and x3-directions from the minimum edge of the grid
  - Counting is done as though entire grid is at the refinement level of the MeshBlock in question
- `x1f`/`x2f`/`x3f`
  - (NBlocks)&times;(nx1/nx2/nx3+1) array of 32-bit floats
  - Values of interface locations along x1/x2/x3-direction
- `x1v`/`x2v`/`x3v`
  - (NBlocks)&times;(nx1/nx2/nx3) array of 32-bit floats
  - Values of cell centers along x1/x2/x3-direction
- Cell-centered datasets
  - One for each entry in `DatasetNames`
  - Each one is an (NVars)&times;(NBlocks)&times;(nx3)&times;(nx2)&times;(nx1) array of 32-bit floats
  - NVars varies between datasets and is given by corresponding entry in `NumVariables`

### Names of quantities
The following names of datasets and variables may be output depending on what is requested via the `variable` argument in the `<output>` block. (Variables and datasets will not be included if they are not included in the simulation due to the selected physics.)

<table>
  <tr> <th>Output<br>Variable</th> <th>Dataset<br>Names</th> <th>Variable<br>Names</th> </tr>
  <tr> <td rowspan="2">prim</td> <td>prim</td> <td>rho, press, vel1, vel2, vel3</td> </tr>
  <tr> <td>B</td> <td>Bcc1, Bcc2, Bcc3</td> </tr>
  <tr> <td rowspan="2">cons</td> <td>cons</td> <td>dens, Etot, mom1, mom2, mom3</td> </tr>
  <tr> <td>B</td> <td>Bcc1, Bcc2, Bcc3</td> </tr>
  <tr> <td>d</td> <td>hydro</td> <td>rho</td> </tr>
  <tr> <td>p</td> <td>hydro</td> <td>press</td> </tr>
  <tr> <td>v</td> <td>hydro</td> <td>vel1, vel2, vel3</td> </tr>
  <tr> <td>D</td> <td>hydro</td> <td>dens</td> </tr>
  <tr> <td>E</td> <td>hydro</td> <td>Etot</td> </tr>
  <tr> <td>m</td> <td>hydro</td> <td>mom1, mom2, mom3</td> </tr>
  <tr> <td>bcc</td> <td>B</td> <td>Bcc1, Bcc2, Bcc3</td> </tr>
  <tr> <td>uov</td> <td>uov</td> <td>user_out_var0, user_out_var1, ...</td> </tr>
</table>

### HDF5 library compatibility
When [[Configuring]] Athena+++ with both the `-mpi` and `-hdf5` flags, it is important that the HDF5 library was built with MPI support. The solver will not compile if only the serial HDF5 routines are available.

If you are building HDF5 from source, the HDF5 library source code contains [both serial HDF5 and Parallel HDF5 (PHDF5)](https://portal.hdfgroup.org/pages/viewpage.action?pageId=48809596). Simply specify the MPI compiler when configuring the installer and pass the following flag:
```
CC=mpicc ./configure --enable-parallel ...
```
A POSIX compliant file system and an MPI library with MPI-I/O are requirements for Parallel HDF5. A parallel file system is necessary for good performance. To check if an installed HDF5 library was linked with the MPI library, check the output of `h5cc -showconfig`.  

On the Princeton University [Research Computing (RC) clusters](https://researchcomputing.princeton.edu/systems-and-services/available-systems), the following compilers and libraries are currently recommended for compiling Athena++ with MPI and Parallel HDF5:
```
intel/18.0/64/18.0.2.199
intel-mpi/intel/2018.2/64
hdf5/intel-17.0/intel-mpi/1.10.0
```
They are loaded using the [Environment Modules](http://modules.sourceforge.net/) package via `module load intel-mpi/intel/2018.2/64`, for example.

Note, there is no need to load an HDF5 module for serial HDF5 on the RC clusters; the library is already loaded in the compiler and linker search paths.

*Last updated on 6/9/18*
