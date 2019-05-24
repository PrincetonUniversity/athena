### `<output[n]>` Block in Input File
Data output in Athena++ is controlled by the `<output[n]>` blocks in the input file, where `[n]` is an integer. There is no limit on the total number of outputs. Here is an example:
```
    <output1>
    file_type   = vtk       # VTK data dump
    variable    = prim      # variables to be output
    id          = primitive # file identifier
    dt          = 0.01      # time interval between outputs
    x3_slice    = 0.0       # slicing at z=0.0
    ghost_zones = 1         # include ghost zones
```

#### Description of Output Block Parameters
* `file_type` (mandatory) : output format (see below)
* `variable` (depends; mandatory, except for `.rst` and `.hst`): output variable (see below)
* `dt` (mandatory) : time interval between outputs
* `data_format` (optional; used only for `.tab` and `.hst`) : format specifier string used for writing data (e.g. `%12.5e`)
* `id` (optional) : output ID string used in file names (default = `out[n]`)
* `x[1,2,3]_slice` (optional) : sliced output at specified position
* `x[1,2,3]_sum` (optional) : summed output along axis (default = `false`)

* `ghost_zones` (optional) : include ghost zones in output (default = `false`)
* `cartesian_vector` (optional) : add vector variables converted into Cartesian coordinates (default = `false`). This feature currently only works for spherical-polar and cylindrical coordinate systems, and it is able to transform the `v`, `m`, and `bcc` vector variables below. Face-centered magnetic field components are treated as separate scalars when written as outputs, because each component is represented at separate spatial locations (longitudinal cell faces). 
* `next_time` (optional) : override the time for the end first output interval `dt`

#### Output Variables
Either a single variable or a set (`cons` vs. `prim` vs. ...) of variables can be specified in the variable field. However, only one of the following (case-sensitive) values can be used in the `variable` field for a single output block:
* `prim` : primitive variables (density, velocities, pressure, and cell-centered magnetic fields)
* `cons` : conserved variables (density, momenta, total energy, and cell-centered magnetic fields)
* `d` : primitive (rest-frame) gas density
* `D` : conservative (lab-frame) gas density (the same as `d` in non-relativistic simulations)
* `p` : gas pressure
* `E` : total energy
* `v` : gas velocity vector
  * `v1`, `v2`, `v3` : gas velocity components
  * `vx`, `vy`, `vz` : gas velocity components
* `m` : gas momentum vector 
  * `m1`, `m2`, `m3` : gas momenta components
* `bcc` : (cell-centered) magnetic field vector
  * `bcc1`, `bcc2`, `bcc3` : (cell-centered) magnetic field components
* `b` : (face-centered) magnetic field vector
  * `b1`, `b2`, `b3` : (face-centered) magnetic field components
* `phi` : gravitational potential (see [[Self-Gravity with FFT]])
  * Also automatically included in `prim` and `cons` outputs if self-gravity is enabled
* `rN` (integer `0<=N<NSCALARS`) : dimensionless passive scalar concentration (in [0, 1]) of the `N`th species (see [[Passive Scalars]])
  * Also automatically included in `prim` outputs if `NSCALARS > 0` 
* `sN` (integer `0<=N<NSCALARS`) : passive scalar density of the `N`th species
  * Also automatically included in `cons` outputs if `NSCALARS > 0` 
* `uov` (user_out_var) : user output variables (see below)

Note, these variables names often do not match the `name` field that is assigned to the variable when writing to file. For example, writing an HDF5 file with `cons` variables would include an HDF5 dataset named `dens` containing the lab-frame density; if `variable=prim` the rest-frame density would be written with the name `rho`. Neither dataset name matches the corresponding individual `output[n]/variable` option, `D` and `d`. Similarly, if `cartesian_vector=true` in a spherical-polar or cylindrical coordinates simulation, the output of the cell-centered magnetic field components will be labeled as `Bcc_xyz1`, `Bcc_xyz2`, `Bcc_xyz3`.

### Output Formats
Currently the following output formats are supported.
* History File (`file_type = hst`)*^
* Formatted Table (`tab`)^
* VTK (`vtk`)
* HDF5 (`hdf5`)*
* Restart (`rst`)*

\* = parallel output

\^ = Plain text with ASCII encoding

#### History File
A history file (`.hst`) is a special file that contains global sums of basic physical variables such as the mass, momenta, and energy. The format is a simple text table. This file is useful for roughly monitoring the simulation as well as a sanity check of the code. By default, it includes the time, time step, total mass, momenta, kinetic energies in three directions, total energy, and magnetic energies in three directions. Users can add user-defined variables (see below).

#### Formatted Table
This is a simple text output, and is useful for relatively small 1D and 2D simulations. The data can be visualized using [[gnuplot|http://www.gnuplot.info/]] or similar plotting software. Because the data size tends to become large, this format is not recommended for large 3D simulations. The file extension is `.tab`.

#### VTK
VTK (Visualization Tool Kit) is a standard data format commonly used in numerical simulations. It can be easily visualized using various visualization software, such as [[VisIt|https://visit.llnl.gov/]] and [[ParaView|http://www.paraview.org/]]. The file extension is `.vtk`.

#### HDF5
Athena++ can output files formatted according to the HDF5 (Hierarchical Data Format) standard. This format is the most suitable one for simulations with mesh refinement. It also has the advantage of working with parallel IO, meaning more than one process can write into a single file, scaling well when used with a parallel file system like GPFS.

This output generates two files per output timestep: `.athdf` and `.athdf.xdmf`. The `.athdf` file is HDF5 and contains the data. The `.athdf.xdmf` (eXtensible Data Model and Format) file is an auxiliary file containing a description of the data structure stored in the `.athdf` file. This file is used when visualizing the results using VisIt or ParaView. Also, [[yt|http://yt-project.org/]] supports the Athena++ HDF5 output.

In order to enable this option, the HDF5 library is needed and the code must be configured with the `-hdf5` option. Also note that a parallel version of the HDF5 library is required for parallel simulations with MPI. By default, the floating point output data is written in 32-bit single precision, corresponding to type `H5T_NATIVE_FLOAT`. To create `H5T_NATIVE_DOUBLE` double precision output, use `-h5double` when configuring. 

HDF5 files are self-describing, rather than adhering to a single, exact layout. One can view the structure using a tool like [[HDFView|https://portal.hdfgroup.org/display/support/Download+HDFView]] provided by the HDF Group. Athena++ has its own particular [[HDF5 format|HDF5-Format]].

#### Restart File
A restarting file (`.rst`) is used literally for restarting simulations. For example, on many supercomputing systems, the time allowed for one run is limited, and you may want to continue your simulations over multiple runs. This is a parallel format using MPI-IO, so only one file per output timestep is generated but no external library other than MPI is required. Because this file necessarily contains all the information in the simulation, the `variable` field is ignored. 

When the code is launched with the `-t hh:mm:ss` option and a `rst` formatted output specified, a restarting file named `*.final.rst` that contains the data at the end of the simulation is automatically created if/when the solver reaches this time limit. Allow some time for the final output as the code makes the final output _after_ the specified time limit is reached. This option can be useful in conjunction with a time limit imposed by job schedulers on shared clusters. E.g. if Athena++ is started with a Slurm time limit of `sbatch --t=6:00:00 ...`, it may be valuable to set a separate safety time limit of `./athena -t 5:50:00`. If the solver is at risk of being killed by the job scheduler that does not provide any "grace time" beyond the given allocation, this ensures that a final restart file is generated and no progress is lost. However, for some very large simulations 10 minutes may be insufficient to write the final restart file.

To restart a simulation, simply specify a restarting file using the `-r` option instead of `-i`
```
    > athena -r example.out1.00010.rst
```
Note that you can change the number of MPI processes when you restart, as long as the number of MPI processes is fewer than that of MeshBlocks. This is useful for AMR simulations, when many refinement levels are created as simulations proceed.

### File Names
For serial output formats (Formatted table and VTK), one file is created for each MeshBlock per every output timestep. Output filenames follow the convention
```
    [problem_id].block[block_id].[output_id].?????.[type]
```
where problem_id is the Problem ID specified in the input file, `[block_id]` is the global ID of MeshBlock, `[output_id]` is the output ID, `?????` is the output step, and `[type]` is the file extension for the specified output format.

For parallel output formats (HDF5), only one file is created per every output timestep, no matter how many MeshBlocks exist. This is convenient for massively parallel simulations. Output filenames follow the convention
```
    [problem_id].[output_id].?????.[type]
```
And the restarting output uses the following file names:
```
    [problem_id].?????.rst
```
Please note that these file names are different from those in the old Athena code. In particular, now MeshBlock ID is used instead of process ID in Athena.

### User Output Variables
User Output Variables (`user_out_var` in the code) are user-defined variables associated with each cell. To enable these variables, specify the number of fields in the `MeshBlock::InitUserMeshBlockData` function in your problem generator file:
```c++
    void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
    {
      AllocateUserOutputVariables(2);
      return;
    }
```
This function sets the class data member `MeshBlock::nuser_out_var` to be 2, and allocates a 4-dimensional `AthenaArray<Real>`user_out_var`. Then, specify an output block in your input file:
```
    <output2>
    file_type   = hdf5
    dt          = 0.1
    variable    = uov  # (or user_out_var)
```
Store the data you want to output in the `MeshBlock::UserWorkBeforeOutput` function in your problem generator file (note: formerly `UserWorkInLoop` was used for this purpose, but now this new dedicated interface is provided because calculation of `uov` can be expensive). For example, in order to calculate and save the gas temperature and plasma beta:
```c++
    void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
    {
      for(int k=ks; k<=ke; k++) {
        for(int j=js; j<=je; j++) {
          for(int i=is; i<=ie; i++) {
            Real pmag = 0.5*(SQR(pfield->bcc(IB1,k,j,i))
                            +SQR(pfield->bcc(IB2,k,j,i))
                            +SQR(pfield->bcc(IB3,k,j,i)));
            user_out_var(0,k,j,i) = phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i);
            user_out_var(1,k,j,i) = phydro->w(IPR,k,j,i)/pmag;
          }
        }
      }
    }
```
Note that if you want to output the data in the boundary ghost cells, the loop limits should be from `ie-NGHOST` to `ie+NGHOST`, etc. The `user_out_var` array will be automatically destroyed at the end of the simulation, so you do not have to delete it.

Also, this function is called before output, except the history output because the history output is intended for more frequent and quick analysis. To add user-defined variables in the history file, please see the next section.

### User-Defined History Output
Users can add any variable in the history output. First, define a function to calculate a sum of user-defined variables within a MeshBlock. The history output automatically calculates a global sum using this function. The example below calculates divergence of magnetic fields (note: Athena++'s constrained transport scheme should conserve &nabla;&#8901;**B** = 0 to machine precision as long as the code and the initial condition are correct).

```c++
Real DivergenceB(MeshBlock *pmb, int iout)
{
  Real divb=0;
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> face1, face2p, face2m, face3p, face3m;
  FaceField &b = pmb->pfield->b;

  face1.NewAthenaArray((ie-is)+2*NGHOST+2);
  face2p.NewAthenaArray((ie-is)+2*NGHOST+1);
  face2m.NewAthenaArray((ie-is)+2*NGHOST+1);
  face3p.NewAthenaArray((ie-is)+2*NGHOST+1);
  face3m.NewAthenaArray((ie-is)+2*NGHOST+1);

  for(int k=ks; k<=ke; k++) {
    for(int j=js; j<=je; j++) {
      pmb->pcoord->Face1Area(k,   j,   is, ie+1, face1);
      pmb->pcoord->Face2Area(k,   j+1, is, ie,   face2p);
      pmb->pcoord->Face2Area(k,   j,   is, ie,   face2m);
      pmb->pcoord->Face3Area(k+1, j,   is, ie,   face3p);
      pmb->pcoord->Face3Area(k,   j,   is, ie,   face3m);
      for(int i=is; i<=ie; i++) {
        divb+=(face1(i+1)*b.x1f(k,j,i+1)-face1(i)*b.x1f(k,j,i)
              +face2p(i)*b.x2f(k,j+1,i)-face2m(i)*b.x2f(k,j,i)
              +face3p(i)*b.x3f(k+1,j,i)-face3m(i)*b.x3f(k,j,i));
      }
    }
  }

  return divb;
}
```

Then, in the `InitUserMeshData` function, set the number of the additional variables you want using the `AllocateUserHistoryOutput` function, and enroll the function defined above.
```c++
    void Mesh::InitUserMeshData(ParameterInput *pin)
    {
      AllocateUserHistoryOutput(1);
      EnrollUserHistoryOutput(0, DivergenceB, "divB");
      return;
    }
```
The function signature of the enrolling function is:
```c++
  void EnrollUserHistoryOutput(int i, HistoryOutputFunc my_func, const char *name,
                               UserHistoryOperation op=UserHistoryOperation::sum);
```
The first parameter of `EnrollUserHistoryOutput` is the index of the user-defined variable, and the third is its name. The fourth and final function parameter is optional and allows the user to choose a reduction operation (across MeshBlocks and across MPI ranks) other than the default choice, summation. Currently, `UserHistoryOption::max` and `UserHistoryOption::min` are the other valid enumerators/choices. If you need more history outputs, simply put a larger number in `AllocateUserHistoryOutput(n)` and enroll the output function with a larger index. The index should be from `0` to `n-1` where `n` is the number of the allocated output variables.

A user-defined history output function receives the index of the output in the second parameter (e.g., `iout = 0` means it is the first user-defined output, and `iout = 1` means it is the second). Users can reuse the same output function using this parameter to control its behavior. For example, one can calculate the enclosed masses within different radii using the same function (e.g., r < 1 for `iout = 0`, r < 2 for `iout = 1`, ...).

### Adding a New Output at Restart
When restarting a simulation from a restart file, one can also update the run-time options by specifying an input file:
```
    > athena -r example.00010.rst -i athinput.new
```
An additional output file (e.g. `<output3>`) can be specified in this input file. The default for this case, is to start the file index at zero and to have the zeroth file be written to disk after the first hydro time step (i.e. one hydro step after the restart file). All subsequent files will be written when the mesh time passes an integer multiple of the output `dt`.

### Visualization Software
See [[Analysis Tools]].
