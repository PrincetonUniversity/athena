### Problem Generator File
When configuring the code, a problem generator must be specified using the `--prob` option. A list of available problem generators is given in the help message of the configuration script (use the `-h` option).

Problem generators stored in `src/pgen` contain problem-specific functions and variables. These include:
* Initial conditions
* User-defined boundary conditions
* User-defined physical source terms
* User-defined Mesh spacing functions
* User-defined mesh refinement criteria
* User-defined analysis functions
* User-defined Mesh data field
* User-defined MeshBlock data field
* etc...

It is probably necessary to read the Programmer Guide in order to understand the data structures in Athena++ well enough to write a new problem generator. The existing files in the `/src/prob` directory can be used as starting templates. In this section, we briefly explain the structure of problem generators.

The following interfaces (and any other functions needed) can be defined in a problem generator file:
* `void Mesh::InitUserMeshData(ParameterInput *pin)`
* `void Mesh::UserWorkInLoop(void)`
* `void Mesh::UserWorkAfterLoop(void)`
* `void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)`
* `void MeshBlock::ProblemGenerator(ParameterInput *pin)`
* `void MeshBlock::UserWorkInLoop(void)`
* `void MeshBlock::UserWorkBeforeOutput(Parameter Input *pin)`

These functions are optional (although in practice a `MeshBlock::ProblemGenerator` method must be always defined). When these functions are not defined, default ones that do nothing are used. The functions are now described in detail:

#### `Mesh::InitUserMeshData`
This function is called at the beginning of a simulation (for both the initial  and any restarted runs) for allocating and initializing global variables shared with all the MeshBlocks, and for enrolling user-defined functions. The following functions can be enrolled here:
* Boundary conditions (see [[Boundary Conditions]])
* Mesh generator (mesh spacing) (see [[Coordinate Systems]])
* Source terms (see below)
* Time step (see below)
* Mesh refinement criteria (see [[Adaptive Mesh Refinement]])
* User-defined general relativity coordinates (see [[Arbitrary Coordinates]])

Also, the user-defined data field in the Mesh class can be declared and allocated in this function. This feature enables users to store data that are automatically saved and loaded when the code is restarted. See below for details.

#### `Mesh::UserWorkInLoop`
This function is called at the end of every cycle (i.e. after the whole task list). This function can be used for analysis over MeshBlocks. As this function is called after `MeshBlock::UserWorkInLoop`, this function can be used to collect the results from MeshBlocks. Also, since this is called only once per node, MPI all-to-all communications can be used here.


#### `Mesh::UserWorkAfterLoop`
This function is called at the end of the simulation for cleaning up the resources allocated in `MeshBlock::InitUserMeshData`.

Users do not need to delete the user-defined data arrays allocated in `Mesh:InitUserMeshData` using `Allocate(Real/Int)UserMeshBlockDataField` as they are automatically deleted.

#### `MeshBlock::InitUserMeshBlockData`
This is similar to `Mesh::InitUserMeshData`, and is called at the beginning of a simulation (for both the initial and any restarted runs) for allocating and initializing local variables belonging to each MeshBlock. This is mainly intended for  User MeshBlock Data. This function is not for initializing the simulation setup, which should be done in `MeshBlock::ProblemGenerator`.

#### `MeshBlock::ProblemGenerator`
This function sets the initial conditions for the problem. For hydrodynamics, the cell-centered *conservative* variables (`phydro->u`) must be set here. The face-centered magnetic fields (`pfield->bx1f`, `bx2f`, and `bx3f`) also must be set for MHD. When [[Passive Scalars]] are enabled, the density of each species (`pscalars->s`) must be initialized within this function. Other variables such as primitive variables and cell-centered magnetic fields are automatically derived by the code (though with GR the primitives must also be set). For non-trivial (i.e. non-uniform) magnetic fields, the magnetic fields should be initialized using the vector potential in order to satisfy the &nabla;&#8901;**B** = 0 constraint exactly.

For details, see the Programmer Guide and sample files in the `src/pgen` directory.

#### `MeshBlock::UserWorkInLoop`
This function is called at the end of every timestep (note: it is not called at the half timestep). A user can do analysis in this function. This is intended only for analysis, and not for manipulating the data. In such a case, one should use a user-defined source term function and/or the boundary conditions.

#### `MeshBlock::UserWorkBeforeOutput`
This function is called at the end of a time step only when output files (except history output) are about to be generated. This function is dedicated for calculating the user-defined output variables (see [[Outputs]]).

### Reading Parameters from the Input File
Users can read parameters from the input file in `Mesh::InitUserMeshData(ParameterInput *pin)` and `MeshBlock::ProblemGenerator(ParameterInput *pin)`. For this purpose, the following functions are provided:
* `int ParameterInput::GetInteger(std::string block, std::string name)`
* `Real ParameterInput::GetReal(std::string block, std::string name)`
* `std::string GetString(std::string block, std::string name)`
* `int ParameterInput::GetOrAddInteger(std::string block, std::string name, int value)`
* `Real ParameterInput::GetOrAddReal(std::string block, std::string name, Real value)`
* `std::string ParameterInput::GetOrAddString(std::string block, std::string name, std::string value)`

The first three functions return a parameter with `name` in `<block>`. These functions fail when that parameter is not defined in the input file. The latter three functions are similar but with a default value. For example, in order to read a plasma beta parameter named "beta" in the `<problem>` block,
```c++
    Real beta = pin->GetReal("problem", "beta");
```
And if you want to set a default value if no value is specified in the input file,
```c++
    Real beta = pin->GetReal("problem", "beta", 1000.0);
```
For details, again, see the Programmer Guide and sample files in the `src/pgen` directory.

### User-Defined Physical Source Term
A user-Defined physical source term such as a cooling/heating function can be implemented in the problem generator file. Note that this feature is designed for time-explicit integration, and does not work well with any time-implicit or sub-cycling method. If you need these integrators, you have to implement the integrator for yourself in the `MeshBlock::UserWorkInLoop` function.

First, define a source term function with a certain prototype:
```c++
    void MySource(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
```
Then enroll it in the `Mesh::InitUserMeshData` function:
```c++
    EnrollUserExplicitSourceFunction(MySource);
```
Within the source function, the conservative variable should be updated (cons). This source function is called after the MHD update but before conservative to primitive conversion. Therefore, the MHD update is already reflected in the conservative variables but the primitive variables (prim) and cell-centered magnetic fields are not updated yet. The conservative variables should be updated using these primitive variables and cell-centered magnetic fields. For example, a simple cooling function can be implemented like:
```c++
    void cooling(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                 const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
    {
      Real g = pmb->peos->GetGamma();
      Real tau = 0.01;
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real temp = (g-1.0)*prim(IEN,k,j,i)/prim(IDN,k,j,i);
            cons(IEN,k,j,i) -= dt*prim(IDN,k,j,i)*(temp - 10.0)/tau/(g-1.0);
          }
        }
      }
      return;
    }
```
If your source term function requires a time-step constraint, it should be implemented and enrolled as explained below. And if the time step is too restrictive, consider writing a time-implicit integrator.

### User-Defined Time Step
In order to impose additional time-step constraint, write and enroll a time step function that returns the minimum required time step within a MeshBlock. This time step is compared with the minimum time step for the MHD part, and the smaller one is used. Note that no additional safety factor (like the CFL number) is multiplied for the user-defined time step.
```c++
    Real MyTimeStep(MeshBlock *pmb)
    {
      Real min_dt=FLT_MAX;
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real dt;
            dt = ... // calculate your own time step here
            min_dt = std::min(min_dt, dt);
          }
        }
      }
      return min_dt;
    }
```
Then enroll it in the `Mesh::InitUserMeshData` function:
```c+++
      EnrollUserTimeStepFunction(MyTimeStep);
```

### User-Defined Data Field
While any variables can be defined in the problem generator file, they are not saved or transferred when the simulation is restarted. Athena++ provides interfaces for users to allocate variables that can be automatically saved and loaded when the simulation is restarted. The data can be arrays of any size, and can be Real or integer. The data field can be stored in Mesh and MeshBlock. User Mesh Data can be used to store the data shared among nodes (i.e. global, common data), while User MeshBlock Data are accessible only within a MeshBlock.

In order to use this feature, first declare the number of the data field:
* `void Mesh::AllocateRealUserMeshDataField(int n)`
* `void Mesh::AllocateIntUserMeshDataField(int n)`
* `void MeshBlock::AllocateRealUserMeshBlockDataField(int n)`
* `void MeshBlock::AllocateIntUserMeshBlockDataField(int n)`

Then, allocate the data arrays such as
* `AthenaArray<int> Mesh::iuser_mesh_data[i]`
* `AthenaArray<Real> Mesh::ruser_mesh_data[i]`
* `AthenaArray<int> MeshBlock::iuser_meshblock_data[i]`
* `AthenaArray<Real> MeshBlock::ruser_meshblock_data[i]`

These functions should be called in `Mesh::InitUserMeshData` or `MeshBlock::InitUserMeshBlockData`.
For example, in order to allocate two Real arrays, one is just one variable and the other is a 16x16 2D array, in Mesh, write the following code:
```c++
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  AllocateRealUserMeshDataField(2);
  ruser_mesh_data[0].NewAthenaArray(1);
  ruser_mesh_data[1].NewAthenaArray(16,16);
  for(int j=0; j<16; j++) {
    for(int i=0; i<16; i++) {
      // initialization
      ruser_mesh_data[1](j,i)=0.0;
    }
  }
  return;
}
```
Then these data can be used from other functions that can see the Mesh class. These data are automatically saved and loaded when the simulation is restarted. However, note that only the data on the master node (rank 0) is written in the restarting file, and it is users' responsibility to maintain consistency between nodes.

Similarly, in order to allocate a 3D integer array in MeshBlocks, write:
```c++
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateIntUserMeshBlockDataField(1);
  iuser_meshblock_data[0].NewAthenaArray(16,16,16);
  for(int k=0; k<16; k++)
    for(int j=0; j<16; j++) {
      for(int i=0; i<16; i++) {
        // initialization
        iuser_meshblock_data[0](k,j,i)=0;
      }
    }
  }
  return;
}
```
Because these data are local data within MeshBlocks, each MeshBlock outputs its own data into the restarting file.

Note that the data fields allocated using this feature are automatically cleaned up at the end of the simulation, so users do not need to delete them.

### Initializing with Preexisting Data

The `from_array` problem generator (`src/pgen/from_array.cpp`) and associated input file (`inputs/mhd/athinput.from_array`) are designed to facilitate starting a simulation from initial values contained in a data file. The problem generator supports Newtonian and special-relativistic hydro and MHD.

The initial values must be stored in an HDF5 data file with the following specifications:
  - The conserved values must be in a 5D dataset with the following axes:
    - Physical quantity: conserved density and momentum (and energy if applicable), with ordering described by the input file.
    - Global MeshBlock ID number (`gid` in the code). The data must be organized by MeshBlock.
    - x3-index of active zone within the MeshBlock, starting at 0 and not including any ghost zones.
    - x2-index of active zone within the MeshBlock, starting at 0 and not including any ghost zones.
    - x1-index of active zone within the MeshBlock, starting at 0 and not including any ghost zones.
  - If magnetic fields are enabled, each of the 3 **face-centered** fields must be in **separate** 4D datasets with the following axes:
    - Global MeshBlock ID number (`gid` in the code). The data must be organized by MeshBlock.
    - x3-index of active zone within the MeshBlock, starting at 0 and not including any ghost zones.
    - x2-index of active zone within the MeshBlock, starting at 0 and not including any ghost zones.
    - x1-index of active zone within the MeshBlock, starting at 0 and not including any ghost zones.

The conserved values can generally be taken from HDF5 outputs made by the code. However the standard HDF5 output yields cell-centered magnetic fields, which are not sufficient for initializing a problem.