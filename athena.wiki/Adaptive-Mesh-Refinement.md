### Enabling Adaptive Mesh Refinement
Like [[Static Mesh Refinement]], no configuration option is needed to enable AMR. You need to set the refinement parameters in the input file. For example,

    <mesh>
    ...
    refinement     = adaptive
    numlevel       = 4
    deref_count    = 5

This means the finer levels are created up to 4 levels (= the root level + 3 finer levels).
The `deref_count` parameter means that at least 5 timesteps are required to be flagged for derefinement before a MeshBlock is actually derefined. This suppresses destruction of newly-refined MeshBlocks immediately after refinement.

### The AMR Condition, `RefinementCondition`
AMR requires a function to check whether a MeshBlock should be refined or derefined. This function must be defined and enrolled in the [[Problem Generator|Problem Generators]]. For example, in `src/pgen/dmr.cpp` (2D Double Mach Reflection):
```c++
int RefinementCondition(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps=0.0;
  int k=pmb->ks;
  for(int j=pmb->js; j<=pmb->je; j++) {
    for(int i=pmb->is; i<=pmb->ie; i++) {
      Real epsr= (std::abs(w(IDN,k,j,i+1)-2.0*w(IDN,k,j,i)+w(IDN,k,j,i-1))
                 +std::abs(w(IDN,k,j+1,i)-2.0*w(IDN,k,j,i)+w(IDN,k,j-1,i)))/w(IDN,k,j,i);
      Real epsp= (std::abs(w(IEN,k,j,i+1)-2.0*w(IEN,k,j,i)+w(IEN,k,j,i-1))
                 +std::abs(w(IEN,k,j+1,i)-2.0*w(IEN,k,j,i)+w(IEN,k,j-1,i)))/w(IEN,k,j,i);
      Real eps = std::max(epsr, epsp);
      maxeps = std::max(maxeps, eps);
    }
  }
  if(maxeps > 0.01) return 1;
  if(maxeps < 0.005) return -1;
  return 0;
}
```

Then this function has to be enrolled in the `Mesh::InitUserMeshData` function in your problem generator along with other boundary conditions:
```c++
  void Mesh::InitUserMeshData(ParameterInput *pin)
  {
    EnrollUserBoundaryFunction(INNER_X1, DMRInnerX1);
    EnrollUserBoundaryFunction(INNER_X2, DMRInnerX2);
    EnrollUserBoundaryFunction(OUTER_X2, DMROuterX2);
    if(adaptive==true)
      EnrollUserRefinementCondition(RefinementCondition);
  }
```

This RefinementCondition function checks the curvature (second derivative) of the density and pressure in a MeshBlock given as an argument, which is used as a simple shock detector. When at least one cell exceeds a certain threshold (0.01), it flags this block to be refined (`return 1`). When all the cells have curvature smaller than 0.005, then flag it to be derefined (`return -1`). Otherwise, it does nothing (`return 0`). While the primitive variables (`phydro->w`) are used in this example, the conservative variables (`phydro->u`) or cell-centered magnetic fields (`pfield->bcc`) can be used as well.

When a MeshBlock is flagged to be refined, it is always refined at the end of the timestep as long as the number of levels does not reach the limit. Also, a MeshBlock can be refined even if it is not flagged to be refined; this occurs when neighboring MeshBlocks are refined and the level difference between neighboring MeshBlocks becomes larger than one. In such a situation, the code automatically refines MeshBlocks and constructs a consistent tree of MeshBlocks.

On the other hand, even if a MeshBlock is flagged to be derefined, it is not necessarily derefined. The following steps are considered:

  1. It must be flagged to be derefined for `deref_count` consecutive time steps.
  2. All the eight (in 3D, four in 2D) MeshBlocks to be merged into one coarser MeshBlock must be flagged to be derefined.
  3. The MeshBlock after derefinement must only be in contact with MeshBlocks on the same level or one level different.

A MeshBlock is derefined only when all of these conditions are satisfied.

The RefinementCondition function should satisfy some requirements. Most importantly, it must be consistent - when a MeshBlock is flagged to be refined, the newly created MeshBlock should not satisfy the derefinement condition, and vice versa. Also it should not be too expensive.

### Running AMR Simulations
As discussed in [[Using MPI and OpenMP]], the size of MeshBlocks should be determined considering the balance between flexibility and performance. For AMR, it is also important for load balancing. If a process has only a few MeshBlocks, the difference of computational load can be significant. If one process has about 10 MeshBlocks, the load difference is only about 10% (assuming all the MeshBlocks have the same costs), which is probably acceptable. A user must choose the MeshBlock size carefully.

In the initialization process, the code checks the refinement criteria on all the MeshBlocks. If refinement is needed, it creates finer MeshBlocks and calls the Problem Generator. It repeats this process until the refinement criteria are satisfied. This process ensures that the initial condition is resolved well enough, but it is difficult to predict how many MeshBlocks are created in advance. The code will give a warning when many MeshBlocks are refined during the initialization. Unfortunately, the `-m` option [[for checking MPI load balancing|Using-MPI-and-OpenMP#mpi-parallelization]] does not work here because this mode does not run the problem generator and therefore refinement is not performed.

In order to assign an adequate number of processes for an AMR simulation, "two step initialization" is needed. First, run the simulation just one step with a [[restarting output|Outputs#restart-file]], just for creating the initial conditions. This can be done as a serial run, but if it requires a large amount memory or you want to save time, more processes can be used.

    > mpiexec -n 16 athena -i athinput.amrexample time/nlim=0

The created restarting file contains all the refined MeshBlocks. Then use the `-m` option for this restarting file to check the number of MeshBlocks, and determine how many processes are needed for the simulation.

    > athena -r amrexample.out1.00000.rst -m 64

Then run the simulation using the restarting file.

    > mpiexec -n 64 athena -r amrexample.out1.00000.rst

Note that this technique can be used also for adjusting the number of processes during the simulation. If the number of MeshBlocks significantly changes, you can increase or decrease the number of processes when you restart the simulation.

### Data Output
With adaptive mesh refinement, the number, level and location of MeshBlocks change as the simulation runs. Serial output formats like VTK are not very useful with AMR; the file ID (= MeshBlock ID) changes when mesh refinement occurs, and for now the `join_vtk++` tool (in `vis/vtk`) does not support data with mesh refinement. Therefore we strongly recommend the HDF5 output for AMR simulations. This format creates only two files per output, regardless of the number of processes or MeshBlocks. For details, please read [[Outputs]] and [[Analysis Tools]].

*Note*: [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) does not update the number of `MeshBlocks` stored as HDF5 data containers when you change the timeslice in an open database. This may result in missing/empty patches of the domain when AMR derefines a `MeshBlock`; reopening the database will fix the missing data. See [[Resampling HDF5 Outputs]] for an alternative solution to this issue. 