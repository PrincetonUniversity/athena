### Domain Decomposition: MeshBlock
For parallel simulations with MPI, the computing domain is decomposed into small units. In Athena++, this decomposition unit is called a MeshBlock, and all the MeshBlocks have the same logical size (i.e., the number of cells). These MeshBlocks are stored on a tree structure, and have unique integer IDs numbered by Z-ordering.

The MeshBlock size is specified by `<meshblock>` parameters in an input file. The following example is decomposing a Mesh with 256<sup>3</sup> into MeshBlocks with 64<sup>3</sup> cells, resulting in 64 MeshBlocks. Obviously, the size of the Mesh must be divisible by that of the MeshBlocks.

    <mesh>
    nx1     =    256
    ...
    nx2     =    256
    ...
    nx3     =    256
    ...
    <meshblock>
    nx1     =    64
    nx2     =    64
    nx3     =    64

With output for non-parallelized formats (e.g. VTK), one file is generated per MeshBlock regardless of the actual number of processes. We recommend the HDF5 output because it combines all the MeshBlocks and outputs only two files per output timestep. For details, see [[Outputs]].

### MPI Parallelization
Athena++ is parallelized using standard Message Passing Interface (MPI-2). To enable this, configure the code with `-mpi` option. Each MPI process owns one or more MeshBlocks. The number of MeshBlocks per process may differ, but of course the best load balance is achieved when the computation load is distributed evenly. For MHD, the computational cost of a process is proportional to the number of MeshBlocks, but Athena++ supports weighting based on actual computational costs when additional physics causes load-inbalance (see [[Load Balancing]]). To check the load balance, use the `-m [nproc]` option before running the simulation.

    > athena -i athinput.example -m 64

This can be done with a single process. This tells you how many MeshBlocks are assigned to each process (see also [[Static Mesh Refinement]]).

In the previous example, up to 64 processes can be launched. To start simulations, simply launch the code using the mpiexec/mpirun commands, etc.

    > mpiexec -n 64 athena -i athinput.example

Please consult the document of your system for details.

### OpenMP Parallelization
OpenMP is a standard shared-memory parallelization within a node. OpenMP parallelizes calculations within each MeshBlock. To enable this, configure the code with the `-omp` option and set `num_threads` (**per MPI task**) in the `<mesh>` block in your input file.

    <mesh>
    ...
    num_threads = 4

Also, you may need to set an environment parameter to specify the number of threads available to the solver. Generally this environment variable is `OMP_NUM_THREADS`, but please check the documentation for your system.

Unlike common implementation of OpenMP parallelization over loops, Athena++ uses more coarse-grained parallelization over MeshBlocks. Therefore, to use OpenMP parallelization, you have to assign MeshBlocks more than or equal to `num_threads` per MPI rank. In other words, each OpenMP thread must be assigned at least one MeshBlock, similar to the flat MPI case. Compared to the common loop-based threading, our implementation performs better for this application as it requires less synchronization. Moreover, the simple design allows additional physics to be automatically parallelized as long as the code fits in the MeshBlock design. The result is more maintainable code that is less susceptible to performance and/or correctness bugs. 

### Hybrid Parallelization
Athena++ supports hybrid parallelization using both MPI and OpenMP. For this, a thread-safe MPI library supporting the [`MPI_THREAD_MULTIPLE`](http://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-2.0/node165.htm) option is required. Sometimes this feature can be enabled through an environmental variable. Please consult the documentation of the system you are using. Also, note that it is users' responsibility to keep the code such as the problem generator and user-defined functions (source, boundary, etc.) thread-safe. Access to the `ParameterInput` object is thread-safe. Lastly, please read the following note on performance to find the best configuration.

### Note on Performance
Generally speaking, larger MeshBlocks (especially with large `meshblock/nx1`) are better for performance, but it is a matter of balance between performance and time to get the solution. Using Haswell Xeon E5-2690v3 and Flat MPI parallelization using 24 processes / node, Athena++ can achieve about 7&times;10<sup>5</sup> cells per second per process for MHD (and almost twice for hydrodynamics) and its weak-scaling is almost perfect when 64<sup>3</sup> cells per process are used. In other words, one timestep takes less than 0.4 second in this situation. If the performance is significantly lower than these values, it means something is wrong.

Generally, we recommend to use MeshBlocks with at least 32<sup>3</sup> cells, preferably 64<sup>3</sup> cells per process (thread). For [[Adaptive Mesh Refinement]], smaller MeshBlocks like 16<sup>3</sup> are useful allowing more flexibility, but if it is too small the overhead will become significant. However, these numbers depend on computers. Therefore we recommend you test the performance of your system before you start production runs.

OpenMP parallelization is not very scalable in general. Generally speaking, we recommend to use flat-MPI or hybrid parallelization with a few (2-4) threads per physical core, except when you are using only one shared-memory node on which full OpenMP parallelization works fine. For CPUs with rich cores, such as Intel Xeon processors, flat-MPI parallelization often gives the best performance. Although modern CPUs support HyperThreading which allows you to launch more than one thread per core, usually you get no performance gain. On the other hand, for architectures with simpler cores but with fast hardware threads, including Intel Xeon Phi (Knights Landing; KNL) and IBM Blue Gene/Q (possibly newer POWER8/9 as well but we have not tested them), it is probably the best to launch a few threads per core (usually 4 for KNL and BG/Q) to keep the cores busy. Because OpenMP threads can share some data, especially the MeshBlock tree, it saves some memory. When you are running extremely large parallel simulations, this will be helpful. In short, hybrid parallelization is recommended only when one or more of the following conditions are met:
* You are using CPUs with fast hardware threads (e.g. KNL and BG/Q).
* You need to save memory foot print by sharing some resources.
* Your system vendor requires it.


#### Changes from the old Athena
Because we completely redesigned the data structure, the domain decomposition is considerably different from Athena. In Athena, the number of processes in each direction was specified, or automatic decomposition was available. In Athena++, now a user must specify the MeshBlock size explicitly. Also, each MeshBlock is numbered using Z-ordering (Morton ordering). These changes are mainly for implementation of Adaptive Mesh Refinement.

Also, previously the data hierarchy was organized as Mesh (the whole computing domain) → Domain (SMR level) → Grid (a decomposition unit). Athena++ does not have Domain, and the level is only a property of each MeshBlock.
