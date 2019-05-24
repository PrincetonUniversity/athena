## Welcome to the Athena++ Code Project.

Athena++ is a complete re-write of the Athena astrophysical magnetohydrodynamics (MHD) code into C++.  Compared to earlier versions, the Athena++ code has (1) much more flexible coordinate and grid options including adaptive mesh refinement, (2) new physics including general relativity, (3) significantly improved performance and scalability, and (4) improved source code clarity and modularity.

The code is freely available to the community.  New features are always under development, and will be made public once they are thoroughly tested.

The current public version supports the following physics:
* compressible hydrodynamics and MHD in 1D, 2D, and 3D;
* special and general relativistic hydrodynamics and MHD.

In addition, this version supports the following grid and algorithm options:
* Cartesian, cylindrical, or spherical polar coordinates;
* static or adaptive mesh refinement in any coordinate system;
* mixed parallelization with both OpenMP and MPI;
* a task-based execution model for improved load balancing, scalability and modularity.

### Learn More

Documentation is provided using the wiki pages.

See [[the Athena++ website|http://princetonuniversity.github.io/athena/]] for additional supporting material (like test suites, performance metrics, etc.).
