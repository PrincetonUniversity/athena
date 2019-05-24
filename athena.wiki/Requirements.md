### Mandatory requirements:
* This is a Git repository, and therefore Git is required to clone a local version. The GitHub UI also offers the option to "Download ZIP" for the latest version of the repository on the main page, and packaged `.zip` and `.tar.gz` version snapshots are regularly made available under Releases. 
* Athena++ is written in C++ and requires a C++ compiler such as `g++`, `clang++`, or `icc`.
* The code is configured with a Python script; therefore Python is required to effectively use the code. The configure script and auxiliary visualization and testing scripts support Python 2.7+ and Python 3.4+.

Note: The latest Intel C++ compiler is recommended for use on Intel processors because it results in significantly better performance through the use of OpenMP 4.0 SIMD vectorization extensions.  The IBM XLC and Cray CC compiler are also supported.

### Optional requirements (for use with certain features):
* An OpenMP enabled compiler (e.g. `gcc4.2` or above) is required to use shared-memory parallelism.
* An MPI-2 (Message Passing Interface version 2) library is required to use distributed memory parallelism.
* For hybrid parallelization using both MPI and OpenMP, a thread-safe MPI library supporting `MPI_THREAD_MULTIPLE` is needed. Sometimes this feature can be enabled through an environmental variable.
* HDF5 outputs (which are essential for AMR) require an HDF library.
  * The library must be configured for Parallel HDF5 (PHD5) if the solver uses MPI. See [[HDF5 Format]] for details.
* FFTW library is required to use FFT for the self-gravity solver and turbulence driving for serial and parallel simulations. 
  * Note, Athena++ is packaged with Plimpton's library for parallel FFT when running the solver with MPI, but it is configured to use FFTW's serial FFT for the underlying operations.

<!-- Bump MPI-2 to MPI-3? -->
