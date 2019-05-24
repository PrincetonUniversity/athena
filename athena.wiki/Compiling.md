### Makefile
The configuration script creates a `Makefile` in the top-level directory, and you can build the code using standard make tools. Before building the code, it is strongly recommended that any stale or temporary files created by previous compilations be cleaned up.
```
    > make clean
```
Then simply,
```
    > make
```
Depending on the compiler and level of optimization, the compilation may take a while. You can speed it up using parallel make
```
    > make -j
```
This process creates the executable `bin/athena`.


### Compatibility
Because Athena++ is written in standard C++11 (although we do not use most of its new features), it should be able to be compiled using any compiler supporting the standard. This includes the GNU C++ compiler (`g++`), Intel C++ compiler (`icc`/`icpc`), Clang, Cray C++ Compiler, PGI C++ Compiler, IBM XL C++ Compiler, etc. Athena++ uses `#pragma omp simd` defined in OpenMP 4.0 and later, but the pragma statements should not cause compatibility problems. Our configuration script supports GNU, Intel, Cray, and IBM. If you want to use other compilers, configure the code for `g++` and modify `Makefile` to specify appropriate options, flags, and libraries (see below).

Athena++ uses OpenMP for shared-memory parallelization and MPI-2 for distributed-memory parallelization. Most of the latest compilers support OpenMP used in Athena++, but if you are using Clang, we recommend the latest version as the older versions did not support OpenMP. For distributed-memory parallelization, any MPI library supporting the MPI-2 standard should work, including OpenMPI, MPICH, MVAPICH, Intel, Cray, IBM, etc.

### For best performance
On Intel-based systems, we strongly recommend use of the Intel C++ Compiler suite (`icc`, version 15 or later) because it results in a significantly faster executable compared to other compilers. In addition, static linking (`-static`) may improve the performance of highly parallel simulations. Note that `Makefile` must be modified whenever the code is (re)configured, as the configuration script overwrites it.

### Known issues:
- Versions of the OpenMPI library around and before v3.0.0 may have issues when writing many output files (>200 per simulation). This issue was patched in late 2017/early 2018 and should be fixed in newer issues of OpenMPI. See [open-mpi/ompi Issue #4336](https://github.com/open-mpi/ompi/issues/4336)