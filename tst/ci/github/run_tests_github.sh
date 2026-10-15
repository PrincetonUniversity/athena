#!/bin/bash

# SCRIPT: run_tests_github.sh
# PURPOSE: Run the subset of the regression test suite used by GitHub Actions on Ubuntu
#          runners. Assumes compilers, MPI, FFTW, serial HDF5, and the Python
#          dependencies of run_tests.py are already installed (see .github/workflows/ci.yml)
# USAGE: CXX_CHOICE=g++ MPI_CHOICE=openmpi ./tst/ci/github/run_tests_github.sh
#        (run from the top-level athena/ directory)

# Terminate script at first cmd w/ non-zero exit status & echo commands before executing
set -ev

CXX_CHOICE=${CXX_CHOICE:-g++}      # g++, clang++
MPI_CHOICE=${MPI_CHOICE:-openmpi}  # openmpi, mpich

# Ubuntu installs the serial HDF5 library outside of the default search paths
export CPATH=$(pkg-config --cflags-only-I hdf5 | sed 's/-I//g; s/ /:/g')${CPATH:+:$CPATH}
export LIBRARY_PATH=$(pkg-config --libs-only-L hdf5 | sed 's/-L//g; s/ /:/g')${LIBRARY_PATH:+:$LIBRARY_PATH}

# Use the selected serial compiler inside the MPI compiler wrappers
export OMPI_CXX=$CXX_CHOICE
export MPICH_CXX=$CXX_CHOICE
if [ "$MPI_CHOICE" == "openmpi" ]; then
    MPI_OPTS=--oversubscribe
    # Disable OpenMPI vader CMA due to namespace permission issues in containers
    export OMPI_MCA_btl_vader_single_copy_mechanism=none
fi

$CXX_CHOICE --version
mpicxx --version
mpirun --version

cd tst/regression/
CONFIG="--config=--cxx=$CXX_CHOICE --config=--ccmd=$CXX_CHOICE"

# Do not use --silent with pgen_compile, so that compiler warnings appear in the log
time python3 run_tests.py pgen/pgen_compile $CONFIG \
     --config=--cflag="$(../ci/set_warning_cflag.sh $CXX_CHOICE)" -v
time python3 run_tests.py pgen/hdf5_reader_serial $CONFIG -v
time python3 run_tests.py mpi $CONFIG --mpirun_opts=$MPI_OPTS --silent
time python3 run_tests.py hybrid $CONFIG --mpirun_opts=$MPI_OPTS --silent
time python3 run_tests.py omp $CONFIG --silent
time python3 run_tests.py grav $CONFIG --mpirun_opts=$MPI_OPTS --silent  # requires FFTW
time python3 run_tests.py gr $CONFIG --silent
time python3 run_tests.py amr $CONFIG --silent
time python3 run_tests.py hydro $CONFIG -c=--flux=hllc --silent
time python3 run_tests.py outputs $CONFIG --silent
time python3 run_tests.py curvilinear $CONFIG --silent
time python3 run_tests.py symmetry $CONFIG --silent
time python3 run_tests.py eos $CONFIG --silent
time python3 run_tests.py scalars $CONFIG --silent
