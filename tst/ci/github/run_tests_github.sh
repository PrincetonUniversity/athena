#!/bin/bash

# SCRIPT: run_tests_github.sh
# PURPOSE: Run a subset of the regression test suite on GitHub Actions runners. Assumes
#          compilers, MPI, FFTW, serial HDF5, and the Python dependencies of run_tests.py
#          are already installed (see .github/workflows/ci.yml)
# USAGE: CXX_CHOICE=g++ MPI_CHOICE=openmpi ./tst/ci/github/run_tests_github.sh
#        (run from the top-level athena/ directory)
# ENVIRONMENT:
#   CXX_CHOICE  configure.py --cxx choice: g++, clang++, icpx
#   CXX_CMD     C++ compiler command (default: $CXX_CHOICE), e.g. g++-15
#   MPI_CHOICE  openmpi, mpich, intelmpi
#   OMP_CHOICE  configure.py --cxx choice for the omp/ and hybrid/ tests (default:
#               $CXX_CHOICE); Apple Clang has no OpenMP, so macOS uses Homebrew GCC
#   OMP_CMD     C++ compiler command for the omp/ and hybrid/ tests (default: $CXX_CMD)
#   TEST_SET    default: the sets formerly run on Travis CI
#               extended: the longer sets, run in a separate job

# Terminate script at first cmd w/ non-zero exit status & echo commands before executing
set -ev

CXX_CHOICE=${CXX_CHOICE:-g++}
CXX_CMD=${CXX_CMD:-$CXX_CHOICE}
MPI_CHOICE=${MPI_CHOICE:-openmpi}
OMP_CHOICE=${OMP_CHOICE:-$CXX_CHOICE}
OMP_CMD=${OMP_CMD:-$CXX_CMD}
TEST_SET=${TEST_SET:-default}

# HDF5 (Ubuntu) and HDF5 and FFTW (Homebrew) are installed outside of the default search
# paths of the compilers
add_pkg_paths () {
    for flag in $(pkg-config --cflags-only-I "$1"); do
        CPATH=${flag#-I}${CPATH:+:$CPATH}
    done
    for flag in $(pkg-config --libs-only-L "$1"); do
        LIBRARY_PATH=${flag#-L}${LIBRARY_PATH:+:$LIBRARY_PATH}
    done
}
add_pkg_paths hdf5
add_pkg_paths fftw3
export CPATH LIBRARY_PATH

# Use the selected serial compiler inside the MPI compiler wrappers
export OMPI_CXX=$CXX_CMD MPICH_CXX=$CXX_CMD I_MPI_CXX=$CXX_CMD
MPIRUN_OPTS=
if [ "$MPI_CHOICE" == "openmpi" ]; then
    # Runners may report fewer cores than the number of ranks in the tests
    MPIRUN_OPTS=--mpirun_opts=--oversubscribe
    # Disable OpenMPI vader CMA due to namespace permission issues in containers
    export OMPI_MCA_btl_vader_single_copy_mechanism=none
elif [ "$MPI_CHOICE" == "intelmpi" ]; then
    # No network fabric on the runners; use shared memory only
    export I_MPI_FABRICS=shm
fi

$CXX_CMD --version
mpicxx --version
mpirun --version

cd tst/regression/
CONFIG="--config=--cxx=$CXX_CHOICE --config=--ccmd=$CXX_CMD"
OMP_CONFIG="--config=--cxx=$OMP_CHOICE --config=--ccmd=$OMP_CMD"

if [ "$TEST_SET" == "extended" ]; then
    time python3 run_tests.py sr $CONFIG --silent
    time python3 run_tests.py mhd $CONFIG --silent
    time python3 run_tests.py shearingbox $CONFIG $MPIRUN_OPTS --silent
    time python3 run_tests.py turb $CONFIG $MPIRUN_OPTS --silent  # requires FFTW
    time python3 run_tests.py diffusion $CONFIG --silent
    time python3 run_tests.py cr $CONFIG --silent
    time python3 run_tests.py nr_radiation $CONFIG $MPIRUN_OPTS --silent
    time python3 run_tests.py implicit_radiation $CONFIG $MPIRUN_OPTS --silent
    time python3 run_tests.py multi_group $CONFIG --silent
    time python3 run_tests.py chemistry $CONFIG --silent
    exit 0
fi

# Do not use --silent with pgen_compile, so that compiler warnings appear in the log
time python3 run_tests.py pgen/pgen_compile $CONFIG \
     --config=--cflag="$(../ci/set_warning_cflag.sh $CXX_CHOICE)" -v
time python3 run_tests.py pgen/hdf5_reader_serial $CONFIG -v
time python3 run_tests.py mpi $CONFIG $MPIRUN_OPTS --silent
time python3 run_tests.py grav $CONFIG $MPIRUN_OPTS --silent  # requires FFTW
(export OMPI_CXX=$OMP_CMD MPICH_CXX=$OMP_CMD I_MPI_CXX=$OMP_CMD
 time python3 run_tests.py hybrid $OMP_CONFIG $MPIRUN_OPTS --silent)
time python3 run_tests.py omp $OMP_CONFIG --silent
time python3 run_tests.py gr $CONFIG --silent
time python3 run_tests.py amr $CONFIG --silent
time python3 run_tests.py hydro $CONFIG -c=--flux=hllc --silent
time python3 run_tests.py outputs $CONFIG --silent
time python3 run_tests.py curvilinear $CONFIG --silent
# icpx defaults to -fp-model=fast, which breaks the bit-exact symmetry checked by
# symmetry/mhd_blast_mirror (see tst/ci/jenkins/run_jenkins_stellar.sh)
SYMMETRY_CFLAG=()
if [ "$CXX_CHOICE" == "icpx" ]; then
    SYMMETRY_CFLAG=(--config=--cflag="-fp-model=precise -ffp-contract=off")
fi
time python3 run_tests.py symmetry $CONFIG "${SYMMETRY_CFLAG[@]}" --silent
time python3 run_tests.py eos $CONFIG --silent
time python3 run_tests.py scalars $CONFIG --silent
