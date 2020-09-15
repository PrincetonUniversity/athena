#!/bin/bash

# Terminate script at first cmd w/ non-zero exit status & echo commands before executing (-v) as lines are read
set -ev # -x # for tracing/debugging commands after variable expansion
# Do not use "python run_tests.py" to run all tests, so that:
# - Each set/directory of tests are timed separately (although entire script is timed as one unit in Travis CI)
# - Script fails after first broken set of tests
# (Could alternatively group sets of tests with && operator)

# Assume script is called from the top-level athena/ directory
cd tst/regression/

if [ "$MPI_CHOICE" == "openmpi" ]; then
    PATH=$TRAVIS_BUILD_DIR/openmpi/bin/:$PATH
    MPI_OPTS=--oversubscribe
    # Disable OpenMPI 3.1 vader CMA due to namespace permission issues on Travis CI / Docker containers
    export OMPI_MCA_btl_vader_single_copy_mechanism=none
else
    PATH=$TRAVIS_BUILD_DIR/mpich/bin/:$PATH
fi

# --silent option refers only to stdout of Makefile calls for condensed build logs. Don't use with pgen_compile.py
time python3 run_tests.py pgen/pgen_compile --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --config=--cflag="$(../ci/set_warning_cflag.sh $TEMP_CXX)" -v
# Only building serial HDF5 library on Travis CI (skip "pgen/hdf5_reader_parallel"):
time python3 run_tests.py pgen/hdf5_reader_serial --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD -v

# need to switch serial compiler to Homebrew's GCC instead of /usr/bin/gcc -> Apple Clang for OpenMP functionality
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    time python3 run_tests.py mpi --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --mpirun_opts=$MPI_OPTS --silent
    # TODO(felker): improve selection of 'gcc-9' so when 'brew install gcc' formula instead installs gcc-10, this won't break
    export OMPI_CC=/usr/local/bin/gcc-9
    export OMPI_CXX=/usr/local/bin/g++-9
    export MPICH_CC=/usr/local/bin/gcc-9
    export MPICH_CXX=/usr/local/bin/g++-9
    time python3 run_tests.py hybrid --config=--cxx=g++ --config=--ccmd=/usr/local/bin/g++-9 \
	 --config=--mpiccmd='mpicxx -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX' --mpirun_opts=$MPI_OPTS --silent
    time python3 run_tests.py omp --config=--cxx=g++ --config=--ccmd=/usr/local/bin/g++-9 --silent
    time python3 run_tests.py grav --config=--cxx=g++ --config=--ccmd=/usr/local/bin/g++-9 \
	 --config=--mpiccmd='mpicxx -DMPICH_SKIP_MPICXX -DOMPI_SKIP_MPICXX' --mpirun_opts=$MPI_OPTS --silent # requires FFTW library
else
    export OMPI_CC=$TEMP_CCMD
    export OMPI_CXX=$TEMP_CCMD
    export MPICH_CC=$TEMP_CCMD
    export MPICH_CXX=$TEMP_CCMD
    time python3 run_tests.py mpi --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --mpirun_opts=$MPI_OPTS --silent
    # Fix for broken libomp.h with Travis CI's clang installation on Ubuntu images:
    export LD_LIBRARY_PATH=/usr/local/clang/lib:$LD_LIBRARY_PATH
    time python3 run_tests.py hybrid --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --mpirun_opts=$MPI_OPTS --silent
    time python3 run_tests.py omp --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --silent
    time python3 run_tests.py grav --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --mpirun_opts=$MPI_OPTS --silent # requires FFTW library
fi
time python3 run_tests.py gr --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --silent
#time python3 run_tests.py diffusion --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --silent
time python3 run_tests.py amr --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --silent
time python3 run_tests.py hydro --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD -c=--flux=hllc --silent
time python3 run_tests.py outputs --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --silent
time python3 run_tests.py curvilinear --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --silent
#time python3 run_tests.py sr --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --silent  # 9x tests take about 11-15m on Travis CI
# ~10 min runtime for 2x shearingbox/ tests (2D and 3D MRI)
#time python3 run_tests.py shearingbox --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --silent
time python3 run_tests.py symmetry --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --silent
time python3 run_tests.py eos --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD  --silent
time python3 run_tests.py scalars --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD  --silent

# mhd/ currently contains the longest set of tests. The following command often times-out after 10 m on Travis CI
# time python3 run_tests.py mhd --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --silent

# hydro4/ fourth-order tests currently take >30 min on Travis CI
# time python3 run_tests.py hydro4 --config=--cxx=$TEMP_CXX -c=--ccmd=$TEMP_CCMD --silent
