#!/usr/bin/env bash

# SCRIPT: run_jenkins_perseus.sh
# AUTHOR: Kyle Gerard Felker - kfelker@princeton.edu
# DATE: 4/10/2018
# PURPOSE: Run style and regression test suites using PICSciE's Jenkins server for continuous
# integration (CI) Current workers include 4x Intel Broadwell nodes on Perseus cluster.

# USAGE: salloc -N1 -n4 --time=0:60:00 ./run_jenkins_perseus.sh
# or similar command in the Jenkins build "Execute shell" step (run from athena/ root dir)

set -e # terminate script at first error/non-zero exit status
# Install Python dependencies
pip install -q --user h5py # outputs/all_outputs.py uses athena_read.athdf() reader
pip install -q --user flake8

# Build step #0: Test source code style consistency
# step #0a: lint Python files
python -m flake8
echo "Finished linting Python files with flake8"

# step #0b: lint C++ files
cd tst/style/; ./cpplint_athena.sh
cd ../regression/

# Build step #1: regression tests using GNU compiler and OpenMPI library
module purge
module load rh # latest GNU compiler
module load openmpi/gcc/1.10.2/64 # openmpi/gcc/3.0.0/64 does not work right now
# Do NOT "module load hdf5" = hdf5/intel-17.0/openmpi-1.10.2/1.10.0
# output/all_outputs.py regression test uses non-MPI HDF5 writer
# (Perseus will error w/ missing mpi.h header if MPI HDF5 is loaded w/o mpicxx)
# grav/ regression tests require MPI and FFTW
module load hdf5/gcc/1.10.0
module load fftw/gcc/3.3.4
module list

# Run regression test sets. Need to specify Slurm mpirun wrapper, srun
# --silent option refers only to stdout of Makefile calls for condensed build logs. Don't use with pgen_compile.py
time python ./run_tests.py pgen/pgen_compile --config=--cflag="$(../ci/set_warning_cflag.sh g++)"
time python ./run_tests.py pgen/hdf5_reader_serial --silent
time python ./run_tests.py grav --mpirun=srun --silent
time python ./run_tests.py mpi --mpirun=srun --silent
time python ./run_tests.py hydro --silent
# MHD is currenlty the longest regression test set:
time python ./run_tests.py mhd --silent
time python ./run_tests.py amr --silent
time python ./run_tests.py outputs --silent
time python ./run_tests.py sr --silent
time python ./run_tests.py gr --silent
time python ./run_tests.py curvilinear --silent
time python ./run_tests.py shearingbox --silent
time python ./run_tests.py diffusion --silent
time python ./run_tests.py symmetry --silent

# High-order solver regression tests w/ GCC
time python ./run_tests.py hydro4 --silent

# Swap serial HDF5 library module for parallel HDF5 library:
module unload hdf5/gcc/1.10.0
module load hdf5/gcc/openmpi-1.10.2/1.10.0
# Workaround issue with parallel HDF5 modules compiled with OpenMPI on Perseus--- linker still takes serial HDF5 library in /usr/lib64/
# due to presence of -L flag in mpicxx wrapper that overrides LIBRARY_PATH environment variable
time python ./run_tests.py pgen/hdf5_reader_parallel --mpirun=srun --config=--lib=/usr/local/hdf5/gcc/openmpi-1.10.2/1.10.0/lib64--silent

# Build step #2: regression tests using Intel compiler and MPI library
module purge
module load intel
module load intel-mpi
module load fftw/gcc/3.3.4
module load hdf5/intel-17.0/1.10.0 # hdf5/intel-17.0/intel-mpi/1.10.0
module load rh
module list

time python ./run_tests.py pgen/pgen_compile --config=--cxx=icc --config=--cflag="$(../ci/set_warning_cflag.sh icc)"
time python ./run_tests.py pgen/hdf5_reader_serial --silent
time python ./run_tests.py grav --config=--cxx=icc --mpirun=srun --silent
time python ./run_tests.py mpi --config=--cxx=icc --mpirun=srun --silent
time python ./run_tests.py hydro --config=--cxx=icc --silent
time python ./run_tests.py mhd --config=--cxx=icc --silent
time python ./run_tests.py amr --config=--cxx=icc --silent
time python ./run_tests.py outputs --config=--cxx=icc --silent
time python ./run_tests.py sr --config=--cxx=icc --silent
time python ./run_tests.py gr --config=--cxx=icc --silent
time python ./run_tests.py curvilinear --config=--cxx=icc --silent
time python ./run_tests.py shearingbox --config=--cxx=icc --silent
time python ./run_tests.py diffusion --config=--cxx=icc --silent
time python ./run_tests.py symmetry --config=--cxx=icc --silent

# High-order solver regression tests w/ Intel compiler
time python ./run_tests.py hydro4 --config=--cxx=icc --silent

# Swap serial HDF5 library module for parallel HDF5 library:
module unload  hdf5/intel-17.0/1.10.0
module load hdf5/intel-17.0/intel-mpi/1.10.0
time python ./run_tests.py pgen/hdf5_reader_parallel --config=--cxx=icc --mpirun=srun --silent

# Test OpenMP 4.5 SIMD-enabled function correctness by disabling IPO and forced inlining w/ Intel compiler flags
# Check subset of regression test sets to try most EOS functions (which heavily depend on vectorization) that are called in rsolvers
time python ./run_tests.py pgen/pgen_compile --config=--cxx=icc-debug --config=--cflag="$(../ci/set_warning_cflag.sh icc)"
time python ./run_tests.py hydro --config=--cxx=icc-debug --silent
time python ./run_tests.py mhd --config=--cxx=icc-debug --silent
time python ./run_tests.py sr --config=--cxx=icc-debug --silent
time python ./run_tests.py gr --config=--cxx=icc-debug --silent

set +e
# end regression tests

# Codecov coverage analysis
# Pipe to bash (Jenkins)
curl -s https://codecov.io/bash | bash -s - -t ccdc959e-e2c3-4811-95c6-512151b39471
