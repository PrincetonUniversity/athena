#!/usr/bin/env bash

# SCRIPT: run_jenkins_perseus.sh
# AUTHOR: Kyle Gerard Felker - kfelker@princeton.edu
# DATE: 4/10/2018
# PURPOSE: Run regression test suite with PICSciE's Jenkins server, 4x Perseus
# computer Intel Broadwell worker nodes for continuous integration (CI)

# USAGE: salloc -N1 -n4 --time=0:60:00 ./run_jenkins_perseus.sh
# or similar command in the Jenkins build "Execute shell" step

set -e # quit at first error
cd tst/regression

# Build step #0: Test source code style consistency
cd ../style/; ./cpplint_athena.sh; cd ../regression/

# Build step #1: GNU compiler and OpenMPI library
module purge
module load rh # latest GNU compiler
module load openmpi/gcc # /1.10.2/64
# output/all_outputs.py regression test uses non-MPI HDF5 writer
# (Perseus will error w/ missing mpi.h header if MPI HDF5 is loaded w/o mpicxx)
module load hdf5/gcc/1.10.0  # openmpi-1.10.2/1.8.16
# grav/ regression tests require MPI and FFTW
module load fftw/gcc/3.3.4
module list

# Run regression test sets. Need to specify Slurm mpirun wrapper, srun
time python ./run_tests.py pgen --config=--cflag="$(../ci/set_warning_cflag.sh g++)"
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

# Build step #2: Intel compiler and MPI library
module purge
module load intel
module load intel-mpi
module load fftw/gcc/3.3.4
module load hdf5/intel-17.0/1.10.0 # hdf5/intel-17.0/intel-mpi/1.10.0
module load rh
module list

time python ./run_tests.py pgen --config=--cxx=icc --config=--cflag="$(../ci/set_warning_cflag.sh icc)"
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

set +e
# end regression tests

# Codecov coverage analysis
# Pipe to bash (Jenkins)
curl -s https://codecov.io/bash | bash -s - -t ccdc959e-e2c3-4811-95c6-512151b39471
