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

# Build step #1: GNU compiler and OpenMPI library
module purge
module load rh # latest GNU compiler
module load openmpi/gcc # /1.10.2/64
# grav/ regression tests require MPI and FFTW
module load fftw/gcc/3.3.4
module list

# Run regression test sets. Need to specify Slurm mpirun wrapper, srun
time python ./run_tests.py pgen --mpirun=srun --cflag="-Wall -Wextra -Wno-unused-private-field -Wno-unused-variable -Wno-unknown-pragmas -Wno-unused-function -Werror"
time python ./run_tests.py grav --mpirun=srun
time python ./run_tests.py mpi --mpirun=srun
time python ./run_tests.py hydro --mpirun=srun
# MHD is currenlty the longest regression test set:
time python ./run_tests.py mhd --mpirun=srun
time python ./run_tests.py amr --mpirun=srun
time python ./run_tests.py outputs --mpirun=srun
time python ./run_tests.py sr --mpirun=srun
time python ./run_tests.py gr --mpirun=srun
time python ./run_tests.py curvilinear --mpirun=srun
time python ./run_tests.py shearingbox --mpirun=srun

# Build step #2: Intel compiler and MPI library
module purge
module load intel
module load intel-mpi
module load fftw/gcc/3.3.4
module load rh
module list

time python ./run_tests.py pgen --cxx=icc --mpirun=srun --cflag="-Wall -Wextra -Wno-unused-private-field -Wno-unused-variable -Wno-unknown-pragmas -Wno-unused-function -Werror"
time python ./run_tests.py grav --cxx=icc --mpirun=srun
time python ./run_tests.py mpi --cxx=icc --mpirun=srun
time python ./run_tests.py hydro --cxx=icc --mpirun=srun
time python ./run_tests.py mhd --cxx=icc --mpirun=srun
time python ./run_tests.py amr --cxx=icc --mpirun=srun
time python ./run_tests.py outputs --cxx=icc --mpirun=srun
time python ./run_tests.py sr --cxx=icc --mpirun=srun
time python ./run_tests.py gr --cxx=icc --mpirun=srun
time python ./run_tests.py curvilinear --cxx=icc --mpirun=srun
time python ./run_tests.py shearingbox --cxx=icc --mpirun=srun

set +e
# end regression tests

# Codecov coverage analysis
# Pipe to bash (Jenkins)
curl -s https://codecov.io/bash | bash -s - -t ccdc959e-e2c3-4811-95c6-512151b39471
