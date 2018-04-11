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
module load rh # latest GNU compiler
module load openmpi/gcc # /1.10.2/64
# grav/ regression tests require MPI and FFTW
module load fftw/gcc/3.3.4
module list

# Run regression test sets. Need to specify Slurm mpirun wrapper, srun
python ./run_tests.py grav --mpirun=srun
python ./run_tests.py mpi --mpirun=srun
python ./run_tests.py hydro --mpirun=srun
# MHD is currenlty the longest regression test set:
python ./run_tests.py mhd --mpirun=srun
python ./run_tests.py amr --mpirun=srun
python ./run_tests.py outputs --mpirun=srun
python ./run_tests.py pgen --mpirun=srun
python ./run_tests.py sr --mpirun=srun
python ./run_tests.py gr --mpirun=srun
python ./run_tests.py curvilinear --mpirun=srun
#python ./run_tests.py shearingbox


# Build step #2: Intel compiler and MPI library
module purge
module load intel
module load intel-mpi
module load fftw/gcc/3.3.4
module list

python ./run_tests.py grav --cxx=icc --mpirun=srun
python ./run_tests.py mpi --cxx=icc --mpirun=srun
python ./run_tests.py hydro --cxx=icc --mpirun=srun
python ./run_tests.py mhd --cxx=icc --mpirun=srun
python ./run_tests.py amr --cxx=icc --mpirun=srun
python ./run_tests.py outputs --cxx=icc --mpirun=srun
python ./run_tests.py pgen --cxx=icc --mpirun=srun
python ./run_tests.py sr --cxx=icc --mpirun=srun
python ./run_tests.py gr --cxx=icc --mpirun=srun
python ./run_tests.py curvilinear --cxx=icc --mpirun=srun
#python ./run_tests.py shearingbox --cxx=icc

set +e
# end regression tests
