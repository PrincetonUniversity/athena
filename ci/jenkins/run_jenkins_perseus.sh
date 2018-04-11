#!/usr/bin/env bash

# SCRIPT: run_jenkins_perseus.sh
# AUTHOR: Kyle Gerard Felker - kfelker@princeton.edu
# DATE: 4/10/2018
# PURPOSE: Run regression test suite in Jenkins continuous integration (CI)

# USAGE: salloc -N1 -n4 --time=0:60:00 ./run_jenkins_perseus.sh
# or similar command in the Jenkins build "Execute shell" step

set -e # quit at first error

# Only MPI compiler on PICSciE Perseus that passes MPI regression test without
# additional --cxx flag to ./configure.py. Need to make "#pragma omp simd"
# usage safe in code regardless of inlining of functions
module load openmpi/gcc # /1.10.2/64

cd tst/regression

# Run regression test sets
# Recall, need to modify tst/regression/scripts/utils.athena.py to use Slurm
# and correct Intel flags
python ./run_tests.py mpi
python ./run_tests.py hydro
python ./run_tests.py mhd # Longest regression test set
python ./run_tests.py amr
python ./run_tests.py outputs
python ./run_tests.py pgen
python ./run_tests.py sr
python ./run_tests.py gr
python ./run_tests.py curvilinear
#python ./run_tests.py shearingbox

module load fftw/gcc/3.3.4
python ./run_tests.py grav # requires MPI and FFTW

set +e
# end regression tests
