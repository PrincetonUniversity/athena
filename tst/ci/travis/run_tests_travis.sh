#!/bin/bash

# Fail at first error & echo commands before executing
set -ev
# Do not use "python run_tests.py" to run all tests, so that:
# - Each set/directory of tests are timed separately
# - Script fails after first broken set of tests
# (Could alternatively group sets of tests with && operator)
cd regression/
# python3 run_tests.py pgen --config=--cxx=$TEMP_CXX --config=--cflag="$(./ci/set_warning_cflag.sh $TEMP_CXX)"

python3 run_tests.py mpi --config=--cxx=$TEMP_CXX --mpirun='mpirun --oversubscribe' # --silent
python3 run_tests.py grav/jeans_3d --config=--cxx=$TEMP_CXX --mpirun='mpirun --oversubscribe' # requires FFTW library

# python3 run_tests.py grav/unstable_jeans_3d_fft grav/unstable_jeans_3d_mg --config=--cxx=$TEMP_CXX --silent # requires FFTW library
# python3 run_tests.py amr --config=--cxx=$TEMP_CXX --silent

# mhd/ contains the longest set of tests. Timeout after 10 m on Travis CI
# python3 run_tests.py mhd

# python3 run_tests.py hydro --config=--cxx=$TEMP_CXX --silent
# python3 run_tests.py outputs --config=--cxx=$TEMP_CXX --silent
# python3 run_tests.py curvilinear --config=--cxx=$TEMP_CXX --silent
# python3 run_tests.py gr --config=--cxx=$TEMP_CXX --silent
# python3 run_tests.py sr --config=--cxx=$TEMP_CXX --silent
# python3 run_tests.py shearingbox --config=--cxx=$TEMP_CXX --silent
# python3 run_tests.py diffusion --config=--cxx=$TEMP_CXX --silent
