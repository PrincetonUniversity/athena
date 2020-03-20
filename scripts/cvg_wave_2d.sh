#!/bin/sh
###############################################################################

###############################################################################
# Short script for taking care of (optional) compilation and running a problem
###############################################################################

###############################################################################
# configure here
export RUN_NAME=002_test
export BIN_NAME=wave_2d
export REL_OUTPUT=outputs/wave
export REL_INPUT=scripts/problems
export INPUT_NAME=cvg_wave_2d_trig.inp

# if compilation is chosen
export DIR_INST=$soft/usr  # correct for location of installed libraries
export COMPILE_STR="--prob=cvg_wave_2d_trig -w --nghost=2
                    --cxx=g++ -omp -debug
                    -hdf5 -h5double --hdf5_path=$DIR_INST"
###############################################################################

###############################################################################
# ensure paths are adjusted and directory structure exists
. utils/provide_paths.sh

###############################################################################
# compile
. utils/compile_run.sh
###############################################################################


# >:D
