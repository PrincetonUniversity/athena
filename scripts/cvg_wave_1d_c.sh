#!/bin/sh
###############################################################################

###############################################################################
# Short script for taking care of (optional) compilation and running a problem
###############################################################################

###############################################################################
# configure here
export RUN_NAME=cvg_wave_1d
export BIN_NAME=wave
export REL_OUTPUT=outputs/wave_c
export REL_INPUT=scripts/problems
export INPUT_NAME=cvg_wave_1d.inp

# if compilation is chosen
export DIR_INST=$soft/usr  # correct for location of installed libraries

export COMPILE_STR="--prob=wave_1d_cvg_trig -w --nghost=2
                    --cxx g++ -omp -debug"

# add hdf5 support
export COMPILE_STR="${COMPILE_STR} -hdf5 -h5double --hdf5_path=$DIR_INST"

###############################################################################

###############################################################################
# ensure paths are adjusted and directory structure exists
. utils/provide_paths.sh

###############################################################################
# compile
. utils/compile_force.sh
###############################################################################

# >:D
