#!/bin/sh
###############################################################################

###############################################################################
# Short script for taking care of (optional) compilation and running a problem
###############################################################################

###############################################################################
# configure here
export RUN_NAME=one_puncture
export BIN_NAME=z4c
export REL_OUTPUT=outputs/z4c_v
export REL_INPUT=scripts/problems
export INPUT_NAME=z4c_one_puncture.inp

# if compilation is chosen
export DIR_INST=$soft/usr  # correct for location of installed libraries
export COMPILE_STR="--prob=z4c_one_puncture -z
                    --cxx g++ -omp -debug -vertex
                    --nghost=2
                    --ncghost=2
                    --nextrapolate=4"

# add hdf5 support
export COMPILE_STR="${COMPILE_STR} -hdf5 -h5double --hdf5_path=$DIR_INST"

###############################################################################

###############################################################################
# ensure paths are adjusted and directory structure exists
. utils/provide_paths.sh

###############################################################################
# submit (compile if required)
. utils/compile_force.sh
###############################################################################

# >:D
