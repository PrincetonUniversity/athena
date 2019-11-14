#!/bin/sh
###############################################################################

###############################################################################
# Short script for taking care of (optional) compilation and running a problem
###############################################################################

###############################################################################
# configure here
export RUN_NAME=advection_1d_cvg_Gaussian_P
export BIN_NAME=advection
export REL_OUTPUT=outputs/advection
export REL_INPUT=scripts/problems
export INPUT_NAME=advection_1d_cvg_Gaussian_P.inp

# if compilation is chosen
export DIR_INST=$soft/usr  # correct for location of installed libraries

export COMPILE_STR="--prob=advection_1d_cvg_Gaussian_P -a --nghost=2
                    --cxx g++ -omp -debug
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
