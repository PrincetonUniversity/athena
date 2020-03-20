#!/bin/sh
###############################################################################

###############################################################################
# Short script for taking care of (optional) compilation and running a problem
###############################################################################

###############################################################################
# configure here
export RUN_NAME=cvg_wave_1d_Gaussian
export BIN_NAME=wave
export REL_OUTPUT=outputs/wave
export REL_INPUT=scripts/problems
export INPUT_NAME=cvg_wave_1d_Gaussian.inp

# if compilation is chosen
export DIR_INST=$soft/usr  # correct for location of installed libraries

export COMPILE_STR="--prob=wave_1d_cvg_Gaussian_P -w --nghost=2
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
