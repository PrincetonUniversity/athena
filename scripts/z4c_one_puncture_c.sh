#!/bin/bash
###############################################################################

###############################################################################
# Short script for taking care of (optional) compilation and running a problem
###############################################################################
export FN=$(readlink -f "$0"); export DIR_SCRIPTS=$(dirname "${FN}")

###############################################################################
# configure here
export RUN_NAME=one_puncture
export BIN_NAME=z4c
export REL_OUTPUT=outputs/z4c_c
export REL_INPUT=scripts/problems
export INPUT_NAME=z4c_one_puncture.inp

# if compilation is chosen
export DIR_USR=${soft}/usr                            # local lib. installation
export DIR_HDF5=${DIR_USR}/hdf5_serial/1.10.5
export COMPILE_STR="--prob=z4c_one_puncture -z -z_tracker
                    --cxx g++ -omp -debug
                    --nghost=2"

# apply caching compiler together with gold linker
export COMPILE_STR="${COMPILE_STR} -ccache -link_gold"

# hdf5 compile str
export COMPILE_STR="${COMPILE_STR} -hdf5 -h5double --hdf5_path=${DIR_HDF5}"
###############################################################################


###############################################################################
# ensure paths are adjusted and directory structure exists
source ${DIR_SCRIPTS}/utils/provide_library_paths.sh ${DIR_HDF5}
source ${DIR_SCRIPTS}/utils/provide_compile_paths.sh
###############################################################################

###############################################################################
# prepare external
# ...
###############################################################################

###############################################################################
# compile
source ${DIR_SCRIPTS}/utils/compile_athena.sh
###############################################################################

###############################################################################
# dump information
source ${DIR_SCRIPTS}/utils/dump_info.sh
###############################################################################

###############################################################################
# execute
source utils/exec.sh
###############################################################################


# >:D