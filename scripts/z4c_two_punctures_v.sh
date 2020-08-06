#!/bin/bash
###############################################################################

###############################################################################
# Short script for taking care of (optional) compilation and running a problem
###############################################################################
export FN=$(readlink -f "$0"); export DIR_SCRIPTS=$(dirname "${FN}")

###############################################################################
# configure here
export RUN_NAME=two_puncture
export BIN_NAME=z4c
export REL_OUTPUT=outputs/z4c_c_tp
export REL_INPUT=scripts/problems

# Will be populated with defaults instead.
export INPUT_NAME=z4c_two_punctures.inp

# if compilation is chosen
export DIR_USR=${soft}/usr                            # local lib. installation
export DIR_HDF5=${DIR_USR}/hdf5_serial/1.10.5
export DIR_GSL=${DIR_USR}/gsl/2.6

export COMPILE_STR="--prob=z4c_two_punctures -z
                    -z_tracker
                    -z_eta_conf
                    -z_assert_is_finite
                    --cxx g++ -omp -vertex
                    --nghost=4 --ncghost=4 --nextrapolate=5"

# add wave extraction
export COMPILE_STR="${COMPILE_STR} -z_wext --nrad 2"

# apply caching compiler together with gold linker
export COMPILE_STR="${COMPILE_STR} -ccache -link_gold"

# hdf5 compile str
export COMPILE_STR="${COMPILE_STR} -hdf5 -h5double --hdf5_path=${DIR_HDF5}"

# gsl compile str [required for twopunctures]
export COMPILE_STR="${COMPILE_STR} -gsl --gsl_path=${DIR_GSL}"

# two punctures relative
export COMPILE_STR="${COMPILE_STR}
  --two_punctures_path=../../twopuncturesc/master"

echo "COMPILE_STR"
echo ${COMPILE_STR}
###############################################################################


###############################################################################
# ensure paths are adjusted and directory structure exists
source ${DIR_SCRIPTS}/utils/provide_library_paths.sh ${DIR_HDF5}
source ${DIR_SCRIPTS}/utils/provide_library_paths.sh ${DIR_GSL}
source ${DIR_SCRIPTS}/utils/provide_compile_paths.sh
###############################################################################

###############################################################################
# prepare external [two punctures; must be pre-compiled]
# source ${DIR_SCRIPTS}/utils/provide_extern.sh \
#   initial_data                                \
#   two_punctures                               \
#   ../../../../twopuncturesc/master
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

# stupid check out outer bc extrap
# export PYTHONPATH=$PYTHONPATH:/mnt/grottoop/_Repositories/NR/athena/development/vis/python
# export data_file="/mnt/grottoop/_Repositories/NR/athena/development/outputs/z4c_c_tp/two_puncture/z4c.out3.00002.athdf"
# export pycmd="import athena_read as ar; import numpy as np"
# export pycmd="${pycmd}; data = ar.athdf('${data_file}', raw=True)"
# export pycmd="${pycmd}; print(data['z4c.alpha'][-1,-1,-1])"

# echo ${pycmd} > pyout
# python pyout

# >:D
