#!/bin/sh
###############################################################################

# ensure environmental variables are correctly adjusted and provide
# directory structures

# add paths to local libs
export PATH=${DIR_INST}/bin:$PATH
export LD_LIBRARY_PATH=${DIR_INST}/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=${DIR_INST}/lib:$LIBRARY_PATH                   # for make
export C_INCLUDE_PATH=${DIR_INST}/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=${DIR_INST}/include:$CPLUS_INCLUDE_PATH

###############################################################################
# determine location and set up paths

export LOC_SCRIPT=$(readlink -f "$0")
export DIR_SCRIPTS=$(dirname "$LOC_SCRIPT")

# make structure
cd ${DIR_SCRIPTS}
cd ..
mkdir -p ${REL_OUTPUT}/${RUN_NAME} > /dev/null 2>&1

export DIR_ATHENA=$PWD
export DIR_OUTPUT=${DIR_ATHENA}/${REL_OUTPUT}/${RUN_NAME}

export EXEC_NAME=${BIN_NAME}"_"${RUN_NAME}


cd ${DIR_SCRIPTS}

# >:D
