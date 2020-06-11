#!/bin/bash
###############################################################################

###############################################################################
# determine location and set up (relative) paths

# export LOC_SCRIPT=$(readlink -f "$0")
# export DIR_SCRIPTS=$(dirname "$LOC_SCRIPT")

# make directory structure for compilation
cd ${DIR_SCRIPTS}
cd ..
mkdir -p ${REL_OUTPUT}/${RUN_NAME} > /dev/null 2>&1

export DIR_ATHENA=$PWD
export DIR_OUTPUT=${DIR_ATHENA}/${REL_OUTPUT}/${RUN_NAME}

export EXEC_NAME=${BIN_NAME}"_"${RUN_NAME}

# cd ${DIR_SCRIPTS}
###############################################################################

# >:D
