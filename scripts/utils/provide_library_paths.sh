#!/bin/bash
###############################################################################

export PATH=${1}/bin:$PATH
export LD_LIBRARY_PATH=${1}/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=${1}/lib:$LIBRARY_PATH
export C_INCLUDE_PATH=${1}/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=${1}/include:$CPLUS_INCLUDE_PATH

###############################################################################

# >:D
