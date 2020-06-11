#!/bin/bash
###############################################################################

###############################################################################
# provide compilation as code-block
function compile {
    echo "Compiling..."
    cd ${DIR_ATHENA}
    rm -f bin/${BIN_NAME} > /dev/null 2>&1
    # remove executable at target if extant
    rm -f ${DIR_ATHENA}/${REL_OUTPUT}/${RUN_NAME}/${EXEC_NAME}.x \
       > /dev/null 2>&1
    make clean
    # extern_relink
    python configure.py $COMPILE_STR
    make -j8 &> ${DIR_OUTPUT}/compile_info
    mv bin/athena bin/${BIN_NAME}
} # &> /dev/null
###############################################################################

compile

###############################################################################
# copy to target
cp -f ${DIR_ATHENA}/bin/${BIN_NAME} \
   ${DIR_OUTPUT}/${EXEC_NAME}.x  > /dev/null 2>&1
rm -rf ${DIR_ATHENA}/bin/  > /dev/null 2>&1

cd ${DIR_SCRIPTS}
###############################################################################

# >:D
