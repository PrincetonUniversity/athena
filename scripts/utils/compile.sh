#!/bin/sh
###############################################################################

# relink externally included libraries
function extern_relink {
    echo "Linking external libraries to obj:"
    find -L $DIR_ATHENA/extern -name \*.o \
         -exec ln -s "{}" $DIR_ATHENA/obj \; -exec ls "{}" \;
}


# compile if required..
function compile {
    echo "compile..."
    cd ${DIR_ATHENA}
    rm -f bin/${BIN_NAME} > /dev/null 2>&1
    # remove executable at target if extant
    rm -f ${DIR_ATHENA}/${REL_OUTPUT}/${RUN_NAME}/${EXEC_NAME}.x \
       > /dev/null 2>&1
    make clean
    # extern_relink
    python configure.py $COMPILE_STR
    make -j8
    mv bin/athena bin/${BIN_NAME}
} # &> /dev/null

compile

###############################################################################
# copy to target
cp -f ${DIR_ATHENA}/bin/${BIN_NAME} \
   ${DIR_OUTPUT}/${EXEC_NAME}.x  > /dev/null 2>&1
rm -rf ${DIR_ATHENA}/bin/  > /dev/null 2>&1
###############################################################################


###############################################################################
# execute

echo "> Executing: ${EXEC_NAME} in ${REL_OUTPUT}/${RUN_NAME} ..."
echo "> Using input: ${REL_INPUT}/${INPUT_NAME} ..."

echo "> Prepared:"
echo "  ${EXEC_NAME} in ${REL_OUTPUT}/${RUN_NAME}"
echo "> Switching to directory:"
echo "  ${DIR_ATHENA}/${REL_OUTPUT}/${RUN_NAME}"
echo "> Execution command:"
echo "  ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME}"

echo "Done >:D"
cd ${DIR_SCRIPTS}
###############################################################################

# >:D
