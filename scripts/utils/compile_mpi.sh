#!/bin/sh
###############################################################################


# Optionally compile
read -p "Compilation required [Yy]?"
TO_COMPILE=$REPLY
echo $TO_COMPILE

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


if [[ $TO_COMPILE =~ ^[Yy]$ ]]
then
    compile
fi

###############################################################################
# copy to target
cp -f ${DIR_ATHENA}/bin/${BIN_NAME} \
   ${DIR_OUTPUT}/${EXEC_NAME}.x  > /dev/null 2>&1
rm -rf ${DIR_ATHENA}/bin/  > /dev/null 2>&1

# copy .inp and .sub
cp -f ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME} \
   ${DIR_OUTPUT}/${INPUT_NAME} > /dev/null 2>&1

cp -f ${DIR_ATHENA}/${REL_INPUT}/${SUBMIT_NAME} \
   ${DIR_OUTPUT}/${SUBMIT_NAME} > /dev/null 2>&1
###############################################################################


###############################################################################
# execute

echo "> Submitting: ${EXEC_NAME} in ${REL_OUTPUT}/${RUN_NAME} ..."
echo "> Using input: ${REL_INPUT}/${INPUT_NAME} ..."


cd ${DIR_ATHENA}/${REL_OUTPUT}/${RUN_NAME}

# use as arg to .sub
export EXEC_CMD="./${EXEC_NAME}.x -i ${INPUT_NAME}"
# ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME} -m 1
# ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME}

echo $PWD
echo $EXEC_CMD
# sbatch $SUBMIT_NAME $EXEC_CMD

echo "Done >:D"
cd ${DIR_SCRIPTS}
###############################################################################

# >:D
