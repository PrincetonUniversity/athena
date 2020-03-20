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
#else
    # only recompile what is required
#    cd ${DIR_ATHENA}
#    make -j6
#    mv bin/athena bin/${BIN_NAME}
fi

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


cd ${DIR_ATHENA}/${REL_OUTPUT}/${RUN_NAME}
# ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME} -m 1
# ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME}

# gdb -ex 'break wave_1d_cvg_trig.cpp:70' \
#     -ex 'info b' \
#     -ex 'set print pretty on' \
#     -ex=r --args ./$EXEC_NAME.x -i $DIR_ATHENA/$REL_INPUT/$INPUT_NAME


# gdb -ex 'break wave_1d_cvg_trig.cpp:104' \
gdb -q \
    -ex 'break bvals_refine.cpp:227' \
    -ex 'info b' \
    -ex 'set print pretty on' \
    -ex=r --args ./$EXEC_NAME.x -i $DIR_ATHENA/$REL_INPUT/$INPUT_NAME


# make call log [inspect with qcachegrind]
# valgrind --tool=callgrind --log-file=callgrind.log ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME}
#valgrind --leak-check=full -s --show-leak-kinds=all ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME}
# valgrind --leak-check=full -s ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME}

echo "Done >:D"
cd ${DIR_SCRIPTS}
###############################################################################

# >:D
