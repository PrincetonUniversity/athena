#!/bin/bash
###############################################################################

###############################################################################
# execute
# cd ${DIR_ATHENA}/${REL_OUTPUT}/${RUN_NAME}
cd ${DIR_OUTPUT}

# ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME}

# echo "> Executing: ${EXEC_NAME} in ${REL_OUTPUT}/${RUN_NAME} ..."
# echo "> Using input: ${REL_INPUT}/${INPUT_NAME} ..."


# cd ${DIR_ATHENA}/${REL_OUTPUT}/${RUN_NAME}
# ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME} -m 1
./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME}

# gdb -ex 'info b' \
#     -ex 'set print pretty on' \
#     -ex=r --args ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME}


# gdb -ex 'break wave_1d_cvg_trig.cpp:104' \
# gdb -q \
#     -ex 'break bvals_refine.cpp:227' \
#     -ex 'info b' \
#     -ex 'set print pretty on' \
#     -ex=r --args ./$EXEC_NAME.x -i $DIR_ATHENA/$REL_INPUT/$INPUT_NAME


# make call log [inspect with qcachegrind]
#valgrind --tool=callgrind --log-file=callgrind.log ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME}
# valgrind --leak-check=full -s --show-leak-kinds=all ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME}
# valgrind --leak-check=full -s ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME} -m 2

echo "Done >:D"
cd ${DIR_SCRIPTS}
###############################################################################

# >:D
