#!/bin/bash
###############################################################################

###############################################################################
# execute (assemble and provide info.)

echo ""
echo "> tail compile_info:"
tail -n 5 ${DIR_OUTPUT}/compile_info
echo ""
echo "> Prepared:"
echo "  ${EXEC_NAME} in ${REL_OUTPUT}/${RUN_NAME}"
echo ""
echo "> Switch to directory [abs]:"
echo "  ${DIR_ATHENA}/${REL_OUTPUT}/${RUN_NAME}"
echo ""
echo "> Input/Par:"
echo "  ${REL_INPUT}/${INPUT_NAME}"
echo ""
echo "> Execute command:"
echo "  ./${EXEC_NAME}.x -i ${DIR_ATHENA}/${REL_INPUT}/${INPUT_NAME}"
echo ""

echo "Done >:D"
###############################################################################

# >:D
