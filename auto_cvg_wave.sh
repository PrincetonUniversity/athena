#!/bin/sh

#################################################################
# Automatic compilation on change of files (recursive monitoring)
#
# Requires:
#  inotify-tools
#
# Example: execute 1d wave equation for cell-centered nodes
#  bash auto_cvg_wave.sh -d 1 -n c
#
# Note:
#  root execution and compilation directory specified below
#################################################################
# defaults
dimension=1
node_type=c
MONITOR_DIRECTORY=$PWD                              # just use the current dir
COMPILATION_DIRECTORY=/dev/shm/athena/vertex_test

function usage {
    echo "usage: [[[-d dimension (1-3) ] [-n node_type (c or v)]] | [-h]]"
}

while [ "$1" != "" ]; do
    case $1 in
        -d | --dimension )      shift
                                dimension=$1
                                ;;
        -n | --node_type )      shift
                                node_type=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

# script that will be run
BASE_FILENAME=cvg_wave_${dimension}d_${node_type}.sh
SCRIPTS_DIRECTORY=${COMPILATION_DIRECTORY}/cmp_${node_type}_${dimension}d/scripts

function runscript {
    echo "> Script: ${BASE_FILENAME}"
    rm -rf ${COMPILATION_DIRECTORY}/cmp_${node_type}_${dimension}d/src
    rm -rf ${COMPILATION_DIRECTORY}/cmp_${node_type}_${dimension}d/scripts

    rm ${COMPILATION_DIRECTORY}/cmp_${node_type}_${dimension}d/Make*
    cp ${MONITOR_DIRECTORY}/Make* ${COMPILATION_DIRECTORY}/cmp_${node_type}_${dimension}d/
    cp ${MONITOR_DIRECTORY}/*.py ${COMPILATION_DIRECTORY}/cmp_${node_type}_${dimension}d/

    cp -R src/ ${COMPILATION_DIRECTORY}/cmp_${node_type}_${dimension}d/
    cp -R scripts/ ${COMPILATION_DIRECTORY}/cmp_${node_type}_${dimension}d/
    cd ${SCRIPTS_DIRECTORY}

    output_dump=${COMPILATION_DIRECTORY}/cmp_${node_type}_${dimension}d/out
    bash ${BASE_FILENAME} > ${output_dump}
    tail -n 10000 ${output_dump}
    python ${SCRIPTS_DIRECTORY}/testing/parse_output.py -p ${output_dump}
    # {
    # } &> /dev/null
    cd ${MONITOR_DIRECTORY}
}

function remove {
    echo "> kill exec."
    {
        rm -rf ${COMPILATION_DIRECTORY}/cmp_${node_type}_${dimension}d/obj
        rm -rf ${COMPILATION_DIRECTORY}/cmp_${node_type}_${dimension}d/outputs
    } &> /dev/null

}

# prepare /dev/shm
function prepare_env {
    echo "> prepare compilation environment"
    {
      #rm -rf ${COMPILATION_DIRECTORY}
      mkdir -p ${COMPILATION_DIRECTORY}
      mkdir -p ${COMPILATION_DIRECTORY}/cmp_${node_type}_${dimension}d
    }
}


# main loop here
while true
do
	inotifywait --timefmt '%d/%m/%y %H:%M:%S' --format '%T %w %f' \
	--quiet \
	--recursive \
        --exclude '(.*#.*|.*pyg|.*log|.*fls|.*athdf)' \
	-r -e close_write ${MONITOR_DIRECTORY} \
		| while read date time dir file; do
		STATUS="Modification observed @ ${time} on ${date}..."
		echo "$STATUS"

    sleep 1  # allow oper. completion
    prepare_env
		runscript
    remove
	done
done
