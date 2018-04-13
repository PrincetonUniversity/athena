#!/usr/bin/env bash

# SCRIPT: join_all_vtk.sh
# AUTHOR: Kyle Gerard Felker - kfelker@princeton.edu
# DATE:   4/9/2018
# PURPOSE:  Wrapper script to join_vtk++.c application to build list of per MeshBlock
# output files at each output step

# USAGE: ./join_all_vtk.sh PROBLEM_ID OUTPUT_ID MAX_BLOCKID MAX_OUTSTEP
# e.g. ./join_all_vtk.sh OrszagTang 2 3 101
# will call ./join_vtk++ -o OrszagTang.joined.out2.00000.vtk \
#           OrszagTang.block0.out2.00000.vtk OrszagTang.block1.out2.00000.vtk \
#           OrszagTang.block2.out2.00000.vtk OrszagTang.block3.out2.00000.vtk
# for all output steps from 00000 to 00101
# ========================================

# job/problem_id in athinput file
PROBLEM_ID=$1
# integer index in <output[N]> vtk athinput block
OUTPUT_ID=$2
# seq FIRST LAST command includes LAST in return, so make this NBLOCKS-1
MAX_BLOCKID=$3
MAX_OUTSTEP=$4

for i in $(seq 0 $MAX_OUTSTEP)
do
    # Format Output Step as fixed width intger of length 5
    printf -v output_step "%05d" $i
    output_file="${PROBLEM_ID}.joined.out${OUTPUT_ID}.${output_step}.vtk"
    input_files=()
    # Build list of each MPI rank's output files at this step
    for blockid in $(seq 0 $MAX_BLOCKID)
    do
	input_files+=("${PROBLEM_ID}.block${blockid}.out${OUTPUT_ID}.${output_step}.vtk")
    done
    #echo "./join_vtk++ -o ${output_file} ${input_files[*]}"
    ./join_vtk++ -o ${output_file} ${input_files[*]}
done
