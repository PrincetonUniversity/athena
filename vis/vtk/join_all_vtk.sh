#!/usr/bin/env bash

# SCRIPT: join_all_vtk.sh
# AUTHOR: Kyle Gerard Felker - kfelker@princeton.edu
# DATE:   4/9/2018
# PURPOSE:  Wrapper script to join_vtk++.c application to build list of per MeshBlock
# output files at each output step

# USAGE: ./join_all_vtk.sh problem_id output_id max_blockid max_outstep
# e.g. ./join_all_vtk.sh OrszagTang 2 3 101
# will call ./join_vtk++ -o OrszagTang.joined.out2.00000.vtk \
#           OrszagTang.block0.out2.00000.vtk OrszagTang.block1.out2.00000.vtk \
#           OrszagTang.block2.out2.00000.vtk OrszagTang.block3.out2.00000.vtk
# for all output steps from 00000 to 00101
# ========================================

# job/problem_id in athinput file
problem_id=$1
# integer index in <output[N]> vtk athinput block
output_id=$2
# seq FIRST LAST command includes LAST in return, so make this NBLOCKS-1
max_blockid=$3
max_outstep=$4

for i in $(seq 0 $max_outstep)
do
    # Format "output_step" as fixed width intger of length 5
    printf -v output_step "%05d" $i
    output_file="${problem_id}.joined.out${output_id}.${output_step}.vtk"
    input_files=()
    # Build list of each MPI rank's output files at this step
    for blockid in $(seq 0 $max_blockid)
    do
	input_files+=("${problem_id}.block${blockid}.out${output_id}.${output_step}.vtk")
    done
    #echo "./join_vtk++ -o ${output_file} ${input_files[*]}"
    ./join_vtk++ -o ${output_file} ${input_files[*]}
done
