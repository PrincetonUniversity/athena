#!/bin/sh
###############################################################################

###############################################################################
# Short script for taking care of (optional) compilation and running a problem
###############################################################################

###############################################################################
# configure here
export RUN_NAME=cvg_wave_1d
export BIN_NAME=wave
export REL_OUTPUT=outputs/wave_c
export REL_INPUT=scripts/problems
export INPUT_NAME=cvg_wave_1d.inp

# if compilation is chosen
export DIR_INST=$soft/usr  # correct for location of installed libraries

# export COMPILE_STR="--prob=wave_1d_cvg_trig -w --nghost=2
#                     --cxx g++ -omp -debug
#                     -fill_wave_bnd_sl
#                     -fill_wave_bnd_frc
#                     -fill_wave_bnd_frf"
#                    -fill_wave_interior
#                    -hdf5 -h5double --hdf5_path=$DIR_INST"

export COMPILE_STR="--prob=wave_1d_cvg_trig -w --nghost=2
                    --cxx g++ -omp -debug"

# fill MeshBlock interior with exact solution
# export COMPILE_STR="${COMPILE_STR} -fill_wave_interior"

# fill MeshBlock boundaries [same level] with exact solution
# export COMPILE_STR="${COMPILE_STR} -fill_wave_bnd_sl"

# fill MeshBlock boundaries [from finer] with exact solution
# export COMPILE_STR="${COMPILE_STR} -fill_wave_bnd_frf"

# fill MeshBlock boundaries [from coarser] with exact solution
# export COMPILE_STR="${COMPILE_STR} -fill_wave_bnd_frc"

# fill MeshBlock coarse buffer with exact solution (prior to prolongation)
# export COMPILE_STR="${COMPILE_STR} -fill_wave_coarse_p"




###############################################################################

###############################################################################
# ensure paths are adjusted and directory structure exists
. utils/provide_paths.sh

###############################################################################
# compile
. utils/compile_force.sh
###############################################################################

# >:D
