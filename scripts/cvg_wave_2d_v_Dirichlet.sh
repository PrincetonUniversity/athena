#!/bin/sh
###############################################################################

###############################################################################
# Short script for taking care of (optional) compilation and running a problem
###############################################################################

###############################################################################
# configure here
export RUN_NAME=cvg_wave_2d_Dirichlet
export BIN_NAME=wave
export REL_OUTPUT=outputs/wave_v
export REL_INPUT=scripts/problems
export INPUT_NAME=cvg_wave_2d_Dirichlet.inp

# if compilation is chosen
export DIR_INST=$soft/usr  # correct for location of installed libraries

export COMPILE_STR="--prob=wave_2d_cvg_Dirichlet -w
                    --cxx g++ -omp -debug -vertex
                    --nghost=3
                    --ncghost=4
                    --nextrapolate=7"

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

# debug vertex consistency
# export COMPILE_STR="${COMPILE_STR} -dbg_vc_consistency"

# add hdf5 support
export COMPILE_STR="${COMPILE_STR} -hdf5 -h5double --hdf5_path=$DIR_INST"

###############################################################################

###############################################################################
# ensure paths are adjusted and directory structure exists
. utils/provide_paths.sh

###############################################################################
# compile
. utils/compile_force.sh
###############################################################################

# >:D
