#!/bin/bash
###############################################################################

###############################################################################
# Short script for taking care of (optional) compilation and running a problem
###############################################################################

function usage {
  echo "-d 'dimension' -ng 'number of ghosts'"
  exit 1
}

while [ "$1" != "" ]; do
  case $1 in
    -d | --dimension )      shift
                            dim=$1
                            ;;
    -ng | --nghost )        shift
                            nghost=$1
                            ;;
    -ncg | --ncghost )      shift
                            ncghost=$1
                            ;;
    -nex | --nextrapolate ) shift
                            nextrapolate=$1
                            ;;
    -v | --vertex )         vertex="-vertex"
                            ;;
    -D | --Dirichlet )      problem_subtype="Dirichlet"
                            ;;
    -h | --help )           usage
                            exit
                            ;;
    * )                     usage
                            exit 0
  esac
  shift
done

if [ "${dim}" == "" ]; then
  usage
  exit 0
fi

if [ "${nghost}" == "" ]; then
  usage
  exit 0
fi

if [ "${ncghost}" == "" ]; then
  ncghost=$(( ${nghost} + 1 ))
fi

if [ "${nextrapolate}" == "" ]; then
  nextrapolate=$(( 2 * ${nghost} + 1 ))
fi

if [ "${problem_subtype}" == "" ]; then
  problem_subtype=trig
fi


###############################################################################
# configure here
export RUN_NAME=cvg_wave_${dim}d_${problem_subtype}
export BIN_NAME=wave
export REL_OUTPUT=outputs/wave_c
export REL_INPUT=scripts/problems
export INPUT_NAME=cvg_wave_${dim}d.inp

# if compilation is chosen
export DIR_INST=$soft/usr  # correct for location of installed libraries

export COMPILE_STR="--prob=wave_${dim}d_cvg_${problem_subtype} -w ${vertex}
                    --nghost=${nghost}
                    --ncghost=${ncghost}
                    --nextrapolate=${nextrapolate}
                    --cxx g++ -omp -debug -ccache -link_gold"

# add hdf5 support
export COMPILE_STR="${COMPILE_STR} -hdf5 -h5double --hdf5_path=$DIR_INST"
###############################################################################

echo "Preparing ${dim}d wave equation"

###############################################################################
# ensure paths are adjusted and directory structure exists
source utils/provide_paths.sh

###############################################################################
# compile
source utils/compile_force.sh
###############################################################################

# >:D
