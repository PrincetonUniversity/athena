#!/usr/bin/env bash

# SCRIPT: run_jenkins_stellar.sh
# AUTHOR: Kyle Gerard Felker - felker@anl.gov
# AUTHOR: Chang-Goo Kim - changgoo@princeton.edu
# DATE: 6/4/2023
# PURPOSE: Run style and regression test suites using PICSciE's Jenkins server for continuous
# integration (CI) Current workers include 4x Intel Cascadelake processors on Stellar cluster.

# USAGE: salloc -N1 -n4 -c2 --kill-command=SIGTERM --time=9:00:00 \
#            --job-name=PrincetonUniversity_athena_jenkins_PR_$BUILD_NUMBER \
#            ./tst/ci/jenkins/run_jenkins_stellar.sh
# or similar command in the Jenkins build "Execute shell" step (run from athena/ root dir)

# Slurm diagnostics: see all timing info when job first exits the queue and actually starts
# Jenkins build time may be misleading, since it includes time sitting in Slurm queue.
sacct --jobs=$SLURM_JOB_ID --format=JobID,JobName%50,Submit,Start,Elapsed,Timelimit  #--noheader

git clean -f -d

set -e # terminate script at first error/non-zero exit status
# Store absolute path of project's root directory for Lcov (realpath is GNU coreutils, not macOS)
athena_rel_path='./'
athena_abs_path=$(realpath $athena_rel_path)

# Install Python dependencies
pip install -q --user --upgrade setuptools # required for matplotlib
pip install -q --user flake8 colorama termcolor matplotlib
pip install -q --user h5py    # needed for outputs/all_outputs.py, pgen/hdf5*, eos/eos_hdf5_table.py tests
pip install -q --user scipy   # needed in scripts/utils/ for eos/ tests

# Build step #0: Test source code style consistency
# step #0a: lint Python files
python -m flake8
echo "Finished linting Python files with flake8"

# step #0b: lint C++ files
cd tst/style/; ./check_athena_cpp_style.sh
cd ../regression/

# Build step #1a: pgen_compile test using Intel compiler and MPI library
module purge
module load anaconda3/2023.3 intel/2022.2.0 intel-mpi/intel/2021.7.0 hdf5/intel-2021.1/1.10.6 fftw/intel-2021.1/3.3.9
#module load anaconda3/2020.11 intel/2021.1.2 intel-mpi/intel/2021.1.1 hdf5/intel-2021.1/1.10.6 fftw/intel-2021.1/3.3.9
module list

time python -u ./run_tests.py pgen/pgen_compile --config=--cxx=icpx --config=--cflag="$(../ci/set_warning_cflag.sh icpx)"

# Build step #1b: pgen_compile test using GNU compiler and OpenMPI library
module purge
module load anaconda3/2020.11 gcc-toolset/10 openmpi/gcc/4.1.0 fftw/gcc/3.3.9 hdf5/gcc/1.10.6
module list

# Build step #2: regression tests using GNU compiler and OpenMPI library
# Run regression test sets. Need to specify Slurm mpirun wrapper, srun
# In order to condense the build log, --silent option suppresses only the stdout of Makefile calls. Don't use with pgen_compile.py:
time python -u ./run_tests.py pgen/pgen_compile --config=--cflag="$(../ci/set_warning_cflag.sh g++)"
time python -u ./run_tests.py pgen/hdf5_reader_serial --silent
time python -u ./run_tests.py grav --mpirun=srun --mpirun_opts=--job-name='GCC grav/jeans_3d' --silent
time python -u ./run_tests.py turb --mpirun=srun --mpirun_opts=--job-name='GCC turb/' --silent
time python -u ./run_tests.py mpi --mpirun=srun --mpirun_opts=--job-name='GCC mpi/mpi_linwave' --silent
time python -u ./run_tests.py omp --silent
timeout --signal=TERM 60m time python -u ./run_tests.py hybrid --mpirun=srun \
	--mpirun_opts=--job-name='GCC hybrid/hybrid_linwave' \
	--silent
time python -u ./run_tests.py hydro --silent
time python -u ./run_tests.py amr --silent
time python -u ./run_tests.py outputs --silent
time python -u ./run_tests.py curvilinear --silent
time python -u ./run_tests.py symmetry --silent
time python -u ./run_tests.py eos --silent
time python -u ./run_tests.py scalars/mignone_radial_1d --silent
# Exclude gr/compile*.py regression tests from code coverage analysis (nothing is executed in these tests):
time python -u ./run_tests.py sr --silent
time python -u ./run_tests.py gr/compile_kerr-schild gr/compile_minkowski gr/compile_schwarzschild --silent
time python -u ./run_tests.py gr/mhd_shocks_hlld gr/mhd_shocks_hlle gr/mhd_shocks_llf --silent
time python -u ./run_tests.py gr/hydro_shocks_hllc gr/hydro_shocks_hlle gr/hydro_shocks_llf --silent
time python -u ./run_tests.py gr/hydro_shocks_hlle_no_transform gr/hydro_shocks_llf_no_transform --silent
#
time python -u ./run_tests.py mhd --silent  # (mhd/mhd_linwave.py is currenlty the slowest regression test):
time python -u ./run_tests.py shearingbox --silent
time python -u ./run_tests.py diffusion --silent
time python -u ./run_tests.py cr --silent
time python -u ./run_tests.py nr_radiation --silent
time python -u ./run_tests.py implicit_radiation --silent
time python -u ./run_tests.py multi_group --silent

# ----------------
# Install CVODE
wget https://github.com/LLNL/sundials/releases/download/v6.2.0/cvode-6.2.0.tar.gz
tar -xvf cvode-6.2.0.tar.gz
mkdir cvode_build || true
mkdir cvode_install || true
export CVODE_PATH=$PWD/cvode_install
cd cvode_build

cmake -D CMAKE_INSTALL_LIBDIR=lib -D CMAKE_INSTALL_PREFIX=${CVODE_PATH} -D ENABLE_MPI=ON -D ENABLE_OPENMP=ON -D CMAKE_C_COMPILER=gcc -D MPI_C_COMPILER=mpicc -D MPIEXEC_EXECUTABLE=srun -D CMAKE_BUILD_TYPE=RelWithDebInfo ../cvode-6.2.0

# cmake -MAKE_INSTALL_LIBDIR=${CVODE_PATH} -ENABLE_MPI=ON -ENABLE_OPENMP=ON -MAKE_C_COMPILER=icx -MPI_C_COMPILER=mpicc -MPIEXEC_EXECUTABLE=srun -MAKE_BUILD_TYPE=RelWithDebInfo ../cvode-6.2.0
make -j16
make install
cd ..
# ----------------

time python -u ./run_tests.py chemistry --silent

# High-order solver regression tests w/ GCC
time python -u ./run_tests.py hydro4 --silent

# Swap serial HDF5 library module for parallel HDF5 library:
module unload hdf5/gcc/1.10.6
module load hdf5/gcc/openmpi-4.1.0/1.10.6
mpi_hdf5_library_path='/usr/local/hdf5/gcc/openmpi-4.1.0/1.10.6/lib64'
module list

# Workaround issue with parallel HDF5 modules compiled with OpenMPI on Perseus--- linker still chooses serial HDF5 library in /usr/lib64/
# due to presence of -L flag in mpicxx wrapper that overrides LIBRARY_PATH environment variable
time python -u ./run_tests.py pgen/hdf5_reader_parallel \ #--coverage="${lcov_capture_cmd}" \
     --mpirun=srun --mpirun_opts=--job-name='GCC pgen/hdf5_reader_parallel' \
     --config=--lib_path=${mpi_hdf5_library_path} --silent

# Build step #2: regression tests using Intel compiler and MPI library
module purge
module load anaconda3/2023.3 intel/2022.2.0 intel-mpi/intel/2021.7.0 hdf5/intel-2021.1/1.10.6 fftw/intel-2021.1/3.3.9
#module load anaconda3/2020.11 intel/2021.1.2 intel-mpi/intel/2021.1.1 hdf5/intel-2021.1/1.10.6 fftw/intel-2021.1/3.3.9
module list

time python -u ./run_tests.py pgen/hdf5_reader_serial --silent --config=--cxx=icpx
time python -u ./run_tests.py grav --config=--cxx=icpx --mpirun=srun --mpirun_opts=--job-name='ICC grav/jeans_3d' --silent
#(changgoo) use icpc instead of icpx to make mpi and serial solutions identical
#(KGF)      stellar intel-mpi/ modules are built against icpc, not icpx
time python -u ./run_tests.py turb --config=--cxx=icpc --mpirun=srun --mpirun_opts=--job-name='ICC turb/' --silent
time python -u ./run_tests.py mpi --config=--cxx=icpc --mpirun=srun --mpirun_opts=--job-name='ICC mpi/mpi_linwave' --silent
time python -u ./run_tests.py omp --config=--cxx=icpx --silent
timeout --signal=TERM 60m time python -u ./run_tests.py hybrid --config=--cxx=icpc \
	--mpirun=srun --mpirun_opts=--job-name='ICC hybrid/hybrid_linwave' --silent
time python -u ./run_tests.py hydro --config=--cxx=icpx --silent
time python -u ./run_tests.py mhd --config=--cxx=icpx --silent
time python -u ./run_tests.py amr --config=--cxx=icpx --silent
time python -u ./run_tests.py outputs --config=--cxx=icpx --silent
#(changgoo) suppress gr/sr tests
time python -u ./run_tests.py sr --config=--cxx=icpx --silent
time python -u ./run_tests.py gr --config=--cxx=icpx --silent
time python -u ./run_tests.py curvilinear --config=--cxx=icpx --silent
time python -u ./run_tests.py shearingbox --config=--cxx=icpx --silent
#(changgoo) suppress diffusion test
# time python -u ./run_tests.py diffusion --config=--cxx=icpc --silent
time python -u ./run_tests.py symmetry --config=--cxx=icpx --silent
time python -u ./run_tests.py eos --config=--cxx=icpx --silent
time python -u ./run_tests.py scalars/mignone_radial_1d --config=--cxx=icpx --silent

time python -u ./run_tests.py cr --config=--cxx=icpx --silent
time python -u ./run_tests.py nr_radiation --config=--cxx=icpx --silent
time python -u ./run_tests.py implicit_radiation --config=--cxx=icpx --silent
time python -u ./run_tests.py multi_group --config=--cxx=icpx --silent

time python -u ./run_tests.py chemistry --config=--cxx=icpx --silent

# High-order solver regression tests w/ Intel compiler
time python -u ./run_tests.py hydro4 --config=--cxx=icpx --silent

# Swap serial HDF5 library module for parallel HDF5 library:
module unload hdf5/intel-2021.1/1.10.6
module load hdf5/intel-2021.1/intel-mpi/1.10.6
mpi_hdf5_library_path='/usr/local/hdf5/intel-2021.1/intel-mpi/1.10.6/lib64'
# TODO: get newer serial and parallel HDF5 libraries compiled with ICX 2022 or 2023
# module unload hdf5/intel-2022.2/1.10.6
# module load hdf5/intel-2022.2/intel-mpi/1.10.6
# mpi_hdf5_library_path='/usr/local/hdf5/intel-2022.2/intel-mpi/1.10.6/lib64'

module list
# Workaround issue with parallel HDF5 modules compiled with OpenMPI on Perseus--- linker still takes serial HDF5 library in /usr/lib64/
# due to presence of -L flag in mpicxx wrapper that overrides LIBRARY_PATH environment variable
time python -u ./run_tests.py pgen/hdf5_reader_parallel --config=--cxx=icpc \
     --mpirun=srun --mpirun_opts=--job-name='ICC pgen/hdf5_reader_parallel' \
     --config=--lib_path=${mpi_hdf5_library_path} --silent

set +e
# end regression tests

# Slurm diagnostics: see all timing info when build script finishes
# (should run in Jenkins "Execute shell" build step when Slurm allocation is released)
sacct --jobs=$SLURM_JOB_ID --format=JobID,JobName%50,Submit,Start,Elapsed,Timelimit,End  #--noheader
