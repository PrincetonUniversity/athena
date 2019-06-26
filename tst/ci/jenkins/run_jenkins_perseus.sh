#!/usr/bin/env bash

# SCRIPT: run_jenkins_perseus.sh
# AUTHOR: Kyle Gerard Felker - kfelker@princeton.edu
# DATE: 4/10/2018
# PURPOSE: Run style and regression test suites using PICSciE's Jenkins server for continuous
# integration (CI) Current workers include 4x Intel Broadwell nodes on Perseus cluster.

# USAGE: salloc -N1 -n4 -c2 --kill-command=SIGTERM --time=9:00:00 \
#            --job-name=PrincetonUniversity_athena_jenkins_PR_$BUILD_NUMBER \
#            ./tst/ci/jenkins/run_jenkins_perseus.sh
# or similar command in the Jenkins build "Execute shell" step (run from athena/ root dir)

# Slurm diagnostics: see all timing info when job first exits the queue and actually starts
# Jenkins build time may be misleading, since it includes time sitting in Slurm queue.
sacct --jobs=$SLURM_JOB_ID --format=JobID,JobName%50,Submit,Start,Elapsed,Timelimit  #--noheader

set -e # terminate script at first error/non-zero exit status
# Store absolute path of project's root directory for Lcov (realpath is GNU coreutils, not macOS)
athena_rel_path='./'
athena_abs_path=$(realpath $athena_rel_path)

# Install Python dependencies
pip install -q --user --upgrade setuptools # required for matplotlib
pip install -q --user flake8 colorama termcolor matplotlib
pip install -q --user h5py    # needed for outputs/all_outputs.py, pgen/hdf5*, eos/eos_hdf5_table.py tests
pip install -q --user scipy   # needed in scripts/utils/ for eos/ tests

# module load anaconda  # Python 2
# conda create -n athena
# conda activate athena
# conda install termcolor flake8
# included with conda distro's base environment:
#conda install h5py scipy colorama matplotlib

# Build step #0: Test source code style consistency
# step #0a: lint Python files
python -m flake8
echo "Finished linting Python files with flake8"

# step #0b: lint C++ files
cd tst/style/; ./check_athena_cpp_style.sh
cd ../regression/

# Build step #1: regression tests using GNU compiler and OpenMPI library
module purge
# Load latest GCC built by Red Hat Developer Toolset:
# (vs. /usr/bin/gcc v4.8.5 (released 2015-06-23)
module load rh/devtoolset/7  # GCC 7.3.1 (v7.3 released on 2018-01-25)
#module load openmpi/gcc/1.10.2/64  # OpenMPI v1.10.2 released on 2016-01-21
module load openmpi/gcc/3.0.3/64
# OpenMPI v3.0.3 was released on 2018-10-29
# OpenMPI v3.0.0 was released on 2017-09-12. Originally, was only installed on Perseus
# without development files (mpicc, etc.) as a VisIt 2.13.1 dependency

# Do NOT "module load hdf5" = hdf5/intel-17.0/openmpi-1.10.2/1.10.0
# output/all_outputs.py regression test uses non-MPI HDF5 writer
# (Perseus will error w/ missing mpi.h header if MPI HDF5 is loaded w/o mpicxx)
# grav/ regression tests require MPI and FFTW
module load hdf5/gcc/1.10.0
module load fftw/gcc/3.3.4
module list

# Lcov command stub used for capturing tracefile and for combining multiple tracefiles:
lcov_cmd="lcov --rc lcov_branch_coverage=1 --no-external --gcov-tool=gcov"
regression_abs_path=$(pwd)
lcov_capture_cmd="${lcov_cmd} --directory=${regression_abs_path}/obj/ --capture --base-directory=${athena_abs_path}"

# Run regression test sets. Need to specify Slurm mpirun wrapper, srun
# In order to condense the build log, --silent option suppresses only the stdout of Makefile calls. Don't use with pgen_compile.py:
time python -u ./run_tests.py pgen/pgen_compile --config=--cflag="$(../ci/set_warning_cflag.sh g++)"
# For (most) regression tests compiled with GCC, perform Gcov code coverage analysis via Lcov front end:
time python -u ./run_tests.py pgen/hdf5_reader_serial --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py grav --mpirun=srun --mpirun_opts=--job-name='GCC grav/jeans_3d' \
     --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py turb --mpirun=srun --mpirun_opts=--job-name='GCC turb/' \
     --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py mpi --mpirun=srun --mpirun_opts=--job-name='GCC mpi/mpi_linwave' \
     --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py omp --coverage="${lcov_capture_cmd}" --silent
timeout --signal=TERM 60m time python -u ./run_tests.py hybrid --mpirun=srun \
	--mpirun_opts=--job-name='GCC hybrid/hybrid_linwave' \
	--coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py hydro --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py amr --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py outputs --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py sr --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py curvilinear --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py symmetry --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py eos --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py scalars/mignone_radial_1d --coverage="${lcov_capture_cmd}" --silent
# Exclude gr/compile*.py regression tests from code coverage analysis (nothing is executed in these tests):
time python -u ./run_tests.py gr/compile_kerr-schild gr/compile_minkowski gr/compile_schwarzschild --silent
time python -u ./run_tests.py gr/mhd_shocks_hlld gr/mhd_shocks_hlle gr/mhd_shocks_llf \
     --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py gr/hydro_shocks_hllc gr/hydro_shocks_hlle gr/hydro_shocks_llf \
     --coverage="${lcov_capture_cmd}" --silent
time python -u ./run_tests.py gr/hydro_shocks_hlle_no_transform gr/hydro_shocks_llf_no_transform \
     --coverage="${lcov_capture_cmd}" --silent
# For regression tests with unacceptably long runtimes with -O0 optimization, "sample" the code coverage by running each test twice:
# - 1x normally (-O3) without --coverage=CMD to check correctness
# - 1x with --coverage=CMD (and hence -O0) and small cycle limit, ignoring failure in subsequent test.analyze() step
time python -u ./run_tests.py mhd --coverage="${lcov_capture_cmd}" -r="time/nlim=10" --silent || true
time python -u ./run_tests.py mhd --silent  # (mhd/mhd_linwave.py is currenlty the slowest regression test):

time python -u ./run_tests.py shearingbox --coverage="${lcov_capture_cmd}" -r="time/nlim=10" --silent || true
time python -u ./run_tests.py shearingbox --silent

time python -u ./run_tests.py diffusion --coverage="${lcov_capture_cmd}" -r="time/nlim=10" --silent || true
time python -u ./run_tests.py diffusion --silent

# High-order solver regression tests w/ GCC
time python -u ./run_tests.py hydro4 --coverage="${lcov_capture_cmd}" -r="time/nlim=10" --silent || true
time python -u ./run_tests.py hydro4 --silent

# Swap serial HDF5 library module for parallel HDF5 library:
module unload hdf5/gcc/1.10.0
module load hdf5/gcc/openmpi-3.0.0/1.10.0
mpi_hdf5_library_path='/usr/local/hdf5/gcc/openmpi-3.0.0/1.10.0/lib64'
module list
# This P-HDF5 library, built with OpenMPI 1.10.2, is incompatible with OpenMPI 3.0.0 on Perseus:
#module load hdf5/gcc/openmpi-1.10.2/1.10.0

# Workaround issue with parallel HDF5 modules compiled with OpenMPI on Perseus--- linker still chooses serial HDF5 library in /usr/lib64/
# due to presence of -L flag in mpicxx wrapper that overrides LIBRARY_PATH environment variable
time python -u ./run_tests.py pgen/hdf5_reader_parallel --coverage="${lcov_capture_cmd}" \
     --mpirun=srun --mpirun_opts=--job-name='GCC pgen/hdf5_reader_parallel' \
     --config=--lib_path=${mpi_hdf5_library_path} --silent

# Combine Lcov tracefiles from individaul regression tests:
# All .info files in current working directory tst/regression/ -> lcov.info
# (remove '-maxdepth 1' to recursively search subfolders for more .info)
lcov_counter=0
set +e  # Don't quit on errors during Lcov processing / don't let the build fail here
while read filename; do
    # Accumulate string variable containing all tracefiles joined by '-a '
    lcov_input_files="$lcov_input_files -a \"$filename\""
    # Alternative to uploading single unified "lcov.info" tracefile: attempt to upload each Lcov
    # test_name.info tracefile separately with a Codecov Flag matching test_name (or test_set/ group?)
    codecov_flag=$(basename ${filename} .info) # "flags must match pattern ^[\w\,]+$"
    # basename command is in GNU coreutils, but here is Bash Parameter Expansion alternative for stripping extension and path:
    #codecov_flag=${${filename%.info}##*/}
    curl -s https://codecov.io/bash | bash -s - -X gcov -t ccdc959e-e2c3-4811-95c6-512151b39471 \
	-F ${codecov_flag} -f "${filename}" || echo "Codecov did not collect coverage reports"
    lcov_counter=$((lcov_counter + 1))
done < <( find . -maxdepth 1 -name '*.info' )
eval "${lcov_cmd}" "${lcov_input_files}" -o lcov.info
# Explicitly return count of individual Lcov tracefiles, and monitor any changes to this number (53 expected as of 2018-12-04):
# (most Lcov failures will be silent and hidden in build log;, missing reports will be hard to notice in Lcov HTML and Codecov reports)
echo "Detected ${lcov_counter} individual tracefiles and combined them -> lcov.info"

# Generate Lcov HTML report and backup to home directory on Perseus (never used by Codecov):
gendesc scripts/tests/test_descriptions.txt --output-filename ./regression_tests.desc
lcov_dir_name="${SLURM_JOB_NAME}_lcov_html"
# TODO(felker): Address "lcov: ERROR: no valid records found in tracefile ./eos_eos_comparison_eos_hllc.info"
genhtml --legend --show-details --keep-descriptions --description-file=regression_tests.desc \
	--branch-coverage -o ${lcov_dir_name} lcov.info
mv lcov.info ${lcov_dir_name}
# GNU (but not BSD) tar supports --remove-files option for cleaning up files (and directories) after adding them to the archive:
tar --remove-files -cvzf "${lcov_dir_name}.tar.gz" ${lcov_dir_name}
mv "${lcov_dir_name}.tar.gz" $HOME  # ~2 MB. Manually rm HTML databases from $HOME on a reg. basis
# genhtml requires that src/ is unmoved since compilation; works from $HOME on Perseus,
# but lcov.info tracefile is not portable across sytems (without --to-package, etc.)
#cp lcov.info $HOME  # ~30 MB --- tracefile is too large to store long-term

# Ensure that no stale tracefiles are kept in Jenkins cached workspace
rm -rf *.info
set -e

# Build step #2: regression tests using Intel compiler and MPI library
module purge
# Delete version info from module names to automatically use latest default version of these libraries as Princeton Research Computing updates them:
# (Was using pinned Intel 17.0 Release 5 versions as of November 2018 due to bugs on Perseus installation of ICC 19.0.
# Intel's MPI Library 2019 version was never installed on Perseus since it is much slower than 2018 version on Mellanox Infiniband)
module load intel/18.0/64/18.0.3.222 # intel/17.0/64/17.0.5.239 # intel ---intel/19.0/64/19.0.3.199 latest version as of 2019-05-04
module load intel-mpi/intel/2017.5/64 # intel-mpi --- intel-mpi/intel/2018.3/64
# Always pinning these modules to a specific version, since new library versions are rarely compiled:
module load fftw/gcc/3.3.4
module load hdf5/intel-17.0/1.10.0 # hdf5/intel-17.0/intel-mpi/1.10.0
# Note, do not mix w/ "module load rh" to ensure that Intel shared libraries are used by the loader (especially OpenMP?)

module list

# Use of --config=--cflag=-gxx-name=/opt/rh/devtoolset-7/root/usr/bin/g++ (For GCC >6 Enumerator Attributes) possibly causing:
# src/eos/general/general_hydro.cpp(163): (col. 10) remark: function was not vectorized: condition too complex
# ": internal error: ** The compiler has encountered an unexpected problem. ** Segmentation violation signal raised. **
# Access violation or stack overflow. Please contact Intel Support for assistance.

time python -u ./run_tests.py pgen/pgen_compile --config=--cxx=icpc --config=--cflag="$(../ci/set_warning_cflag.sh icpc)"
time python -u ./run_tests.py pgen/hdf5_reader_serial --silent
time python -u ./run_tests.py grav --config=--cxx=icpc --mpirun=srun --mpirun_opts=--job-name='ICC grav/jeans_3d' --silent
time python -u ./run_tests.py turb --config=--cxx=icpc --mpirun=srun --mpirun_opts=--job-name='ICC turb/' --silent
time python -u ./run_tests.py mpi --config=--cxx=icpc --mpirun=srun --mpirun_opts=--job-name='ICC mpi/mpi_linwave' --silent
time python -u ./run_tests.py omp --config=--cxx=icpc --silent
timeout --signal=TERM 60m time python -u ./run_tests.py hybrid --config=--cxx=icpc \
	--mpirun=srun --mpirun_opts=--job-name='ICC hybrid/hybrid_linwave' --silent
time python -u ./run_tests.py hydro --config=--cxx=icpc --silent
time python -u ./run_tests.py mhd --config=--cxx=icpc --silent
time python -u ./run_tests.py amr --config=--cxx=icpc --silent
time python -u ./run_tests.py outputs --config=--cxx=icpc --silent
time python -u ./run_tests.py sr --config=--cxx=icpc --silent
time python -u ./run_tests.py gr --config=--cxx=icpc --silent
time python -u ./run_tests.py curvilinear --config=--cxx=icpc --silent
time python -u ./run_tests.py shearingbox --config=--cxx=icpc --silent
time python -u ./run_tests.py diffusion --config=--cxx=icpc --silent
time python -u ./run_tests.py symmetry --config=--cxx=icpc --silent
time python -u ./run_tests.py eos --config=--cxx=icpc --silent
time python -u ./run_tests.py scalars/mignone_radial_1d --config=--cxx=icpc --silent

# High-order solver regression tests w/ Intel compiler
time python -u ./run_tests.py hydro4 --config=--cxx=icpc --silent

# Swap serial HDF5 library module for parallel HDF5 library:
module unload hdf5/intel-17.0/1.10.0
module load hdf5/intel-17.0/intel-mpi/1.10.0
mpi_hdf5_library_path='/usr/local/hdf5/intel-17.0/intel-mpi/1.10.0/lib64'
module list
# Workaround issue with parallel HDF5 modules compiled with OpenMPI on Perseus--- linker still takes serial HDF5 library in /usr/lib64/
# due to presence of -L flag in mpicxx wrapper that overrides LIBRARY_PATH environment variable
time python -u ./run_tests.py pgen/hdf5_reader_parallel --config=--cxx=icpc \
     --mpirun=srun --mpirun_opts=--job-name='ICC pgen/hdf5_reader_parallel' \
     --config=--lib_path=${mpi_hdf5_library_path} --silent

# Test OpenMP 4.5 SIMD-enabled function correctness by disabling IPO and forced inlining w/ Intel compiler flags
# Check subset of regression test sets to try most EOS functions (which heavily depend on vectorization) that are called in rsolvers
time python -u ./run_tests.py pgen/pgen_compile --config=--cxx=icpc-debug --config=--cflag="$(../ci/set_warning_cflag.sh icpc)"
time python -u ./run_tests.py hydro --config=--cxx=icpc-debug --silent
time python -u ./run_tests.py mhd --config=--cxx=icpc-debug --silent
time python -u ./run_tests.py sr --config=--cxx=icpc-debug --silent
time python -u ./run_tests.py gr --config=--cxx=icpc-debug --silent

set +e
# end regression tests

# Alternative: Upload single combined tracefile for Codecov analysis of test coverage reports:
# Use curl-pipe to Codecov Bash Uploader (recommended approach for Jenkins).
# If using the default options (no -f PATTERN), any Lcov tracefile must be named "lcov.info".
#curl -s https://codecov.io/bash | bash -s - -X gcov -t ccdc959e-e2c3-4811-95c6-512151b39471 || echo "Codecov did not collect coverage reports"
# (Default will always exit with 0. Use -Z to exit with 1 if not successful.)
# "exit 0" in Codecov Bash uploader script is not fool-proof. Preventing build failures with catch-all echo statement to ensure exit status=0

# Slurm diagnostics: see all timing info when build script finishes
# (should run in Jenkins "Execute shell" build step when Slurm allocation is released)
sacct --jobs=$SLURM_JOB_ID --format=JobID,JobName%50,Submit,Start,Elapsed,Timelimit,End  #--noheader
