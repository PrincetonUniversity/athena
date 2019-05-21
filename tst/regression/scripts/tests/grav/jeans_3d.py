# Regression test based on Newtonian hydro linear wave convergence problem
#
# Runs a linear wave convergence test in 3D including SMR and checks L1 errors (which
# are computed by the executable automatically and stored in the temporary file
# linearwave_errors.dat)

# Modules
import logging
import os
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('mpi', 'fft',
                     prob='jeans',
                     grav='fft', **kwargs)
    athena.make()
    os.system('mv bin/athena bin/athena_mpi_fft')
    os.system('mv obj obj_mpi_fft')

    athena.configure('mpi',
                     prob='jeans',
                     grav='mg', **kwargs)
    athena.make()
    os.system('mv bin/athena bin/athena_mpi_mg')
    os.system('mv obj obj_mpi_mg')

    athena.configure('fft',
                     prob='jeans',
                     grav='fft',
                     **kwargs)
    athena.make()
    os.system('mv bin/athena bin/athena_fft')
    os.system('mv obj obj_fft')

    athena.configure(prob='jeans',
                     grav='mg',
                     **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    arguments = ['time/ncycle_out=10',
                 'mesh/nx1=64', 'mesh/nx2=32', 'mesh/nx3=32',
                 'meshblock/nx1=16',
                 'meshblock/nx2=16',
                 'meshblock/nx3=16',
                 'output2/dt=-1', 'time/tlim=1.0', 'problem/compute_error=true']
    athena.run('hydro/athinput.jeans_3d', arguments, lcov_test_suffix='mg')

    os.system('rm -rf obj')
    os.system('mv obj_fft obj')
    os.system('mv bin/athena_fft bin/athena')
    athena.run('hydro/athinput.jeans_3d', arguments, lcov_test_suffix='fft')

    os.system('rm -rf obj')
    os.system('mv obj_mpi_mg obj')
    os.system('mv bin/athena_mpi_mg bin/athena')
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  1, 'hydro/athinput.jeans_3d', arguments)
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  2, 'hydro/athinput.jeans_3d', arguments)
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  4, 'hydro/athinput.jeans_3d', arguments,
                  lcov_test_suffix='mpi_mg')

    os.system('rm -rf obj')
    os.system('mv obj_mpi_fft obj')
    os.system('mv bin/athena_mpi_fft bin/athena')
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  1, 'hydro/athinput.jeans_3d', arguments)
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  2, 'hydro/athinput.jeans_3d', arguments)
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  4, 'hydro/athinput.jeans_3d', arguments,
                  lcov_test_suffix='mpi_fft')
    return 'skip_lcov'


# Analyze outputs
def analyze():
    analyze_status = True
    # read data from error file
    filename = 'bin/jeans-errors.dat'
    data = athena_read.error_dat(filename)

    logger.info("%g %g", data[0][4], data[1][4])
    logger.info("%g %g %g", data[2][4], data[3][4], data[4][4])
    logger.info("%g %g %g", data[5][4], data[6][4], data[7][4])

    # check errors between runs w/wo MPI and different numbers of cores
    if data[0][4] != data[2][4]:
        logger.warning(
            "Linear wave error with one core w/wo MPI not identical for MG gravity %g %g",
            data[0][4],
            data[2][4])
        analyze_status = False
    if data[0][4] > 1.e-7:
        logger.warning("Linear wave error is too large for MG gravity %g", data[0][4])
        analyze_status = False
    if data[1][4] > 1.e-7:
        logger.warning("Linear wave error is too large for FFT gravity %g", data[1][4])
        analyze_status = False
    if data[1][4] != data[5][4]:
        logger.warning(
            "Linear wave error with one core w/wo MPI not identical for FFT gravity"
            + " %g %g", data[1][4], data[5][4])
        analyze_status = False
    if abs(data[3][4]-data[2][4]) > 5.0e-4:
        logger.warning(
            "Linear wave error between 2 and 1 cores too large for MG gravity %g %g",
            data[3][4], data[2][4])
        analyze_status = False
    if abs(data[4][4]-data[2][4]) > 5.0e-4:
        logger.warning(
            "Linear wave error between 4 and 1 cores too large for MG gravity %g %g",
            data[4][4], data[2][4])
        analyze_status = False
    if abs(data[6][4]-data[5][4]) > 5.0e-4:
        logger.warning(
            "Linear wave error between 2 and 1 cores too large for FFT gravity %g %g",
            data[6][4], data[5][4])
        analyze_status = False
    if abs(data[7][4]-data[5][4]) > 5.0e-4:
        logger.warning(
            "Linear wave error between 4 and 1 cores too large for FFT gravity %g %g",
            data[7][4], data[5][4])
        analyze_status = False

    return analyze_status
