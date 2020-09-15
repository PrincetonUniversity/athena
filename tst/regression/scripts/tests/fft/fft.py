# Regression test based on intercomparison between serial and MPI runs.
#
# Runs a driven turbulence test in 3D and checks L1 errors between
# serial and MPI runs using kinetic energy history.
#
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
                     prob='fft',
                     **kwargs)
    athena.make()
    os.system('mv bin/athena bin/athena_mpi_fft')
    os.system('mv obj obj_mpi_fft')

    athena.configure('fft',
                     prob='fft',
                     **kwargs)
    athena.make()
    os.system('mv bin/athena bin/athena_fft')
    os.system('mv obj obj_fft')


# Run Athena++
def run(**kwargs):
    os.system('mv obj_fft obj')
    os.system('mv bin/athena_fft bin/athena')
    athena.run('hydro/athinput.fft', ['job/problem_id=fft_serial'],
               lcov_test_suffix='fft')

    os.system('mv obj_mpi_fft obj')
    os.system('mv bin/athena_mpi_fft bin/athena')
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  1, 'hydro/athinput.fft', ['job/problem_id=fft_mpi1'])
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  2, 'hydro/athinput.fft', ['job/problem_id=fft_mpi2'])
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  4, 'hydro/athinput.fft', ['job/problem_id=fft_mpi4'],
                  lcov_test_suffix='mpi_fft')

    return 'skip_lcov'


# Analyze outputs
def analyze():
    analyze_status = True
    # read data from error file
    filename = 'bin/fft-errors.dat'
    data = athena_read.error_dat(filename)

    logger.info("%g %g %g", data[0][4], data[0][5], data[0][6])

    # check errors between runs w/wo MPI and different numbers of cores
    if data[0][5] != data[1][5]:
        logger.warning("FFT runs with serial and 1 MPI rank are different %g %g",
                       data[0][5], data[1][5])
        analyze_status = False
    if (data[0][5] > 1.e-10) or (data[0][6] > 1.e-10):
        logger.warning("FFT error is too large for a serial run real=%g, imag=%g",
                       data[0][5], data[0][6])
        analyze_status = False
    if (data[2][5] > 1.e-10) or (data[3][5] > 1.e-10):
        logger.warning("FFT error is too large for MPI runs np2=%g np4=%g",
                       data[2][5], data[3][5])
        analyze_status = False
    if (data[0][4] < 1.e6):
        logger.warning("FFT is too slow for a serial run zcs=%g",
                       data[0][4])
    if (data[3][4] < 1.e6):
        logger.warning("FFT is too slow for a MPI run (np=4) zcs=%g",
                       data[3][4])

    return analyze_status
