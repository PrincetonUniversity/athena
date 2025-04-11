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
import numpy as np
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('mpi', 'fft',
                     prob='turb',
                     **kwargs)
    athena.make()
    os.system('mv bin/athena bin/athena_mpi_fft')
    os.system('mv obj obj_mpi_fft')

    athena.configure('fft',
                     prob='turb',
                     **kwargs)
    athena.make()
    os.system('mv bin/athena bin/athena_fft')
    os.system('mv obj obj_fft')


# Run Athena++
def run(**kwargs):
    os.system('mv obj_fft obj')
    os.system('mv bin/athena_fft bin/athena')
    arguments = ['time/ncycle_out=10',
                 'mesh/nx1=64', 'mesh/nx2=32', 'mesh/nx3=32',
                 'meshblock/nx1=16',
                 'meshblock/nx2=16',
                 'meshblock/nx3=16',
                 'turbulence/turb_flag=3',
                 'turbulence/rseed=1',
                 'output2/dt=-1', 'time/tlim=0.3']
    athena.run('hydro/athinput.turb', arguments + ['job/problem_id=turb_serial'],
               lcov_test_suffix='fft')

    os.system('mv obj_mpi_fft obj')
    os.system('mv bin/athena_mpi_fft bin/athena')
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  1, 'hydro/athinput.turb', arguments + ['job/problem_id=turb_mpi1'])
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  2, 'hydro/athinput.turb', arguments + ['job/problem_id=turb_mpi2'])
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'],
                  4, 'hydro/athinput.turb', arguments + ['job/problem_id=turb_mpi4'],
                  lcov_test_suffix='mpi_fft')

    return 'skip_lcov'


# Analyze outputs
def analyze():
    analyze_status = True
    # read data from error file
    filenames = ['bin/turb_serial.hst',
                 'bin/turb_mpi1.hst',
                 'bin/turb_mpi2.hst',
                 'bin/turb_mpi4.hst']

    hst = athena_read.hst(filenames[0])
    KE0 = hst['1-KE']+hst['2-KE']+hst['3-KE']

    KE_final = [KE0[-1]]
    diff = []
    for hstfile in filenames[1:]:
        hst = athena_read.hst(hstfile)
        KE = hst['1-KE']+hst['2-KE']+hst['3-KE']
        KE_final.append(KE[-1])
        diff.append(np.sum(np.abs(KE-KE0)))

    logger.info("%g %g %g %g", KE0[-1], KE_final[1], KE_final[2], KE_final[3])

    # check errors between runs w/wo MPI and different numbers of cores
    if diff[0] > 1.e-7:
        logger.warning("Turb runs with serial and 1 MPI rank are different %g", diff[0])
        analyze_status = False
    if diff[1] > 1.e-7:
        logger.warning("Turb runs with serial and 2 MPI ranks are different %g", diff[1])
        analyze_status = False
    if diff[2] > 1.e-7:
        logger.warning("Turb runs with serial and 4 MPI ranks are different %g", diff[2])
        analyze_status = False

    return analyze_status
