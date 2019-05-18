# Regression test based on Newtonian MHD linear wave convergence problem with MPI+OpenMP
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


# Prepare Athena++ w/wo MPI+OpenMP
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b', 'mpi', 'omp', prob='linear_wave', coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()
    os.system('mv bin/athena bin/athena_hybrid')
    os.system('mv obj obj_hybrid')

    athena.configure('b', prob='linear_wave', coord='cartesian', flux='hlld', **kwargs)
    athena.make()


# Run Athena++: serial, then hybrid MPI+OpenMP
def run(**kwargs):
    # L-going fast wave
    arguments = ['time/ncycle_out=0',
                 'problem/wave_flag=0', 'problem/vflow=0.0', 'mesh/refinement=static',
                 'mesh/nx1=32', 'mesh/nx2=16', 'mesh/nx3=16',
                 'meshblock/nx1=8',
                 'meshblock/nx2=8',
                 'meshblock/nx3=8',
                 'output2/dt=-1', 'time/tlim=2.0', 'problem/compute_error=true']
    athena.run('mhd/athinput.linear_wave3d', arguments, lcov_test_suffix='serial')

    os.system('rm -rf obj')
    os.system('mv obj_hybrid obj')
    os.system('mv bin/athena_hybrid bin/athena')
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'], 1,
                  'mhd/athinput.linear_wave3d', arguments + ['mesh/num_threads=1'])
    # 4 total threads = 2 MPI ranks x 2 OpenMP threads / rank
    athena.mpirun(kwargs['mpirun_cmd'], kwargs['mpirun_opts'], 2,
                  'mhd/athinput.linear_wave3d', arguments + ['mesh/num_threads=2'],
                  lcov_test_suffix='hybrid')
    return 'skip_lcov'


# Analyze outputs
def analyze():
    analyze_status = True
    # read data from error file
    filename = 'bin/linearwave-errors.dat'
    data = athena_read.error_dat(filename)

    logger.info("%g %g %g", data[0][4], data[1][4], data[2][4])

    # check errors between runs: serial vs. hybrid w/ 1 thread vs. hybrid w/ 4 threads
    if data[0][4] != data[1][4]:
        logger.warning(
            "Linear wave error from serial calculation vs."
            "single thread not identical %g %g", data[0][4], data[1][4])
        analyze_status = False
    if abs(data[2][4] - data[0][4]) > 5.0e-4:
        logger.warning("Linear wave error differences between 4 threads vs."
                       "serial is too large %g %g", data[2][4], data[0][4])
        analyze_status = False

    return analyze_status
