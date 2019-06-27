# Regression test based on MHD circularly polarized Alfven wave convergence problem
#
# Runs an isothermal cpaw convergence test in 2D including SMR and checks L1 errors (which
# are computed by the executable automatically and stored in the temporary file
# linearwave_errors.dat)

# Modules
import logging
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b', prob='cpaw', eos='isothermal', flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    # run R-going wave at two resolutions
    for i in (128, 256):
        arguments = ['time/ncycle_out=0',
                     'mesh/refinement=static',
                     'mesh/nx1=' + repr(i), 'mesh/nx2=' + repr(i/2),
                     'meshblock/nx1=' + repr(i/4), 'meshblock/nx2=' + repr(i/8),
                     'output2/dt=-1', 'time/tlim=1.0', 'problem/compute_error=true']
        athena.run('mhd/athinput.cpaw2d', arguments)
    # run L-going wave
    arguments = [
        'time/ncycle_out=0',
        'mesh/refinement=static',
        'mesh/nx1=256',
        'mesh/nx2=128',
        'meshblock/nx1=64',
        'meshblock/nx2=32',
        'output2/dt=-1',
        'time/tlim=1.0',
        'problem/compute_error=true',
        'problem/dir=2']
    athena.run('mhd/athinput.cpaw2d', arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    # read data from error file
    filename = 'bin/cpaw-errors.dat'
    data = athena_read.error_dat(filename)

    logger.info("%g", data[0][4])
    logger.info("%g", data[1][4])
    logger.info("%g", data[2][4])

    # check absolute error and convergence
    if data[1][4] > 2.0e-4:
        logger.warning("error in L-going fast wave too large %g", data[1][4])
        analyze_status = False
    if data[1][4]/data[0][4] > 0.3:
        logger.warning("not converging for L-going fast wave %g %g",
                       data[0][4], data[1][4])
        analyze_status = False

    # check error identical for waves in each direction
    if abs(data[2][4] - data[1][4]) > 2.0e-6:
        logger.warning("error in L/R-going Alfven waves not equal %g %g",
                       data[2][4], data[0][4])
        analyze_status = False

    return analyze_status
