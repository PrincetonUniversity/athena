# Regression test based on Newtonian 2D MHD linear wave test problem with AMR
#
# Runs a 2D linear wave test with AMR, using a refinement condition that tracks the
# velocity maxima.  Then checks L1 and L_infty (max) error.  This test is very sensitive
# to finding errors in AMR prolongation/restriction/boundaries

# Modules
import logging
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     prob='linear_wave',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    # L-going fast wave (set by default in input)
    arguments = ['time/ncycle_out=10',
                 'time/cfl_number=0.3',  # default =0.4, but tolerances measured w/ 0.3
                 'output1/dt=-1',
                 'output2/file_type=vtk',
                 'output2/dt=-1',
                 ]
    athena.run('mhd/athinput.linear_wave2d_amr', arguments)


# Analyze outputs
def analyze():
    # read data from error file
    filename = 'bin/linearwave-errors.dat'
    data = athena_read.error_dat(filename)

    analyze_status = True
    if data[0][4] > 2.0e-8:
        logger.warning("RMS error in L-going fast wave too large %g", data[0][4])
        analyze_status = False
    if data[0][13] > 5.5:
        logger.warning("maximum relative error in L-going fast wave too large %g",
                       data[0][13])
        analyze_status = False

    return analyze_status
