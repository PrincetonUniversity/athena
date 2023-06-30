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
    athena.configure('implicit_radiation',
                     prob='rad_linearwave',
                     coord='cartesian',
                     flux='hllc', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    # L-going fast wave (set by default in input)
    arguments = ['time/cfl_number=0.3',  # default =0.4, but tolerances measured w/ 0.3
                 'time/ncycle_out=100'
                 ]
    athena.run('radiation/athinput.rad_linearwave_amr', arguments)


# Analyze outputs
def analyze():
    # read data from error file
    filename = 'bin/linearwave-errors.dat'
    data = []
    with open(filename, 'r') as f:
        raw_data = f.readlines()
        for line in raw_data:
            if line.split()[0][0] == '#':
                continue
            data.append([float(val) for val in line.split()])

    if data[0][4] > 1.05e-8:
        print("error in regime 8: ", data[0][4])
        return False

    return True
