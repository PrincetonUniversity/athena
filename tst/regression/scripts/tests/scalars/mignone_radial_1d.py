# Regression test based on Mignone (2014) section 5.1.1 tests: radial advection of passive
# scalar profiles in cylindrical and spherical-polar curvilinear coordinate systems. Tests
# the importance of curvilinear corrections to PLM and PPM4 limiters and PPM4 stencil
# weights for interpolation of face-averages from cell-avergages.


# Modules
import logging
import os
import scripts.utils.athena as athena
from termcolor import colored
import numpy as np
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

# List of time/integrator and time/xorder combinations to test:
solvers = [('rk3', 2), ('rk3', 3), ]
coords = ['cylindrical', 'spherical_polar']
diff_threshold = 0.005  # tolerate up to 0.5% difference of L1 errors vs. Table 1 values
iprob = 1  # 1 = radial 1D tests, section 5.1.1 from Mignone (2014)
cases = {
    # [a, b] parameters for each case A, B
    'A': [10.0, 0.0],
    'B': [16.0, 0.5],
}

mignone_tbl1 = {
    'cylindrical': {
        'A': {
            # PLM monotonized central (MC) limiter:
            2: [3.36e-4, 1.02e-4, 2.07e-5, 4.61e-6, 1.08e-6, 2.63e-7, 6.48e-8],
            # PLM van Leer (VL) limiter (unpublished--- from Athena++):
            # 2: [],
            3: [1.41e-4, 1.42e-5, 1.13e-6, 7.18e-8, 4.50e-9, 2.83e-10, 1.77e-11],
        },
        'B': {
            # PLM MC
            2: [1.48e-2, 5.93e-3, 2.15e-3, 6.35e-4, 1.82e-4, 4.92e-5, 1.28e-5],
            # PLM VL
            # 2: [],
            3: [1.23e-2, 3.78e-3, 9.72e-4, 2.09e-4, 4.54e-5, 9.82e-6, 2.11e-6],
        },
    },
    'spherical_polar': {
        'A': {
            # PLM MC
            2: [3.04e-5, 1.16e-5, 2.34e-6, 5.13e-7, 1.20e-7, 2.90e-8, 7.13e-9],
            # PLM VL
            # 2: [],
            3: [1.24e-5, 1.23e-6, 8.65e-8, 5.37e-9, 3.32e-10, 2.05e-11, 1.27e-12],
        },
        'B': {
            # PLM MC
            2: [5.58e-3, 2.27e-3, 8.62e-4, 2.48e-4, 6.99e-5, 1.88e-5, 4.86e-6],
            # PLM VL
            # 2: [],
            3: [4.68e-3, 1.41e-3, 3.36e-4, 7.56e-5, 1.64e-5, 3.57e-6, 7.75e-7],
        },
    }
}

resolution_range = [32, 64, 128, 256, 512, 1024, 2048]
num_nx1 = len(resolution_range)
# Number of times Athena++ is run for each above configuration:
nrows_per_solver = 1*num_nx1 + 2

# TODO(felker): conditionally import "termcolors" (if available) and use of:
percentage_formatter = lambda x: colored(
    "{:.2%}".format(x), None if np.abs(x) < diff_threshold
    else 'green' if x >= 0 else 'red')


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    for coord_ in coords:
        athena.configure(
            nghost=3,  # required for PPM
            prob='mignone_advection',
            eos='isothermal',
            nscalars=1,
            coord=coord_,
            **kwargs)
        athena.make()
        os.system('mv bin/athena bin/athena_{}'.format(coord_))
        os.system('mv obj obj_{}'.format(coord_))


# Run Athena++
def run(**kwargs):
    for coord_ in coords:
        os.system('mv obj_{} obj'.format(coord_))
        os.system('mv bin/athena_{} bin/athena'.format(coord_))
        for (case_, case_params) in cases.items():
            for (torder, xorder) in solvers:
                for i in resolution_range:
                    arguments = ['time/ncycle_out=0',
                                 'time/xorder={}'.format(xorder),
                                 'time/integrator={}'.format(torder),
                                 'mesh/nx1={}'.format(i),
                                 # suppress .hst and .tab output, by default
                                 'output1/dt=-1',
                                 'output2/dt=-1',
                                 'problem/iprob={}'.format(iprob),
                                 'problem/a_width={}'.format(case_params[0]),
                                 'problem/b_center={}'.format(case_params[1]),
                                 'problem/compute_error=true']
                    athena.run('hydro/athinput.mignone_radial', arguments,
                               lcov_test_suffix=coord_)
                default_error_file = os.path.join('bin', 'mignone_radial-errors.dat')
                error_file = os.path.join('bin', 'errors_{}_case_{}_{}_xorder_{}'.format(
                    coord_, case_, torder, xorder))
                os.system('mv {} {}'.format(default_error_file, error_file))
    return 'skip_lcov'


# Analyze outputs
def analyze():
    analyze_status = True
    for coord_ in coords:
        logger.info("Testing --coord={}".format(coord_))
        for (case_, case_params) in cases.items():
            logger.info("Mignone 1D radial test problem, case {}".format(case_))
            for (torder, xorder) in solvers:
                logger.info('{} + {}'.format(torder.upper(), xorder))
                error_file = os.path.join('bin', 'errors_{}_case_{}_{}_xorder_{}'.format(
                    coord_, case_, torder, xorder))
                # read Athena++ data from error file
                data = athena_read.error_dat(error_file)
                # compare with above hard-coded values from Mignone (2014)
                mignone_errs = mignone_tbl1[coord_][case_][xorder]
                percent_diff = (mignone_errs - data[:, 4])/mignone_errs
                logger.info("Reference L1 errors: Mignone (2014) Table 1, N=[32, 2048]")
                logger.info(mignone_errs)
                logger.info("Observed Athena++ L1 errors, Nx1=[32, 2048]")
                logger.info(data[:, 4])
                logger.info("Percentage differences relative to reference values")
                with np.printoptions(formatter={'float_kind': percentage_formatter}):
                    logger.info(percent_diff)
                # perhaps allow Athena++ L1 errors to be much LOWER than Mignone errors?
                if np.max(np.abs(percent_diff)) > diff_threshold:
                    logger.warning(
                        "Observed Athena++ L1 error disagrees with Mignone Table 1 value"
                        " by > {:2%}".format(diff_threshold))
                    analyze_status = False
    return analyze_status
