# Regression test based on Mignone (2014) section 5.1.1 tests: radial advection of passive
# scalar profiles in cylindrical and spherical-polar curvilinear coordinate systems. Tests
# the importance of curvilinear corrections to PLM and PPM4 limiters and PPM4 stencil
# weights for interpolation of face-averages from cell-avergages.


# Modules
import logging
import os
import scripts.utils.athena as athena
import numpy as np
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa

athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module
# Switches to automatically use utilities in plot_mignone/ subfolder:
plot_profiles = False
nx1_profile = 64
plot_convergence = False
if plot_convergence or plot_profiles:
    from scripts.utils.plot_mignone.section_5_1_1 import \
        figure2_profiles, figure3_convergence # noqa

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
            # 2: [3.36e-4, 1.02e-4, 2.07e-5, 4.61e-6, 1.08e-6, 2.63e-7, 6.48e-8],
            # PLM van Leer (VL) limiter (unpublished--- from Athena++ 79bc2b8751):
            2: [2.18e-4, 5.93e-5, 1.54e-5, 3.97e-6, 1.01e-6, 2.54e-7, 6.36e-8],
            3: [1.41e-4, 1.42e-5, 1.13e-6, 7.18e-8, 4.50e-9, 2.83e-10, 1.77e-11],
        },
        'B': {
            # PLM MC
            # 2: [1.48e-2, 5.93e-3, 2.15e-3, 6.35e-4, 1.82e-4, 4.92e-5, 1.28e-5],
            # PLM VL
            2: [1.65e-2, 7.48e-3, 2.63e-3, 7.92e-4, 2.25e-4, 5.92e-5, 1.49e-5],
            3: [1.23e-2, 3.78e-3, 9.72e-4, 2.09e-4, 4.54e-5, 9.82e-6, 2.11e-6],
        },
    },
    'spherical_polar': {
        'A': {
            # PLM MC
            # 2: [3.04e-5, 1.16e-5, 2.34e-6, 5.13e-7, 1.20e-7, 2.90e-8, 7.13e-9],
            # PLM VL
            2: [1.57e-5, 5.56e-6, 1.61e-6, 4.26e-7, 1.09e-7, 2.77e-08, 6.97e-09],
            3: [1.24e-5, 1.23e-6, 8.65e-8, 5.37e-9, 3.32e-10, 2.05e-11, 1.27e-12],
        },
        'B': {
            # PLM MC
            # 2: [5.58e-3, 2.27e-3, 8.62e-4, 2.48e-4, 6.99e-5, 1.88e-5, 4.86e-6],
            # PLM VL
            2: [6.17e-3, 2.83e-3, 1.04e-3, 3.05e-4, 8.59e-5, 2.25e-5, 5.67e-6],
            3: [4.68e-3, 1.41e-3, 3.36e-4, 7.56e-5, 1.64e-5, 3.57e-6, 7.75e-7],
        },
    }
}

resolution_range = [32, 64, 128, 256, 512, 1024, 2048]
base_x1rat = 1.0
num_nx1 = len(resolution_range)
# Number of times Athena++ is run for each above configuration:
nrows_per_solver = 1*num_nx1 + 2

try:
    # Filter ANSI escape character sequences from text sent to stdout and replace with
    # Win32 calls for adding color to text
    import colorama
    colorama.init()  # does nothing if not on Windows
except ImportError:
    pass
try:
    # convenience package for colorizing and stylizing text using terminal ANSI escape
    # codes
    from termcolor import colored
    percentage_formatter = lambda x: colored(
        "{:.2%}".format(x), None if np.abs(x) < diff_threshold
        else 'green' if x >= 0 else 'red')
except ImportError:
    # otherwise, do not use color annotations at all
    percentage_formatter = lambda x: "{:.2%}".format(x)

error_formatter = lambda x: "{:.2e}".format(x)


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    for coord_ in coords:
        athena.configure(
            nghost=3,  # required for PPM
            prob='mignone_advection',
            eos='isothermal',
            # flux='upwind',
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
                for i, nx_ in enumerate(resolution_range):
                    if i > 0:
                        # take sqrt of geometric mesh spacing factor every time the radial
                        # resolution is doubled:
                        x1rat = np.power(base_x1rat, 0.5**i)
                    else:
                        x1rat = base_x1rat
                    dt_tab = 1.0 if plot_profiles and nx_ == nx1_profile else -1
                    arguments = ['time/ncycle_out=0',
                                 'time/xorder={}'.format(xorder),
                                 'time/integrator={}'.format(torder),
                                 'mesh/nx1={}'.format(nx_),
                                 'mesh/x1rat={}'.format(x1rat),
                                 # suppress .hst output, by default
                                 'output1/dt=-1',
                                 # .tab output:
                                 'output2/dt={}'.format(dt_tab),
                                 'problem/iprob={}'.format(iprob),
                                 'problem/a_width={}'.format(case_params[0]),
                                 'problem/b_center={}'.format(case_params[1]),
                                 'problem/compute_error=true']
                    athena.run('hydro/athinput.mignone_radial', arguments,
                               lcov_test_suffix=coord_)
                    if plot_profiles and nx_ == nx1_profile:
                        tab_file = os.path.join('bin',
                                                'MignoneRadial.block0.out2.00001.tab')
                        new_tab_file = os.path.join(
                            'bin',
                            '{}_case_{}_{}_xorder_{}_nx1_{}.tab'.format(
                                coord_, case_, torder, xorder, nx_))
                        os.system('mv {} {}'.format(tab_file, new_tab_file))
                default_error_file = os.path.join('bin', 'mignone_radial-errors.dat')
                error_file = os.path.join(
                    'bin', 'errors_{}_case_{}_{}_xorder_{}.dat'.format(coord_, case_,
                                                                       torder, xorder))
                os.system('mv {} {}'.format(default_error_file, error_file))
    return 'skip_lcov'


# Analyze outputs
def analyze():
    analyze_status = True
    for coord_ in coords:
        logger.info("\n\nTesting --coord={}\n".format(coord_))
        for (case_, case_params) in cases.items():
            logger.info("\nMignone 1D radial test problem, case {}".format(case_))
            for (torder, xorder) in solvers:
                logger.info('{} + xorder={}'.format(torder.upper(), xorder))
                error_file = os.path.join(
                    'bin', 'errors_{}_case_{}_{}_xorder_{}.dat'.format(coord_, case_,
                                                                       torder, xorder))
                # read Athena++ data from error file
                data = athena_read.error_dat(error_file)
                # compare with above hard-coded values from Mignone (2014) (or Athena++
                # for VL limiter for PLM)
                mignone_errs = mignone_tbl1[coord_][case_][xorder]
                percent_diff = (mignone_errs - data[:, 4])/mignone_errs
                logger.info("Reference L1 errors: Mignone (2014) Table 1, N=[32, 2048]")
                with np.printoptions(formatter={'float_kind': error_formatter}):
                    logger.info(mignone_errs)
                    logger.info("Observed Athena++ L1 errors, Nx1=[32, 2048]")
                    logger.info(data[:, 4])
                logger.info("Percentage differences relative to reference values")
                with np.printoptions(formatter={'float_kind': percentage_formatter}):
                    logger.info(percent_diff)
                # currently allowing Athena++ L1 errors to be more than threshold % LOWER
                # than Mignone errors
                if np.min(percent_diff) < -diff_threshold:
                    # if np.max(np.abs(percent_diff)) > diff_threshold:
                    logger.warning(
                        "Observed Athena++ L1 error is higher than Mignone Table 1 value"
                        " by > {:2%}".format(diff_threshold))
                    analyze_status = False

    # (optional) produce plots of results:
    if plot_profiles:
        # section_5_1_1.figure2_profiles()
        figure2_profiles()
    if plot_convergence:
        figure3_convergence()
    return analyze_status
