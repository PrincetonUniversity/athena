# Regression test based on Mignone (2014) section 5.1.2 tests: meridional advection of
# passive scalar profiles in spherical-polar curvilinear coordinate system. Tests the
# importance of curvilinear corrections to PLM and PPM4 limiters and PPM4 stencil weights
# for interpolation of face-averages from cell-avergages.


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
nx2_profile = 64
plot_convergence = False
if plot_convergence or plot_profiles:
    from scripts.utils.plot_mignone.section_5_1_2 import \
        plot_profiles, figure4_convergence # noqa

# List of time/integrator and time/xorder combinations to test:
solvers = [('rk3', 2), ('rk3', 3), ]
diff_threshold = 0.005  # tolerate up to 0.5% difference of L1 errors vs. Table 1 values
iprob = 2  # 1 = meridional 1D tests, section 5.1.2 from Mignone (2014)
cases = {
    # [a, b] parameters for each case A, B
    'A': [10.0, 0.0],
    'B': [16.0, 0.196349540849362077],  # = pi/16.0
}

mignone_tbl4 = {
    'A': {
        # PLM monotonized central (MC) limiter:
        # 2: [3.95e-4, 1.31e-4, 2.95e-5, 6.92e-6, 1.67e-6, 4.12e-7, 1.02e-7],
        # PLM van Leer (VL) limiter (unpublished--- from Athena++ ???):
        2: [2.17e-4, 7.84e-5, 2.26e-5, 6.07e-6, 1.57e-6, 3.99e-7, 1.01e-7],
        3: [8.76e-5, 9.97e-6, 8.44e-7, 6.41e-8, 4.66e-9,
            # Unknown error convergence floor in 2D Athena++ variant of 1D Mignone problem
            # that affects monotonic case's PPM4 results at nx2=1024, 2048:
            3.45e-10, 8.17e-11],
        # If analytic x2-face area-average v_theta(theta, r) is used to upwind scalars:
        # (nearly the same for all dx1, unlike general technique's O(dx1^2) flux error)
        # 3.45e-10, 8.72e-11],

        # Actual Mignone (2014) Table 4 errors:
        # 3.33e-10, 2.42e-11],
    },
    'B': {
        # PLM MC
        # 2: [3.97e-3, 1.30e-3, 3.90e-4, 9.74e-5, 2.43e-5, 6.05e-6, 1.50e-6],
        # PLM VL
        2: [4.73e-3, 1.47e-3, 4.35e-4, 1.15e-4, 2.88e-5, 7.14e-6, 1.75e-6],
        3: [2.57e-3, 5.87e-4, 1.30e-4, 2.81e-5, 6.10e-6, 1.33e-6,
            # See above comments:
            2.84e-07]
        # 2.82e-7],
    },
}

resolution_range = [32, 64, 128, 256, 512, 1024, 2048]
base_x2rat = 1.0
num_nx2 = len(resolution_range)
# Number of times Athena++ is run for each above configuration:
nrows_per_solver = 1*num_nx2 + 2

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
    athena.configure(
        nghost=3,  # required for PPM
        prob='mignone_advection',
        eos='isothermal',
        # flux='upwind',
        nscalars=1,
        coord='spherical_polar',
        **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    for (case_, case_params) in cases.items():
        for (torder, xorder) in solvers:
            for i, nx_ in enumerate(resolution_range):
                if i > 0:
                    # take sqrt of geometric mesh spacing factor every time the x2
                    # resolution is doubled:
                    x2rat = np.power(base_x2rat, 0.5**i)
                else:
                    x2rat = base_x2rat
                dt_tab = 1.0 if plot_profiles and nx_ == nx2_profile else -1
                arguments = ['time/ncycle_out=10',
                             'time/xorder={}'.format(xorder),
                             'time/integrator={}'.format(torder),
                             'mesh/nx2={}'.format(nx_),
                             'mesh/x2rat={}'.format(x2rat),
                             # suppress .hst output, by default
                             'output1/dt=-1',
                             # .tab output:
                             'output2/dt={}'.format(dt_tab),
                             'problem/iprob={}'.format(iprob),
                             'problem/a_width={}'.format(case_params[0]),
                             'problem/b_center={}'.format(case_params[1]),
                             'problem/compute_error=true']
                athena.run('hydro/athinput.mignone_meridional', arguments)
                if plot_profiles and nx_ == nx2_profile:
                    tab_file = os.path.join('bin',
                                            'MignoneMeridional.block0.out2.00001.tab')
                    new_tab_file = os.path.join(
                        'bin',
                        'case_{}_{}_xorder_{}_nx2_{}.tab'.format(case_, torder,
                                                                 xorder, nx_))
                    os.system('mv {} {}'.format(tab_file, new_tab_file))
            default_error_file = os.path.join('bin', 'mignone_meridional-errors.dat')
            error_file = os.path.join(
                'bin', 'errors_case_{}_{}_xorder_{}.dat'.format(case_, torder, xorder))
            os.system('mv {} {}'.format(default_error_file, error_file))
    return


# Analyze outputs
def analyze():
    analyze_status = True
    for (case_, case_params) in cases.items():
        logger.info("\nMignone 1D meridional test problem, case {}".format(case_))
        for (torder, xorder) in solvers:
            logger.info('{} + xorder={}'.format(torder.upper(), xorder))
            error_file = os.path.join(
                'bin', 'errors_case_{}_{}_xorder_{}.dat'.format(case_, torder, xorder))
            # read Athena++ data from error file
            data = athena_read.error_dat(error_file)
            # compare with above hard-coded values from Mignone (2014) (or Athena++
            # for VL limiter for PLM)
            mignone_errs = mignone_tbl4[case_][xorder]
            percent_diff = (mignone_errs - data[:, 4])/mignone_errs
            logger.info("Reference L1 errors: Mignone (2014) Table 4, N=[32, 2048]")
            with np.printoptions(formatter={'float_kind': error_formatter}):
                logger.info(mignone_errs)
                logger.info("Observed Athena++ L1 errors, Nx2=[32, 2048]")
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
        plot_profiles()
    if plot_convergence:
        figure4_convergence()
    return analyze_status
