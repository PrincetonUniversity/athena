# Regression test based on the Newtonian MHD linear wave convergence problem, 1D

# Confirms fourth-order convergence rates of the RK3 + PPM + Laplacian flux correction
# + UCT4 (Felker & Stone 2018) solver configuration for each MHD linear wave mode on a
# 1D uniform Cartesian grid. The exact edge-averaged vector potential initial condition
# and cell-averaged analytic errors are used (correct_ic, correct_err).

# Modules
import logging
import scripts.utils.athena as athena
from math import log
import numpy as np
import sys
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

# List of time/integrator and time/xorder combinations to test:
solvers = [('rk3', '4')]

# Matching above list of solver configurations, provide bounds on error metrics:
# for each tested resolution (excluding lowest Nx1=16) and wave_flag.
# Tolerances are ~2x the errors measured on 2026-07 (macOS/AppleClang, -O3):
#   fast:    1.60e-9, 1.99e-10, 2.49e-11
#   Alfven:  3.06e-10, 2.26e-11, 2.10e-12
#   slow:    3.04e-10, 1.96e-11, 3.91e-12
#   entropy: 2.93e-10, 1.90e-11, 1.29e-12
# Upper bound on RMS-L1 errors:
error_tols = [
    # RK3 + MHD4
    ((3.2e-9, 4.0e-10, 5.0e-11), (6.2e-10, 4.6e-11, 4.3e-12),
     (6.1e-10, 4.0e-11, 8.0e-12), (5.9e-10, 3.9e-11, 2.7e-12)),
]
# for each wave_flag, lower bound on convergence rate at Nx1=128. The fast wave rate is
# limited to ~3.0 by the O(dt^3) RK3 truncation error; the slow wave error approaches
# its floor (~4e-12) by Nx1=128
rate_tols = [
    # RK3 + MHD4
    (2.9, 3.2, 2.1, 3.5),
]

resolution_range = [16, 32, 64, 128]
num_nx1 = len(resolution_range)
# Number of times Athena++ is run for each above configuration:
# L-going wave for each mode is run num_nx1 times, then L and R going fast waves are run
# at a single resolution for a symmetry test
nrows_per_solver = 4*num_nx1 + 2
wave_mode_names = ['fast', 'Alfven', 'slow', 'entropy']


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     nghost=6,  # required for fourth-order MHD configurations
                     prob='linear_wave',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    for (torder, xorder) in solvers:
        # L-going fast/Alfven/slow waves
        for w in (0, 1, 2):
            tlim = max(0.5, w)
            for i in resolution_range:
                arguments = ['time/ncycle_out=0',
                             'time/xorder=' + xorder, 'time/integrator=' + torder,
                             'problem/wave_flag={}'.format(w), 'problem/vflow=0.0',
                             'mesh/nx1={}'.format(i),
                             'output2/dt=-1', 'time/tlim={}'.format(tlim),
                             'time/correct_err=true',
                             'time/correct_ic=true',
                             'problem/compute_error=true']
                athena.run('mhd/athinput.linear_wave1d', arguments)
        # L-going entropy wave
        for i in resolution_range:
            arguments = ['time/ncycle_out=0',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         'problem/wave_flag=3', 'problem/vflow=1.0',
                         'mesh/nx1={}'.format(i),
                         'output2/dt=-1', 'time/tlim=1.0',
                         'time/correct_err=true',
                         'time/correct_ic=true',
                         'problem/compute_error=true']
            athena.run('mhd/athinput.linear_wave1d', arguments)
        # L/R-going fast wave for symmetry comparison
        for w in (0, 6):
            arguments = ['time/ncycle_out=0',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         'problem/wave_flag={}'.format(w),
                         'output2/dt=-1', 'time/tlim=0.5',
                         'time/correct_err=true',
                         'time/correct_ic=true',
                         'problem/compute_error=true']
            athena.run('mhd/athinput.linear_wave1d', arguments)


# formatting rules for numpy arrays of floats:
def error_formatter(x):
    return "{:.2e}".format(x)


def rate_formatter(x):
    return "{:.2f}".format(x)


# Analyze outputs
def analyze():
    analyze_status = True
    # read data from error file
    filename = 'bin/linearwave-errors.dat'
    data = athena_read.error_dat(filename)

    for ((torder, xorder), err_tol, rate_tol) in zip(solvers, error_tols, rate_tols):
        # effectively list.pop() range of rows for this solver configuration
        solver_results = np.array(data[0:nrows_per_solver])
        data = np.delete(data, np.s_[0:nrows_per_solver], 0)

        # Compute error convergence rates with Richardson extrapolation per wave flag
        logger.info('------------------------------')
        logger.info('{} + time/xorder={}'.format(torder.upper(), xorder))
        logger.info('------------------------------')
        logger.info('Solver wave mode error tolerances at each resolution:')
        logger.info('nx1=' + repr(resolution_range[1:]))
        with np.printoptions(formatter={'float_kind': error_formatter}):
            logger.info(err_tol)
        logger.info('Wave mode convergence rate tolerances (at nx1=128)')
        with np.printoptions(formatter={'float_kind': rate_formatter}):
            logger.info(rate_tol)
        for w, wave_mode in enumerate(wave_mode_names):
            # L-going wave of each mode
            logger.info("{} wave error convergence:".format(wave_mode.capitalize()))
            logger.info("nx1   |   rate   |   RMS-L1")
            rms_errs = solver_results[0:num_nx1, 4]
            nx1_range = solver_results[0:num_nx1, 0]
            solver_results = np.delete(solver_results, np.s_[0:num_nx1], 0)
            for i in range(1, num_nx1):
                rate = log(rms_errs[i-1]/rms_errs[i])/log(nx1_range[i]/nx1_range[i-1])
                logger.info("%d %g %g", int(nx1_range[i]), rate, rms_errs[i])
                if (nx1_range[i] == 128 and rate < rate_tol[w]):
                    logger.warning(
                        "L-going {} wave converging at rate {} slower than {}".format(
                            wave_mode, rate, rate_tol[w]))
                    analyze_status = False
                if (rms_errs[i] > err_tol[w][i-1]):
                    logger.warning(
                        "L-going {} wave error {} is larger than tolerance {}".format(
                            wave_mode, rms_errs[i], err_tol[w][i-1]))
                    analyze_status = False

        # Check that errors are identical for fast waves in each direction. The RMS-L1
        # errors here are ~2.5e-11 (near the 1D error floor), so allow an absolute
        # difference of a few ULP-accumulations (measured L/R difference: ~8e-16)
        if not np.allclose(solver_results[-2, 4], solver_results[-1, 4],
                           atol=5e-15, rtol=1e-5):
            msg = "L/R-going fast wave errors, {} and {}"
            msg += ", have a difference that is not close to round-off"
            logger.warning(msg.format(solver_results[-2, 4], solver_results[-1, 4]))
            analyze_status = False

    return analyze_status
