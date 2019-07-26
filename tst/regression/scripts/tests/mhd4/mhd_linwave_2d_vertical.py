# Regression test based on Newtonian MHD linear wave convergence problem

# Check errors of reconstruction and time integrator options other than default VL2+PLM
# primitive reconstruction. In particular, confirm fourth-order convergence rate for
# semi-discrete integration with RK4 + PPM + Laplacian flux correction terms + UCT4.
# 2D uniform square grid, no SMR--- midpoint assumption used in init. and in error
# calculations (like the hydro4/ tests)

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
solvers = [
    # ('vl2', '2'),
    ('rk3', '4')]  # , ('rk4', '4c'), ('ssprk5_4', '4')]
# [('vl2', '2c'), ('vl2', '3')]  # , ('rk2', '3c')]


# Matching above list of solver configurations, provide bounds on error metrics:
# for each tested resolution (excluding lowest Nx1=16) and wave_flag
# Upper bound on RMS-L1 errors:
error_tols = [
    # VL2 + PLMc + UCT2
    ((1.6e-7, 5.2e-8, 1.3e-8, 2.9e-9), (1.4e-7, 5.1e-8, 1.3e-8, 2.9e-9),
     (2.0e-7, 6.8e-8, 1.8e-8, 4.3e-9), (1.1e-7, 3.9e-8, 9.9e-9, 2.4e-9)),
    # VL2 + PPM + UCT2
    ((1.3e-7, 3.0e-8, 7.3e-9, 1.9e-9), (7.6e-8, 1.9e-8, 4.5e-9, 1.2e-9),
     (7.5e-8, 1.8e-8, 4.4e-9, 1.1e-9), (3.4e-8, 7.6e-9, 1.9e-9, 4.7e-10)),
    # --------
    # RK2 + PPMc + UCT2 --- unstable? branch mhd4_3d @ f7e0f6a has ~25% higher Alfven wave
    # errors compared to master @ 9c53b6b
    # ((3.9e-8, 1.1e-8, 2.9e-9, 7.2e-10), (5.0e-9, 1.6e-9, 4.5e-10, 2.8e-9),
    # all wave modes have strangely high convergence at nx1=32
    # alfven wave stops converging (error increases by 5x between nx1=128,256)
    # slow wave stops converging (error inceases gradually after nx=64
    # entropy wave error starts anomalously low, increases at nx1=64, then converges at
    # sub-second-order
    # (1.75e-8, 3.1e-9, 3.9e-9, 4.6e-9), (3.0e-10, 7.5e-10, 2.4e-10, 6.1e-11)),
    # --------

    # RK3 + MHD4
    # RK4 + MHD4c
    # SSPRK(5,4) + MHD4
]
# for each wave_flag, lower bound on convergence rate at Nx1=128 asymptotic convergence
# regime. Linear hydro waves stop converging around RMS-L1 error 1e-11 to 1e-12
rate_tols = [
    # VL2 + PLMc + UCT2
    (2.0, 2.0, 1.9, 1.9),
    # VL2 + PPM + UCT2
    (2.0, 2.0, 2.0, 2.0),
    # --------
    # RK2 + PPMc + UCT2--- see above comments
    # (1.9, 1.8, -1.0, 1.0),
    # --------
    # RK3 + MHD4
    # RK4 + MHD4c
    # SSPRK(5,4) + MHD4
]
# this metric is redundant with above error_tols, but it is simpler...

resolution_range = [16, 32, 64, 128]  # , 256]  # , 512]
num_nx2 = len(resolution_range)
# Number of times Athena++ is run for each above configuration:
# L-going wave for each mode is run num_nx2 times, then L and R going fast waves are run
# at single resolution for symmetry test
nrows_per_solver = 4*num_nx2 + 2
wave_mode_names = ['fast', 'Alfven', 'slow', 'entropy']


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     nghost=6,  # required for fourth-order configurations
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
                             'problem/ang_3_vert=true',
                             'mesh/x1max=0.5', 'mesh/x2max=1.0',
                             'problem/wave_flag={}'.format(w), 'problem/vflow=0.0',
                             'mesh/nx1={}'.format(i/2), 'mesh/nx2={}'.format(i),
                             'output2/dt=-1', 'time/tlim={}'.format(tlim),
                             'time/correct_err=true',
                             'time/correct_ic=true',
                             'problem/compute_error=true']
                athena.run('mhd/athinput.linear_wave2d_aligned', arguments)
        # L-going entropy wave
        for i in resolution_range:
            arguments = ['time/ncycle_out=0',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         'problem/ang_3_vert=true',
                         'mesh/x1max=0.5', 'mesh/x2max=1.0',
                         'problem/wave_flag=3', 'problem/vflow=1.0',
                         'mesh/nx1={}'.format(i/2), 'mesh/nx2={}'.format(i),
                         'output2/dt=-1', 'time/tlim=1.0',
                         'time/correct_err=true',
                         'time/correct_ic=true',
                         'problem/compute_error=true']
            athena.run('mhd/athinput.linear_wave2d_aligned', arguments)
        # L/R-going fast wave for symmetry comparison
        for w in (0, 6):
            arguments = ['time/ncycle_out=0',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         'problem/ang_3_vert=true',
                         'mesh/x1max=0.5', 'mesh/x2max=1.0',
                         'problem/wave_flag={}'.format(w),
                         'mesh/nx1={}'.format(i/2), 'mesh/nx2={}'.format(i),
                         'output2/dt=-1', 'time/tlim=0.5',
                         'time/correct_err=true',
                         'time/correct_ic=true',
                         'problem/compute_error=true']
            athena.run('mhd/athinput.linear_wave2d_aligned', arguments)


# formatting rules for numpy arrays of floats:
error_formatter = lambda x: "{:.2e}".format(x)
rate_formatter = lambda x: "{:.2f}".format(x)


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

        # Compute error convergence rates with Richardson extrapolation for each wave flag
        # --------------------
        logger.info('------------------------------')
        logger.info('{} + time/xorder={}'.format(torder.upper(), xorder))
        logger.info('------------------------------')
        logger.info('Solver wave mode error tolerances at each resolution:')
        logger.info('nx2=' + repr(resolution_range[1:]))
        with np.printoptions(formatter={'float_kind': error_formatter}):
            logger.info(err_tol)
        logger.info('Wave mode convergence rate tolerances (at nx2=128)')
        with np.printoptions(formatter={'float_kind': rate_formatter}):
            logger.info(rate_tol)
        for w, wave_mode in enumerate(wave_mode_names):
            # L-going wave of each mode
            logger.info("{} wave error convergence:".format(wave_mode.capitalize()))
            logger.info("nx2   |   rate   |   RMS-L1")
            rms_errs = solver_results[0:num_nx2, 4]
            nx2_range = solver_results[0:num_nx2, 1]
            # print(solver_results)
            # solver_results = np.delete(solver_results, np.s_[0:len(wave_mode_names)], 0)
            solver_results = np.delete(solver_results, np.s_[0:num_nx2], 0)
            for i in range(1, num_nx2):
                rate = log(rms_errs[i-1]/rms_errs[i])/log(nx2_range[i]/nx2_range[i-1])
                logger.info("%d %g %g", int(nx2_range[i]), rate, rms_errs[i])
                if (nx2_range[i] == 128 and rate < rate_tol[w]):
                    logger.warning(
                        "L-going {} wave converging at rate {} slower than {}".format(
                            wave_mode, rate, rate_tol[w]))
                    analyze_status = False
                if (rms_errs[i] > err_tol[w][i-1]):
                    logger.warning(
                        "L-going {} wave error {} is larger than tolerance {}".format(
                            wave_mode, rms_errs[i], err_tol[w][i-1]))
                    analyze_status = False

        # Check that errors are identical for fast waves in each direction at Nx2=256
        if (not np.allclose(solver_results[-2, 4], solver_results[-1, 4],
                            atol=5e-16, rtol=1e-5) and xorder != '3c'):
            msg = "L/R-going fast wave errors, {} and {}"
            msg += ", have a difference that is not close to round-off"
            logger.warning(msg.format(solver_results[-2, 4], solver_results[-1, 4]))
            analyze_status = False

    return analyze_status
