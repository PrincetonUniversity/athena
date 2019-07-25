# Regression test based on Newtonian hydro linear wave convergence problem

# Check errors of reconstruction and time integrator options other than default VL2+PLM
# primitive reconstruction. In particular, confirm fourth-order convergence rate for
# semidiscrete integration with RK4 + PPM + Laplacian flux correction terms. 2D uniform
# square grid, no SMR--- midpoint assumption used in init. and in error calculations

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
solvers = [('vl2', '2c'), ('vl2', '3'), ('rk2', '3c'),
           ('rk3', '4'), ('rk4', '4c'), ('ssprk5_4', '4')]
# Matching above list of solver configurations, provide bounds on error metrics:
# for each tested resolution (excluding lowest Nx1=16) and wave_flag
# Upper bound on RMS-L1 errors:
error_tols = [((1.4e-7, 4.6e-8, 1.1e-8, 2.5e-9), (1.1e-7, 3.7e-8, 9.3e-9, 2.2e-9)),
              ((9.6e-8, 2.4e-8, 5.8e-9, 1.5e-9), (4.5e-8, 1.1e-8, 2.6e-9, 6.4e-10)),
              ((3.7e-8, 1.1e-8, 2.7e-9, 6.7e-10), (4.8e-9, 2.0e-9, 5.3e-10, 1.4e-10)),
              ((5.5e-9, 4.0e-10, 3.6e-11, 6.2e-12), (3.7e-9, 2.5e-10, 1.6e-11, 1.1e-12)),
              ((5.2e-9, 3.4e-10, 2.2e-11, 5.6e-12), (3.8e-9, 2.4e-10, 1.6e-11, 1.7e-12)),
              ((5.2e-9, 3.4e-10, 2.1e-11, 5.6e-12), (3.8e-9, 2.4e-10, 1.6e-11, 1.1e-12))
              ]
# for each wave_flag, lower bound on convergence rate at Nx1=128 asymptotic convergence
# regime. Linear hydro waves stop converging around RMS-L1 error 1e-11 to 1e-12
rate_tols = [(2.0, 1.9), (2.0, 2.0), (1.95, 1.85),
             (3.4, 3.95), (3.95, 3.95),  (3.95, 3.95)]
# this metric is redundant with above error_tols, but it is simpler...

resolution_range = [16, 32, 64, 128, 256]  # , 512]
num_nx1 = len(resolution_range)
# Number of times Athena++ is run for each above configuration:
nrows_per_solver = 2*num_nx1 + 2


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(
        nghost=4,  # required for fourth-order configurations
        prob='linear_wave',
        coord='cartesian',
        flux='hllc', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    for (torder, xorder) in solvers:
        # L-going sound wave
        for i in resolution_range:
            arguments = ['time/ncycle_out=0',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         'problem/wave_flag=0', 'problem/vflow=0.0',
                         'mesh/nx1=' + repr(i), 'mesh/nx2=' + repr(i/2),
                         'output2/dt=-1', 'time/tlim=1.0',
                         'problem/compute_error=true']
            athena.run('hydro/athinput.linear_wave2d', arguments)
        # L-going entropy wave
        for i in resolution_range:
            arguments = ['time/ncycle_out=0',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         'problem/wave_flag=3', 'problem/vflow=1.0',
                         'mesh/nx1=' + repr(i), 'mesh/nx2=' + repr(i/2),
                         'output2/dt=-1', 'time/tlim=1.0',
                         'problem/compute_error=true']
            athena.run('hydro/athinput.linear_wave2d', arguments)
        # L/R-going sound wave, no SMR
        for w in (0, 4):
            arguments = ['time/ncycle_out=0',
                         'time/xorder=' + xorder, 'time/integrator=' + torder,
                         'problem/wave_flag=' + repr(w),
                         'output2/dt=-1', 'time/tlim=1.0',
                         'problem/compute_error=true']
            athena.run('hydro/athinput.linear_wave2d', arguments)


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
        logger.info('Solver sound/entropy wave error tolerances at each resolution:')
        logger.info('nx1=' + repr(resolution_range[1:]))
        with np.printoptions(formatter={'float_kind': error_formatter}):
            logger.info(err_tol)
        logger.info('Sound/entropy wave convergence rate tolerances (at nx1=128)')
        with np.printoptions(formatter={'float_kind': rate_formatter}):
            logger.info(rate_tol)
        # L-going sound wave
        logger.info("Sound wave error convergence:")
        logger.info("nx1   |   rate   |   RMS-L1")
        rms_errs = solver_results[0:num_nx1, 4]
        nx1_range = solver_results[0:num_nx1, 0]
        for i in range(1, num_nx1):
            rate = log(rms_errs[i-1]/rms_errs[i])/log(nx1_range[i]/nx1_range[i-1])
            logger.info("%d %g %g", int(nx1_range[i]), rate, rms_errs[i])
            # old rate calculation from hydro/hydro_linwave.py:
            # logger.info(rms_errs[i]/rms_errs[i-1])
            if (nx1_range[i] == 128 and rate < rate_tol[0]):
                logger.warning(
                    "L-going sound wave converging at rate {} slower than {}".format(
                        rate, rate_tol[0]))
                analyze_status = False
            if (rms_errs[i] > err_tol[0][i-1]):
                logger.warning(
                    "L-going sound wave error {} is larger than tolerance {}".format(
                        rms_errs[i], err_tol[0][i-1]))
                analyze_status = False

        # L-going entropy wave
        logger.warning("Entropy wave error convergence:")
        logger.warning("nx1   |   rate   |   RMS-L1")
        rms_errs = solver_results[num_nx1:2*num_nx1, 4]
        nx1_range = solver_results[num_nx1:2*num_nx1, 0]
        for i in range(1, num_nx1):
            rate = log(rms_errs[i-1]/rms_errs[i])/log(nx1_range[i]/nx1_range[i-1])
            logger.info("%d %g %g", int(nx1_range[i]), rate, rms_errs[i])
            # old rate calculation from hydro/hydro_linwave.py:
            # logger.warning(rms_errs[i]/rms_errs[i-1])
            if (nx1_range[i] == 128 and rate < rate_tol[1]):
                logger.warning(
                    "L-going entropy wave converging at rate {} slower than {}".format(
                        rate, rate_tol[1]))
                analyze_status = False
            if (rms_errs[i] > err_tol[1][i-1]):
                logger.warning(
                    "L-going entropy wave error {} is larger than tolerance {}".format(
                        rms_errs[i], err_tol[1][i-1]))
                analyze_status = False

        # logger.info(solver_results[-2, 4] - solver_results[-1, 4])
        # logger.info(solver_results[-2, 4], solver_results[-1, 4])
        # logger.info("numpy bound = {}".format(5e-16 + 1e-5*abs(solver_results[-1, 4])))

        # Check that errors are identical for sound waves in each direction at Nx1=256
        if (not np.allclose(solver_results[-2, 4], solver_results[-1, 4],
                            atol=5e-16, rtol=1e-5) and xorder != '3c'):
            msg = "L/R-going sound wave errors, {} and {}"
            msg += ", have a difference that is not close to round-off"
            logger.warning(msg.format(solver_results[-2, 4], solver_results[-1, 4]))
            analyze_status = False
        # Unlike comparison for VL2+PLM, the errors become small enough in the high-order
        # solver that round-off needs to be accounted for. Assuming true solution A is
        # O(1), then the difference betwee L/R-going RMS-L1 errors approach:
        # Diff = sum(A - sol1) - sum(A - sol2) ---> 1 + epsilon

        # Currently, RK2 + 3c produces the largest directional asymmetry w/
        # diff=-2.22e-14. Characteristic projection procedure may be introducing asymmetry
    return analyze_status
