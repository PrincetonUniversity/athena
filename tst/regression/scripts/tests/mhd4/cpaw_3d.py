# Regression test based on the 3D circularly polarized Alfven wave (CPAW) problem

# Confirms fourth-order convergence of the RK3 + PPM + Laplacian flux correction + UCT4
# (Felker & Stone 2018) solver configuration for a finite-amplitude, exactly nonlinear
# smooth MHD solution propagating along the diagonal of a 3D uniform square-cell
# Cartesian grid. This is the strongest validator of the 3D corner-EMF (E1, E2)
# reconstruction and upwinding machinery. Multiple MeshBlocks per run exercise the
# fourth-order ghost-zone exchanges.

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

# Tolerances are ~2x the errors measured on 2026-07 (macOS/AppleClang, -O3, single
# MeshBlock): 1.39e-2 (16), 9.28e-4 (32), 6.72e-5 (64); measured rates 3.91, 3.79
resolution_range = [16, 32, 64]
error_tols = (2.8e-2, 1.9e-3, 1.4e-4)
rate_tol = 3.4  # at nx1=64


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     nghost=6,  # required for fourth-order MHD configurations
                     prob='cpaw',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    for i in resolution_range:
        arguments = ['time/ncycle_out=0',
                     'time/xorder=4', 'time/integrator=rk3',
                     'time/correct_ic=true', 'time/correct_err=true',
                     'mesh/nx1={}'.format(i), 'mesh/nx2={}'.format(i//2),
                     'mesh/nx3={}'.format(i//2),
                     'meshblock/nx1={}'.format(i//2),
                     'meshblock/nx2={}'.format(i//2),
                     'meshblock/nx3={}'.format(i//2),
                     'output1/dt=-1', 'output2/dt=-1', 'time/tlim=1.0',
                     'problem/compute_error=1']
        athena.run('mhd/athinput.cpaw3d', arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    # read data from error file
    filename = 'bin/cpaw-errors.dat'
    data = athena_read.error_dat(filename)
    data = np.array(data)

    logger.info('RK3 + time/xorder=4 (UCT4) 3D CPAW convergence:')
    logger.info("nx1   |   rate   |   RMS-L1")
    rms_errs = data[:, 4]
    nx1_range = data[:, 0]
    for i in range(len(resolution_range)):
        rate = None
        if i > 0:
            rate = log(rms_errs[i-1]/rms_errs[i])/log(nx1_range[i]/nx1_range[i-1])
        logger.info("%d %s %g", int(nx1_range[i]),
                    "{:.2f}".format(rate) if rate is not None else "----", rms_errs[i])
        if rms_errs[i] > error_tols[i]:
            logger.warning("CPAW 3D error {} at nx1={} exceeds tolerance {}".format(
                rms_errs[i], int(nx1_range[i]), error_tols[i]))
            analyze_status = False
        if i == len(resolution_range) - 1 and rate is not None and rate < rate_tol:
            logger.warning("CPAW 3D convergence rate {} at nx1={} is below {}".format(
                rate, int(nx1_range[i]), rate_tol))
            analyze_status = False

    return analyze_status
