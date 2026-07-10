# Regression test based on the Newtonian MHD RJ2a shock tube problem, 1D, xorder=4

# Runs the Ryu & Jones "2a" shock tube in the x1 direction with the RK3 + PPM +
# Laplacian flux correction + UCT4 (Felker & Stone 2018) solver configuration, and
# checks L1 errors against the reference solution (computed by the executable and
# stored in the temporary file shock-errors.dat). This exercises the limiting and
# shock-capturing robustness of the fourth-order MHD configuration (JCP paper sec. 5.5).
# NOTE: outflow boundary conditions are permitted for 4th-order MHD in 1D only.

# Modules
import logging
import scripts.utils.athena as athena
import sys
import numpy as np
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

_nxs = [256, 512]  # resolutions to test (ascending)
# Upper bounds on L1 errors at each resolution, and lower bound on the L1 convergence
# rate between them (shocks limit convergence to ~1st order). Characteristic
# reconstruction (xorder=4c) is used, as in the JCP paper's shock tube results.
# Tolerances are ~1.5x the errors measured on 2026-07: 1.073e-2 (256), 5.959e-3 (512)
_error_tols = [0.016, 0.009]
_rate_tol = 0.7


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     nghost=6,  # required for fourth-order MHD configurations
                     prob='shock_tube',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    for nx in _nxs:
        arguments = ['time/ncycle_out=0',
                     'time/xorder=4c', 'time/integrator=rk3',
                     'time/cfl_number=0.3',
                     'mesh/nx1={}'.format(nx),
                     'problem/shock_dir=1',
                     'problem/compute_error=true']
        athena.run('mhd/athinput.rj2a', arguments)


# Analyze outputs
def analyze():
    analyze_status = True
    # read data from error file
    filename = 'bin/shock-errors.dat'
    data = athena_read.error_dat(filename)
    data = np.array(data)

    logger.info('RK3 + time/xorder=4c (UCT4) 1D RJ2a shock tube errors:')
    logger.info("nx1   |   rate   |   L1")
    errs = data[:, 4]
    nxs = data[:, 0]
    for i in range(len(_nxs)):
        rate = None
        if i > 0:
            rate = np.log(errs[i-1]/errs[i])/np.log(nxs[i]/nxs[i-1])
        logger.info("%d %s %g", int(nxs[i]),
                    "{:.2f}".format(rate) if rate is not None else "----", errs[i])
        if errs[i] > _error_tols[i]:
            logger.warning("RJ2a L1 error {} at nx1={} exceeds tolerance {}".format(
                errs[i], int(nxs[i]), _error_tols[i]))
            analyze_status = False
        if i > 0 and rate < _rate_tol:
            logger.warning("RJ2a convergence rate {} is below {}".format(
                rate, _rate_tol))
            analyze_status = False

    return analyze_status
