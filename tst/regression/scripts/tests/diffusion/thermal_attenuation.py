# Regression test based on the decaying linear wave due to thermal
# conduction. The decay rate is fit and then compared with analytic
# solution.

# Modules
import logging
import numpy as np
from numpy.polynomial import Polynomial
import sys
import scripts.utils.athena as athena
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

_kappa = 0.04

resolution_range = [32, 64]
method = 'Explicit'
# Upper bound on relative L1 error for each above nx1:
error_rel_tols = [0.38, 0.10]


def prepare(*args, **kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(*args,
                     prob='linear_wave',
                     flux='hllc',
                     eos='adiabatic', **kwargs)
    athena.make()


def run(**kwargs):
    for i in resolution_range:
        arguments = ['output1/dt=0.03',
                     'output2/dt=-1',  # disable .vtk outputs
                     'time/tlim=3.0',
                     'time/ncycle_out=0',
                     # L-going sound wave
                     'problem/wave_flag=0',
                     'problem/amp=1.0e-4',
                     'problem/vflow=0.0',
                     'problem/kappa_iso=0.04',
                     'mesh/nx1=' + repr(i),
                     'mesh/nx2=' + repr(i/2),
                     'mesh/nx3=' + repr(i/2),
                     'meshblock/nx1=' + repr(i),
                     'meshblock/nx2=' + repr(i/2),
                     'meshblock/nx3=' + repr(i/2),
                     'job/problem_id=DecayLinWave-{}'.format(i)]
        athena.run('hydro/athinput.linear_wave3d', arguments)


def analyze():
    # Lambda=1 for Athena++'s linear wave setups in 1D, 2D, and 3D:
    L = 1.0
    ksqr = (2.0*np.pi/L)**2
    # The decay rate for a sound wave with a thermal conduction term is given by
    # decay rate = ((\gamma-1)^2*kappa/gamma/2)*k^2
    decay_rate = 2.0*_kappa/15.0*ksqr
    analyze_status = True
    errors_abs = []

    for (nx, err_tol) in zip(resolution_range, error_rel_tols):
        logger.info('[Decaying 3D Linear Wave {}]: '
                    'Mesh size {} x {} x {}'.format(method, nx, nx/2, nx/2))
        filename = 'bin/DecayLinWave-{}.hst'.format(nx)
        hst_data = athena_read.hst(filename)
        tt = hst_data['time']
        max_vy = hst_data['max-v2']

        # estimate the decay rate from simulation, using weighted least-squares (WLS)
        yy = np.log(np.abs(max_vy))
        p, [resid, rank, sv, rcond] = Polynomial.fit(tt, yy, 1, w=np.sqrt(max_vy),
                                                     full=True)
        resid_normal = np.sum((yy - p(tt)) ** 2)
        r2 = 1 - resid_normal/(yy.size*yy.var())
        pnormal = p.convert(domain=(-1, 1))
        fit_rate = -pnormal.coef[-1]

        error_abs = np.fabs(decay_rate - fit_rate)
        errors_abs += [error_abs]
        error_rel = np.fabs(decay_rate/fit_rate - 1.0)
        err_rel_tol_percent = err_tol*100.

        logger.info('[Decaying 3D Linear Wave {}]: R-squared of WLS regression = {}'.
                    format(method, r2))
        logger.info('[Decaying 3D Linear Wave {}]: Analytic decay rate = {}'.format(
            method, decay_rate))
        logger.info('[Decaying 3D Linear Wave {}]: Measured decay rate = {}'.format(
            method, fit_rate))
        logger.info('[Decaying 3D Linear Wave {}]: Decay rate absolute error = {}'.format(
            method, error_abs))
        logger.info('[Decaying 3D Linear Wave {}]: Decay rate relative error = {}'.format(
            method, error_rel))

        if error_rel > err_tol:
            logger.warning('[Decaying 3D Linear Wave {}]: decay rate disagrees'
                           ' with prediction by >{}%'.format(method, err_rel_tol_percent))
            analyze_status = False
        else:
            logger.info('[Decaying 3D Linear Wave {}]: decay rate is within '
                        '{}% of analytic value'.format(method, err_rel_tol_percent))
            logger.info('')

    return analyze_status
