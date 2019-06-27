# Regression test based on the decaying linear wave due to viscosity,
# Ohmic resistivity and thermal conduction. The decay rate is fit and
# then compared with analytic solution, given by:

# Ryu, D., Jones, T. W., & Frank, A. (1995). Numerical Magnetohydrodynamics in
# Astrophysics: Algorithm and Tests for Multidimensional Flow. The Astrophysical Journal,
# 452, 785. doi:10.1086/176347


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

_nu = 0.01
_kappa = _nu*2.0
_eta = _kappa
_c_s = 0.5  # slow mode wave speed of Athena++ linear wave configuration

resolution_range = [32, 64]

method = 'Explicit'
# Upper bound on relative L1 error for each above nx1:
error_rel_tols = [0.22, 0.05]
# lower bound on convergence rate at final (Nx1=64) asymptotic convergence regime
rate_tols = [2.0]  # convergence rate > 3.0 for this particular resolution, sovler

# NOTE: the linear wave convergence test is currently insufficiently discriminatory for
# both STS and non-STS handling of the diffusion terms. pgen/linear_wave.cpp does NOT
# initialize the correct diffusive eigenmodes, nor are we checking the exact analytic
# solution (unlike simpler 1D diffusion/ tests like viscous_diffusion.py) but instead are
# checking a decay rate predicted by linear analysis that becomes invalid at large
# diffusion coefficients (where we would be able to observe the different convergence
# rates of the STS and non-STS solvers).


def prepare(*args, **kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b', *args,
                     prob='linear_wave',
                     flux='hlld',
                     eos='adiabatic', **kwargs)
    athena.make()


def run(**kwargs):
    for i in resolution_range:
        arguments = ['output1/dt=0.03',
                     'output2/dt=-1',  # disable .vtk outputs
                     'time/tlim=3.0',
                     'time/ncycle_out=0',
                     # L-going slow wave
                     'problem/wave_flag=2',
                     'problem/amp=1.0e-4',
                     'problem/vflow=0.0',
                     'problem/nu_iso={}'.format(_nu),
                     'problem/eta_ohm={}'.format(_eta),
                     'problem/kappa_iso={}'.format(_kappa),
                     'mesh/nx1=' + repr(i),
                     'mesh/nx2=' + repr(i/2),
                     'mesh/nx3=' + repr(i/2),
                     'meshblock/nx1=' + repr(i),
                     'meshblock/nx2=' + repr(i/2),
                     'meshblock/nx3=' + repr(i/2),
                     'job/problem_id=DecayLinWave-{}'.format(i)]
        athena.run('mhd/athinput.linear_wave3d', arguments)


def analyze():
    # Lambda=1 for Athena++'s linear wave setups in 1D, 2D, and 3D:
    L = 1.0
    ksqr = (2.0*np.pi/L)**2
    # Equation 3.13 from Ryu, et al. (modified to add thermal conduction term)
    # fast mode decay rate = (19\nu/4 + 3\eta + 3(\gamma-1)^2*kappa/gamma/4)*(2/15)*k^2
    # Equation 3.14 from Ryu, et al. (modified to add thermal conduction term)
    # slow mode decay rate = (4\nu + 3\eta/4 + 3(\gamma-1)^2*kappa/gamma)*(2/15)*k^2
    slow_mode_rate = (4.0*_nu + 3.0*_eta/4.0 + _kappa*4.0/5.0)*(2.0/15.0)*ksqr

    # Equation 3.16
    re_num = (4.0*np.pi**2 * _c_s)/(L*slow_mode_rate)
    analyze_status = True
    errors_abs = []

    for (nx, err_tol) in zip(resolution_range, error_rel_tols):
        logging.info('[Decaying 3D Linear Wave {}]: '
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

        error_abs = np.fabs(slow_mode_rate - fit_rate)
        errors_abs += [error_abs]
        error_rel = np.fabs(slow_mode_rate/fit_rate - 1.0)
        err_rel_tol_percent = err_tol*100.

        logger.info('[Decaying 3D Linear Wave {}]: Reynolds number of slow mode: {}'.
                    format(method, re_num))
        logger.info('[Decaying 3D Linear Wave {}]: R-squared of WLS regression = {}'.
                    format(method, r2))
        logger.info('[Decaying 3D Linear Wave {}]: Analytic decay rate = {}'.format(
            method, slow_mode_rate))
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

    # Check 2nd order convergence rate of solver (STS should only converge at 1st order)
    # SEE ABOVE NOTE
    rate = np.log(errors_abs[-2]/errors_abs[-1]) / (
        np.log(resolution_range[-1]/resolution_range[-2]))
    logger.info('[Decaying 3D Linear Wave]: convergence rate of decay rate error = {}'.
                format(rate))
    if rate < rate_tols[-1]:
        logger.warning('[Decaying 3D Linear Wave]: convergence of decay rate absolute '
                       'error is slower than {}'.format(rate_tols[-1]))
        analyze_status = False
    else:
        logger.info('[Decaying 3D Linear Wave]: convergence of decay rate absolute error '
                    'is at least {}'.format(rate_tols[-1]))

    return analyze_status
