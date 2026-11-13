# Regression test: decay of a 3D oblique sound wave by fourth-order thermal conduction
# (HydroDiffusion::ThermalFluxIsoFourth) with RK3 + xorder=4, at two background
# densities. The measured decay rate is compared with the exact root of the linear
# dispersion relation for a sound wave damped by the conductive flux
# F = -kappa rho grad T.
#
# Note: the (ideal) linear_wave eigenvector is not the conductive eigenmode, so the IC
# also excites the counter-propagating sound wave at O(kappa k / c_s) relative amplitude,
# and the fitted rate is accurate to ~0.5% at nx1 >= 32 independent of resolution. This
# test therefore bounds the magnitude of the conductive flux (a missing density factor
# changes the d0=2 rate by 2x) but cannot resolve the convergence order of the operator,
# which differs from the second-order operator by only (k dx)^2/12 ~ 0.3% at nx1=32.

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
_gamma = 5.0/3.0
_p0 = 1.0/_gamma  # linear_wave background pressure
resolution_range = [32, 64]
density_range = [1.0, 2.0]
# Upper bound on the relative decay-rate error (measured 2026-10: 0.5-0.7%)
error_rel_tol = 0.03


def prepare(*args, **kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(*args,
                     nghost=4,  # required for fourth-order configurations
                     prob='linear_wave',
                     flux='hllc',
                     eos='adiabatic', **kwargs)
    athena.make()


def run(**kwargs):
    for d0 in density_range:
        for nx in resolution_range:
            arguments = ['output1/dt=0.03',
                         'output2/dt=-1',  # disable .vtk outputs
                         'time/tlim=3.0',
                         'time/ncycle_out=0',
                         'time/xorder=4', 'time/integrator=rk3',
                         'time/correct_ic=true',
                         # L-going sound wave
                         'problem/wave_flag=0',
                         'problem/amp=1.0e-4',
                         'problem/vflow=0.0',
                         'problem/d0={}'.format(d0),
                         'problem/kappa_iso={}'.format(_kappa),
                         'mesh/nx1={}'.format(nx),
                         'mesh/nx2={}'.format(nx//2),
                         'mesh/nx3={}'.format(nx//2),
                         'meshblock/nx1={}'.format(nx),
                         'meshblock/nx2={}'.format(nx//2),
                         'meshblock/nx3={}'.format(nx//2),
                         'job/problem_id=DecayLinWave-d{}-{}'.format(d0, nx)]
            athena.run('hydro/athinput.linear_wave3d', arguments)


def exact_decay_rate(d0):
    # Linearized perturbations ~ exp(i(k x - omega t)) of a gas with p = rho T and
    # conductive flux -kappa rho grad T; eliminating rho1 and p1 from the continuity,
    # momentum, and energy equations gives the cubic
    #   -i w^3 + K w^2 + i k^2 cs^2 w - K k^2 p0/rho0 = 0,  K = (gamma-1) kappa k^2,
    # whose two oscillatory roots are the damped sound waves (the third is the
    # conductively damped entropy mode). To leading order in K the damping rate is
    # (gamma-1)^2 kappa k^2/(2 gamma), as in diffusion/thermal_attenuation.py.
    k = 2.0*np.pi  # wavelength 1 for Athena++'s linear wave setups
    bigk = (_gamma - 1.0)*_kappa*k**2
    cs2 = _gamma*_p0/d0
    roots = np.roots([-1j, bigk, 1j*k**2*cs2, -bigk*k**2*_p0/d0])
    sound = [w for w in roots if abs(w.real) > 1e-8]
    return -sound[0].imag


def analyze():
    analyze_status = True
    for d0 in density_range:
        decay_rate = exact_decay_rate(d0)
        for nx in resolution_range:
            label = '[Decaying 3D Linear Wave, xorder=4, d0={}, nx1={}]'.format(d0, nx)
            filename = 'bin/DecayLinWave-d{}-{}.hst'.format(d0, nx)
            hst_data = athena_read.hst(filename)
            tt = hst_data['time']
            max_vy = hst_data['max-v2']

            # estimate the decay rate from simulation, using weighted least-squares (WLS)
            yy = np.log(np.abs(max_vy))
            p = Polynomial.fit(tt, yy, 1, w=np.sqrt(max_vy))
            fit_rate = -p.convert(domain=(-1, 1)).coef[-1]
            error_rel = np.fabs(decay_rate/fit_rate - 1.0)

            logger.info('{}: exact decay rate = {}, measured = {}, relative error = {}'
                        .format(label, decay_rate, fit_rate, error_rel))
            if error_rel > error_rel_tol:
                logger.warning('{}: decay rate disagrees with the dispersion relation'
                               ' by >{}%'.format(label, error_rel_tol*100.))
                analyze_status = False

    return analyze_status
