# Regression test based on the decaying linear wave due to viscosity,
# Ohmic resistivity and thermal conduction. The decay rate is fit and
# then compared with analytic solution, given by:

# Ryu, D., Jones, T. W., & Frank, A. (1995). Numerical Magnetohydrodynamics in
# Astrophysics: Algorithm and Tests for Multidimensional Flow. The Astrophysical Journal,
# 452, 785. doi:10.1086/176347


# Modules
import numpy as np
from numpy.polynomial import Polynomial
import sys
import scripts.utils.athena as athena
sys.path.insert(0, '../../vis/python')
import athena_read                             # noqa
athena_read.check_nan_flag = True

_nu = 0.01
_kappa = _nu*2.0
_eta = _kappa
_c_s = 0.5  # slow mode wave speed of Athena++ linear wave configuration


def prepare(**kwargs):
    athena.configure('b',
                     prob='linear_wave',
                     flux='hlld',
                     eos='adiabatic')
    athena.make()


def run(**kwargs):
    arguments0 = ['output2/dt=0.03',
                  'time/tlim=3.0',
                  'time/ncycle_out=0',
                  # L-going slow wave
                  'problem/wave_flag=2',
                  'problem/amp=1.0e-4',
                  'problem/vflow=0.0',
                  'problem/nu_iso={}'.format(_nu),
                  'problem/eta_ohm={}'.format(_eta),
                  'problem/kappa_iso={}'.format(_kappa)]
    arguments = arguments0 + ['job/problem_id=DecayLinWave']
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
    print("Reynolds number of slow mode linear wave in 3D: {}".format(re_num))

    basename = 'bin/DecayLinWave.block0.out2.'
    nframe = 100
    dumprate = 0.03
    max_vy = np.zeros(nframe)
    tt = np.zeros(nframe)
    # TODO(#237): Replace full 3D dataset output with user-defined .hst output
    for i in range(nframe):
        x1f, x2f, x3f, data = athena_read.vtk(basename + '{:05d}.vtk'.format(i))
        max_vy[i] = np.max(data['vel'][..., 1])
        tt[i] = i*dumprate

    # estimate the decay rate from simulation, using weighted least-squares (WLS)
    yy = np.log(np.abs(max_vy))
    p, [resid, rank, sv, rcond] = Polynomial.fit(tt, yy, 1, w=np.sqrt(max_vy), full=True)
    resid_normal = np.sum((yy - p(tt)) ** 2)
    r2 = 1 - resid_normal/(yy.size*yy.var())
    pnormal = p.convert(domain=(-1, 1))
    fit_rate = -pnormal.coef[-1]

    print('[Decaying 3D Linear Wave]: R-squared of WLS regression = {}'.format(r2))
    print('[Decaying 3D Linear Wave]: analytic decay rate = {}'.format(slow_mode_rate))
    print('[Decaying 3D Linear Wave]: measured decay rate = {}'.format(fit_rate))

    flag = True
    error_rel = np.fabs(slow_mode_rate/fit_rate - 1.0)
    if error_rel > 0.1:
        print('[Decaying 3D Linear Wave]: decay rate disagrees with prediction by >10%')
        flag = False
    else:
        print('[Decaying 3D Linear Wave]: decay rate is within 10% of analytic value')

    return flag
