# Regression test: fourth-order convergence of the 3D UCT4 edge EMFs on a steady state
# whose EMF is a nonlinear product of two varying fields.

# The rotating_column problem is an exact, smooth, stationary solution of ideal MHD:
# a differentially rotating column (v = Omega(r) r phi_hat) threaded by an axial field
# B_z(r), with v NOT parallel to B. The EMF E = -v x B = -Omega r B_z r_hat is curl-free,
# so the error at time t is simply u(t) - u(0). The problem is 2D in the plane normal to
# the axis. With one cell along the axis the 2D CT path is used, in which the UCT4 corner
# states are point values and the edge EMFs are the 4th-order face-averaged fluxes of
# B_z. With several (uniform) cells along the axis the 3D CT path is used, in which the
# corner states are averaged along the edge, and the edge-averaged EMF must be formed
# from them with 4th-order accuracy even though <v_t B_z> != <v_t><B_z> along the edge.
# (CPAW, linear waves and the uniform-velocity field loop cannot detect this: their EMF
# perturbation is linear in the varying quantity.)
#
# Checks: 4th-order convergence of the 2D path (axis 3), 4th-order convergence of the
# axial-field error of the 3D path (axis 3), and that the 3D path's axial-field error is
# within a factor of 2 of the 2D path's at the same resolution, for all three axes.
# Measured 2026-10 (macOS/AppleClang -O3), B_z L1 error at tlim = 0.5:
#   2D path: 2.64e-4 (32), 1.06e-5 (64), 4.71e-7 (128); rates 4.64, 4.49
#   3D path with edge-averaged products: 5.46e-4, 5.04e-5, 1.08e-5; rates 3.44, 2.23

# Modules
import logging
import numpy as np
from math import log
import sys
import scripts.utils.athena as athena
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

_tlim = 0.5
_nx3 = 8  # cells along the axis for the 3D CT path (>= NGHOST = 6)
# resolutions in the plane of the column
_res_2d = [32, 64, 128]
_res_3d = [32, 64]
# upper bounds on the RMS-L1 error (all variables) of the 2D path, ~2x measured
_rms_tols_2d = [1.1e-3, 4.0e-5, 2.2e-6]
_rate_tol = 3.8          # on the RMS-L1 (2D) and axial-field L1 (2D and 3D) errors
_ratio_tol = 2.0         # 3D-path axial-field error / 2D-path error at the same N


def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     nghost=6,  # required for fourth-order MHD configurations
                     prob='rotating_column',
                     coord='cartesian',
                     flux='hlld', **kwargs)
    athena.make()


def _args(n, nx3, axis):
    # square cells: the axis extent is nx3 cells of width 1/n; the mesh is (n, n, nx3)
    # permuted so that the nx3 cells lie along the chosen axis
    size = {1: n, 2: n, 3: n}
    half = {1: 0.5, 2: 0.5, 3: 0.5}
    size[axis] = nx3
    half[axis] = 0.5*nx3/n
    args = ['output1/dt=-1', 'output2/dt={}'.format(_tlim),
            'time/tlim={}'.format(_tlim), 'time/ncycle_out=0',
            'time/xorder=4', 'time/integrator=rk3', 'time/correct_ic=true',
            'problem/axis={}'.format(axis),
            'job/problem_id=RotColumn-{}-{}-{}'.format(axis, nx3, n)]
    for d in (1, 2, 3):
        args += ['mesh/nx{}={}'.format(d, size[d]),
                 'meshblock/nx{}={}'.format(d, size[d]),
                 'mesh/x{}min=-{}'.format(d, half[d]),
                 'mesh/x{}max={}'.format(d, half[d])]
    return args


def run(**kwargs):
    for n in _res_2d:
        athena.run('mhd/athinput.rotating_column', _args(n, 1, 3))
    for n in _res_3d:
        athena.run('mhd/athinput.rotating_column', _args(n, _nx3, 3))
    for axis in (1, 2):
        athena.run('mhd/athinput.rotating_column', _args(_res_3d[-1], _nx3, axis))


def _errors(axis, nx3, n):
    # L1 errors of u(tlim) - u(0) from the two vtk dumps; returns (RMS-L1 over all
    # variables, L1 of the axial field component)
    base = 'bin/RotColumn-{}-{}-{}.block0.out2.0000'.format(axis, nx3, n)
    _, _, _, d0 = athena_read.vtk(base + '0.vtk')
    _, _, _, d1 = athena_read.vtk(base + '1.vtk')
    l1 = [np.mean(np.abs(d1['rho'] - d0['rho'])),
          np.mean(np.abs(d1['press'] - d0['press']))]
    for c in range(3):
        l1.append(np.mean(np.abs(d1['vel'][..., c] - d0['vel'][..., c])))
        l1.append(np.mean(np.abs(d1['Bcc'][..., c] - d0['Bcc'][..., c])))
    baxial = l1[2*axis + 1]  # Bcc component axis-1 is at index 3 + 2*(axis-1)
    return np.sqrt(np.mean(np.square(l1))), baxial


def analyze():
    analyze_status = True
    rms_2d, b_2d = {}, {}
    logger.info('---- rotating column, 2D CT path (axis 3, nx3 = 1) ----')
    logger.info('N  |  RMS-L1  |  rate  |  B_axial L1  |  rate')
    for i, n in enumerate(_res_2d):
        rms_2d[n], b_2d[n] = _errors(3, 1, n)
        rate = brate = float('nan')
        if i > 0:
            m = _res_2d[i-1]
            rate = log(rms_2d[m]/rms_2d[n])/log(n/m)
            brate = log(b_2d[m]/b_2d[n])/log(n/m)
        logger.info('%d %g %.2f %g %.2f', n, rms_2d[n], rate, b_2d[n], brate)
        if rms_2d[n] > _rms_tols_2d[i]:
            logger.warning('2D path RMS-L1 error %g at N=%d exceeds %g',
                           rms_2d[n], n, _rms_tols_2d[i])
            analyze_status = False
        if i == len(_res_2d) - 1 and min(rate, brate) < _rate_tol:
            logger.warning('2D path converging at rate %.2f (RMS-L1), %.2f (B_axial) '
                           'slower than %g', rate, brate, _rate_tol)
            analyze_status = False

    logger.info('---- rotating column, 3D CT path (axis 3, nx3 = %d) ----', _nx3)
    logger.info('N  |  RMS-L1  |  B_axial L1  |  rate  |  ratio to 2D path')
    b_3d = {}
    for i, n in enumerate(_res_3d):
        rms, b_3d[n] = _errors(3, _nx3, n)
        brate = float('nan')
        if i > 0:
            m = _res_3d[i-1]
            brate = log(b_3d[m]/b_3d[n])/log(n/m)
        ratio = b_3d[n]/b_2d[n]
        logger.info('%d %g %g %.2f %.2f', n, rms, b_3d[n], brate, ratio)
        if i == len(_res_3d) - 1 and brate < _rate_tol:
            logger.warning('3D path B_axial error converging at rate %.2f slower than %g',
                           brate, _rate_tol)
            analyze_status = False
        if i == len(_res_3d) - 1 and ratio > _ratio_tol:
            logger.warning('3D path B_axial error is %.2f x the 2D path error at N=%d '
                           '(tolerance %g)', ratio, n, _ratio_tol)
            analyze_status = False

    n = _res_3d[-1]
    for axis in (1, 2):
        rms, b = _errors(axis, _nx3, n)
        ratio = b/b_2d[n]
        logger.info('3D CT path, axis %d, N=%d: RMS-L1 %g, B_axial L1 %g, '
                    '%.2f x the 2D path error', axis, n, rms, b, ratio)
        if ratio > _ratio_tol:
            logger.warning('3D path (axis %d) B_axial error is %.2f x the 2D path error '
                           'at N=%d (tolerance %g)', axis, ratio, n, _ratio_tol)
            analyze_status = False

    return analyze_status
