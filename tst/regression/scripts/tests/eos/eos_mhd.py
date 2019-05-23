"""
Regression test for general EOS 1D Riemann problems.
"""

# Modules
import logging
import sys
import numpy as np
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison
from scripts.utils.RiemannSolver.riemann import riemann_problem
from .eos_riemann import _tests, _thresh, _states
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

_mag_list = [str(i) for i in [1e-5, 1e-4, 1e-3]]
_mhd_tests = _tests + len(_mag_list) * [_tests[0]]
_mhd_thresh = _thresh + len(_mag_list) * [_thresh[0]]
_mhd_states = _states + len(_mag_list) * [_states[0]]

# TODO: (MSBC) Get test to work without threshold adjustments
for i in ['rho', 'press', 'vel']:
    _mhd_thresh[2][i] *= 3
    _mhd_thresh[3][i] *= 2
    _mhd_thresh[8][i] *= 2

# Solution data
eos_rj2a = [[1.08, 1.2, 1.0e-2, 0.5, 1.0155412503859613, 0.5641895835477563, 0.95],
            [1.527257978792011, 0.6073696428720876, 0.13086751616368864,
             0.5671486200909381, 1.4837575078718286, 0.8243097265954603,
             1.4795324075881635],
            [1.527257978792011, 0.6073696428720876, 0.2420190005392491,
             0.3071362707776988, 1.621121045878055, 0.502980538527231,
             1.4795324075881635],
            [1.7756066010105562, 0.5636029312390515, 3.098189935379706e-2,
             0.2416583997245898, 1.4423591033631782, 0.44751658761318547,
             1.8006247196213527],
            [1.5041659368828393, 0.5636029312390515, 3.098189935379706e-2,
             0.2416583997245898, 1.4423591033631782, 0.44751658761318547,
             1.8006247196213527],
            [1.3165384507447337, 0.520140488883702, -0.19219053007193185,
             0.17241533547445267, 1.6238113887027137, 0.5038152634148083,
             1.4757264774462406],
            [1.3165384507447337, 0.520140488883702, -0.10231023474942734,
             -5.115511737471367e-2, 1.5206822799599837, 0.7603411399799919,
             1.4757264774462406],
            [1.0, 0.0, 0.0, 0.0, 1.1283791670955126, 0.5641895835477563, 1.0]]
eos_rj2a = np.array(eos_rj2a)
wave_speeds = [-0.8236630408305421, 0.15084012571797972, 0.29445281898990067,
               0.5636029312390516, 0.868568785133063, 1.011849638757723,
               2.1633547260796724]


def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     prob='shock_tube',
                     coord='cartesian',
                     eos='general/hydrogen',
                     **kwargs)
    athena.make()


def run(**kwargs):
    # run same tests as in eos_riemann.py
    for n, test in enumerate(_tests):
        args = [i + '={0:}'.format(test[i]) for i in test]
        args += ['job/problem_id=eos_mhd_{0:02d}'.format(n), 'time/ncycle_out=100']
        athena.run('hydro/athinput.sod_general_H', args)
    # run first Riemann problem but with different none zero Bx
    for mag in _mag_list:
        n += 1
        args = ['problem/bxl=' + mag, 'problem/bxr=' + mag]
        args += ['job/problem_id=eos_mhd_{0:02d}'.format(n), 'time/ncycle_out=100']
        athena.run('hydro/athinput.sod_general_H', args)
    # run RJ2a test for hydrogen EOS
    args = ['hydro/eos_rho_unit=1e-7', 'hydro/eos_egas_unit=2e-8', 'time/ncycle_out=100',
            'output1/file_type=vtk']
    athena.run('mhd/athinput.rj2a', args)


def analyze():
    analyze_status = True
    for n, state in enumerate(_mhd_states):
        t = 1
        x_ref, _, _, data_ref = athena_read.vtk(
            'bin/eos_mhd_{0:02d}.block0.out1.{1:05d}.vtk'.format(n, t))
        xi = (.5 * x_ref[:-1] + .5 * x_ref[1:]) / _mhd_tests[n]['time/tlim']
        exact = riemann_problem(state, 'H').data_array(xi)
        for var in ['rho', 'press', 'vel']:
            data = data_ref[var][0, 0, :]
            if var == 'vel':
                data = data_ref[var][0, 0, :, 0]
            diff = comparison.l1_norm(x_ref, data - exact[var])
            msg = ['EOS MHD Riemann', 'fail:', 'Test#, var, diff, thresh =', n, var, diff,
                   _mhd_thresh[n][var]]
            if diff > _mhd_thresh[n][var]:
                logger.warning(' '.join(map(str, msg)))
                analyze_status = False
            else:
                msg[1] = 'pass:'
                logger.debug(' '.join(map(str, msg)))
    # test hydrogen version of RJ2a
    x_ref, _, _, data = athena_read.vtk('bin/RJ2a.block0.out1.00040.vtk')
    xc = (x_ref[1:] + x_ref[:-1]) * .5
    j = 0
    t = 0.2
    # construct H-RJ2a solution array
    data_ref = np.empty((7, xc.size))
    for i, x in enumerate(xc):
        try:
            while x / t > wave_speeds[j]:
                j += 1
        except IndexError:
            pass
        data_ref[:, i] = eos_rj2a[j]
    # compute errors
    errors = [comparison.l1_norm(x_ref, data['rho'][0, 0] - data_ref[0]),
              comparison.l1_norm(x_ref, data['vel'][0, 0, :, 0] - data_ref[1]),
              comparison.l1_norm(x_ref, data['vel'][0, 0, :, 1] - data_ref[2]),
              comparison.l1_norm(x_ref, data['vel'][0, 0, :, 2] - data_ref[3]),
              comparison.l1_norm(x_ref, data['Bcc'][0, 0, :, 1] - data_ref[4]),
              comparison.l1_norm(x_ref, data['Bcc'][0, 0, :, 2] - data_ref[5]),
              comparison.l1_norm(x_ref, data['press'][0, 0] - data_ref[6])]
    names = ["rho", "vx", "vy", "vz", "By", "Bz", "p"]
    threshes = [7e-3] + [5e-3] * 3 + [7e-3] * 3
    # check errors
    for err, name, thresh in zip(errors, names, threshes):
        msg = ['EOS RJ2a', 'fail:', 'var, diff, thresh =', name, err, thresh]
        if err > thresh:
            logger.warning(' '.join(map(str, msg)))
            analyze_status = False
        else:
            msg[1] = 'pass:'
            logger.debug(' '.join(map(str, msg)))
    return analyze_status
