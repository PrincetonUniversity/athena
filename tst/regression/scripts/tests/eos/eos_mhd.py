"""
Regression test for general EOS 1D Riemann problems.
"""

# Modules
import logging
import sys
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


def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('b',
                     prob='shock_tube',
                     coord='cartesian',
                     flux='hlld',
                     eos='general/hydrogen',
                     **kwargs)
    athena.make()


def run(**kwargs):
    # run same tests as in eos_riemann.py
    for n, test in enumerate(_tests):
        args = [i + '={0:}'.format(test[i]) for i in test]
        args += ['job/problem_id=eos_mhd_{0:02d}'.format(n), 'time/ncycle_out=0']
        athena.run('hydro/athinput.sod_general_H', args)
    # run first Riemann problem but with different none zero Bx
    for mag in _mag_list:
        n += 1
        args = ['problem/bxl=' + mag, 'problem/bxr=' + mag]
        args += ['job/problem_id=eos_mhd_{0:02d}'.format(n), 'time/ncycle_out=0']
        athena.run('hydro/athinput.sod_general_H', args)
    # run RJ2a test for hydrogen EOS
    args = ['hydro/eos_rho_unit=1e-7', 'hydro/eos_egas_unit=1e-8', 'time/ncycle_out=0']
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
            msg = ['EOS Riemann', 'fail:', 'Test#, var, diff, thresh =', n, var, diff,
                   _mhd_thresh[n][var]]
            if diff > _mhd_thresh[n][var]:
                logger.warning(' '.join(map(str, msg)))
                analyze_status = False
            else:
                msg[1] = 'pass:'
                logger.debug(' '.join(map(str, msg)))
    return analyze_status
