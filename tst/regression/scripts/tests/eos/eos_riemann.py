"""
Regression test for general EOS 1D Riemann problems.
"""

# Modules
import logging
import sys
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison
from scripts.utils.RiemannSolver.riemann import riemann_problem
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

_tests = [[1e-07, 0.00, 0.150, 1.25e-8, 0., 0.062, .25],
          [4e-06, 0.00, 0.120, 4e-08, 0.00, 0.019, 0.3],
          [8e-07, 1.10, 0.006, 4e-07, -1.7, 0.006, 1.5],
          [8e-05, -0.8, 0.095, 8e-05, 0.80, 0.095, .25],
          [5e-09, 1.50, 0.006, 4e-09, -1.8, 0.006, 1.5],
          [8e-05, -0.5, 0.095, 8e-05, 0.90, 0.095, .25]
          ]
_thresh = [[5.1e-10, 6.5e-11, 0.0039],
           [8.0e-09, 1.4e-09, 0.0069],
           [2.0e-07, 2.0e-08, 0.027],
           [5.7e-07, 5.7e-08, 0.0091],
           [4.4e-09, 6.8e-10, 0.077],
           [5.1e-07, 4.7e-08, 0.0062]
           ]
tmp = ['dl', 'ul', 'Tl', 'dr', 'ur', 'Tr']
_states = [dict(zip(tmp, i[:-1])) for i in _tests]
tmp = ['problem/' + i for i in tmp] + ['time/tlim']
_tests = [dict(zip(tmp, i)) for i in _tests]
_thresh = [dict(zip(['rho', 'press', 'vel'], i)) for i in _thresh]


def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(
                     prob='shock_tube',
                     coord='cartesian',
                     flux='hllc',
                     eos='general/hydrogen',
                     **kwargs)
    athena.make()


def run(**kwargs):
    for n, test in enumerate(_tests):
        args = [i + '={0:}'.format(test[i]) for i in test]
        args += ['job/problem_id=eos_riemann_{0:02d}'.format(n), 'time/ncycle_out=0']
        athena.run('hydro/athinput.sod_general_H', args)


def analyze():
    analyze_status = True
    for n, state in enumerate(_states):
        t = 1
        x_ref, _, _, data_ref = athena_read.vtk(
            'bin/eos_riemann_{0:02d}.block0.out1.{1:05d}.vtk'.format(n, t))
        xi = (.5 * x_ref[:-1] + .5 * x_ref[1:]) / _tests[n]['time/tlim']
        exact = riemann_problem(state, 'H').data_array(xi)
        for var in ['rho', 'press', 'vel']:
            data = data_ref[var][0, 0, :]
            if var == 'vel':
                data = data_ref[var][0, 0, :, 0]
            diff = comparison.l1_norm(x_ref, data - exact[var])
            msg = ['EOS Riemann', 'fail:', 'Test#, var, diff, thresh =', n, var, diff,
                   _thresh[n][var]]
            if diff > _thresh[n][var]:
                logger.warning(' '.join(map(str, msg)))
                analyze_status = False
            else:
                msg[1] = 'pass:'
                logger.debug(' '.join(map(str, msg)))
    return analyze_status
