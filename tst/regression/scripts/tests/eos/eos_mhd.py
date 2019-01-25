"""
Regression test for general EOS 1D Riemann problems.
"""

# Modules
import sys                                     # standard Python module to change path
import scripts.utils.athena as athena          # utilities for running Athena++
import scripts.utils.comparison as comparison  # more utilities explicitly for testing
from scripts.utils.RiemannSolver.riemann import riemann_problem
from .eos_riemann import _tests, _thresh, _states
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read  # noqa                     # utilities for reading Athena++ data


def prepare(**kwargs):
    athena.configure('b',
                     prob='shock_tube',
                     coord='cartesian',
                     flux='hlld',
                     eos='general/hydrogen',
                     **kwargs)
    athena.make()


def run(**kwargs):
    for n, test in enumerate(_tests):
        args = [i + '={0:}'.format(test[i]) for i in test]
        args += ['job/problem_id=eos_riemann_{0:02d}'.format(n), 'time/ncycle_out=100']
        athena.run('mhd/athinput.sod_general_H', args)


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
            if diff > _thresh[n][var]:
                print(' '.join(
                    map(str, ['EOS Riemann fail. Test#, var, diff, thresh =', n,
                        var, diff, _thresh[n][var]])))
                analyze_status = False
    return analyze_status
