"""
Regression test for general EOS 1D Riemann problems.
"""

# Modules
import logging
import os
import sys
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison
from scripts.utils.RiemannSolver.riemann import riemann_problem
from shutil import move
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module
_fluxes = ['hllc']
_exec = os.path.join('bin', 'athena')

_tests = [[1e-07, 0.00, 0.150, 1.25e-8, 0., 0.062, .25],
          [4e-06, 0.00, 0.120, 4e-08, 0.00, 0.019, 0.3],
          [8e-07, 1.10, 0.006, 4e-07, -1.7, 0.006, 1.5],
          [5e-07, 1.50, 0.006, 4e-07, -1.8, 0.006, 1.5],
          [8e-05, -0.8, 0.095, 8e-05, 0.80, 0.095, .25],
          [6e-05, -0.5, 0.095, 8e-05, 0.90, 0.095, .25]
          ]
_thresh = [[5.1e-10, 6.5e-11, 0.0039],
           [8.0e-09, 1.4e-09, 0.0069],
           [2.0e-07, 3.0e-08, 0.0350],
           [4.0e-07, 8.0e-08, 0.1000],
           [5.7e-07, 5.7e-08, 0.0091],
           [4.4e-07, 4.1e-08, 0.0062]
           ]
tmp = ['dl', 'ul', 'Tl', 'dr', 'ur', 'Tr']
_states = [dict(zip(tmp, i[:-1])) for i in _tests]
tmp = ['problem/' + i for i in tmp] + ['time/tlim']
_tests = [dict(zip(tmp, i)) for i in _tests]
_thresh = [dict(zip(['rho', 'press', 'vel'], i)) for i in _thresh]


def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    global _fluxes
    for i in athena.global_config_args:
        tmp = i.split('=')
        if tmp[0] == '--flux' and len(tmp) == 2:
            _fluxes = [tmp[1]]
    for flux in _fluxes:
        athena.configure(
                         prob='shock_tube',
                         coord='cartesian',
                         flux=flux,
                         eos='general/hydrogen',
                         **kwargs)
        # to save time, reuse compiled .o files for all executables created in this test:
        athena.make(clean_first=False)
        move(_exec, _exec + '_' + flux)
        os.system('cp -r obj obj_' + flux)
    os.system('rm -rf obj')


def run(**kwargs):
    for flux in _fluxes:
        move(_exec + '_' + flux, _exec)
        os.system('mv obj_' + flux + ' obj')
        for n, test in enumerate(_tests):
            args = [i + '={0:}'.format(test[i]) for i in test]
            args += ['job/problem_id=eos_riemann_{0:}_{1:02d}'.format(flux, n),
                     'time/ncycle_out=0']
            athena.run('hydro/athinput.sod_general_H', args)


def analyze():
    analyze_status = True
    for flux in _fluxes:
        for n, state in enumerate(_states):
            # the double shock tests are too hard for hlle
            t = 1
            fn = 'bin/eos_riemann_{0:}_{1:02d}.block0.out1.{2:05d}.vtk'
            x_ref, _, _, data_ref = athena_read.vtk(fn.format(flux, n, t))
            xi = (.5 * x_ref[:-1] + .5 * x_ref[1:]) / _tests[n]['time/tlim']
            exact = riemann_problem(state, 'H').data_array(xi)
            for var in ['rho', 'press', 'vel']:
                data = data_ref[var][0, 0, :]
                if var == 'vel':
                    data = data_ref[var][0, 0, :, 0]
                diff = comparison.l1_norm(x_ref, data - exact[var])
                msg = 'Test#, var, diff, thresh = ' + '{:d}, {:>5}' + ', {:.03e}' * 2
                msg = msg.format(n + 1, var, diff, _thresh[n][var])
                msg = ['EOS Riemann ({0:})'.format(flux), 'FAIL:', msg]
                if diff > _thresh[n][var]:
                    logger.warning(' '.join(msg))
                    analyze_status = False
                else:
                    msg[1] = 'pass:'
                    logger.debug(' '.join(msg))
    return analyze_status
