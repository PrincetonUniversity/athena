"""
Regression test for 1D Sod shock tube comparing different general EOS implementations.
"""

# Modules
import logging
import numpy as np
import sys
import os
from shutil import move
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison
from scripts.utils.EquationOfState.writeEOS import mk_ideal, write_H
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module
_gammas = [1.1, 1.4, 5./3.]


def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(
                     prob='shock_tube',
                     coord='cartesian',
                     flux='hllc',
                     eos='general/eos_table',
                     **kwargs)
    athena.make()
    src = os.path.join('bin', 'athena')
    dst = os.path.join('bin', 'athena_eos_hllc')
    move(src, dst)
    os.system('mv obj obj_eos_hllc')

    athena.configure(
                     prob='shock_tube',
                     coord='cartesian',
                     flux='hllc',
                     eos='general/hydrogen',
                     **kwargs)
    athena.make()
    src = os.path.join('bin', 'athena')
    dst = os.path.join('bin', 'athena_H')
    move(src, dst)
    os.system('mv obj obj_H')

    athena.configure(
                     prob='shock_tube',
                     coord='cartesian',
                     flux='hllc',
                     eos='general/ideal',
                     **kwargs)
    athena.make()
    src = os.path.join('bin', 'athena')
    dst = os.path.join('bin', 'athena_ideal')
    move(src, dst)
    os.system('mv obj obj_ideal')

    athena.configure(
                     prob='shock_tube',
                     coord='cartesian',
                     flux='hllc',
                     eos='adiabatic',
                     **kwargs)
    athena.make()

    write_H()
    for g in _gammas:
        mk_ideal(g)
        mk_ideal(g, out_type='ascii')


def run(**kwargs):
    arguments = {0: ['hydro/gamma={0:}', 'job/problem_id=Sod_adiabatic_{1:}',
                     'time/ncycle_out=0', 'output1/file_type=vtk']}
    arguments['adiabatic'] = arguments[0][:]
    for i, g in enumerate(_gammas):
        args = [j.format(g, i) for j in arguments['adiabatic']]
        athena.run('hydro/athinput.sod', args, lcov_test_suffix='adiabatic')

    os.system('rm -rf obj')
    os.system('mv obj_ideal obj')
    src = os.path.join('bin', 'athena_ideal')
    dst = os.path.join('bin', 'athena')
    move(src, dst)

    arguments['ideal'] = arguments['adiabatic'][:]
    arguments['ideal'][1] = 'job/problem_id=Sod_ideal_{1:}'
    for i, g in enumerate(_gammas):
        args = [j.format(g, i) for j in arguments['ideal']]
        athena.run('hydro/athinput.sod', args, lcov_test_suffix='adiabatic')

    os.system('rm -rf obj')
    os.system('mv obj_eos_hllc obj')
    src = os.path.join('bin', 'athena_eos_hllc')
    dst = os.path.join('bin', 'athena')
    move(src, dst)

    arguments['binary'] = arguments[0][:]
    arguments['binary'][1] = 'job/problem_id=Sod_eos_hllc_{1:}'
    arguments['binary'].append('hydro/eos_file_name=gamma_is_{0:.3f}.data')
    arguments['ascii'] = arguments[0][:]
    arguments['ascii'][1] = 'job/problem_id=Sod_eos_hllc_ascii_{1:}'
    arguments['ascii'].append('hydro/eos_file_name=gamma_is_{0:.3f}.tab')
    for i, g in enumerate(_gammas):
        arg = [j.format(g, i) for j in arguments['binary']]
        athena.run('hydro/athinput.sod', arg, lcov_test_suffix='eos_hllc')
        arg = [j.format(g, i) for j in arguments['ascii']]
        athena.run('hydro/athinput.sod', arg, lcov_test_suffix='eos_hllc')
    # now run with simple H table

    arguments['binary'][1] = 'job/problem_id=Sod_eos_H_binary'
    arguments['binary'][-1] = 'hydro/eos_file_name=SimpleHydrogen.data'
    arguments['binary'].append('mesh/nx1=512')
    arguments['ascii'][1] = 'job/problem_id=Sod_eos_H_ascii'
    arguments['ascii'][-1] = 'hydro/eos_file_name=SimpleHydrogen.tab'
    arguments['ascii'].append('mesh/nx1=512')

    tmp = ['dl', 'ul', 'pl', 'dr', 'ur', 'pr']
    tmp = ['problem/' + i + '={0:}'for i in tmp] + ['time/tlim={0:}']
    tmp = zip(tmp, [1e-07, 0.00, 3e-8, 1.25e-8, 0.00, 1e-9, .25])
    ic = [i[0].format(i[1]) for i in tmp]
    athena.run('hydro/athinput.sod', ic + arguments['binary'],
               lcov_test_suffix='eos_hllc')
    athena.run('hydro/athinput.sod', ic + arguments['ascii'], lcov_test_suffix='eos_hllc')

    os.system('rm -rf obj')
    os.system('mv obj_H obj')
    src = os.path.join('bin', 'athena_H')
    dst = os.path.join('bin', 'athena')
    move(src, dst)
    arguments['H'] = arguments[0][:] + ['mesh/nx1=512']
    arguments['H'][1] = 'job/problem_id=Sod_eos_H'
    athena.run('hydro/athinput.sod', ic + arguments['H'], lcov_test_suffix='H')
    return 'skip_lcov'


def analyze():
    analyze_status = True
    lbls = ['ideal', 'binary', 'ascii']
    ids = ['ideal', 'eos_hllc', 'eos_hllc_ascii']
    id = 'bin/Sod_{0:}_{1:}.block0.out1.{2:05d}.vtk'
    tolerances = [0.0, 0.0, 1e-6]
    for i, g in enumerate(_gammas):
        for t in [10, 25]:
            x_ref, _, _, data_ref = athena_read.vtk(id.format('adiabatic', i, t))
            loc = tuple([0, 0, slice(None)])
            for var in ['rho', 'press']:
                for j in range(len(lbls)):
                    x_new, _, _, data_new = athena_read.vtk(id.format(ids[j], i, t))
                    norm = comparison.l1_norm(x_ref, data_ref[var][loc])
                    diff = comparison.l1_diff(
                        x_ref, data_ref[var][loc], x_new, data_new[var][loc]) / norm
                    msg = '({0:}). var, err, gamma ='.format(lbls[j])
                    if diff > tolerances[j] or np.isnan(diff):
                        line = 'Eos ideal table test fail '
                        logger.warning(' '.join(map(str, [line, msg, var, diff, g])))
                        analyze_status = False
                    else:
                        line = 'Eos ideal table test pass '
                        logger.debug(' '.join(map(str, [line, msg, var, diff, g])))

    tol = .004
    for t in [10, 25]:
        x_ref, _, _, data_ref = athena_read.vtk(
            'bin/Sod_eos_H.block0.out1.{:05d}.vtk'.format(t))
        x_new, _, _, data_new = athena_read.vtk(
            'bin/Sod_eos_H_binary.block0.out1.{:05d}.vtk'.format(t))
        x_ascii, _, _, data_ascii = athena_read.vtk(
            'bin/Sod_eos_H_ascii.block0.out1.{:05d}.vtk'.format(t))
        loc = tuple([0, 0, slice(None)])
        for var in ['rho', 'press']:
            norm = comparison.l1_norm(x_ref, data_ref[var][loc])
            diff = comparison.l1_diff(
                x_ref, data_ref[var][loc], x_new, data_new[var][loc]) / norm
            msg = ['Eos H table test', 'fail', '(binary). var, err =', var, diff]
            if diff > tol or np.isnan(diff):
                logger.warning(' '.join(map(str, msg)))
                analyze_status = False
            else:
                msg[1] = 'pass'
                logger.debug(' '.join(map(str, msg)))

            diff = comparison.l1_diff(
                x_ref, data_ref[var][loc], x_ascii, data_ascii[var][loc]) / norm
            msg = ['Eos H table test', 'fail', '(ascii). var, err =', var, diff]
            if diff > tol or np.isnan(diff):
                logger.warning(' '.join(map(str, msg)))
                analyze_status = False
            else:
                msg[1] = 'pass'
                logger.debug(' '.join(map(str, msg)))

    return analyze_status
