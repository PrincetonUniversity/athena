"""
Regression test for general EOS 1D Sod shock tube with HDF5 tables.
"""

# Modules
import logging
import numpy as np
import sys
import os
from shutil import move
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison
from .eos_comparison import mk_ideal, write_H
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module

_gammas = [1.1, 1.4, 5./3.]


def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('hdf5',
                     prob='shock_tube',
                     coord='cartesian',
                     flux='hllc',
                     eos='general/eos_table',
                     **kwargs)
    athena.make()
    src = os.path.join('bin', 'athena')
    dst = os.path.join('bin', 'athena_eos_hllc_hdf5')
    move(src, dst)
    os.system('mv obj obj_eos_hllc_hdf5')

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
                     eos='adiabatic',
                     **kwargs)
    athena.make()

    write_H(binary=False, ascii=False, hdf5=True)
    for g in _gammas:
        mk_ideal(g, out_type='hdf5')


def run(**kwargs):
    arguments0 = ['hydro/gamma={0:}', 'job/problem_id=Sod_ideal_{1:}',
                  'time/ncycle_out=0', 'output1/file_type=vtk']
    for i, g in enumerate(_gammas):
        arguments = [j.format(g, i) for j in arguments0]
        athena.run('hydro/athinput.sod', arguments, lcov_test_suffix='adiabatic')

    os.system('rm -rf obj')
    os.system('mv obj_eos_hllc_hdf5 obj')
    src = os.path.join('bin', 'athena_eos_hllc_hdf5')
    dst = os.path.join('bin', 'athena')
    move(src, dst)
    arguments0[1] = 'job/problem_id=Sod_eos_hllc_hdf5_{1:}'
    arguments0.extend(['hydro/eos_file_name=gamma_is_{0:.3f}.hdf5'])
    for i, g in enumerate(_gammas):
        arguments = [j.format(g, i) for j in arguments0]
        logger.debug(' '.join(arguments))
        athena.run('hydro/athinput.sod', arguments, lcov_test_suffix='eos_hllc_hdf5')
    # now run with simple H table
    arguments0[0] = 'hydro/gamma=1.6667'
    arguments0[1] = 'job/problem_id=Sod_eos_H_hdf5'
    arguments0[-1] = 'hydro/eos_file_name=SimpleHydrogen.hdf5'
    tmp = ['dl', 'ul', 'pl', 'dr', 'ur', 'pr']
    tmp = ['problem/' + i + '={0:}'for i in tmp] + ['time/tlim={0:}']
    tmp = zip(tmp, [1e-07, 0.00, 3e-8, 1.25e-8, 0.00, 1e-9, .25])
    ic = [i[0].format(i[1]) for i in tmp]
    athena.run('hydro/athinput.sod', ic + arguments0, lcov_test_suffix='eos_hllc_hdf5')

    os.system('rm -rf obj')
    os.system('mv obj_H obj')
    src = os.path.join('bin', 'athena_H')
    dst = os.path.join('bin', 'athena')
    move(src, dst)
    arguments0[1] = 'job/problem_id=Sod_eos_H'
    logger.debug(ic + arguments0[:-1])
    athena.run('hydro/athinput.sod', ic + arguments0[:-1])
    return 'skip_lcov'


def analyze():
    analyze_status = True
    for i, g in enumerate(_gammas):
        for t in [10, 25]:
            x_ref, _, _, data_ref = athena_read.vtk(
                'bin/Sod_ideal_{0:}.block0.out1.{1:05d}.vtk'.format(i, t))
            x_new, _, _, data_new = athena_read.vtk(
                'bin/Sod_eos_hllc_hdf5_{0:}.block0.out1.{1:05d}.vtk'.format(i, t))
            loc = tuple([0, 0, slice(None)])
            for var in ['rho', 'press']:
                diff = comparison.l1_diff(
                    x_ref, data_ref[var][loc], x_new, data_new[var][loc])
                diff /= comparison.l1_norm(x_ref, data_ref[var][loc])
                msg = ['Eos hdf5 table', 'failed.', 'var, diff, gamma =', var, diff, g]
                if diff > 1e-8 or np.isnan(diff):
                    logger.warning(' '.join(map(str, msg)))
                    analyze_status = False
                else:
                    msg[1] = 'passed.'
                    logger.debug(' '.join(map(str, msg)))

    tol = [3e-3, 7e-4]
    for i, t in enumerate([10, 25]):
        x_ref, _, _, data_ref = athena_read.vtk(
            'bin/Sod_eos_H.block0.out1.{:05d}.vtk'.format(t))
        x_new, _, _, data_new = athena_read.vtk(
            'bin/Sod_eos_H_hdf5.block0.out1.{:05d}.vtk'.format(t))
        loc = tuple([0, 0, slice(None)])
        for var in ['rho', 'press']:
            norm = comparison.l1_norm(x_ref, data_ref[var][loc])
            diff = comparison.l1_diff(
                x_ref, data_ref[var][loc], x_new, data_new[var][loc]) / norm
            msg = ['Eos H table test', 'failed.', 'var, err =', var, diff]
            if diff > tol[i] or np.isnan(diff):
                logger.warning(' '.join(map(str, msg)))
                analyze_status = False
            else:
                msg[1] = 'passed.'
                logger.debug(' '.join(map(str, msg)))

    return analyze_status
