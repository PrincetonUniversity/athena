"""
Regression test for general EOS 1D Sod shock tube with ASCII plain-text tables.
"""

# Modules
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import os
from shutil import move                        # moves/renames files
import scripts.utils.athena as athena          # utilities for running Athena++
import scripts.utils.comparison as comparison  # more utilities explicitly for testing
from scripts.utils.EquationOfState.writeEOS import mk_ideal, write_H
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read  # noqa                     # utilities for reading Athena++ data

_gammas = [1.1, 1.4, 5./3.]


def prepare(**kwargs):
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
                     eos='adiabatic',
                     **kwargs)
    athena.make()

    write_H()
    for g in _gammas:
        mk_ideal(g)
        mk_ideal(g, out_type='ascii')


def run(**kwargs):
    arguments0 = ['hydro/gamma={0:}', 'job/problem_id=Sod_ideal_{1:}',
                  'time/ncycle_out=0', 'output1/file_type=vtk']
    for i, g in enumerate(_gammas):
        arguments = [j.format(g, i) for j in arguments0]
        athena.run('hydro/athinput.sod', arguments, lcov_test_suffix='adiabatic')

    os.system('rm -rf obj')
    os.system('mv obj_eos_hllc obj')
    src = os.path.join('bin', 'athena_eos_hllc')
    dst = os.path.join('bin', 'athena')
    move(src, dst)
    arguments0[1] = 'job/problem_id=Sod_eos_hllc_{1:}'
    arguments1 = arguments0[:]
    arguments0.extend(
        ['hydro/eos_file_name=gamma_is_{0:.3f}.data'])
    arguments1[1] = 'job/problem_id=Sod_eos_hllc_ascii_{1:}'
    arguments1.extend(
        ['hydro/eos_file_name=gamma_is_{0:.3f}.tab'])
    for i, g in enumerate(_gammas):
        arguments = [j.format(g, i) for j in arguments0]
        athena.run('hydro/athinput.sod', arguments, lcov_test_suffix='eos_hllc')
        arguments = [j.format(g, i) for j in arguments1]
        athena.run('hydro/athinput.sod', arguments, lcov_test_suffix='eos_hllc')
    # now run with simple H table
    arguments0[0] = 'hydro/gamma=1.6667'
    arguments0[1] = 'job/problem_id=Sod_eos_H_binary'
    arguments0[-1] = 'hydro/eos_file_name=SimpleHydrogen.data'
    arguments1[0] = 'hydro/gamma=1.6667'
    arguments1[1] = 'job/problem_id=Sod_eos_H_ascii'
    arguments1[-1] = 'hydro/eos_file_name=SimpleHydrogen.tab'
    arguments0.append('mesh/nx1=512')
    arguments1.append('mesh/nx1=512')
    tmp = ['dl', 'ul', 'pl', 'dr', 'ur', 'pr']
    tmp = ['problem/' + i + '={0:}'for i in tmp] + ['time/tlim={0:}']
    tmp = zip(tmp, [1e-07, 0.00, 3e-8, 1.25e-8, 0.00, 1e-9, .25])
    ic = [i[0].format(i[1]) for i in tmp]
    athena.run('hydro/athinput.sod', ic + arguments0, lcov_test_suffix='eos_hllc')
    athena.run('hydro/athinput.sod', ic + arguments1, lcov_test_suffix='eos_hllc')

    os.system('rm -rf obj')
    os.system('mv obj_H obj')
    src = os.path.join('bin', 'athena_H')
    dst = os.path.join('bin', 'athena')
    move(src, dst)
    arguments0[1] = 'job/problem_id=Sod_eos_H'
    athena.run('hydro/athinput.sod', ic + arguments0[:-1], lcov_test_suffix='H')
    return 'skip_lcov'


def analyze():
    analyze_status = True
    for i, g in enumerate(_gammas):
        for t in [10, 25]:
            x_ref, _, _, data_ref = athena_read.vtk(
                'bin/Sod_ideal_{0:}.block0.out1.{1:05d}.vtk'.format(i, t))
            x_new, _, _, data_new = athena_read.vtk(
                'bin/Sod_eos_hllc_{0:}.block0.out1.{1:05d}.vtk'.format(i, t))
            x_ascii, _, _, data_ascii = athena_read.vtk(
                'bin/Sod_eos_hllc_ascii_{0:}.block0.out1.{1:05d}.vtk'.format(i, t))
            loc = tuple([0, 0, slice(None)])
            for var in ['rho', 'press']:
                norm = comparison.l1_norm(x_ref, data_ref[var][loc])
                diff = comparison.l1_diff(
                    x_ref, data_ref[var][loc], x_new, data_new[var][loc]) / norm
                if diff > 0.0 or np.isnan(diff):
                    line = ['Eos ideal table test fail (binary). var, err, gamma =',
                            var, diff, g]
                    print(' '.join(map(str, line)))
                    analyze_status = False
                diff = comparison.l1_diff(
                    x_ref, data_ref[var][loc], x_ascii, data_ascii[var][loc]) / norm
                if diff > 1e-6 or np.isnan(diff):
                    line = ['Eos ideal table test fail (ascii). var, err, gamma =',
                            var, diff, g]
                    print(' '.join(map(str, line)))
                    analyze_status = False
    tol = .005
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
            if diff > tol or np.isnan(diff):
                line = ['Eos H table test fail (binary). var, err =', var, diff]
                print(' '.join(map(str, line)))
                analyze_status = False
            diff = comparison.l1_diff(
                x_ref, data_ref[var][loc], x_ascii, data_ascii[var][loc]) / norm
            if diff > tol or np.isnan(diff):
                line = ['Eos H table test fail (ascii). var, err =', var, diff]
                print(' '.join(map(str, line)))
                analyze_status = False

    return analyze_status
