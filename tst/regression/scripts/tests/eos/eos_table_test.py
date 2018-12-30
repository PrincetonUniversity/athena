"""
Regression test for general EOS 1D Sod shock tube.
"""

# Modules
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import os
from shutil import move                        # moves/renames files
import scripts.utils.athena as athena          # utilities for running Athena++
import scripts.utils.comparison as comparison  # more utilities explicitly for testing
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read                             # utilities for reading Athena++ data

_gammas = [1.1, 1.4, 5./3.]
_names = ['e/p(p/rho,rho)', 'asq*rho/p(p/rho,rho)', 'T*rho/p(p/rho,rho)',
          'p/e(e/rho,rho)', 'asq*rho/e(e/rho,rho)', 'T*rho/e(e/rho,rho)',
          'e/h(h/rho,rho)', 'asq*rho/h(h/rho,rho)', 'T*rho/h(h/rho,rho)']


def write_varlist(dlim, elim, varlist, fn=None, eOp=1.5,
                  ftype='float64', sdim=0, ascii=False, opt=None, hdf5=False):
    if ascii:
        if opt is None:
            opt = {'sep': ' ', 'format': '%.4e'}
    if opt is None:
        opt = {'sep': ''}
    iopt = {'sep': opt.get('sep', '')}
    if fn is None:
        fn = 'eos_tables.data'
    dlim = np.atleast_1d(dlim)  # .astype(ftype)
    elim = np.atleast_1d(elim)  # .astype(ftype)
    nd = np.array(varlist[0].shape[0], 'int32')
    ne = np.array(varlist[0].shape[1], 'int32')
    eOp = np.array(eOp)  # , ftype)
    if hdf5:
        import h5py
        with h5py.File(fn, 'w') as f:
            f.create_dataset('LogDensLim', data=dlim.astype(ftype))
            f.create_dataset('LogEspecLim', data=elim.astype(ftype))
            f.create_dataset('ratios', data=np.array(
                [1, eOp, eOp / (1 + eOp)], dtype=ftype))
            for i, d in enumerate(varlist):
                f.create_dataset(_names[i], data=np.log10(d).astype(ftype))
    elif ascii:
        with open(fn, 'w') as f:
            nvar = len(varlist)
            f.write('# Entries must be space separated.\n')
            f.write('# n_var, n_rho, n_espec\n')
            np.array([nvar, nd, ne]).tofile(f, **iopt)
            f.write('\n# Log espec lim\n')
            elim.tofile(f, **opt)
            f.write('\n# Log rho lim\n')
            dlim.tofile(f, **opt)
            f.write('\n# 1, eint/pres, eint/h, 0, ..., 0 (length nvar)\n')
            np.array([1., eOp, eOp / (1 + eOp)]
                     + [0 for i in range(len(varlist) - 3)]).tofile(f, **opt)
            f.write('\n')
            for i, d in enumerate(varlist):
                f.write('# ' + _names[i] + '\n')
                np.savetxt(f, np.log10(d), opt['format'], delimiter=opt['sep'])
    else:
        with open(fn, 'wb') as f:
            nd.tofile(f, **iopt)
            dlim.tofile(f, **opt)
            ne.tofile(f, **iopt)
            elim.tofile(f, **opt)
            eOp.tofile(f, **opt)
            np.array(len(varlist), 'int32').tofile(f, **iopt)
            out = np.stack(varlist, axis=sdim).astype(ftype)
            np.log10(out).tofile(f, **opt)
    return


def mk_ideal(gamma=5./3., n=2, fn=None, mu=.6, R=None, ascii=False, hdf5=False):
    dlim = np.linspace(-24., 4., n)
    elim = np.linspace(-10., 20., n)
    if R is None:
        Rinv = mu * 1.660538921e-24 / 1.3807e-16
    else:
        Rinv = 1. / R

    e, d = np.meshgrid(1e1**elim, 1e1**dlim)
    #eint = e * d
    g = gamma
    gm1 = g - 1.

    varlist = [gm1, g * gm1, gm1 * Rinv, 1. / gm1, g, Rinv, 1. / g, gm1, gm1 / g * Rinv]
    varlist = [np.ones(e.shape) * i for i in varlist]

    if fn is None:
        fn = 'bin/gamma_is_{0:.3f}.'.format(g)
        if ascii:
            fn += 'tab'
        elif hdf5:
            fn += 'hdf5'
        else:
            fn += 'data'
    opt = dict(fn=fn, eOp=1. / gm1, ascii=ascii, hdf5=hdf5)
    if hdf5:
        opt['ftype'] = 'float32'
    write_varlist(dlim[[0, -1]], elim[[0, -1]], varlist, **opt)
    return


def prepare(**kwargs):
    """
    Configure and make the executable.

    This function is called first. It is responsible for calling the configure script and
    make to create an executable. It takes no inputs and produces no outputs.
    """

    athena.configure(
                     prob='shock_tube',
                     coord='cartesian',
                     flux='hllc',
                     eos='eos_table',
                     **kwargs)
    athena.make()
    src = os.path.join('bin', 'athena')
    dst = os.path.join('bin', 'athena_eos_hllc')
    move(src, dst)

    athena.configure(
                     prob='shock_tube',
                     coord='cartesian',
                     flux='hllc',
                     eos='adiabatic',
                     **kwargs)
    athena.make()

    for g in _gammas:
        mk_ideal(g)
        mk_ideal(g, ascii=True)


def run(**kwargs):
    """``
    Run the executable.

    This function is called second. It is responsible for calling the Athena++ binary in
    such a way as to produce testable output. It takes no inputs and produces no outputs.
    """

    arguments0 = ['hydro/gamma={0:}', 'job/problem_id=Sod_ideal_{1:}', 'time/ncycle_out=0',
                  'output1/file_type=vtk']
    for i, g in enumerate(_gammas):
        arguments = [j.format(g, i) for j in arguments0]
        athena.run('hydro/athinput.sod', arguments)

    src = os.path.join('bin', 'athena_eos_hllc')
    dst = os.path.join('bin', 'athena')
    move(src, dst)
    arguments0[1] = 'job/problem_id=Sod_eos_hllc_{1:}'
    arguments1 = arguments0[:]
    arguments0.extend(
        ['hydro/EOS_file_name=gamma_is_{0:.3f}.data', 'hydro/EOS_file_type=binary'])
    arguments1[1] = 'job/problem_id=Sod_eos_hllc_ascii_{1:}'
    arguments1.extend(
        ['hydro/EOS_file_name=gamma_is_{0:.3f}.tab', 'hydro/EOS_file_type=ascii'])
    for i, g in enumerate(_gammas):
        arguments = [j.format(g, i) for j in arguments0]
        athena.run('hydro/athinput.sod', arguments)
        arguments = [j.format(g, i) for j in arguments1]
        athena.run('hydro/athinput.sod', arguments)


def analyze():
    """
    Analyze the output and determine if the test passes.

    This function is called third; nothing from this file is called after it. It is
    responsible for reading whatever data it needs and making a judgment about whether or
    not the test passes. It takes no inputs. Output should be True (test passes) or False
    (test fails).
    """

    analyze_status = True
    for i, g in enumerate(_gammas):
        for t in [10, 26]:
            x_ref, _, _, data_ref = athena_read.vtk(
                'bin/Sod_ideal_{0:}.block0.out1.{1:05d}.vtk'.format(i, t))
            x_new, _, _, data_new = athena_read.vtk(
                'bin/Sod_eos_hllc_{0:}.block0.out1.{1:05d}.vtk'.format(i, t))
            x_ascii, _, _, data_ascii = athena_read.vtk(
                'bin/Sod_eos_hllc_ascii_{0:}.block0.out1.{1:05d}.vtk'.format(i, t))
            loc = [0, 0, slice(None)]
            for var in ['rho', 'press']:
                norm = comparison.l1_norm(x_ref, data_ref[var][loc])
                diff = comparison.l1_diff(
                    x_ref, data_ref[var][loc], x_new, data_new[var][loc]) / norm
                if diff > 1e-3 or np.isnan(diff):
                    print(
                        ' '.join(map(str, ['Eos table test fail (binary). var, diff, gamma =', var, diff, g])))
                    analyze_status = False
                diff = comparison.l1_diff(
                    x_ref, data_ref[var][loc], x_ascii, data_ascii[var][loc]) / norm
                if diff > 1e-3 or np.isnan(diff):
                    print(
                        ' '.join(map(str, ['Eos table test fail (ascii). var, diff, gamma =', var, diff, g])))
                    analyze_status = False

    return analyze_status
