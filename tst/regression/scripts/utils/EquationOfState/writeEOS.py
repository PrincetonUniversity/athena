# Modules
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import os

_names = ['p/e(e/rho,rho)','e/p(e/rho,rho)', 'asq*rho/p(p/rho,rho)', 'asq*rho/h(h/rho,rho)']

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
                [1., eOp, eOp, eOp / (1 + eOp)], dtype=ftype))
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
            f.write('\n# 1, eint/pres, eint/pres, eint/h\n')
            np.array([1., eOp, eOp, eOp / (1 + eOp)]).tofile(f, **opt)
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

    varlist = [gm1, 1. / gm1, g, gm1]
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
