# Modules
import numpy as np                             # standard Python module for numerics
from . import eos

default_names = ['p/e(e/rho,rho)', 'e/p(p/rho,rho)', 'asq*rho/p(p/rho,rho)']


def write_varlist(dlim, elim, varlist, fn=None, out_type=None, eOp=1.5,
                  ratios=None, ftype=None, sdim=0, opt=None, var_names=None):
    if var_names is None:
        var_names = default_names
    nvar = len(varlist)
    if out_type is None:
        ext = fn.split('.')[-1]
        if ext in ['tab', 'txt', 'ascii']:
            out_type = 'ascii'
        elif ext in ['hdf5']:
            out_type = 'hdf5'
    if ftype is None:
        ftype = 'float64'
        if out_type == 'hdf5':
            ftype = 'float32'
    if out_type == 'ascii':
        if opt is None:
            opt = {'sep': ' ', 'format': '%.4e'}
    if opt is None:
        opt = {'sep': ''}
    iopt = {'sep': opt.get('sep', '')}
    if fn is None:
        fn = 'eos_tables.data'
    dlim = np.atleast_1d(dlim).astype(ftype)
    elim = np.atleast_1d(elim).astype(ftype)
    nd = np.array(varlist[0].shape[1], 'int32')
    ne = np.array(varlist[0].shape[0], 'int32')
    eOp = np.array(eOp)  # , ftype)
    if ratios is None:
        ratios = [1., eOp, eOp]
        if nvar > 3:
            ratios += [1] * (nvar - 3)
        ratios = np.array(ratios, dtype=ftype)
    if out_type == 'hdf5':
        import h5py
        with h5py.File(fn, 'w') as f:
            f.create_dataset('LogDensLim', data=dlim.astype(ftype))
            f.create_dataset('LogEspecLim', data=elim.astype(ftype))
            f.create_dataset('ratios', data=ratios.astype(ftype))
            for i, d in enumerate(varlist):
                f.create_dataset(var_names[i], data=np.log10(d).astype(ftype))
    elif out_type == 'ascii':
        with open(fn, 'w') as f:
            f.write('# Entries must be space separated.\n')
            f.write('# n_var, n_espec, n_rho\n')
            np.array([nvar, ne, nd]).tofile(f, **iopt)
            f.write('\n# Log espec lim\n')
            elim.tofile(f, **opt)
            f.write('\n# Log rho lim\n')
            dlim.tofile(f, **opt)
            f.write('\n# 1, eint/pres, eint/pres, eint/h\n')
            ratios.tofile(f, **opt)
            f.write('\n')
            for i, d in enumerate(varlist):
                f.write('# ' + var_names[i] + '\n')
                np.savetxt(f, np.log10(d), opt['format'], delimiter=opt['sep'])
    else:  # if binary:
        with open(fn, 'wb') as f:
            np.array([nvar, ne, nd], dtype='int32').tofile(f, **iopt)
            elim.tofile(f, **opt)
            dlim.tofile(f, **opt)
            ratios.tofile(f, **opt)
            out = np.stack(varlist, axis=sdim).astype(ftype)
            np.log10(out).tofile(f, **opt)
    return


def mk_ideal(gamma=5./3., n=2, fn=None, out_type=None):
    dlim = np.linspace(-24., 4., n)
    elim = np.linspace(-10., 20., n)

    e, d = np.meshgrid(1e1**elim, 1e1**dlim)
    g = gamma
    gm1 = g - 1.

    varlist = [gm1, 1. / gm1, g]
    varlist = [np.ones(e.shape) * i for i in varlist]

    if fn is None:
        fn = 'bin/gamma_is_{0:.3f}.'.format(g)
        if out_type == 'ascii':
            fn += 'tab'
        elif out_type == 'hdf5':
            fn += 'hdf5'
        else:
            fn += 'data'
    opt = dict(fn=fn, eOp=1. / gm1, out_type=out_type)
    if out_type == 'hdf5':
        opt['ftype'] = 'float32'
    write_varlist(dlim[[0, -1]], elim[[0, -1]], varlist, **opt)
    return


def write_H(nEspec=256, nRho=64, logEspecLim=None, logRhoLim=None, eOp=1.5,
            binary=True, ascii=True, hdf5=False, ret=False):
    fn = 'bin/SimpleHydrogen'
    if logEspecLim is None:
        logEspecLim = np.array([-2, 1])
    if logRhoLim is None:
        logRhoLim = np.array([-9, -6])
    ld = np.linspace(*logRhoLim, num=nRho)
    le = np.linspace(*logEspecLim, num=nEspec)
    one = np.ones((nEspec, nRho))
    rho = 10**ld
    es = 10**le
    Heos = eos.SimpleHydrogen()
    p = np.empty_like(one)
    e = np.empty_like(one)
    a2p = np.empty_like(one)
    for i in range(nEspec):
        for j in range(nRho):
            p[i, j] = Heos.p_of_rho_es(rho[j], es[i]) / (es[i] * rho[j])
            p0 = rho[j] * es[i] / eOp
            temp = Heos.T_of_rho_p(rho[j], p0)
            e[i, j] = Heos.ei_of_rho_T(rho[j], temp) / p0
            a2p[i, j] = Heos.gamma1(rho[j], temp)
    args = logRhoLim, logEspecLim, [p, e, a2p]
    kwargs = dict(eOp=eOp)
    if binary:
        write_varlist(*args, fn=fn + '.data', **kwargs)
    if ascii:
        write_varlist(*args, fn=fn + '.tab', **kwargs)
    if hdf5:
        write_varlist(*args, fn=fn + '.hdf5', **kwargs)
    return
