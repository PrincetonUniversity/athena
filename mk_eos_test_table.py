#! /usr/bin/env python

import numpy as np
import sys

def write_varlist(dlim, elim, varlist, fn=None, log=True, eOp=1.5, ftype='float32', sdim=0):
  if fn is None:
    fn = 'eos_tables.data'
  dlim = np.atleast_1d(dlim).astype(ftype)
  elim = np.atleast_1d(elim).astype(ftype)
  nd = np.array(varlist[0].shape[0], 'int32')
  ne = np.array(varlist[0].shape[1], 'int32')
  eOp = np.array(eOp, ftype)
  with open(fn, 'wb') as f:
    nd.tofile(f)
    dlim.tofile(f)
    ne.tofile(f)
    elim.tofile(f)
    eOp.tofile(f)
    np.array(len(varlist), 'int32').tofile(f)

    out = np.stack(varlist, axis=sdim).astype(ftype)
    if log:
      (np.log10(out)).tofile(f)
    else:
      out.tofile(f)

    f.close()
  return


def mk_ideal(gamma=5./3., n=2, fn=None, mu=.6, R=None):
  dlim = np.linspace(-24., 0., n)
  elim = np.linspace(8., 16., n)
  if R is None:
    R = 1.3807e-16 / (mu * 1.660538921e-24)
  Rinv = 1. / R

  e, d = np.meshgrid(1e1**elim, 1e1**dlim)
  eint = e * d
  g = gamma
  gm1 = g - 1.

  p = (g - 1.) * eint
  h = eint + p
  asq = g * (g - 1.) * e
  T = p /(d * R)

  varlist = [gm1, g * gm1, gm1 * Rinv, 1. / gm1, g, Rinv, 1. / g, gm1, gm1 / g * Rinv]
  varlist = [np.ones(e.shape) * i for i in varlist]

  if fn is None:
    fn = 'gamma_is_{0:.3f}.data'.format(g)
  write_varlist(dlim[[0,-1]], elim[[0,-1]], varlist, fn=fn, eOp=1. / gm1)
  return


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Compute test EOS tables for Athena++.')
    parser = argparse.ArgumentParser()
    parser.add_argument('--gamma', default=5./3., type=float)
    parser.add_argument('--n', default=2, type=int)
    parser.add_argument('--fn', default='bin/eos_table.data', type=str)
    parser.add_argument('--mu', default=.6, type=float)
    parser.add_argument('--R', default=None, type=float)
    p = parser.parse_args()
    opts = vars(p)

    mk_ideal(**opts)
