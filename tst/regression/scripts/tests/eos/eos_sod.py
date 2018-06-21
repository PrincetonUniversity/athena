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

def write_varlist(dlim, elim, varlist, fn=None, log=True, eOp=1.5, ftype='float64', sdim=0):
  if fn is None:
    fn = 'eos_tables.data'
  dlim = np.atleast_1d(dlim)#.astype(ftype)
  elim = np.atleast_1d(elim)#.astype(ftype)
  nd = np.array(varlist[0].shape[0], 'int32')
  ne = np.array(varlist[0].shape[1], 'int32')
  eOp = np.array(eOp)#, ftype)
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
    fn = 'bin/gamma_is_{0:.3f}.data'.format(g)
  write_varlist(dlim[[0,-1]], elim[[0,-1]], varlist, fn=fn, eOp=1. / gm1)
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
  arguments0.append('hydro/EosFn=gamma_is_{0:.3f}.data')
  for i, g in enumerate(_gammas):
      arguments = [j.format(g, i) for j in arguments0]
      athena.run('hydro/athinput.sod', arguments)

def analyze():
  """
  Analyze the output and determine if the test passes.

  This function is called third; nothing from this file is called after it. It is
  responsible for reading whatever data it needs and making a judgment about whether or
  not the test passes. It takes no inputs. Output should be True (test passes) or False
  (test fails).
  """

  for i, g in enumerate(_gammas):
    for t in [10,26]:
        x_ref,_,_,data_ref = athena_read.vtk('bin/Sod_ideal_{0:}.block0.out1.{1:05d}.vtk'.format(i, t))
        x_new,_,_,data_new = athena_read.vtk('bin/Sod_eos_hllc_{0:}.block0.out1.{1:05d}.vtk'.format(i, t))
        loc = [0,0,slice(None)]
        for var in ['rho', 'press']:
            diff = comparison.l1_diff(x_ref, data_ref[var][loc], x_new, data_new[var][loc])
            diff /= comparison.l1_norm(x_ref, data_ref[var][loc])
            if diff > 1e-8 or np.isnan(diff):
              return False

  return True
