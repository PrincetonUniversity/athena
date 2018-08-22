"""
Regression test for general EOS 1D Riemann problems.
"""

# Modules
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import os
from shutil import move                        # moves/renames files
import scripts.utils.athena as athena          # utilities for running Athena++
import scripts.utils.comparison as comparison  # more utilities explicitly for testing
from scripts.utils.RiemannSolver.riemann import riemann_problem
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read                             # utilities for reading Athena++ data


_tests = [[1e-07, 0.00, 0.150, 1.25e-8, 0., 0.062, .25],
         [4e-06, 0.00, 0.120, 4e-08, 0.00, 0.019, 0.3],
         [8e-07, 1.10, 0.006, 4e-07, -1.7, 0.006, 1.5],
         [8e-05, -0.8, 0.095, 8e-05, 0.80, 0.095, .25],
         [5e-09, 1.50, 0.006, 4e-09, -1.8, 0.006, 1.5],
         [8e-05, -0.5, 0.095, 8e-05, 0.90, 0.095, .25]
         ]
_thresh = [dict(rho=6e-10, press=1.2e-10, vel=.006),
          dict(rho=1.1e-8, press=5e-9, vel=.02),
          dict(rho=1.7e-7, press=2.2e-8, vel=.03),
          dict(rho=5e-7, press=6e-8, vel=.008),
          dict(rho=5e-9, press=7e-10, vel=.09),
          dict(rho=5e-7, press=5e-8, vel=.006),
          ]
tmp = ['dl', 'ul', 'Tl', 'dr', 'ur', 'Tr']
_states = [dict(zip(tmp, i[:-1])) for i in _tests]
tmp = ['problem/' + i for i in tmp] + ['time/tlim']
_tests = [dict(zip(tmp, i)) for i in _tests]

def prepare(**kwargs):
  """
  Configure and make the executable.

  This function is called first. It is responsible for calling the configure script and
  make to create an executable. It takes no inputs and produces no outputs.
  """

  athena.configure(
                   prob='hydrogen_shock_tube',
                   coord='cartesian',
                   flux='hllc',
                   eos='hydrogen',
                   **kwargs)
  athena.make()

def run(**kwargs):
  """
  Run the executable.

  This function is called second. It is responsible for calling the Athena++ binary in
  such a way as to produce testable output. It takes no inputs and produces no outputs.
  """

  for n, test in enumerate(_tests):
    args = [i + '={0:}'.format(test[i]) for i in test]
    args += ['job/problem_id=eos_riemann_{0:02d}'.format(n)]
    athena.run('hydro/athinput.cc18', args)

def analyze():
  """
  Analyze the output and determine if the test passes.

  This function is called third; nothing from this file is called after it. It is
  responsible for reading whatever data it needs and making a judgment about whether or
  not the test passes. It takes no inputs. Output should be True (test passes) or False
  (test fails).
  """

  for n, state in enumerate(_states):
    t = 1
    x_ref,_,_,data_ref = athena_read.vtk('bin/eos_riemann_{0:02d}.block0.out1.{1:05d}.vtk'.format(n, t))
    xi = (.5 * x_ref[:-1] + .5 * x_ref[1:]) / _tests[n]['time/tlim']
    exact = riemann_problem(state, 'H').data_array(xi)
    for var in ['rho', 'press', 'vel']:
      data = data_ref[var][0,0,:]
      if var == 'vel':
          data = data_ref[var][0,0,:,0]
      diff = comparison.l1_norm(x_ref, data - exact[var])
      if diff > _thresh[n][var]:
        print('FAIL', n, var, diff, _thresh[n][var])
        return False
  return True
