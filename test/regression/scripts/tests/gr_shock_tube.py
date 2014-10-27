# Test script for relativistic shock tube using GR framework

# Modules
import numpy as np
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison

# Primary function
def run_test():
  athena.configure('g',
      prob='shock_tube_gr',
      coord='minkowski_cartesian')
  athena.make()
  arguments = [
      'job/problem_id=gr_shock_tube',
      'output1/file_type=tab',
      'output1/variable=cons',
      'output1/data_format=%24.16e',
      'output1/dt=0.4',
      'time/cfl_number=0.4',
      'time/tlim=0.4',
      'mesh/nx1=400']
  athena.run('hydro_sr/athinput.mb_1', arguments)
  headings = ['x', 'D', 'E', 'M1', 'M2', 'M3']
  data_ref = athena.read('data/gr_shock_tube.tab', headings)
  data_new = athena.read('bin/gr_shock_tube.out1.0001.tab', headings)
  faces_ref = np.linspace(-0.5, 0.5, len(data_ref['x'])+1)
  faces_new = np.linspace(-0.5, 0.5, len(data_new['x'])+1)
  eps_d = comparison.l1_diff(faces_ref, data_ref['D'], faces_new, data_new['D'])
  eps_d /= comparison.l1_norm(faces_ref, data_ref['D'])
  if eps_d < 0.02:
    return True
  else:
    return False
