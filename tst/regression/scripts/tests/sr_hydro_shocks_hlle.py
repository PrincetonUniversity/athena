# Test script for relativistic hydro shock tubes with HLLE

# Modules
import numpy as np
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison

# Prepare Athena++
def prepare():
  athena.configure('s',
      prob='shock_tube_sr',
      coord='cartesian')
  athena.make()

# Run Athena++
def run():
  arguments = [
      '',
      'output1/file_type=vtk',
      'output1/variable=cons',
      'output1/dt=0.4',
      'time/cfl_number=0.4',
      'time/tlim=0.4',
      'mesh/nx1=400']
  for i in range(1,5):
    arguments[0] = 'job/problem_id=sr_hydro_shock' + repr(i)
    athena.run('hydro_sr/athinput.mb_'+repr(i), arguments)

# Analyze outputs
def analyze():
  for i in range(1,5):
    x_ref,_,_,data_ref = athena.read_vtk('data/sr_hydro_shock{0}_hlle.vtk'.format(i))
    x_new,_,_,data_new = \
        athena.read_vtk('bin/sr_hydro_shock{0}.block0.out1.00001.vtk'.format(i))
    for key in data_ref.iterkeys():
      print(key)
    exit()
  eps_d = comparison.l1_diff(faces_ref, data_ref['D'], faces_new, data_new['D'])
  eps_d /= comparison.l1_norm(faces_ref, data_ref['D'])
  if eps_d < 0.02:
    return True
  else:
    return False
