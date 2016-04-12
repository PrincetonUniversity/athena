# Regression test based on Newtonian hydro linear wave convergence problem
#
# Runs a linear wave convergence test in 3D including SMR and checks L1 errors (which
# are computed by the executable automatically and stored in the temporary file
# linearwave_errors.dat)

# Modules
import numpy as np
import math
import sys
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison
sys.path.insert(0, '../../vis/python')

# Prepare Athena++
def prepare():
  athena.configure(
      prob='linear_wave',
      coord='cartesian',
      flux='hllc')
  athena.make()

# Run Athena++
def run():
  # L-going sound wave
  for i in (48,96):
    arguments = [
      'problem/wave_flag=0','problem/vflow=0.0',
      'mesh/nx1=' + repr(i), 'mesh/nx2=' + repr(i/2), 'mesh/nx3=' + repr(i/2),
      'meshblock/nx1=' + repr(i/4),
      'meshblock/nx2=' + repr(i/4),
      'meshblock/nx3=' + repr(i/4),
      'output2/dt=10', 'time/tlim=1.0', 'problem/compute_error=1']
    athena.run('hydro/athinput.linear_wave3d', arguments)
  # entropy wave
  for i in (48,96):
    arguments = [
      'problem/wave_flag=3','problem/vflow=1.0',
      'mesh/nx1=' + repr(i), 'mesh/nx2=' + repr(i/2), 'mesh/nx3=' + repr(i/2),
      'meshblock/nx1=' + repr(i/4),
      'meshblock/nx2=' + repr(i/4),
      'meshblock/nx3=' + repr(i/4),
      'output2/dt=10', 'time/tlim=1.0', 'problem/compute_error=1']
    athena.run('hydro/athinput.linear_wave3d', arguments)
  # R-going sound wave
  for i in (48,96):
    arguments = [
      'problem/wave_flag=4','problem/vflow=0.0',
      'mesh/nx1=' + repr(i), 'mesh/nx2=' + repr(i/2), 'mesh/nx3=' + repr(i/2),
      'meshblock/nx1=' + repr(i/4),
      'meshblock/nx2=' + repr(i/4),
      'meshblock/nx3=' + repr(i/4),
      'output2/dt=10', 'time/tlim=1.0', 'problem/compute_error=1']
    athena.run('hydro/athinput.linear_wave3d', arguments)

# Analyze outputs
def analyze():
  # read data from error file
  filename = 'bin/linearwave-errors.dat'
  data = []
  with open(filename, 'r') as f:
    raw_data = f.readlines()
    for line in raw_data:
      if line.split()[0][0] == '#':
        continue
      data.append([float(val) for val in line.split()])

  # check absolute error and convergence of all three waves
  if data[1][4] > 2.5e-8:
    print "error in L-going sound wave too large",data[1][4]
    return False
  if data[1][4]/data[0][4] > 0.3:
    print "not converging for L-going sound wave",data[0][4],data[1][4]
    return False

  if data[3][4] > 2.0e-8:
    print "error in entropy wave too large",data[3][4]
    return False
  if data[3][4]/data[2][4] > 0.3:
    print "not converging for entropy wave",data[2][4],data[3][4]
    return False

  if data[5][4] > 2.5e-8:
    print "error in R-going sound wave too large",data[5][4]
    return False
  if data[5][4]/data[4][4] > 0.3:
    print "not converging for R-going sound wave",data[4][4],data[5][4]
    return False

  # check error identical for waves in each direction
  if data[1][4] != data[5][4]:
    print "error in L/R-going sound waves not equal",data[1][4],data[5][4]
    return False

  return True
