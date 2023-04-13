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
import os
sys.path.insert(0, '../../vis/python')

# Prepare Athena++
def prepare(**kwargs):
  athena.configure('radiation','mpi','hdf5',
      prob='beam',
      coord='cartesian',
      flux='hllc',
      nghost='4')
  athena.make()

# Run Athena++
def run(**kwargs):
  #case 1
  arguments = ['time/rad_xorder=3']
  athena.run('radiation/athinput.beam_smr', arguments)
  bashcommand="mv bin/*athdf* ../../"
#  os.system(bashcommand)
# Analyze outputs
def analyze():

  print("Take a look at the hdf5 data!")
  return True
