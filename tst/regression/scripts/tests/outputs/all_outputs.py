# Regression test for all output types
#
# Runs Orszag Tang vortex test, restarting the job twice and making history (hst),
# formatted table (.tab), VTK, and HDF5 (if available) outputs.  Then reads last
# version of each file to make sure output data is correct

# Modules
import numpy as np
import math
import sys
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison
sys.path.insert(0, '../../vis/python')
import athena_read

# Prepare Athena++
def prepare(**kwargs):
  athena.configure('b',
      prob='orszag_tang',
      flux='hlld')
  athena.make()

# Run Athena++
def run():
  arguments = ['time/nlim=80']
  athena.run('mhd/athinput.test_outputs', arguments)
  arguments = ['time/nlim=330']
  athena.restart('TestOutputs.00001.rst', arguments)
  arguments = ['time/nlim=-1']
  athena.restart('TestOutputs.00004.rst', arguments)

# Analyze outputs
def analyze():
  # check density max and Vz and Bz components in tab slice
  slice_data = athena_read.tab(filename='bin/TestOutputs.block0.out2.00010.tab', raw=True,
      dimensions=1)
  if max(slice_data[1,:]) < 0.25:
    return False
  if max(slice_data[5,:]) != 0.0:
    return False
  if max(slice_data[8,:]) != 0.0:
    return False

  # check density max and Vz and Bz components in tab sum
  sum_data = athena_read.tab(filename='bin/TestOutputs.block0.out3.00010.tab', raw=True,
      dimensions=1)
  if max(sum_data[1,:]) < 15.0 and max(sum_data[:,1]) > 20.0:
    return False
  if max(sum_data[5,:]) != 0.0:
    return False
  if max(sum_data[8,:]) != 0.0:
    return False

  # check data in VTK dump
#  xf,yf,_,vtk_data = athena_read.vtk(filename='bin/TestOutputs.block0.out4.00010.vtk')
#  if max(xf) != 1.0 and min(xf) != 0.0:
#    return False
#  if max(yf) != 1.0 and min(yf) != 0.0:
#    return False
#  print(vtk_data['dens'])

#  print(max(vtk_data[:,:,:,'dens']), min(vtk_data[:,:,:,'dens']))

  return True
