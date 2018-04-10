"""
Regression test of shearing box with 3d MRI
"""

# Modules
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import scripts.utils.athena as athena          # utilities for running Athena++
import scripts.utils.comparison as comparison  # more utilities explicitly for testing
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read                             # utilities for reading Athena++ data

def prepare():
  """
  Configure and make the executable.

  This function is called first. It is responsible for calling the configure script and
  make to create an executable. It takes no inputs and produces no outputs.
  """

  # Configure as though we ran
  #     python configure.py -b -sh --prob=hb3 --flux=hlld
  #     
  # from the athena/ directory. Note that additional -<flag> command-line arguments can
  # be specified as additional '<flag>' arguments before the <key>='<value>' arguments
  # to athena.configure(). Any number of --<key>=<value> command-line arguments can also
  # be supplied. Note athena.configure() expects the values only to be quoted, e.g.
  # --<key>='<value>'.
  athena.configure('b', 'sh',
      prob='hgb',
      flux='hlld',
      eos='isothermal')

  # Call make as though we ran
  #     make clean
  #     make
  # from the athena/ directory.
  athena.make()

def run():
  """
  Run the executable.

  This function is called second. It is responsible for calling the Athena++ binary in
  such a way as to produce testable output. It takes no inputs and produces no outputs.
  """

  # Create list of runtime arguments to override the athinput file. Each element in the
  # list is simply a string of the form '<block>/<field>=<value>', where the contents of
  # the string are exactly what one would type on the command line run running Athena++.
  arguments = [
      'job/problem_id=HGB',
      'output1/file_type=hst','output1/dt=0.062831853',
      'output2/file_type=vtk','output2/variable=prim','output2/dt=31.4616',
      'time/cfl_number=0.3','time/tlim=62.83185','time/nlim=10000',
      'mesh/nx1=32','mesh/x1min=-0.5','mesh/x1max=0.5','mesh/ix1_bc=shear_periodic','mesh/ox1_bc=shear_periodic',
      'mesh/nx2=24','mesh/x2min=-1.57079632679','mesh/x2max=1.57079632679','mesh/ix2_bc=periodic','mesh/ox2_bc=periodic',
      'mesh/nx3=32','mesh/x3min=-0.5','mesh/x3max=0.5','mesh/ix3_bc=periodic','mesh/ox3_bc=periodic',
      'meshblock/nx1=32','meshblock/nx2=24','meshblock/nx3=32',
      'hydro/iso_sound_speed=1.0',
      'problem/beta=100','problem/amp=0.025','problem/ipert=1','problem/ifield=1',
      'problem/Omega0=1.0','problem/qshear=1.5']

  athena.run('mhd/athinput.hgb', arguments)

def analyze():
  """
  Analyze the output and determine if the test passes.

  This function is called third; nothing from this file is called after it. It is
  responsible for reading whatever data it needs and making a judgment about whether or
  not the test passes. It takes no inputs. Output should be True (test passes) or False
  (test fails).
  """

  fname= 'data/mhd_mri_3d.hst'
  omg  = 1.0
  rho0 = 1.0
  cs   = 1.0
  pres = rho0*cs**2
  vol  = 1.0*np.pi*1.0
  index= -500
  a = np.loadtxt(fname, dtype=np.dtype([('me1','f8'),('me2','f8'),('me3','f8'),('stress','f8')]), skiprows=2, usecols=(9,10,11,12))
  me = (a['me1']+a['me2']+a['me3'])
  ref_stress = np.average(a['stress'][index:]/vol/pres)
  ref_me = np.average(me[index:]/vol/pres)
  ref_ratio = ref_me/ref_stress

  fname='bin/HGB.hst'
  b = np.loadtxt(fname, dtype=np.dtype([('me1','f8'),('me2','f8'),('me3','f8'),('stress','f8')]), skiprows=2, usecols=(9,10,11,12))
  me = (b['me1']+b['me2']+b['me3'])
  new_stress = np.average(b['stress'][index:]/vol/pres)
  new_me = np.average(me[index:]/vol/pres)
  new_ratio = new_me/new_stress

  print 'Ref: (stress,ME,ratio) = ',ref_stress,ref_me,ref_ratio
  print 'New: (stress,ME,ratio) = ',new_stress,new_me,new_ratio
  flag = True
  error_rel = np.fabs(new_stress/ref_stress-1.0)
  if error_rel > 0.5:
    print('[MRI-3D]: averaged stress is off by a factor > 2')
    flag = False
  error_rel = np.fabs(new_me/ref_me-1.0)
  if error_rel > 0.5:
    print('[MRI-3D]: averaged magnetic energy is off by a factor > 2')
    flag = False
  error_rel = np.fabs(new_ratio-ref_ratio)
  if error_rel > 1.0:
    print('[MRI-3D]: energy-to-stress ratio is off by an amount > 1.0')
    flag = False

  return flag
