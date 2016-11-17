"""
Example test script.

This is a complete, working example that can be run as part of the test suite. It does a
simple test of a relativistic shock tube using the GR framework. There are many comments
in order to make this file self-explanatory, but the actual working code is only 28
lines long.

There are three functions defined here:
    prepare()
    run()
    analyze()
All three must be defined with the same names and no inputs in order to make a working
script. They are called in sequence from the main test script run_tests.py. Additional
support functions can be defined here, to be called by the three primary functions.

Heavy use is made of support utilities defined in scripts/utils/athena.py. These are
general-purpose Python scripts that interact with Athena++. They should be used whenever
possible, since they work together to compile and run Athena++ and read the output data.
In particular, proper use of them will result in all files outside tst/regression/ being
in the same state after the test as they were before (including whatever configured
version of Athena++ existed in athena/bin/), as well as cleaning up any new files
produced by the test.
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
      prob='hb3',
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
      'job/problem_id=HB3',
      'output1/file_type=hst','output1/dt=62.831853',
      'output2/file_type=vtk','output2/variable=prim','output2/dt=62831.853',
      'time/cfl_number=0.4','time/tlim=50265.482','time/nlim=500000',
      'mesh/nx1=64','mesh/x1min=-0.5','mesh/x1max=0.5','mesh/ix1_bc=user','mesh/ox1_bc=user',
      'mesh/nx2=64','mesh/x2min=-0.5','mesh/x2max=0.5','mesh/ix2_bc=periodic','mesh/ox2_bc=periodic',
      'mesh/nx3=1','mesh/x3min=-0.5','mesh/x3max=0.5',
      'hydro/iso_sound_speed=0.00408',
      'problem/beta=4000','problem/amp=0.01','problem/ipert=1','problem/ifield=1',
      'problem/Omega0=1.0e-3','problem/qshear=1.5','problem/shboxcoord=2']

  # Run Athena++ as though we called
  #     ./athena -i ../inputs/hydro_sr/athinput.mb_1 job/problem_id=gr_shock_tube <...>
  # from the bin/ directory. Note we omit the leading '../inputs/' below when specifying
  # the athinput file.)
  athena.run('mhd/athinput.hb3', arguments)

def analyze():
  """
  Analyze the output and determine if the test passes.

  This function is called third; nothing from this file is called after it. It is
  responsible for reading whatever data it needs and making a judgment about whether or
  not the test passes. It takes no inputs. Output should be True (test passes) or False
  (test fails).
  """

  # Read in reference data. The tst/regression/data/ directory has reference runs for
  # comparing future output of the code. We only need to specify file names starting
  # with "data/". Now athena_read.vtk() returns four objects: the x-interface locations,
  # the y-interface locations, the z-interface locations, and the values of the
  # variables themselves. In a 1D problem we ignore the second and third returned
  # values, assigning them to the _ variable as is typical Python style.
  #x_ref,_,_,data_ref = athena_read.vtk('data/sr_hydro_shock1_hlle.vtk')
  fname= 'data/mhd_mri_2d.hst'
  omg  = 1e-3
  rho0 = 1.0
  cs   = 0.00408
  pres = rho0*cs**2
  vol  = 1.0*1.0
  index= -400
  a = np.loadtxt(fname, dtype=np.dtype([('me1','f8'),('me2','f8'),('me3','f8'),('stress','f8')]), skiprows=2, usecols=(9,10,11,12))
  me = (a['me1']+a['me2']+a['me3'])
  ref_stress = np.average(a['stress'][index:]/vol/pres)
  ref_me = np.average(me[index:]/vol/pres)
  ref_ratio = ref_me/ref_stress
  # Read in the data produced during this test. This will usually be stored in the
  # tst/regression/bin/ directory, but again we omit the first part of the path. Note
  # the file name is what we expect based on the job/problem_id field supplied in run().
  #x_new,_,_,data_new = athena_read.vtk('bin/gr_shock_tube.block0.out1.00001.vtk')
  fname='bin/HB3.hst'
  b = np.loadtxt(fname, dtype=np.dtype([('me1','f8'),('me2','f8'),('me3','f8'),('stress','f8')]), skiprows=2, usecols=(9,10,11,12))
  me = (b['me1']+b['me2']+b['me3'])
  new_stress = np.average(b['stress'][index:]/vol/pres)
  new_me = np.average(me[index:]/vol/pres)
  new_ratio = new_me/new_stress

  print '[MRI-2D]: Ref(stress,ME,ratio) = ',ref_stress,ref_me,ref_ratio
  print '[MRI-2D]: New(stress,ME,ratio) = ',new_stress,new_me,new_ratio
  flag = True
  error_rel = np.fabs(new_stress/ref_stress-1.0)
  if error_rel > 0.5:
    print('[MRI-2D]: averaged stress is off by a factor > 2')
    flag = False
  error_rel = np.fabs(new_me/ref_me-1.0)
  if error_rel > 0.5:
    print('[MRI-2D]: averaged magnetic energy is off by a factor > 2')
    flag = False
  error_rel = np.fabs(new_ratio-ref_ratio)
  if error_rel > 1.0:
    print('[MRI-2D]: energy-to-stress ratio is off by an amount > 1.0')
    flag = False

  return flag
