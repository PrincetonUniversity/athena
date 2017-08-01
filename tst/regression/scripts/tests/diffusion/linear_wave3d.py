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
from scipy.optimize import curve_fit

def prepare():
  """
  Configure and make the executable.

  This function is called first. It is responsible for calling the configure script and
  make to create an executable. It takes no inputs and produces no outputs.
  """

  # Configure as though we ran
  #     python configure.py -b --prob=linear_wave3d --flux=hlld 
  #     
  # from the athena/ directory. Note that additional -<flag> command-line arguments can
  # be specified as additional '<flag>' arguments before the <key>='<value>' arguments
  # to athena.configure(). Any number of --<key>=<value> command-line arguments can also
  # be supplied. Note athena.configure() expects the values only to be quoted, e.g.
  # --<key>='<value>'.
  athena.configure('b',
      prob='linear_wave',
      flux='hlld',
      eos='adiabatic') #isothermal')

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
      'job/problem_id=LinWave',
      'output1/file_type=hst','output1/dt=0.01',
      'output2/file_type=vtk','output2/variable=prim','output2/dt=0.05',
      'time/cfl_number=0.3','time/tlim=5.0','time/nlim=-1',
      'mesh/nx1=64','mesh/x1min=0.0','mesh/x1max=3.0','mesh/ix1_bc=periodic','mesh/ox1_bc=periodic',
      'mesh/nx2=32','mesh/x2min=0.0','mesh/x2max=1.5','mesh/ix2_bc=periodic','mesh/ox2_bc=periodic',
      'mesh/nx3=32','mesh/x3min=0.0','mesh/x3max=1.5','mesh/ix3_bc=periodic', 'mesh/ox3_bc=periodic',
      'mesh/num_threads=1','mesh/refinement=none',
      'meshblock/nx1=64','meshblock/nx2=32','meshblock/nx3=32',
      'hydro/iso_sound_speed=1.0',
      'problem/compute_error=false','problem/wave_flag=1',
      'problem/amp=1.0e-4','problem/vflow=0.0',
      'problem/nuiso=0.02','problem/eta_O=0.01']

  # Run Athena++ as though we called
  #     ./athena -i ../inputs/hydro_sr/athinput.mb_1 job/problem_id=gr_shock_tube <...>
  # from the bin/ directory. Note we omit the leading '../inputs/' below when specifying
  # the athinput file.)
  athena.run('mhd/athinput.linear_wave3d', arguments)

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
  
  ksqr = (2.0*np.pi)**2
  rate = 0.5*(0.02+0.01)*ksqr #(nu+eta)*k^2/2 decay rate of Alfven wave

  # Read in the data produced during this test. This will usually be stored in the
  # tst/regression/bin/ directory, but again we omit the first part of the path. Note
  # the file name is what we expect based on the job/problem_id field supplied in run().
  #x_new,_,_,data_new = athena_read.vtk('bin/gr_shock_tube.block0.out1.00001.vtk')
  basename='bin/LinWave.block0.out2.'
  nframe = 101
  dumprate = 0.05
  max_vy = np.zeros(nframe)
  tt     = np.zeros(nframe)
  for i in range(nframe):
    x1f,x2f,x3f,data = athena_read.vtk(basename+str(i).zfill(5)+'.vtk')
    max_vy[i] = np.max(data['vel'][...,1])
    tt[i]     = i*dumprate

  # estimate the decay rate from simulation
  #def func(x,a,b,c):
  #    return a*np.exp(-b*x)+c
  #popt,pcov = curve_fit(func,tt,max_vy)
  #new_rate = popt[1]
  #print '[Decaying Linear Wave-3D]: Ref(decay_rate) = ',rate
  #print '[Decaying Linear Wave-3D]: New(decay_rate) = ',new_rate
  #print 'optimal parameter values: (a,b,c) = ',popt[0],popt[1],popt[2] 
  def func(x,b):
      return max_vy[0]*np.exp(-b*x)
  popt,pcov = curve_fit(func,tt,max_vy)
  new_rate = popt[0]
  print '[Decaying Linear Wave-3D]: Ref(decay_rate) = ',rate
  print '[Decaying Linear Wave-3D]: New(decay_rate) = ',new_rate

  flag = True
  error_rel = np.fabs(rate/new_rate-1.0)
  if error_rel > 0.1:
    print('[Decaying Linear Wave-3D]: decay rate is off by 10 percent')
    flag = False
  else:
    print('[Decaying Linear Wave-3D]: decay rate is within 10 percent precision')

  return flag
