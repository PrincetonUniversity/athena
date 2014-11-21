# Test script for relativistic shock tube using GR framework

# Modules
import numpy as np
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison

# Prepare Athena++
# Inputs: (none)
# Outputs: (none)
# Notes:
#   This function is called first and should configure and make the executable
def prepare():

  # Configure as though we ran
  #     python configure.py -g --prob=shock_tube_gr --coord=minkowski_cartesian
  # from the athena/ directory
  athena.configure('g',
      prob='shock_tube_gr',
      coord='minkowski_cartesian')

  # Call make as though we ran
  #     make clean
  #     make
  # from the athena/ directory
  athena.make()

# Run Athena++
# Inputs: (none)
# Outputs: (none)
# Notes:
#   This function is called second and should run the executable
def run():

  # List of runtime arguments to override athinput file
  arguments = [
      'job/problem_id=gr_shock_tube',
      'output1/file_type=tab',
      'output1/variable=cons',
      'output1/data_format=%24.16e',
      'output1/dt=0.4',
      'time/cfl_number=0.4',
      'time/tlim=0.4',
      'mesh/nx1=400']

  # Run athena as though we ran
  #     ./athena -i ../inputs/hydro_sr/athinput.mb_1 <arguments>
  # from the bin/ directory (note though we omit "../inputs/")
  athena.run('hydro_sr/athinput.mb_1', arguments)

# Analyze outputs
# Inputs: (none)
# Outputs:
#   Returns True if test passes
#   Returns False otherwise
# Notes:
#   This function is called third; nothing from this file is called after it
def analyze():

  # Use the tab file reader utility to read the reference and newly generated data
  #   The returned objects are dictionaries with the supplied headings as keys
  headings = ['x', 'D', 'E', 'M1', 'M2', 'M3']
  data_ref = athena.read_tab('data/gr_shock_tube.tab', headings)
  data_new = athena.read_tab('bin/gr_shock_tube.out1.0001.tab', headings)

  # Generate interface locations
  #   These are not stored in the tab files
  #   These are needed for the L1 error routines, which compute exact integrals
  faces_ref = np.linspace(-0.5, 0.5, len(data_ref['x'])+1)
  faces_new = np.linspace(-0.5, 0.5, len(data_new['x'])+1)

  # Check that the integrated error in conserved density is less than 2% the integrated
  # accepted value
  eps_d = comparison.l1_diff(faces_ref, data_ref['D'], faces_new, data_new['D'])
  eps_d /= comparison.l1_norm(faces_ref, data_ref['D'])
  if eps_d < 0.02:
    return True
  else:
    return False
