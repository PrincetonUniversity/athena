#!/usr/bin/python

# Runs and generates plots for predefined problems.
# Compares to old Athena in some cases

# Modules
import numpy as np
import argparse
import glob
import os
import subprocess
import re

# Matplotlib
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc

# Main function
def main(**kwargs):

  # Plot settings
  rc('text', usetex=True)
  rc('text.latex', preamble='\usepackage{color}')

  # Extract inputs
  plots_needed = not kwargs['computation_only']
  computation_needed = not kwargs['plot_only']
  movie_needed = kwargs['movie']

  # Case out on problem
  problem = kwargs['problem']
  if problem == 'hydro_sr_shockset':  # relativistic hydro shocks
    if computation_needed:
      settings = [['1', '400', '0.4', '0.4', 'prim'],
                  ['2', '400', '0.4', '0.4', 'prim'],
                  ['3', '400', '0.4', '0.4', 'prim'],
                  ['4', '400', '0.4', '0.4', 'prim']]
      run_old_shock('mb', 'hydro_sr_old_', settings)
      run_new_shock('mb_', 'hydro_sr_new_', settings)
    if plots_needed:
      plot_shockset('plots/hydro_sr_shockset')
  elif problem == 'hydro_sr_shockset_gr':  # relativistic hydro shocks in GR framework
    if computation_needed:
      settings = [['1', '400', '0.4', '0.4', 'prim'],
                  ['2', '400', '0.4', '0.4', 'prim'],
                  ['3', '400', '0.4', '0.4', 'prim'],
                  ['4', '400', '0.4', '0.4', 'prim']]
      run_old_shock('mb', 'hydro_sr_gr_old_', settings)
      run_new_shock_gr('mb_', 'hydro_sr_gr_new_', settings)
    if plots_needed:
      plot_shockset('plots/hydro_sr_shockset_gr', gr=True)
  elif problem == 'hydro_schwarzschild_1d' or problem == 'hydro_schwarzschild_2d':  # 1D radial accretion onto Schwarzschild BH
    if computation_needed:
      configure_string = 'python configure.py \
          --with-geometry=schwarzschild \
          --with-coordinates=schwarzschild \
          --with-eos=adiabatic \
          --with-integrator=vl2 \
          --with-rsolver=hlle \
          --with-reconstruct=plm \
          --with-problem=accretion \
          --enable-general-relativity'
      make_string = 'make all'
      dimension_string = '_1d' if problem == 'hydro_schwarzschild_1d' else '_2d'
      name_string = 'hydro_schwarzschild_geodesic' + dimension_string
      if movie_needed:
        run_string = './athena \
            -i ../inputs/hydro_gr/athinput.geodesic{1} \
            job/problem_id={0} \
            time/tlim=100.0 \
            output1/dt=1.0 \
            output2/dt=1.0'.format(name_string, dimension_string)
      else:
        run_string = './athena \
            -i ../inputs/hydro_gr/athinput.geodesic{1} \
            job/problem_id={0} \
            output1/dt=100.0 \
            time/tlim=100.0'.format(name_string, dimension_string)
      run_new(configure_string, make_string, run_string, name_string)
    if plots_needed:
      plot_accretion_1d('hydro_schwarzschild_geodesic'+dimension_string, movie_needed)
  else:
    print('ERROR: problem not recognized')

# Function for plotting shock set
def plot_shockset(filename, gr=False):

  # Read data
  print('Reading data...')
  data_old_1 = read_athena('data/hydro_sr_old_1.0001.tab',
      ['x', 'rho', 'vx', 'vy', 'vz', 'pgas'])
  data_old_2 = read_athena('data/hydro_sr_old_2.0001.tab',
      ['x', 'rho', 'vx', 'vy', 'vz', 'pgas'])
  data_old_3 = read_athena('data/hydro_sr_old_3.0001.tab',
      ['x', 'rho', 'vx', 'vy', 'vz', 'pgas'])
  data_old_4 = read_athena('data/hydro_sr_old_4.0001.tab',
      ['x', 'rho', 'vx', 'vy', 'vz', 'pgas'])
  if gr:
    data_new_1 = read_athena('data/hydro_sr_gr_new_1.0001.tab',
        ['x', 'rho', 'pgas', 'vx', 'vy', 'vz'])
    data_new_2 = read_athena('data/hydro_sr_gr_new_2.0001.tab',
        ['x', 'rho', 'pgas', 'vx', 'vy', 'vz'])
    data_new_3 = read_athena('data/hydro_sr_gr_new_3.0001.tab',
        ['x', 'rho', 'pgas', 'vx', 'vy', 'vz'])
    data_new_4 = read_athena('data/hydro_sr_gr_new_4.0001.tab',
        ['x', 'rho', 'pgas', 'vx', 'vy', 'vz'])
  else:
    data_new_1 = read_athena('data/hydro_sr_new_1.0001.tab',
        ['x', 'rho', 'pgas', 'vx', 'vy', 'vz'])
    data_new_2 = read_athena('data/hydro_sr_new_2.0001.tab',
        ['x', 'rho', 'pgas', 'vx', 'vy', 'vz'])
    data_new_3 = read_athena('data/hydro_sr_new_3.0001.tab',
        ['x', 'rho', 'pgas', 'vx', 'vy', 'vz'])
    data_new_4 = read_athena('data/hydro_sr_new_4.0001.tab',
        ['x', 'rho', 'pgas', 'vx', 'vy', 'vz'])

  # Plot data
  print('Plotting data...')
  plot_shockset_aux(2, 2, 1, data_old_1, data_new_1, [10, 25, 1],
      [-0.1, 1.1], [0.0, 1.0, 6], 4)
  plot_shockset_aux(2, 2, 2, data_old_2, data_new_2, [12, 25, 1],
      [-0.7, 1.0], [-0.5, 1.0, 4], 5)
  plot_shockset_aux(2, 2, 3, data_old_3, data_new_3, [10, 20, 1],
      [-0.1, 1.1], [0.0, 1.0, 6], 4)
  plot_shockset_aux(2, 2, 4, data_old_4, data_new_4, [10, 1000, 1],
      [-0.1, 1.1], [0.0, 1.0, 6], 4)

  # Save figure
  plt.tight_layout()
  plt.savefig(filename)

# Auxiliary function for plotting shock set
def plot_shockset_aux(rows, cols, position, data_old, data_new, divisors, y_limits,
    y_ticks, y_subdivisions):

  # Plot old data
  plt.subplot(rows, cols, position)
  plt.plot(data_old['x'], data_old['rho']/float(divisors[0]), 'r', alpha=0.25, lw=2.0)
  plt.plot(data_old['x'], data_old['pgas']/float(divisors[1]), 'g', alpha=0.25, lw=2.0)
  plt.plot(data_old['x'], data_old['vx']/float(divisors[2]), 'b', alpha=0.25, lw=2.0)

  # Plot new data
  plt.plot(data_new['x'], data_new['rho']/divisors[0],
      'ro', markeredgecolor='r', ms=0.8)
  plt.plot(data_new['x'], data_new['pgas']/divisors[1],
      'go', markeredgecolor='g', ms=0.8)
  plt.plot(data_new['x'], data_new['vx']/divisors[2],
      'bo', markeredgecolor='b', ms=0.8)

  # Set ticks
  plt.xlim([-0.5, 0.5])
  plt.xticks(np.linspace(-0.5, 0.5, 6))
  plt.gca().xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
  plt.ylim(y_limits)
  plt.yticks(np.linspace(y_ticks[0], y_ticks[1], y_ticks[2]))
  plt.gca().yaxis.set_minor_locator(ticker.AutoMinorLocator(y_subdivisions))

  # Set labels
  plt.xlabel(r'$x$')
  divisor_strings = []
  for divisor in divisors:
    if divisor != 1.0:
      divisor_strings.append('/' + repr(divisor))
    else:
      divisor_strings.append('')
  plt.ylabel((r'$\rho{0}\ \mathrm{{(r)}},\ ' +
      r'p_\mathrm{{g}}{1}\ \mathrm{{(g)}},\ ' +
      r'v_x{2}\ \mathrm{{(b)}}$').\
      format(divisor_strings[0], divisor_strings[1], divisor_strings[2]))

# Function for plotting accretion
def plot_accretion_1d(filename, movie_needed):

  # Prepare to make new plots
  print('deleting old plots...')
  os.system('rm -f plots/{0}.png'.format(filename))
  os.system('rm -f plots/{0}_*.png'.format(filename))
  print('generating new plots...')

  # Set parameters
  mass = 1.0
  gamma_adi = 5.0/3.0

  # Prepare list of files
  filenames_actual = glob.glob('data/{0}.*.tab'.format(filename))

  # Create frames
  for filename_actual in filenames_actual:

    # Get frame number
    match = re.match(r'data/{0}.(\d+).tab'.format(filename), filename_actual)
    frame_string = match.group(1)
    frame = int(frame_string)

    # Read and process data
    data = read_athena(filename_actual, ['r', 'rho', 'pgas', 'v1', 'v2', 'v3'])
    if frame == 0:
      alpha = (1.0 - 2.0*mass/data['r'])**0.5
      u0 = 1.0/(1.0 - 2.0*mass/data['r'])
    epsilon = 1.0/(1.0-gamma_adi) * data['pgas']/data['rho']
    d_norm = alpha * data['rho'] * u0
    e_norm = epsilon * d_norm
    s_norm = (d_norm + gamma_adi * e_norm) * (2.0*mass/data['r'])**0.5 \
        / (1.0 - 2.0*mass/data['r'])

    # Calculate expected values
    if frame == 0:
      v_expected = -(2.0*mass/data['r'])**0.5 * (1.0 - 2.0*mass/data['r'])
      d_expected = 1.0/data['r']**2 * (data['r']/(2.0*mass))**0.5 \
          * (1.0 - 2.0*mass/data['r'])**-0.5
      d_expected *= d_norm[-1] / d_expected[-1]
      e_expected = (data['r']**2 * (2.0*mass/data['r'])**0.5)**-gamma_adi \
          * (1.0 - 2.0*mass/data['r'])**(-(gamma_adi+1.0)/4.0)
      # TODO: decide whether to use HSW value (above) or what is probably correct value (below)
      #e_expected = (1.0 - 2.0*mass/data['r'])**-0.5 \
      #    * (data['r']**2 * (2.0*mass/data['r'])**0.5)**-gamma_adi
      e_expected *= e_norm[-1] / e_expected[-1]
      s_expected = (d_expected + gamma_adi * e_expected) * (2.0*mass/data['r'])**0.5 \
          / (1.0 - 2.0*mass/data['r'])
      rho_expected = d_expected / (alpha * u0)
      pgas_expected = -(gamma_adi-1.0) * e_expected / (alpha * u0)

    # Prepare limits
    if frame == 0:
      rho_limits = [5.0e-3, 2.0e-1]
      v_limits = [0.0, 0.5]
      pgas_limits = [2.0e-7, 1.0e-3]
      d_limits = [5.0e-3, 3.0e-1]
      s_limits = [1.0e-3, 1.0e0]
      e_limits = [3.0e-7, 3.0e-3]

    # Plot data
    plt.figure(figsize=(10,6))
    if frame == 0 and not movie_needed:
      continue
    plot_accretion_aux(2, 3, 1, data['r'], rho_expected, data['rho'], [0.0, 20.0],
        rho_limits, r'$\rho$', True)
    plot_accretion_aux(2, 3, 2, data['r'], -v_expected, -data['v1'], [0.0, 20.0],
        v_limits, r'$v_\mathrm{infall} = -v^1$', False)
    plot_accretion_aux(2, 3, 3, data['r'], pgas_expected, data['pgas'], [0.0, 20.0],
        pgas_limits, r'$p_\mathrm{gas}$', True)
    plot_accretion_aux(2, 3, 4, data['r'], d_expected, d_norm, [0.0, 20.0], d_limits,
        r'$D = \alpha \rho u^0$', True)
    plot_accretion_aux(2, 3, 5, data['r'], s_expected, s_norm, [0.0, 20.0], s_limits,
        r'$S_1 = -\alpha T^0{}_1$', True)
    plot_accretion_aux(2, 3, 6, data['r'], -e_expected, -e_norm, [0.0, 20.0], e_limits,
        r'$E = -\alpha \rho \epsilon u^0$', True)
    plt.tight_layout()
    if movie_needed:
      plt.savefig('plots/{0}_{1}.png'.format(filename, frame_string))
    else:
      plt.savefig('plots/{0}.png'.format(filename))
    plt.clf()

  # Make movie
  if movie_needed:
    print('deleting old movie...')
    os.system('rm -f movies/{0}.mp4'.format(filename))
    print('generating movie...')
    num_digits = len(frame_string)
    try:
      os.system('ffmpeg -loglevel warning -r 10.0 -i plots/{0}_%0{1}d.png -pix_fmt \
          yuv420p -y movies/{0}.mp4'.format(filename, num_digits))
    except OSError as err:
      print('OS Error ({0}): {1}'.format(err.errno, err.strerror))
      exit()    

# Auxiliary function for plotting accretion
def plot_accretion_aux(rows, cols, position, r, vals_expected, vals_actual, r_range,
    val_range, val_label, log_plot):
  plt.subplot(rows, cols, position)
  plt.plot(r, vals_expected, c='gray', lw=3)
  plt.plot(r, vals_actual, 'ko', ms=2)
  if log_plot:
    plt.yscale('log')
  plt.xlim(r_range)
  if val_range is not None:
    plt.ylim(val_range)
  plt.xlabel(r'$r$')
  plt.ylabel(val_label)

# Function for running old Athena
def run_old_shock(input_prefix, output_prefix, settings):

  # Prepare strings
  old_athena_env = 'ATHENA_ROOT'
  try:
    old_directory = os.environ['ATHENA_ROOT']
  except KeyError:
    print('ERROR: {0} must be set to the directory containing old Athena'.\
        format(old_athena_env))
    exit()
  old_configure_string = './configure --with-problem=shkset1d --with-gas=hydro \
      --enable-special-relativity --with-integrator=vl --with-order=2p --with-flux=hllc'
  old_run_string = 'bin/athena -i tst/1D-sr-hydro/athinput.' + input_prefix + '{1} \
      -d {0}/data job/problem_id=' + output_prefix + '{1} job/maxout=1 \
      output1/dat_fmt=%24.16e output1/dt={3} output1/out={5} time/cour_no={4} \
      time/nlim=100000 time/tlim={3} domain1/Nx1={2} domain1/x1min=-0.5 \
      domain1/x1max=0.5'

  # Generate data
  print('deleting old Athena data...')
  try:
    data_files = glob.glob('data/{0}*.tab'.format(output_prefix))
    rm_command = 'rm -f'.split()
    rm_command.extend(data_files)
    subprocess.call(rm_command)
  except OSError as err:
    print('OS Error ({0}): {1}'.format(err.errno, err.strerror))
    exit()
  print('running old Athena...')
  try:
    current_directory = os.getcwd()
    os.chdir(old_directory)
    os.system('make clean &> /dev/null')
    os.system(old_configure_string+' &> /dev/null')
    os.system('make all &> /dev/null')
    for case in settings:
      os.system(old_run_string.format(current_directory, case[0], case[1], case[2],
          case[3], case[4]) + ' &> /dev/null')
    os.chdir(current_directory)
  except OSError as err:
    print('OS Error ({0}): {1}'.format(err.errno, err.strerror))
    exit()

# Function for running new Athena
def run_new_shock(input_prefix, output_prefix, settings):

  # Prepare strings
  new_configure_string = 'python configure.py \
      --with-coordinates=cartesian \
      --with-eos=adiabatic \
      --with-integrator=vl2 \
      --with-rsolver=hlle \
      --with-reconstruct=plm \
      --with-problem=shock_tube \
      --enable-special-relativity'
  new_make_string = 'make all'
  # TODO: change when -d option works
  #new_run_string = './athena -i inputs/hydro_sr/athinput.' + input_prefix + '{1} \
  #    -d {0}/data job/problem_id=' + output_prefix + '{1} output1/variable={5} \
  #    output1/data_format=%24.16e output1/dt={3} time/cfl_number={4} time/nlim=-1 \
  #    time/tlim={3} mesh/nx1={2}'
  new_run_string = './athena \
      -i ../inputs/hydro_sr/athinput.' + input_prefix + '{1} \
      job/problem_id=' + output_prefix + '{1} output1/variable={5} \
      output1/data_format=%24.16e output1/dt={3} time/cfl_number={4} time/nlim=-1 \
      time/tlim={3} mesh/nx1={2}'

  # Generate data
  print('deleting new Athena data...')
  try:
    data_files = glob.glob('data/{0}*.tab'.format(output_prefix))
    rm_command = 'rm -f'.split()
    rm_command.extend(data_files)
    subprocess.call(rm_command)
    # TODO: remove when -d option works
    data_files = glob.glob('../bin/{0}*.tab'.format(output_prefix))
    rm_command = 'rm -f'.split()
    rm_command.extend(data_files)
    subprocess.call(rm_command)
  except OSError as err:
    print('OS Error ({0}): {1}'.format(err.errno, err.strerror))
    exit()
  print('running new Athena...')
  try:
    current_directory = os.getcwd()
    os.chdir('..')
    os.system(new_configure_string + ' &> /dev/null')
    os.system('make clean &> /dev/null')
    os.system(new_make_string + ' &> /dev/null')
    os.chdir('bin')
    for case in settings:
      os.system(new_run_string.format(current_directory, case[0], case[1], case[2],
          case[3], case[4]) + ' &> /dev/null')
      # TODO: remove when -d option works
      os.system('mv {0}*.tab {1}/data/.'.format(output_prefix, current_directory))
    os.chdir(current_directory)
  except OSError as err:
    print('OS Error ({0}): {1}'.format(err.errno, err.strerror))
    exit()

# Function for running new Athena with GR
def run_new_shock_gr(input_prefix, output_prefix, settings):

  # Prepare strings
  new_configure_string = 'python configure.py \
      --with-geometry=minkowski \
      --with-coordinates=cartesian \
      --with-eos=adiabatic \
      --with-integrator=vl2 \
      --with-rsolver=hlle \
      --with-reconstruct=plm \
      --with-problem=shock_tube \
      --enable-general-relativity'
  new_make_string = 'make all'
  # TODO: change when -d option works
  #new_run_string = './athena -i inputs/hydro_sr/athinput.' + input_prefix + '{1} \
  #    -d {0}/data job/problem_id=' + output_prefix + '{1} output1/variable={5} \
  #    output1/data_format=%24.16e output1/dt={3} time/cfl_number={4} time/nlim=-1 \
  #    time/tlim={3} mesh/nx1={2}'
  new_run_string = './athena \
      -i ../inputs/hydro_sr/athinput.' + input_prefix + '{1} \
      job/problem_id=' + output_prefix + '{1} output1/variable={5} \
      output1/data_format=%24.16e output1/dt={3} time/cfl_number={4} time/nlim=-1 \
      time/tlim={3} mesh/nx1={2}'

  # Generate data
  print('deleting new Athena data...')
  try:
    data_files = glob.glob('data/{0}*.tab'.format(output_prefix))
    rm_command = 'rm -f'.split()
    rm_command.extend(data_files)
    subprocess.call(rm_command)
    # TODO: remove when -d option works
    data_files = glob.glob('../bin/{0}*.tab'.format(output_prefix))
    rm_command = 'rm -f'.split()
    rm_command.extend(data_files)
    subprocess.call(rm_command)
  except OSError as err:
    print('OS Error ({0}): {1}'.format(err.errno, err.strerror))
    exit()
  print('running new Athena...')
  try:
    current_directory = os.getcwd()
    os.chdir('..')
    os.system(new_configure_string + ' &> /dev/null')
    os.system('make clean &> /dev/null')
    os.system(new_make_string + ' &> /dev/null')
    os.chdir('bin')
    for case in settings:
      os.system(new_run_string.format(current_directory, case[0], case[1], case[2],
          case[3], case[4]) + ' &> /dev/null')
      # TODO: remove when -d option works
      os.system('mv {0}*.tab {1}/data/.'.format(output_prefix, current_directory))
    os.chdir(current_directory)
  except OSError as err:
    print('OS Error ({0}): {1}'.format(err.errno, err.strerror))
    exit()

# Function for running new Athena in general
def run_new(configure_string, make_string, run_string, name_string):

  # Delete old data
  print('deleting old data...')
  try:
    os.system('rm -f data/{0}.*.tab'.format(name_string))
    os.system('rm -f data/{0}.*.vtk'.format(name_string))
    # TODO: remove when -d option works
    os.system('rm -f ../bin/{0}.*.tab'.format(name_string))
    os.system('rm -f ../bin/{0}.*.vtk'.format(name_string))
  except OSError as err:
    print('OS Error ({0}): {1}'.format(err.errno, err.strerror))
    exit()

  # Compile code
  print('compiling code...')
  try:
    current_directory = os.getcwd()
    os.chdir('..')
    os.system(configure_string + ' &> /dev/null')
    os.system('make clean &> /dev/null')
    os.system(make_string + ' &> /dev/null')
    os.chdir('bin')
  except OSError as err:
    print('OS Error ({0}): {1}'.format(err.errno, err.strerror))
    exit()

  # Generate new data
  print('generating new data...')
  try:
    os.system(run_string + ' &> /dev/null')
    # TODO: remove when -d option works
    os.system('mv {0}.*.tab {1}/data/.'.format(name_string, current_directory))
    os.system('mv {0}.*.vtk {1}/data/.'.format(name_string, current_directory))
    os.chdir(current_directory)
  except OSError as err:
    print('OS Error ({0}): {1}'.format(err.errno, err.strerror))
    exit()

# Function for reading Athena data
def read_athena(filename, headings=None):

  # Read raw data
  with open(filename, 'r') as data_file:
    raw_data = data_file.readlines()

  # Create array of data
  data = []
  for line in raw_data:
    if line.split()[0][0] == '#':
      continue
    row = []
    for val in line.split()[1:]:
      row.append(float(val))
    data.append(row)
  data = np.array(data)

  # Create dict if desired
  if headings is not None:
    data_dict = {}
    for i in range(len(headings)):
      data_dict[headings[i]] = data[:,i]
    return data_dict
  else:
    return data

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('problem', type=str, help='name of problem')
  parser.add_argument('-p', '--plot_only', action='store_true', default=False,
      help='flag indicating computations are not to be redone')
  parser.add_argument('-c', '--computation_only', action='store_true', default=False,
      help='flag indicating no plots are to be made')  
  parser.add_argument('-m', '--movie', action='store_true', default=False,
      help='flag indicating movie should be made')
  args = parser.parse_args()
  main(**vars(args))
