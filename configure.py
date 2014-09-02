#!/usr/bin/python

# Modules
import argparse
import re

# Main function
def main(**kwargs):

  # Set filenames
  makefile_input = 'Makefile.in'
  makefile_output = 'Makefile'
  defsfile_input = 'src/defs.hpp.in'
  defsfile_output = 'src/defs.hpp'

  # Extract flags and packages
  flags = {}
  packages = {}
  definitions = {}
  for key,val in kwargs.iteritems():
    match = re.match(r'enable_(.+)', key)
    if match is not None:
      flags[match.group(1)] = val
      continue
    match = re.match(r'with_(.+)', key)
    if match is not None:
      packages[match.group(1)] = val

  # Check that flags have no obvious conflicts
  if flags['special_relativity'] and flags['general_relativity']:
    print('ERROR: must have at most one of special and general relativity')
    exit()
  if packages['geometry'] is None and flags['general_relativity']:
    print('ERROR: general relativity requires geometry specification')
    exit()
  if packages['geometry'] is not None and not flags['general_relativity']:
    print('ERROR: geometry only applies to general relativity')
    exit()

  # Determine coordinates file
  if packages['geometry'] is None:
    if packages['coordinates'] not in ['cartesian', 'cylindrical', 'spherical_polar']:
      print('ERROR: invalid coordinate choice for Newtonian/flat spacetime')
      exit()
  elif packages['geometry'] == 'minkowski':
    if packages['coordinates'] not in ['cartesian', 'cylindrical', 'spherical_polar']:
      print('ERROR: invalid coordinate choice for Minkowski spacetime')
      exit()
    packages['coordinates'] = 'minkowski_' + packages['coordinates']
  elif packages['geometry'] == 'schwarzschild':
    if packages['coordinates'] != 'schwarzschild':
      print('ERROR: Schwarzschild spacetime must use Schwarzschild coordinates')
      exit()
  elif packages['geometry'] == 'kerr':
    print('ERROR: Kerr spacetime not yet supported')
    exit()

  # Determine equation of state file
  packages['eos'] += '_hydro'
  if packages['eos'] == 'adiabatic_hydro':
    if flags['special_relativity']:
      packages['eos'] += '_sr'
    if flags['general_relativity']:
      packages['eos'] += '_gr'
    definitions['NON_BAROTROPIC_EOS'] = '1'
    definitions['NVAR'] = '5'
  elif packages['eos'] == 'isothermal_hydro':
    if flags['special_relativity'] or flags['general_relativity']:
      print('ERROR: isothermal hydrodynamics do not support relativity')
      exit()
    definitions['NON_BAROTROPIC_EOS'] = '0'
    definitions['NVAR'] = '4'

  # Determine integrator file
  if packages['integrator'] == 'vl2':
    packages['integrator'] = 'van_leer2'

  # Determine reconstruction file
  definitions['NGHOST'] = '2'

  # Determine Riemann solver file
  if packages['rsolver'] == 'roe':
    if flags['special_relativity'] or flags['general_relativity']:
      print('ERROR: Roe solver does not support relativity')
      exit()
  elif packages['rsolver'] == 'hlle':
    if flags['special_relativity']:
      packages['rsolver'] += '_sr'
    if flags['general_relativity']:
      packages['rsolver'] += '_gr'
  elif packages['rsolver'] == 'hllc':
    if flags['general_relativity']:
      print('ERROR: HLLC solver does not support general relativity')
      exit()
    if flags['special_relativity']:
      packages['rsolver'] += '_sr'

  # Determine problem generator file
  if packages['problem'] == 'shu_osher':
    if flags['special_relativity'] or flags['general_relativity']:
      print('ERROR: Roe solver does not support relativity')
      exit()
  elif packages['problem'] == 'shock_tube':
    if flags['special_relativity']:
      packages['problem'] += '_sr'
    if flags['general_relativity']:
      packages['problem'] += '_gr'
  elif packages['problem'] == 'accretion':
    if packages['geometry'] != 'schwarzschild':
      print('ERROR: accretion designed for Schwarzschild geometry')
      exit()
    packages['problem'] += '_gr'

  # Determine remaining definitions
  if flags['special_relativity'] or flags['general_relativity']:
    definitions['RELATIVISTIC_DYNAMICS'] = '1'
  else:
    definitions['RELATIVISTIC_DYNAMICS'] = '0'

  # Read input templates
  with open(makefile_input, 'r') as current_file:
    makefile_template = current_file.read()
  with open(defsfile_input, 'r') as current_file:
    defsfile_template = current_file.read()

  # Make substitutions
  for key,val in packages.iteritems():
    if val is not None:
      makefile_template = re.sub(r'@{0}@'.format(key.upper()+'_FILE'), val+'.cpp', makefile_template)
  for key,val in definitions.iteritems():
    defsfile_template = re.sub(r'@{0}@'.format(key), val, defsfile_template)

  # Write output files
  with open(makefile_output, 'w') as current_file:
    current_file.write(makefile_template)
  with open(defsfile_output, 'w') as current_file:
    current_file.write(defsfile_template)

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--with-geometry',
      default=None,
      choices=[None, 'minkowski', 'schwarzschild'],
      help='spacetime geometry (independent of coordinate representation)')
  parser.add_argument('--with-coordinates',
      default='cartesian',
      choices=['cartesian', 'cylindrical', 'spherical_polar', 'schwarzschild'],
      help='type of coordinate system to use')
  parser.add_argument('--with-eos',
      default='adiabatic',
      choices=['adiabatic', 'isothermal'],
      help='fluid equation of state')
  parser.add_argument('--with-integrator',
      default='vl2',
      choices=['vl2'],
      help='integration algorithm')
  parser.add_argument('--with-rsolver',
      default='hlle',
      choices=['roe', 'hlle', 'hllc'],
      help='Riemann solver for finding fluxes at interfaces')
  parser.add_argument('--with-reconstruct',
      default='plm',
      choices=['plm'],
      help='reconstruction method')
  parser.add_argument('--with-problem',
      required=True,
      default='shock_tube',
      choices=['shock_tube', 'accretion'],
      help='problem generator')
  parser.add_argument('-s', '--enable-special-relativity',
      action='store_true',
      default=False,
      help='flag indicating SR is to be used')
  parser.add_argument('-g', '--enable-general-relativity',
      action='store_true',
      default=False,
      help='flag indicating GR is to be used')
  args = parser.parse_args()
  main(**vars(args))
