#!/usr/lib/python2.7/bin/python
#---------------------------------------------------------------------------------------
# configure.py: Athena++ configuration script in python.  Original version by CJW.
#
# When configure.py is run, it uses the command line options and default settings to
# create custom versions of the files Makefile and src/defs.hpp
#
# The following options are implememted:
#   -h  --help        help message
#   --prob=name       use src/pgen/name.cpp as the problem generator
#   --coord=choice    use choice as the coordinate system
#   --eos=choice      use choice as the equation of state
#   --flux=choice     use choice as the Riemann solver
#   -b                enable magnetic fields
#   -s                enable special-relativity
#   -g                enable general-relativity
#   --order=choice    use choice as the spatial reconstruction algorithm
#   --fint=choice     use choice as the fluid time-integration algorithm
#   --cxx=choice      use choice as the C++ compiler
#   --ifov=N          enable N internal fluid output variables 
#   --omp             enable parallelization with OpenMP
#---------------------------------------------------------------------------------------

# Modules
import argparse
import glob
import re

# Set template and output filenames
defsfile_input = 'src/defs.hpp.in'
defsfile_output = 'src/defs.hpp'
makefile_input = 'Makefile.in'
makefile_output = 'Makefile'

# Prepare parser, add each of the arguments 
parser = argparse.ArgumentParser()

# --prob=[name] argument
pgen_directory = 'src/pgen/'
# set choices to list of .cpp files in src/pgen/
pgen_choices = glob.glob(pgen_directory + '*.cpp')
# remove 'src/pgen/' prefix and '.cpp' extension from filenames 
pgen_choices = [choice[len(pgen_directory):-4] for choice in pgen_choices]
parser.add_argument('--prob',
    default='shock_tube',
    choices=pgen_choices,
    help='selects problem generator')

# --coord=[name] argument
parser.add_argument('--coord',
    default='cartesian',
    choices=['cartesian','cylindrical','spherical_polar',\
        'minkowski_cartesian','schwarzschild'],
    help='selects coordinate system')

# --eos=[name] argument
parser.add_argument('--eos',
    default='adiabatic',
    choices=['adiabatic','isothermal'],
    help='selects equation of state')

# --flux=[name] argument
parser.add_argument('--flux',
    default='hlle',
    choices=['hlle','hllc'],
    help='selects Riemann solver')

# -b argument
parser.add_argument('-b',
    action='store_true',
    default=False,
    help='enables magnetic field')

# -s argument
parser.add_argument('-s',
    action='store_true',
    default=False,
    help='enables special relativity')

# -g argument
parser.add_argument('-g',
    action='store_true',
    default=False,
    help='enables general relativity')

# --order=[name] argument
parser.add_argument('--order',
    default='plm',
    choices=['plm'],
    help='selects spatial reconstruction algorithm')

# --fint=[name] argument
parser.add_argument('--fint',
    default='vl2',
    choices=['vl2'],
    help='selects fluid time-integration algorithm')

# --cxx=[name] argument
parser.add_argument('--cxx',
    default='g++',
    choices=['g++','icc'],
    help='selects C++ compiler')

# -omp argument
parser.add_argument('-omp',
    action='store_true',
    default=False,
    help='enable parallelization with OpenMP')

# -ifov=N argument
parser.add_argument('--ifov',
    type=int,
    default=0,
    help='number of internal fluid output variables')

# Parse command-line inputs
args = vars(parser.parse_args())

# Prepare dictionaries of substitutions to be made
definitions = {}
makefile_options = {}

# Set definitions and Makefile options based on above arguments

definitions['PROBLEM'] = args['prob']
makefile_options['PROBLEM_FILE'] = args['prob'] + '.cpp'

definitions['COORDINATE_SYSTEM'] = args['coord']
makefile_options['COORDINATES_FILE'] = args['coord'] + '.cpp'

definitions['NON_BAROTROPIC_EOS'] = '1' if args['eos'] == 'adiabatic' else '0'
makefile_options['EOS_FILE'] = args['eos']
if args['eos'] == 'adiabatic':
  definitions['NFLUID_VARIABLES'] = '5'
if args['eos'] == 'isothermal':
  definitions['NFLUID_VARIABLES'] = '4'

definitions['RSOLVER'] = args['flux']
makefile_options['RSOLVER_FILE'] = args['flux']

definitions['MAGNETIC_FIELDS_ENABLED'] = '1' if args['b'] else '0'
makefile_options['EOS_FILE'] += '_mhd' if args['b'] else '_hydro'
if args['b']:
  definitions['NFIELD_VARIABLES'] = '3'
else:
  definitions['NFIELD_VARIABLES'] = '0'

definitions['RELATIVISTIC_DYNAMICS'] = '1' if args['s'] or args['g'] else '0'
if args['s']:
  makefile_options['EOS_FILE'] += '_sr'
  makefile_options['RSOLVER_FILE'] += '_sr'
if args['g']:
  makefile_options['EOS_FILE'] += '_gr'
  makefile_options['RSOLVER_FILE'] += '_gr'
makefile_options['EOS_FILE'] += '.cpp'
makefile_options['RSOLVER_FILE'] += '.cpp'

definitions['RECONSTRUCT'] = args['order']
makefile_options['RECONSTRUCT_FILE'] = args['order'] + '.cpp'

definitions['FLUID_INTEGRATOR'] = args['fint']
makefile_options['FLUID_INT_FILE'] = args['fint'] + '.cpp'

definitions['COMPILER_CHOICE'] = args['cxx']
makefile_options['COMPILER_CHOICE'] = args['cxx']
if args['cxx'] == 'icc':
  makefile_options['COMPILER_FLAGS'] = '-O3 -xhost -ipo'
  definitions['COMPILER_FLAGS'] = '-O3 -xhost -ipo'
if args['cxx'] == 'g++':
  makefile_options['COMPILER_FLAGS'] = '-O3'
  definitions['COMPILER_FLAGS'] = '-O3'

definitions['OPENMP_OPTION'] = 'OPENMP_PARALLEL' if args['omp'] \
    else 'NOT_OPENMP_PARALLEL'
if args['omp']:
  if args['cxx'] == 'g++':
    makefile_options['COMPILER_FLAGS'] += ' -fopenmp'
    definitions['COMPILER_FLAGS'] += ' -fopenmp'
  if args['cxx'] == 'icc':
    makefile_options['COMPILER_FLAGS'] += ' -openmp'
    definitions['COMPILER_FLAGS'] += ' -openmp'

definitions['NUM_IFOV'] = str(args['ifov'])

# Read templates
with open(defsfile_input, 'r') as current_file:
  defsfile_template = current_file.read()
with open(makefile_input, 'r') as current_file:
  makefile_template = current_file.read()

# Make substitutions
for key,val in definitions.iteritems():
  defsfile_template = re.sub(r'@{0}@'.format(key), val, defsfile_template)
for key,val in makefile_options.iteritems():
  makefile_template = re.sub(r'@{0}@'.format(key), val, makefile_template)

# Write output files
with open(defsfile_output, 'w') as current_file:
  current_file.write(defsfile_template)
with open(makefile_output, 'w') as current_file:
  current_file.write(makefile_template)

# Finish with diagnostic output
print('Your Athena++ distribution has now been configured with the following options:')
print('  Problem generator:       ' + args['prob'])
print('  Coordinate system:       ' + args['coord'])
print('  Equation of state:       ' + args['eos'])
print('  Riemann solver:          ' + args['flux'])
print('  Reconstruction method:   ' + args['order'])
print('  Fluid integrator:        ' + args['fint'])
print('  Compiler and flags:      ' + args['cxx'])
print('  Magnetic fields:         ' + ('enabled' if args['b'] else 'disabled'))
print('  Special relativity:      ' + ('enabled' if args['s'] else 'disabled'))
print('  General relativity:      ' + ('enabled' if args['g'] else 'disabled'))
print('  OpenMP parallelism:      ' + ('enabled' if args['omp'] else 'disabled'))
print('  Internal fluid outvars:  ' + str(args['ifov']))
