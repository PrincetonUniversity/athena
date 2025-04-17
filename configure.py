#!/usr/bin/env python
# ---------------------------------------------------------------------------------------
# configure.py: Athena++ configuration script in python. Original version by CJW.
#
# When configure.py is run, it uses the command line options and default settings to
# create custom versions of the files Makefile and src/defs.hpp from the template files
# Makefile.in and src/defs.hpp.in respectively.
#
# The following options are implememted:
#   -h  --help        help message
#   --prob=name       use src/pgen/name.cpp as the problem generator
#   --coord=xxx       use xxx as the coordinate system
#   --eos=xxx         use xxx as the equation of state
#   --flux=xxx        use xxx as the Riemann solver
#   --nghost=xxx      set NGHOST=xxx
#   --nscalars=xxx    set NSCALARS=xxx
#   --nspecies=xxx    set NSPECIES=xxx
#   -eos_table        enable EOS table
#   -b                enable magnetic fields
#   -s                enable special relativity
#   -g                enable general relativity
#   -t                enable interface frame transformations for GR
#   -debug            enable debug flags (-g -O0); override other compiler options
#   -coverage         enable compiler-dependent code coverage flags
#   -float            enable single precision (default is double)
#   -mpi              enable parallelization with MPI
#   -omp              enable parallelization with OpenMP
#   -hdf5             enable HDF5 output (requires the HDF5 library)
#   --hdf5_path=path  path to HDF5 libraries (requires the HDF5 library)
#   -fft              enable FFT (requires the FFTW library)
#   --fftw_path=path  path to FFTW libraries (requires the FFTW library)
#   --grav=xxx        use xxx as the self-gravity solver
#   --chemistry=xxx   enable chemistry, use xxx as chemical network
#   --kida_rates=xxx  add special rates xxx to kida network
#   --chem_ode_solver=xxx  ode solver xxx for chemistry
#   --cvode_path=path path to CVODE libraries (cvode solver requires the library)
#   --chem_radiation=xxx  enable radiative transfer, use xxx for integrator
#   --cxx=xxx         use xxx as the C++ compiler (works w/ or w/o -mpi)
#   --ccmd=name       use name as the command to call the (non-MPI) C++ compiler
#   --mpiccmd=name    use name as the command to call the MPI C++ compiler
#   --gcovcmd=name    use name as the command to call the gcov utility
#   --cflag=string    append string whenever invoking compiler/linker
#   --include=path    use -Ipath when compiling
#   --lib_path=path   use -Lpath when linking
#   --lib=xxx         use -lxxx when linking
#   -nr_radiation        turn on non-relativistic radiation transport
#   -implicit_radiation  implicit radiation transport module
#   -cr                  enable cosmic ray transport
#   -crdiff              enable cosmic ray diffusion with Multigrid
# ----------------------------------------------------------------------------------------

# Modules
import argparse
import glob
import re


# Set template and output filenames
makefile_input = 'Makefile.in'
makefile_output = 'Makefile'
defsfile_input = 'src/defs.hpp.in'
defsfile_output = 'src/defs.hpp'

# --- Step 1. Prepare parser, add each of the arguments ------------------
athena_description = (
    "Prepare custom Makefile and defs.hpp for compiling Athena++ solver"
)
athena_epilog = (
    "Full documentation of options available at "
    "https://github.com/PrincetonUniversity/athena-public-version/wiki/Configuring"
)
parser = argparse.ArgumentParser(description=athena_description, epilog=athena_epilog)

# --prob=[name] argument
pgen_directory = 'src/pgen/'
# set pgen_choices to list of .cpp files in src/pgen/
pgen_choices = glob.glob(pgen_directory + '*.cpp')
# remove 'src/pgen/' prefix and '.cpp' extension from each filename
pgen_choices = [choice[len(pgen_directory):-4] for choice in pgen_choices]
parser.add_argument('--prob',
                    default='shock_tube',
                    choices=pgen_choices,
                    help='select problem generator')

# --coord=[name] argument
parser.add_argument(
    '--coord',
    default='cartesian',
    choices=[
        'cartesian',
        'cylindrical',
        'spherical_polar',
        'minkowski',
        'schwarzschild',
        'kerr-schild',
        'gr_user'],
    help='select coordinate system')

# --eos=[name] argument
parser.add_argument('--eos',
                    default='adiabatic',
                    choices=['adiabatic', 'isothermal', 'general/eos_table',
                             'general/hydrogen', 'general/ideal'],
                    help='select equation of state')

# --flux=[name] argument
parser.add_argument('--flux',
                    default='default',
                    choices=['default', 'hlle', 'hllc', 'lhllc', 'hlld', 'lhlld', 'roe', 'llf'], # noqa
                    help='select Riemann solver')

# --nghost=[value] argument
parser.add_argument('--nghost',
                    default='2',
                    help='set number of ghost zones')

# --nscalars=[value] argument
parser.add_argument('--nscalars',
                    default='0',
                    help='set number of passive scalars')

# --nspecies=[value] argument
parser.add_argument('--nspecies',
                    default='0',
                    help='set number of chemical species')

# -b argument
parser.add_argument('-b',
                    action='store_true',
                    default=False,
                    help='enable magnetic field')

# -sts argument
parser.add_argument('-sts',
                    action='store_true',
                    default=False,
                    help='enable super-time-stepping')

# -s argument
parser.add_argument('-s',
                    action='store_true',
                    default=False,
                    help='enable special relativity')

# -g argument
parser.add_argument('-g',
                    action='store_true',
                    default=False,
                    help='enable general relativity')

# -t argument
parser.add_argument('-t',
                    action='store_true',
                    default=False,
                    help='enable interface frame transformations for GR')

# -debug argument
parser.add_argument('-debug',
                    action='store_true',
                    default=False,
                    help='enable debug flags; override other compiler options')

# -coverage argument
parser.add_argument('-coverage',
                    action='store_true',
                    default=False,
                    help='enable compiler-dependent code coverage flag')

# -float argument
parser.add_argument('-float',
                    action='store_true',
                    default=False,
                    help='enable single precision')

# -mpi argument
parser.add_argument('-mpi',
                    action='store_true',
                    default=False,
                    help='enable parallelization with MPI')

# -omp argument
parser.add_argument('-omp',
                    action='store_true',
                    default=False,
                    help='enable parallelization with OpenMP')

# --grav=[name] argument
parser.add_argument('--grav',
                    default='none',
                    choices=['none', 'fft', 'mg'],
                    help='select self-gravity solver')

# -fft argument
parser.add_argument('-fft',
                    action='store_true',
                    default=False,
                    help='enable FFT')

# --fftw_path argument
parser.add_argument('--fftw_path',
                    default='',
                    help='path to FFTW libraries')

# --chemistry argument
parser.add_argument('--chemistry',
                    default=None,
                    choices=["gow17", "H2", "kida", "G14Sod"],
                    help='select chemical network')

# --kida_rates argument
parser.add_argument('--kida_rates',
                    default=None,
                    choices=["H2", "gow17", "nitrogen", "nitrogen_gas_Sipila"],
                    help='select special rates for kida network')

# --chem_radiation argument
parser.add_argument('--chem_radiation',
                    default=None,
                    choices=["const", "six_ray"],
                    help='enable and select radiative transfer method for chemistry')

# --chem_ode_solver argument
parser.add_argument('--chem_ode_solver',
                    default=None,
                    choices=["cvode", "forward_euler"],
                    help='ode solver for chemistry')

# --cvode_path argument
parser.add_argument('--cvode_path',
                    type=str,
                    default='',
                    help='path to CVODE libraries')

# -hdf5 argument
parser.add_argument('-hdf5',
                    action='store_true',
                    default=False,
                    help='enable HDF5 Output')

# -h5double argument
parser.add_argument('-h5double',
                    action='store_true',
                    default=False,
                    help='enable double precision HDF5 output')

# --hdf5_path argument
parser.add_argument('--hdf5_path',
                    default='',
                    help='path to HDF5 libraries')

# -nr_radiation argument
parser.add_argument('-nr_radiation',
                    action='store_true',
                    default=False,
                    help='enable non-relativistic radiative transfer')

# -implicit_radiation argument
parser.add_argument('-implicit_radiation',
                    action='store_true',
                    default=False,
                    help='enable radiative transfer')

# -cosmic ray argument
parser.add_argument('-cr',
                    action='store_true',
                    default=False,
                    help='enable cosmic ray transport')

# -cosmic ray diffusion argument
parser.add_argument('-crdiff',
                    action='store_true',
                    default=False,
                    help='enable implicit cosmic ray diffusion')

# The main choices for --cxx flag, using "ctype[-suffix]" formatting, where "ctype" is the
# major family/suite/group of compilers and "suffix" may represent variants of the
# compiler version and/or predefined sets of compiler options. The C++ compiler front ends
# are the main supported/documented options and are invoked on the command line, but the C
# front ends are also acceptable selections and are mapped to the matching C++ front end:
# gcc -> g++, clang -> clang++, icc-> icpc
cxx_choices = [
    'g++',
    'g++-simd',
    'icpx',
    'icpx-old',
    'icpc',
    'icpc-debug',
    'icpc-phi',
    'cray',
    'clang++',
    'clang++-simd',
    'clang++-apple',
    'aocc',
]


def c_to_cpp(arg):
    arg = arg.replace('gcc', 'g++', 1)
    arg = arg.replace('icc', 'icpc', 1)
    arg = arg.replace('icx', 'icpx', 1)
    if arg == 'clang':
        arg = 'clang++'
    else:
        arg = arg.replace('clang-', 'clang++-', 1)
    return arg


# --cxx=[name] argument
parser.add_argument(
    '--cxx',
    default='g++',
    type=c_to_cpp,
    choices=cxx_choices,
    help='select C++ compiler and default set of flags (works w/ or w/o -mpi)')

# --ccmd=[name] argument
parser.add_argument('--ccmd',
                    default=None,
                    help='override for command to use to call (non-MPI) C++ compiler')

# --mpiccmd=[name] argument
parser.add_argument('--mpiccmd',
                    default=None,
                    help='override for command to use to call MPI C++ compiler')

# --gcovcmd=[name] argument
parser.add_argument('--gcovcmd',
                    default=None,
                    help='override for command to use to call Gcov utility in Makefile')

# --cflag=[string] argument
parser.add_argument('--cflag',
                    default=None,
                    help='additional string of flags to append to compiler/linker calls')

# --include=[name] arguments
parser.add_argument(
    '--include',
    default=[],
    action='append',
    help=('extra path for included header files (-I<path>); can be specified multiple '
          'times'))

# --lib_path=[name] arguments
parser.add_argument(
    '--lib_path',
    default=[],
    action='append',
    help=('extra path for linked library files (-L<path>); can be specified multiple '
          'times'))

# --lib=[name] arguments
parser.add_argument(
    '--lib',
    default=[],
    action='append',
    help='name of library to link against (-l<lib>); can be specified multiple times')

# Parse command-line inputs
args = vars(parser.parse_args())

# --- Step 2. Test for incompatible arguments ----------------------------

# Set default flux; HLLD for MHD, HLLC for hydro, HLLE for isothermal hydro or any GR
if args['flux'] == 'default':
    if args['g']:
        args['flux'] = 'hlle'
    elif args['b']:
        args['flux'] = 'hlld'
    elif args['eos'] == 'isothermal':
        args['flux'] = 'hlle'
    else:
        args['flux'] = 'hllc'

# Check Riemann solver compatibility
if args['flux'] == 'hllc' and args['eos'] == 'isothermal':
    raise SystemExit('### CONFIGURE ERROR: HLLC flux cannot be used with isothermal EOS')
if args['flux'] == 'hllc' and args['b']:
    raise SystemExit('### CONFIGURE ERROR: HLLC flux cannot be used with MHD')
if args['flux'] == 'lhllc' and args['eos'] == 'isothermal':
    raise SystemExit('### CONFIGURE ERROR: LHLLC flux cannot be used with isothermal EOS') # noqa
if args['flux'] == 'lhllc' and args['b']:
    raise SystemExit('### CONFIGURE ERROR: LHLLC flux cannot be used with MHD')
if args['flux'] == 'hlld' and not args['b']:
    raise SystemExit('### CONFIGURE ERROR: HLLD flux can only be used with MHD')
if args['flux'] == 'lhlld' and args['eos'] == 'isothermal':
    raise SystemExit('### CONFIGURE ERROR: LHLLD flux cannot be used with isothermal EOS') # noqa
if args['flux'] == 'lhlld' and not args['b']:
    raise SystemExit('### CONFIGURE ERROR: LHLLD flux can only be used with MHD')

# Check relativity
if args['s'] and args['g']:
    raise SystemExit('### CONFIGURE ERROR: '
                     + 'GR implies SR; the -s option is restricted to pure SR')
if args['t'] and not args['g']:
    raise SystemExit('### CONFIGURE ERROR: Frame transformations only apply to GR')
if args['g'] and not args['t'] and args['flux'] not in ('llf', 'hlle'):
    raise SystemExit('### CONFIGURE ERROR: Frame transformations required for {0}'
                     .format(args['flux']))
if args['g'] and args['coord'] in ('cartesian', 'cylindrical', 'spherical_polar'):
    raise SystemExit('### CONFIGURE ERROR: GR cannot be used with {0} coordinates'
                     .format(args['coord']))
if not args['g'] and args['coord'] not in ('cartesian', 'cylindrical', 'spherical_polar'):
    raise SystemExit('### CONFIGURE ERROR: '
                     + args['coord'] + ' coordinates only apply to GR')
if args['eos'] == 'isothermal':
    if args['s'] or args['g']:
        raise SystemExit('### CONFIGURE ERROR: '
                         + 'Isothermal EOS is incompatible with relativity')
if args['eos'][:8] == 'general/':
    if args['s'] or args['g']:
        raise SystemExit('### CONFIGURE ERROR: '
                         + 'General EOS is incompatible with relativity')
    if args['flux'] not in ['hllc', 'hlld']:
        raise SystemExit('### CONFIGURE ERROR: '
                         + 'General EOS is incompatible with flux ' + args['flux'])

if args['chemistry'] is None and args['chem_ode_solver'] is not None:
    raise SystemExit('### CONFIGURE ERROR: must choose chemistry network for ode solver.')

if args['chemistry'] is not None and args['chem_ode_solver'] is None:
    raise SystemExit('### CONFIGURE ERROR: must choose ode solver for chemistry.')

if args['chemistry'] == 'kida' and args['kida_rates'] is None:
    raise SystemExit('### CONFIGURE ERROR: must provide rates for kida chemistry.')

if args['chem_radiation'] == 'six_ray' and (
        args['chemistry'] != 'gow17' and args['chemistry'] != 'kida'):
    raise SystemExit('### CONFIGURE ERROR: six ray radiation'
                     + 'only compatible with gow17 or kida chemistry enabled.')

if args['chem_ode_solver'] == 'cvode' and args['cvode_path'] == '':
    raise SystemExit('### CONFIGURE ERROR: must provide library path to cvode.')

if args['g'] and (args['nr_radiation'] or args['implicit_radiation']):
    raise SystemExit('### CONFIGURE ERROR: '
                     + ' GR is incompatible with nr_radiation or implicit_radiation')

if args['nr_radiation'] and args['implicit_radiation']:
    raise SystemExit('### CONFIGURE ERROR: '
                     + ' nr_radiation and implicit_radiation cannot be used together')

# --- Step 3. Set definitions and Makefile options based on above argument

# Prepare dictionaries of substitutions to be made
definitions = {}
makefile_options = {}
makefile_options['LOADER_FLAGS'] = ''

# --prob=[name] argument
definitions['PROBLEM'] = makefile_options['PROBLEM_FILE'] = args['prob']

# --coord=[name] argument
definitions['COORDINATE_SYSTEM'] = makefile_options['COORDINATES_FILE'] = args['coord']

# --eos=[name] argument
definitions['NON_BAROTROPIC_EOS'] = '0' if args['eos'] == 'isothermal' else '1'
makefile_options['EOS_FILE'] = args['eos']
definitions['EQUATION_OF_STATE'] = args['eos']
# set number of hydro variables for adiabatic/isothermal
definitions['GENERAL_EOS'] = '0'
makefile_options['GENERAL_EOS_FILE'] = 'noop'
definitions['EOS_TABLE_ENABLED'] = '0'
if args['eos'] == 'isothermal':
    definitions['NHYDRO_VARIABLES'] = '4'
elif args['eos'] == 'adiabatic':
    definitions['NHYDRO_VARIABLES'] = '5'
else:
    definitions['GENERAL_EOS'] = '1'
    makefile_options['GENERAL_EOS_FILE'] = 'general'
    definitions['NHYDRO_VARIABLES'] = '5'
    if args['eos'] == 'general/eos_table':
        definitions['EOS_TABLE_ENABLED'] = '1'

# --flux=[name] argument
definitions['RSOLVER'] = makefile_options['RSOLVER_FILE'] = args['flux']

# --nghost=[value] argument
definitions['NUMBER_GHOST_CELLS'] = args['nghost']

# --nscalars=[value] argument
definitions['NUMBER_PASSIVE_SCALARS'] = args['nscalars']

# --nspecies=[value] argument
definitions['NUMBER_CHEMICAL_SPECIES'] = args['nspecies']

# -b argument
# set variety of macros based on whether MHD/hydro or adi/iso are defined
if args['b']:
    definitions['MAGNETIC_FIELDS_ENABLED'] = '1'
    if definitions['GENERAL_EOS'] != '0':
        makefile_options['GENERAL_EOS_FILE'] += '_mhd'
    else:
        makefile_options['EOS_FILE'] += '_mhd'
    definitions['NFIELD_VARIABLES'] = '3'
    makefile_options['RSOLVER_DIR'] = 'mhd/'
    if args['flux'] == 'hlle' or args['flux'] == 'llf' or args['flux'] == 'roe':
        makefile_options['RSOLVER_FILE'] += '_mhd'
    if args['eos'] == 'isothermal':
        definitions['NWAVE_VALUE'] = '6'
        if args['flux'] == 'hlld':
            makefile_options['RSOLVER_FILE'] += '_iso'
    else:
        definitions['NWAVE_VALUE'] = '7'
else:
    definitions['MAGNETIC_FIELDS_ENABLED'] = '0'
    if definitions['GENERAL_EOS'] != '0':
        makefile_options['GENERAL_EOS_FILE'] += '_hydro'
    else:
        makefile_options['EOS_FILE'] += '_hydro'
    definitions['NFIELD_VARIABLES'] = '0'
    makefile_options['RSOLVER_DIR'] = 'hydro/'
    if args['eos'] == 'isothermal':
        definitions['NWAVE_VALUE'] = '4'
    else:
        definitions['NWAVE_VALUE'] = '5'

# -sts argument
if args['sts']:
    definitions['STS_ENABLED'] = '1'
else:
    definitions['STS_ENABLED'] = '0'

# -s, -g, and -t arguments
definitions['RELATIVISTIC_DYNAMICS'] = '1' if args['s'] or args['g'] else '0'
definitions['GENERAL_RELATIVITY'] = '1' if args['g'] else '0'
if args['s']:
    makefile_options['EOS_FILE'] += '_sr'
    if definitions['GENERAL_EOS'] != '0':
        makefile_options['GENERAL_EOS_FILE'] += '_sr'
    makefile_options['RSOLVER_FILE'] += '_rel'
if args['g']:
    makefile_options['EOS_FILE'] += '_gr'
    if definitions['GENERAL_EOS'] != '0':
        makefile_options['GENERAL_EOS_FILE'] += '_gr'
    makefile_options['RSOLVER_FILE'] += '_rel'
    if not args['t']:
        makefile_options['RSOLVER_FILE'] += '_no_transform'


# -radiation argument
definitions['NRAD_VARIABLES'] = '0'

if args['nr_radiation']:
    definitions['NR_RADIATION_ENABLED'] = '1'
    definitions['NRAD_VARIABLES'] = '14'
else:
    definitions['NR_RADIATION_ENABLED'] = '0'

if args['implicit_radiation']:
    definitions['IM_RADIATION_ENABLED'] = '1'
    definitions['NRAD_VARIABLES'] = '14'
else:
    definitions['IM_RADIATION_ENABLED'] = '0'

# -cr argument
definitions['NCR_VARIABLES'] = '0'
if args['cr']:
    definitions['CR_ENABLED'] = '1'
    definitions['NCR_VARIABLES'] = '4'
else:
    definitions['CR_ENABLED'] = '0'

# -crdiff argument
if args['crdiff']:
    definitions['CRDIFFUSION_ENABLED'] = '1'
else:
    definitions['CRDIFFUSION_ENABLED'] = '0'


# --cxx=[name] argument
if args['cxx'] == 'g++':
    # GCC is C++11 feature-complete since v4.8.1 (2013-05-31)
    definitions['COMPILER_CHOICE'] = 'g++'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'g++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11'
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'g++-simd':
    # GCC version >= 4.9, for OpenMP 4.0; version >= 6.1 for OpenMP 4.5 support
    definitions['COMPILER_CHOICE'] = 'g++-simd'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'g++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
        '-O3 -std=c++11 -fopenmp-simd -fwhole-program -flto -ffast-math '
        '-march=native -fprefetch-loop-arrays'
        # -march=skylake-avx512, skylake, core-avx2
        # -mprefer-vector-width=128  # available in gcc-8, but not gcc-7
        # -mtune=native, generic, broadwell
        # -mprefer-avx128
        # -m64 (default)
    )
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'icpx':
    # New Versions of LLVM-based Intel oneAPI DPC++/C++ Compiler
    definitions['COMPILER_CHOICE'] = 'icpx'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'icpx'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    # ICX drivers icx and icpx will accept ICC Classic Compiler options or Clang*/LLVM
    # Compiler options
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xhost -qopenmp-simd '
      '-Wno-tautological-constant-compare -Wno-array-bounds'
    )
    # Currently unsupported, but "options to be supported" according to icpx
    # -qnextgen-diag: '-inline-forceinline'
    # -qopt-prefetch=4 is supported but it crashes with version 2024.0 (2024.2 is OK)
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'icpx-old':
    # Old versions of LLVM-based Intel oneAPI DPC++/C++ Compiler (backward compatibility)
    definitions['COMPILER_CHOICE'] = 'icpx'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'icpx'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xhost -qopenmp-simd '
      '-Wno-tautological-constant-compare -Wno-array-bounds'
    )
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'icpc':
    # ICC is C++11 feature-complete since v15.0 (2014-08-26)
    definitions['COMPILER_CHOICE'] = 'icpc'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xhost -inline-forceinline -qopenmp-simd -qopt-prefetch=4 '
      '-qoverride-limits '  # -qopt-report-phase=ipo (does nothing without -ipo)
      '-diag-disable=10441'  # The Intel(R) C++ Compiler Classic (ICC) is deprecated
    )
    # -qopt-zmm-usage=high'  # typically harms multi-core performance on Skylake Xeon
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'icpc-debug':
    # Disable IPO, forced inlining, and fast math. Enable vectorization reporting.
    # Useful for testing symmetry, SIMD-enabled functions and loops with OpenMP 4.5
    definitions['COMPILER_CHOICE'] = 'icpc'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -xhost -qopenmp-simd -fp-model precise -qopt-prefetch=4 '
      '-qopt-report=5 -qopt-report-phase=openmp,vec -g -qoverride-limits '
      '-diag-disable=10441'
    )
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'icpc-phi':
    # Cross-compile for Intel Xeon Phi x200 KNL series (unique AVX-512ER and AVX-512FP)
    # -xMIC-AVX512: generate AVX-512F, AVX-512CD, AVX-512ER and AVX-512FP
    definitions['COMPILER_CHOICE'] = 'icpc'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xMIC-AVX512 -inline-forceinline -qopenmp-simd '
      '-qopt-prefetch=4 -qoverride-limits'
    )
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'cray':
    # New HPE Cray Compiling Environment based on clang (2019-)
    # e.g. NAOJ's HPE Cray XD2000 with Sapphire Rapids
    definitions['COMPILER_CHOICE'] = 'cray'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'CC'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11 -flto'  # -Ofast
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = '-lm'
if args['cxx'] == 'clang++':
    # Clang is C++11 feature-complete since v3.3 (2013-06-17)
    definitions['COMPILER_CHOICE'] = 'clang++'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11 -flto'
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'clang++-simd':
    # LLVM/Clang version >= 3.9 for most of OpenMP 4.0 and 4.5 (still incomplete; no
    # offloading, target/declare simd directives). OpenMP 3.1 fully supported in LLVM 3.7
    definitions['COMPILER_CHOICE'] = 'clang++-simd'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11 -flto -fopenmp-simd'
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'clang++-apple':
    # Apple LLVM/Clang: forked version of the open-source LLVM project bundled in macOS
    definitions['COMPILER_CHOICE'] = 'clang++-apple'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11'
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'aocc':
    # AMD Optimizing C++ Compiler based on clang
    definitions['COMPILER_CHOICE'] = 'aocc'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11 -flto -zopt'
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''

# --chemistry=[network] argument
makefile_options['CHEMISTRY_FILE'] = \
    'src/chemistry/network_wrapper.cpp src/chemistry/utils/*.cpp'
if args['chemistry'] is not None:
    definitions['CHEMISTRY_ENABLED'] = '1'
    definitions['CHEMNETWORK_HEADER'] = '../chemistry/network/' \
                                        + args['chemistry'] + '.hpp'
    makefile_options['CHEMNET_FILE'] = 'src/chemistry/network/' \
        + args['chemistry'] + '.cpp'
    # specify the number of species for each network
    if args['chemistry'] == "gow17":
        definitions['NUMBER_CHEMICAL_SPECIES'] = '12'
    elif args['chemistry'] == "H2":
        definitions['NUMBER_CHEMICAL_SPECIES'] = '2'
    elif args['chemistry'] == "G14Sod":
        definitions['NUMBER_CHEMICAL_SPECIES'] = '8'
else:
    definitions['CHEMISTRY_ENABLED'] = '0'
    definitions['NUMBER_CHEMICAL_SPECIES'] = '0'
    makefile_options['CHEMNET_FILE'] = 'src/chemistry/network/none.cpp'
    definitions['CHEMNETWORK_HEADER'] = '../chemistry/network/chem_network.hpp'

# check number of species and scalars
if definitions['NUMBER_PASSIVE_SCALARS'] == '0':
    definitions['NUMBER_PASSIVE_SCALARS'] = definitions['NUMBER_CHEMICAL_SPECIES']
elif int(definitions['NUMBER_PASSIVE_SCALARS']) < int(
                                     definitions['NUMBER_CHEMICAL_SPECIES']):
    raise SystemExit(
      '### CONFIGURE ERROR: number of passive scalars ({:s})'.format(
        definitions['NUMBER_PASSIVE_SCALARS'])
      + ' less than the number of chemical species ({:s})!'.format(
        definitions['NUMBER_CHEMICAL_SPECIES']))

# --kida_rates=[rates] argument
if args['kida_rates'] is not None:
    if args['chemistry'] == "kida":
        makefile_options['CHEMNET_FILE'] += (
            ' src/chemistry/network/kida_network_files/'
            + args['kida_rates']
            + '/kida_'
            + args['kida_rates']
            + '.cpp')

# --chem_ode_solver=[solver] argument
if args['chem_ode_solver'] == 'cvode':
    definitions['CVODE_OPTION'] = 'CVODE'
    makefile_options['LIBRARY_FLAGS'] += ' -lsundials_cvode -lsundials_nvecserial'
else:
    definitions['CVODE_OPTION'] = 'NO_CVODE'
if args['chem_ode_solver'] is not None:
    makefile_options['CHEM_ODE_SOLVER_FILE'] = args['chem_ode_solver']+'.cpp'
else:
    makefile_options['CHEM_ODE_SOLVER_FILE'] = 'forward_euler.cpp'

# --cvode_path=[path] argument
if args['cvode_path'] != '':
    makefile_options['PREPROCESSOR_FLAGS'] += '-I%s/include' % args['cvode_path']
    makefile_options['LINKER_FLAGS'] += '-L%s/lib' % args['cvode_path']
    makefile_options['LINKER_FLAGS'] += " -Wl,-rpath," + '%s/lib' % args['cvode_path']

# --chem_radiation=[chem_radiation] argument
if args['chem_radiation'] is not None:
    definitions['CHEMRADIATION_ENABLED'] = '1'
    makefile_options['CHEMRADIATION_FILE'] = args['chem_radiation']+'.cpp'
    definitions['CHEMRADIATION_INTEGRATOR'] = args['chem_radiation']
else:
    definitions['CHEMRADIATION_ENABLED'] = '0'
    makefile_options['CHEMRADIATION_FILE'] = 'const.cpp'
    definitions['CHEMRADIATION_INTEGRATOR'] = 'none'

# -float argument
if args['float']:
    definitions['SINGLE_PRECISION_ENABLED'] = '1'
else:
    definitions['SINGLE_PRECISION_ENABLED'] = '0'

# -debug argument
if args['debug']:
    definitions['DEBUG_OPTION'] = '1'
    # Completely replace the --cxx= sets of default compiler flags, disable optimization,
    # and emit debug symbols in the compiled binaries
    if (args['cxx'] == 'g++' or args['cxx'] == 'g++-simd'
            or args['cxx'] == 'icpx' or args['cxx'] == 'icpx-old'
            or args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug'
            or args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple' or args['cxx'] == 'cray'
            or args['cxx'] == 'aocc'):
        makefile_options['COMPILER_FLAGS'] = '-O0 -std=c++11 -g'  # -Og
    if args['cxx'] == 'icpc-phi':
        makefile_options['COMPILER_FLAGS'] = '-O0 -std=c++11 -g -xMIC-AVX512'
else:
    definitions['DEBUG_OPTION'] = '0'

# -coverage argument
if args['coverage']:
    definitions['EXCEPTION_HANDLING_OPTION'] = 'DISABLE_EXCEPTIONS'
    # For now, append new compiler flags and don't override --cxx set, but set code to be
    # unoptimized (-O0 instead of -O3) to get useful statement annotations. Should we add
    # '-g -fopenmp-simd', by default? Don't combine lines when writing source code!
    if (args['cxx'] == 'g++' or args['cxx'] == 'g++-simd'):
        makefile_options['COMPILER_FLAGS'] += (
            ' -O0 -fprofile-arcs -ftest-coverage'
            ' -fno-inline -fno-exceptions -fno-elide-constructors'
            )
    if (args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug'
            or args['cxx'] == 'icpc-phi'):
        makefile_options['COMPILER_FLAGS'] += ' -O0 -prof-gen=srcpos'
    if (args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple'):
        # Clang's "source-based" code coverage feature to produces .profraw output
        makefile_options['COMPILER_FLAGS'] += (
            ' -O0 -fprofile-instr-generate -fcoverage-mapping'
            )  # use --coverage to produce GCC-compatible .gcno, .gcda output for gcov
    if (args['cxx'] == 'icpx' or args['cxx'] == 'icpx-old'
            or args['cxx'] == 'cray' or args['cxx'] == 'aocc'):
        raise SystemExit(
            '### CONFIGURE ERROR: No code coverage avaialbe for selected compiler!')
else:
    # Enable C++ try/throw/catch exception handling, by default. Disable only when testing
    # code coverage, since it causes Gcov and other tools to report misleadingly low
    # branch coverage statstics due to untested throwing of exceptions from function calls
    definitions['EXCEPTION_HANDLING_OPTION'] = 'ENABLE_EXCEPTIONS'

# --ccmd=[name] argument
if args['ccmd'] is not None:
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = args['ccmd']

# --gcovcmd=[name] argument (only modifies Makefile target)
if args['gcovcmd'] is not None:
    makefile_options['GCOV_COMMAND'] = args['gcovcmd']
else:
    makefile_options['GCOV_COMMAND'] = 'gcov'

# -mpi argument
if args['mpi']:
    definitions['MPI_OPTION'] = 'MPI_PARALLEL'
    if (args['cxx'] == 'g++' or args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug'
            or args['cxx'] == 'icpx' or args['cxx'] == 'icpx-old'
            or args['cxx'] == 'icpc-phi' or args['cxx'] == 'g++-simd'
            or args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple' or args['cxx'] == 'aocc'):
        definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'mpicxx'
    if args['cxx'] == 'cray':
        definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'CC'
    if args['mpiccmd'] is not None:
        definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = args['mpiccmd']  # noqa
else:
    definitions['MPI_OPTION'] = 'NOT_MPI_PARALLEL'

# -omp argument
if args['omp']:
    definitions['OPENMP_OPTION'] = 'OPENMP_PARALLEL'
    if (args['cxx'] == 'g++' or args['cxx'] == 'g++-simd' or args['cxx'] == 'clang++'
            or args['cxx'] == 'clang++-simd' or args['cxx'] == 'cray'
            or args['cxx'] == 'aocc'):
        makefile_options['COMPILER_FLAGS'] += ' -fopenmp'
    if (args['cxx'] == 'clang++-apple'):
        # Apple Clang disables the front end OpenMP driver interface; enable it via the
        # preprocessor. Must install LLVM's OpenMP runtime library libomp beforehand
        makefile_options['COMPILER_FLAGS'] += ' -Xpreprocessor -fopenmp'
        makefile_options['LIBRARY_FLAGS'] += ' -lomp'
    if (args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug' or args['cxx'] == 'icpc-phi'
            or args['cxx'] == 'icpx' or args['cxx'] == 'icpx-old'):
        makefile_options['COMPILER_FLAGS'] += ' -qopenmp'
else:
    definitions['OPENMP_OPTION'] = 'NOT_OPENMP_PARALLEL'
    if (args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug'
            or args['cxx'] == 'icpc-phi'):
        # suppressed messages:
        #   3180: pragma omp not recognized
        makefile_options['COMPILER_FLAGS'] += ' -diag-disable 3180'

# --grav argument
if args['grav'] == "none":
    definitions['SELF_GRAVITY_ENABLED'] = '0'
else:
    if args['grav'] == "fft":
        definitions['SELF_GRAVITY_ENABLED'] = '1'
        if not args['fft']:
            raise SystemExit(
                '### CONFIGURE ERROR: FFT Poisson solver only be used with FFT')
    if args['grav'] == "mg":
        definitions['SELF_GRAVITY_ENABLED'] = '2'


# -fft argument
makefile_options['MPIFFT_FILE'] = ' '
definitions['FFT_OPTION'] = 'NO_FFT'
if args['fft']:
    definitions['FFT_OPTION'] = 'FFT'
    if args['fftw_path'] != '':
        makefile_options['PREPROCESSOR_FLAGS'] += ' -I{0}/include'.format(
            args['fftw_path'])
        makefile_options['LINKER_FLAGS'] += ' -L{0}/lib'.format(args['fftw_path'])
    if args['omp']:
        makefile_options['LIBRARY_FLAGS'] += ' -lfftw3_omp'
    if args['mpi']:
        makefile_options['MPIFFT_FILE'] = ' $(wildcard src/fft/plimpton/*.cpp)'
    makefile_options['LIBRARY_FLAGS'] += ' -lfftw3'

# -hdf5 argument
if args['hdf5']:
    definitions['HDF5_OPTION'] = 'HDF5OUTPUT'

    if args['hdf5_path'] != '':
        makefile_options['PREPROCESSOR_FLAGS'] += ' -I{0}/include'.format(
            args['hdf5_path'])
        makefile_options['LINKER_FLAGS'] += ' -L{0}/lib'.format(args['hdf5_path'])
    if (args['cxx'] == 'g++' or args['cxx'] == 'g++-simd'
            or args['cxx'] == 'cray' or args['cxx'] == 'icpc'
            or args['cxx'] == 'icpx' or args['cxx'] == 'icpx-old'
            or args['cxx'] == 'icpc-debug' or args['cxx'] == 'icpc-phi'
            or args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple' or args['cxx'] == 'aocc'):
        makefile_options['LIBRARY_FLAGS'] += ' -lhdf5'
else:
    definitions['HDF5_OPTION'] = 'NO_HDF5OUTPUT'

# -h5double argument (does nothing if no -hdf5)
if args['h5double']:
    definitions['H5_DOUBLE_PRECISION_ENABLED'] = '1'
else:
    definitions['H5_DOUBLE_PRECISION_ENABLED'] = '0'

# --cflag=[string] argument
if args['cflag'] is not None:
    makefile_options['COMPILER_FLAGS'] += ' '+args['cflag']

# --include=[name] arguments
for include_path in args['include']:
    makefile_options['COMPILER_FLAGS'] += ' -I'+include_path

# --lib_path=[name] arguments
for library_path in args['lib_path']:
    makefile_options['LINKER_FLAGS'] += ' -L'+library_path

# --lib=[name] arguments
for library_name in args['lib']:
    makefile_options['LIBRARY_FLAGS'] += ' -l'+library_name

# Assemble all flags of any sort given to compiler
definitions['COMPILER_FLAGS'] = ' '.join(
    [makefile_options[opt+'_FLAGS'] for opt in
     ['PREPROCESSOR', 'COMPILER', 'LINKER', 'LIBRARY']])

# --- Step 4. Create new files, finish up --------------------------------

# Terminate all filenames with .cpp extension
makefile_options['PROBLEM_FILE'] += '.cpp'
makefile_options['COORDINATES_FILE'] += '.cpp'
makefile_options['EOS_FILE'] += '.cpp'
makefile_options['GENERAL_EOS_FILE'] += '.cpp'
makefile_options['RSOLVER_FILE'] += '.cpp'

# Read templates
with open(defsfile_input, 'r') as current_file:
    defsfile_template = current_file.read()
with open(makefile_input, 'r') as current_file:
    makefile_template = current_file.read()

# Make substitutions
for key, val in definitions.items():
    defsfile_template = re.sub(r'@{0}@'.format(key), val, defsfile_template)
for key, val in makefile_options.items():
    makefile_template = re.sub(r'@{0}@'.format(key), val, makefile_template)

# Write output files
with open(defsfile_output, 'w') as current_file:
    current_file.write(defsfile_template)
with open(makefile_output, 'w') as current_file:
    current_file.write(makefile_template)

# Finish with diagnostic output
# To match show_config.cpp output: use 2 space indent for option, value string starts on
# column 30
self_grav_string = 'OFF'
if args['grav'] == 'fft':
    self_grav_string = 'FFT'
elif args['grav'] == 'mg':
    self_grav_string = 'Multigrid'


def output_config(opt_descr, opt_choice, filehandle=None):
    first_col_width = 32
    first_col_indent = 2
    descr_len = len(opt_descr)
    right_pad_len = first_col_width - (descr_len + first_col_indent + 2)  # include colon
    right_pad = right_pad_len*' ' if right_pad_len >= 0 else ''
    line_str = first_col_indent*' ' + opt_descr + ': ' + right_pad + opt_choice
    print(line_str)
    if (filehandle is not None):
        filehandle.write(line_str + '\n')


# write the configuration optitions into a log file
flog = open('./configure.log', 'w')

output_config('Your Athena++ distribution has now been configured with the following options', '', flog)  # noqa
output_config('Problem generator', args['prob'], flog)
output_config('Coordinate system', args['coord'], flog)
output_config('Equation of state', args['eos'], flog)
output_config('Riemann solver', args['flux'], flog)
output_config('Magnetic fields', ('ON' if args['b'] else 'OFF'), flog)
output_config('Number of scalars', definitions['NUMBER_PASSIVE_SCALARS'], flog)
output_config('Number of chemical species', definitions['NUMBER_CHEMICAL_SPECIES'], flog)
output_config('Special relativity', ('ON' if args['s'] else 'OFF'), flog)
output_config('General relativity', ('ON' if args['g'] else 'OFF'), flog)
output_config('Radiative Transfer', ('ON' if args['nr_radiation'] else 'OFF'), flog)
output_config('Implicit Radiation', ('ON' if args['implicit_radiation'] else 'OFF'), flog)
output_config('Cosmic Ray Transport', ('ON' if args['cr'] else 'OFF'), flog)
output_config('Cosmic Ray Diffusion', ('ON' if args['crdiff'] else 'OFF'), flog)
output_config('Frame transformations', ('ON' if args['t'] else 'OFF'), flog)
output_config('Self-Gravity', self_grav_string, flog)
output_config('Super-Time-Stepping', ('ON' if args['sts'] else 'OFF'), flog)
output_config('Chemistry', (args['chemistry']
                            if args['chemistry'] is not None else 'OFF'), flog)
output_config('KIDA rates', (args['kida_rates']
                             if args['kida_rates'] is not None else 'OFF'), flog)
output_config('ChemRadiation', (args['chem_radiation']
                                if args['chem_radiation'] is not None else 'OFF'), flog)
output_config('chem_ode_solver', (args['chem_ode_solver'] if args['chem_ode_solver']
                                  is not None else 'OFF'), flog)
output_config('Debug flags', ('ON' if args['debug'] else 'OFF'), flog)
output_config('Code coverage flags', ('ON' if args['coverage'] else 'OFF'), flog)
output_config('Linker flags', makefile_options['LINKER_FLAGS'] + ' '
              + makefile_options['LIBRARY_FLAGS'], flog)
output_config('Floating-point precision', ('single' if args['float'] else 'double'), flog)
output_config('Number of ghost cells', args['nghost'], flog)
output_config('MPI parallelism', ('ON' if args['mpi'] else 'OFF'), flog)
output_config('OpenMP parallelism', ('ON' if args['omp'] else 'OFF'), flog)
output_config('FFT', ('ON' if args['fft'] else 'OFF'), flog)
output_config('HDF5 output', ('ON' if args['hdf5'] else 'OFF'), flog)
if args['hdf5']:
    output_config('HDF5 precision', ('double' if args['h5double'] else 'single'), flog)
output_config('Compiler', args['cxx'], flog)
output_config('Compilation command', makefile_options['COMPILER_COMMAND'] + ' '
              + makefile_options['PREPROCESSOR_FLAGS'] + ' '
              + makefile_options['COMPILER_FLAGS'], flog)

flog.close()
