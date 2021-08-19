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
#   --cxx=xxx         use xxx as the C++ compiler (works w/ or w/o -mpi)
#   --ccmd=name       use name as the command to call the (non-MPI) C++ compiler
#   --mpiccmd=name    use name as the command to call the MPI C++ compiler
#   --gcovcmd=name    use name as the command to call the gcov utility
#   --cflag=string    append string whenever invoking compiler/linker
#   --include=path    use -Ipath when compiling
#   --lib_path=path   use -Lpath when linking
#   --lib=xxx         use -lxxx when linking
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
        'sinusoidal',
        'tilted',
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
                    choices=['none', 'fft'],
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
    'icpc',
    'icpc-debug',
    'icpc-phi',
    'cray',
    'bgxlc++',
    'clang++',
    'clang++-simd',
    'clang++-apple',
]


def c_to_cpp(arg):
    arg = arg.replace('gcc', 'g++', 1)
    arg = arg.replace('icc', 'icpc', 1)
    arg = arg.replace('icx', 'icpx', 1)
    if arg == 'bgxl' or arg == 'bgxlc':
        arg = 'bgxlc++'

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
definitions['FRAME_TRANSFORMATIONS'] = '1' if args['t'] else '0'
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
    # Next-gen LLVM-based Intel oneAPI DPC++/C++ Compiler
    definitions['COMPILER_CHOICE'] = 'icpx'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'icpx'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    # ICX drivers icx and icpx will accept ICC Classic Compiler options or Clang*/LLVM
    # Compiler options
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xhost -qopenmp-simd '
    )
    # Currently unsupported, but "options to be supported" according to icpx
    # -qnextgen-diag: '-inline-forceinline -qopt-prefetch=4 '
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'icpc':
    # ICC is C++11 feature-complete since v15.0 (2014-08-26)
    definitions['COMPILER_CHOICE'] = 'icpc'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'icpc'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -std=c++11 -ipo -xhost -inline-forceinline -qopenmp-simd -qopt-prefetch=4 '
      '-qoverride-limits'  # -qopt-report-phase=ipo (does nothing without -ipo)
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
      '-qopt-report=5 -qopt-report-phase=openmp,vec -g -qoverride-limits'
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
    # Cray Compiling Environment 8.4 (2015-09-24) introduces C++11 feature completeness
    # (except "alignas"). v8.6 is C++14 feature-complete
    definitions['COMPILER_CHOICE'] = 'cray'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'CC'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -h std=c++11 -h aggress -h vector3 -hfp3'
    makefile_options['LINKER_FLAGS'] = '-hwp -hpl=obj/lib'
    makefile_options['LIBRARY_FLAGS'] = '-lm'
if args['cxx'] == 'bgxlc++':
    # IBM XL C/C++ for BG/Q is NOT C++11 feature-complete as of v12.1.0.15 (2017-12-22)
    # suppressed messages:
    #   1500-036:  The NOSTRICT option has the potential to alter the program's semantics
    #   1540-1401: An unknown "pragma simd" is specified
    #   1586-083:  ld option ignored by IPA
    #   1586-233:  Duplicate definition of symbol ignored
    #   1586-267:  Inlining of specified subprogram failed due to the presence of a C++
    #                exception handler
    definitions['COMPILER_CHOICE'] = 'bgxlc++'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'bgxlc++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = (
      '-O3 -qhot=level=1:vector -qinline=level=5:auto -qipa=level=1:noobject'
      ' -qstrict=subnormals -qmaxmem=150000 -qlanglvl=extended0x -qsuppress=1500-036'
      ' -qsuppress=1540-1401 -qsuppress=1586-083 -qsuppress=1586-233'
      ' -qsuppress=1586-267'
    )
    makefile_options['LINKER_FLAGS'] = makefile_options['COMPILER_FLAGS']
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'clang++':
    # Clang is C++11 feature-complete since v3.3 (2013-06-17)
    definitions['COMPILER_CHOICE'] = 'clang++'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11'
    makefile_options['LINKER_FLAGS'] = ''
    makefile_options['LIBRARY_FLAGS'] = ''
if args['cxx'] == 'clang++-simd':
    # LLVM/Clang version >= 3.9 for most of OpenMP 4.0 and 4.5 (still incomplete; no
    # offloading, target/declare simd directives). OpenMP 3.1 fully supported in LLVM 3.7
    definitions['COMPILER_CHOICE'] = 'clang++-simd'
    definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'clang++'
    makefile_options['PREPROCESSOR_FLAGS'] = ''
    makefile_options['COMPILER_FLAGS'] = '-O3 -std=c++11 -fopenmp-simd'
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

# -float argument
if args['float']:
    definitions['SINGLE_PRECISION_ENABLED'] = '1'
else:
    definitions['SINGLE_PRECISION_ENABLED'] = '0'

# -debug argument
if args['debug']:
    definitions['DEBUG_OPTION'] = 'DEBUG'
    # Completely replace the --cxx= sets of default compiler flags, disable optimization,
    # and emit debug symbols in the compiled binaries
    if (args['cxx'] == 'g++' or args['cxx'] == 'g++-simd'
            or args['cxx'] == 'icpx'
            or args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug'
            or args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple'):
        makefile_options['COMPILER_FLAGS'] = '-O0 --std=c++11 -g'  # -Og
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] = '-O0 -h std=c++11'
    if args['cxx'] == 'bgxlc++':
        makefile_options['COMPILER_FLAGS'] = '-O0 -g -qlanglvl=extended0x'
    if args['cxx'] == 'icpc-phi':
        makefile_options['COMPILER_FLAGS'] = '-O0 --std=c++11 -g -xMIC-AVX512'
else:
    definitions['DEBUG_OPTION'] = 'NOT_DEBUG'

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
            or args['cxx'] == 'icpx'
            or args['cxx'] == 'icpc-phi'):
        makefile_options['COMPILER_FLAGS'] += ' -O0 -prof-gen=srcpos'
    if (args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple'):
        # Clang's "source-based" code coverage feature to produces .profraw output
        makefile_options['COMPILER_FLAGS'] += (
            ' -O0 -fprofile-instr-generate -fcoverage-mapping'
            )  # use --coverage to produce GCC-compatible .gcno, .gcda output for gcov
    if (args['cxx'] == 'cray' or args['cxx'] == 'bgxlc++'):
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
            or args['cxx'] == 'icpx'
            or args['cxx'] == 'icpc-phi' or args['cxx'] == 'g++-simd'
            or args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple'):
        definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'mpicxx'
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -h mpi1'
    if args['cxx'] == 'bgxlc++':
        definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = 'mpixlcxx'  # noqa
    # --mpiccmd=[name] argument
    if args['mpiccmd'] is not None:
        definitions['COMPILER_COMMAND'] = makefile_options['COMPILER_COMMAND'] = args['mpiccmd']  # noqa
else:
    definitions['MPI_OPTION'] = 'NOT_MPI_PARALLEL'

# -omp argument
if args['omp']:
    definitions['OPENMP_OPTION'] = 'OPENMP_PARALLEL'
    if (args['cxx'] == 'g++' or args['cxx'] == 'g++-simd' or args['cxx'] == 'clang++'
            or args['cxx'] == 'clang++-simd'):
        makefile_options['COMPILER_FLAGS'] += ' -fopenmp'
    if (args['cxx'] == 'clang++-apple'):
        # Apple Clang disables the front end OpenMP driver interface; enable it via the
        # preprocessor. Must install LLVM's OpenMP runtime library libomp beforehand
        makefile_options['COMPILER_FLAGS'] += ' -Xpreprocessor -fopenmp'
        makefile_options['LIBRARY_FLAGS'] += ' -lomp'
    if (args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug' or args['cxx'] == 'icpc-phi'
            or args['cxx'] == 'icpx'):
        makefile_options['COMPILER_FLAGS'] += ' -qopenmp'
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -homp'
    if args['cxx'] == 'bgxlc++':
        # use thread-safe version of compiler
        definitions['COMPILER_COMMAND'] += '_r'
        makefile_options['COMPILER_COMMAND'] += '_r'
        makefile_options['COMPILER_FLAGS'] += ' -qsmp'
else:
    definitions['OPENMP_OPTION'] = 'NOT_OPENMP_PARALLEL'
    if args['cxx'] == 'cray':
        makefile_options['COMPILER_FLAGS'] += ' -hnoomp'
    if (args['cxx'] == 'icpc' or args['cxx'] == 'icpc-debug' or args['cxx'] == 'icpc-phi'
            or args['cxx'] == 'icpx'):
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
            or args['cxx'] == 'icpx'
            or args['cxx'] == 'icpc-debug' or args['cxx'] == 'icpc-phi'
            or args['cxx'] == 'clang++' or args['cxx'] == 'clang++-simd'
            or args['cxx'] == 'clang++-apple'):
        makefile_options['LIBRARY_FLAGS'] += ' -lhdf5'
    if args['cxx'] == 'bgxlc++':
        makefile_options['PREPROCESSOR_FLAGS'] += (
            ' -D_LARGEFILE_SOURCE -D_LARGEFILE64_SOURCE -D_BSD_SOURCE'
            ' -I/soft/libraries/hdf5/1.10.0/cnk-xl/current/include'
            ' -I/bgsys/drivers/ppcfloor/comm/include')
        makefile_options['LINKER_FLAGS'] += (
            ' -L/soft/libraries/hdf5/1.10.0/cnk-xl/current/lib'
            ' -L/soft/libraries/alcf/current/xl/ZLIB/lib')
        makefile_options['LIBRARY_FLAGS'] += ' -lhdf5 -lz -lm'
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

print('Your Athena++ distribution has now been configured with the following options:')
print('  Problem generator:          ' + args['prob'])
print('  Coordinate system:          ' + args['coord'])
print('  Equation of state:          ' + args['eos'])
print('  Riemann solver:             ' + args['flux'])
print('  Magnetic fields:            ' + ('ON' if args['b'] else 'OFF'))
print('  Number of scalars:          ' + args['nscalars'])
print('  Special relativity:         ' + ('ON' if args['s'] else 'OFF'))
print('  General relativity:         ' + ('ON' if args['g'] else 'OFF'))
print('  Frame transformations:      ' + ('ON' if args['t'] else 'OFF'))
print('  Self-Gravity:               ' + self_grav_string)
print('  Super-Time-Stepping:        ' + ('ON' if args['sts'] else 'OFF'))
print('  Debug flags:                ' + ('ON' if args['debug'] else 'OFF'))
print('  Code coverage flags:        ' + ('ON' if args['coverage'] else 'OFF'))
print('  Linker flags:               ' + makefile_options['LINKER_FLAGS'] + ' '
      + makefile_options['LIBRARY_FLAGS'])
print('  Floating-point precision:   ' + ('single' if args['float'] else 'double'))
print('  Number of ghost cells:      ' + args['nghost'])
print('  MPI parallelism:            ' + ('ON' if args['mpi'] else 'OFF'))
print('  OpenMP parallelism:         ' + ('ON' if args['omp'] else 'OFF'))
print('  FFT:                        ' + ('ON' if args['fft'] else 'OFF'))
print('  HDF5 output:                ' + ('ON' if args['hdf5'] else 'OFF'))
if args['hdf5']:
    print('  HDF5 precision:             ' + ('double' if args['h5double'] else 'single'))
print('  Compiler:                   ' + args['cxx'])
print('  Compilation command:        ' + makefile_options['COMPILER_COMMAND'] + ' '
      + makefile_options['PREPROCESSOR_FLAGS'] + ' ' + makefile_options['COMPILER_FLAGS'])
