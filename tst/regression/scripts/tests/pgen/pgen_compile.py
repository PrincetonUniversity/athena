# Test script to make sure all problem generator files in /src/pgen compile

# Modules
import logging
import scripts.utils.athena as athena
import glob
import os

current_dir = os.getcwd()
test = current_dir[0:(len(current_dir) - 14)]
pgen_directory = test + 'src/pgen/'
logger = logging.getLogger('athena' + __name__[7:])

# set pgen_choices to list of .cpp files in src/pgen/
pgen_choices = glob.glob(pgen_directory + '*.cpp')
# remove 'src/pgen/' prefix and '.cpp' extension from each filename
pgen_choices = set([choice[len(pgen_directory):-4] for choice in pgen_choices])

# Skip --prob=default_pgen; weak attr. may not compile with Intel C++ compiler
pgen_choices.remove('default_pgen')
# Skip from_array pgen
pgen_choices.remove('from_array')
# Currently testing GR Hydro, not MHD:
# coord='minkowski', #'-t'
gr_probs = set([pgen for pgen in pgen_choices if pgen[0:3] == 'gr_'])
# MHD-required problems: --rsolver=hlld by default
mhd_probs = set(['hb3', 'hgb', 'jgg', 'cpaw', 'field_loop',
                 'orszag_tang', 'rotor', 'resist', 'magnoh'])
grav_mg_probs = set(['binary_gravity', 'collapse'])
crdiff_probs = set(['cr_diffusion_mg'])
# Chemistry problems
chem_probs = set([pgen for pgen in pgen_choices if pgen[0:5] == 'chem_'])
# Newtonian Hydro-only or MHD-optional problems are all leftover pgen/ files
# --rsolver=hllc by default
hydro_probs = pgen_choices - gr_probs - mhd_probs - grav_mg_probs \
            - chem_probs - crdiff_probs
# Curvilinear problems
# Define configure flags for each set:
gr_args = ['g', 't', '-coord=minkowski']
mhd_args = ['b']
grav_mg_args = ['-grav=mg']
chem_args = ['-chemistry=H2', '-chem_ode_solver=forward_euler']
hydro_args = []
crdiff_args = ['b', 'crdiff']


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    # Check that code compiles all pgen files in single or double precision
    for single_precision in [True, False]:
        for pgen_set, args in zip([gr_probs, mhd_probs, grav_mg_probs,
                                   chem_probs, crdiff_probs, hydro_probs],
                                  [gr_args, mhd_args, grav_mg_args,
                                   chem_args, crdiff_args, hydro_args]):
            args_lcl = list(args)
            if single_precision:
                args_lcl.extend(['float'])
            # "make clean" and link into executable only for the first problem
            # in the set that shares ./configure.py flags (except --pgen)
            first = True
            for pgen in pgen_set:
                if first:
                    athena.configure(*args_lcl, prob=pgen, **kwargs)
                    athena.make(clean_first=True, obj_only=False)
                    first = False
                else:
                    athena.configure(*args_lcl, prob=pgen, **kwargs)
                    athena.make(clean_first=False, obj_only=True)


# Run Athena++
def run(**kwargs):
    pass


# Analyze outputs
def analyze():
    return True
