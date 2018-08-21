# Test script to make sure all problem generator files in /src/pgen compile

# Modules
import scripts.utils.athena as athena
import glob
import os

current_dir = os.getcwd()
test = current_dir[0:(len(current_dir) - 14)]
pgen_directory = test + 'src/pgen/'

# set pgen_choices to list of .cpp files in src/pgen/
pgen_choices = glob.glob(pgen_directory + '*.cpp')
# remove 'src/pgen/' prefix and '.cpp' extension from each filename
pgen_choices = set([choice[len(pgen_directory):-4] for choice in pgen_choices])

# Skip --prob=default_pgen; weak attr. may not compile with Intel C++ compiler
pgen_choices.remove('default_pgen')
# Currently testing GR Hydro, not MHD:
# coord='minkowski', #'-t'
gr_probs = set([pgen for pgen in pgen_choices if pgen[0:3] == 'gr_'])
# MHD-required problems: --rsolver=hlld by default
mhd_probs = set(['cpaw', 'field_loop', 'orszag_tang', 'rotor', 'resist', 'magnoh'])
# Shearing box MHD problems
shear_probs = set(['hb3', 'hgb', 'ssheet', 'strat'])
# Newtonian Hydro-only or MHD-optional problems are all leftover pgen/ files
# --rsolver=hllc by default
hydro_probs = pgen_choices - gr_probs - mhd_probs - shear_probs
# Curvilinear problems

# Define configure flags for each set:
gr_args = ['g', 't', '-coord=minkowski']
mhd_args = ['b']
shear_args = ['b', 'shear']
hydro_args = []


# Prepare Athena++
def prepare(**kwargs):
    # Check that code compiles all pgen files in single or double precision
    for single_precision in [True, False]:
        for pgen_set, args in zip([gr_probs, mhd_probs,
                                   shear_probs, hydro_probs],
                                  [gr_args, mhd_args,
                                   shear_args, hydro_args]):
            if single_precision:
                args.extend(['float'])
            # "make clean" and link into executable only for the first problem
            # in the set that shares ./configure.py flags (except --pgen)
            pgen = pgen_set.pop()
            athena.configure(*args, prob=pgen, **kwargs)
            athena.make(clean_first=True, obj_only=False)
            for pgen in pgen_set:
                athena.configure(*args, prob=pgen, **kwargs)
                athena.make(clean_first=False, obj_only=True)


# Run Athena++
def run(**kwargs):
    pass


# Analyze outputs
def analyze():
    return True
