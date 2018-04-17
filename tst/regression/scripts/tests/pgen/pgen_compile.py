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
pgen_choices = [choice[len(pgen_directory):-4] for choice in pgen_choices]


# Prepare Athena++
def prepare(**kwargs):
    # Check that code compiles all pgen files in single or double precision
    for single_precision in [True, False]:
        for pgen in pgen_choices:
            args = []
            if single_precision:
                args = ['float']
            # Skip GR and default_pgen problems by default
            # (default_pgen may not compile with Intel C++ compiler)
            if (pgen[0:3] == 'gr_' or pgen == 'default_pgen'):
                print(pgen)
                args.extend(['g'])
                # athena.configure(args, coord='minkowski',
                #                  flux='hlle', prob=pgen)
                # athena.make()
            # MHD-required problems
            elif (pgen == 'cpaw' or pgen == 'field_loop' or
                  pgen == 'orszag_tang' or pgen == 'rotor'):
                args.extend(['b'])
                athena.configure(*args, prob=pgen, **kwargs)
                athena.make()
            # Shearinbox MHD problems
            elif (pgen == 'hb3' or pgen == 'hgb' or pgen == 'ssheet' or
                  pgen == 'strat'):
                args.extend(['b', 'shear'])
                athena.configure(*args, prob=pgen, **kwargs)
                athena.make()
            # Hydro-only and MHD-optional problems
            else:
                athena.configure(*args, prob=pgen, **kwargs)
                athena.make()


# Run Athena++
def run(**kwargs):
    pass


# Analyze outputs
def analyze():
    return True
