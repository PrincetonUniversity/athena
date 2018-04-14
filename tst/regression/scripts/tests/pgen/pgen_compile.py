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
    for pgen in pgen_choices:
        if (pgen[0:3] == 'gr_' or pgen == 'default_pgen'):
            print(pgen)
            # athena.configure('g', coord='minkowski', flux='hlle', prob=pgen)
            # athena.make()
        elif (pgen == 'cpaw' or pgen == 'field_loop' or
              pgen == 'orszag_tang' or pgen == 'rotor'):
            athena.configure('b', prob=pgen, **kwargs)
            athena.make()
        elif (pgen == 'hb3' or pgen == 'hgb' or pgen == 'ssheet' or
              pgen == 'strat'):
            athena.configure('b', 'shear', prob=pgen, **kwargs)
            athena.make()
        else:
            athena.configure(prob=pgen, **kwargs)
            athena.make()


# Run Athena++
def run(**kwargs):
    pass


# Analyze outputs
def analyze():
    return True
