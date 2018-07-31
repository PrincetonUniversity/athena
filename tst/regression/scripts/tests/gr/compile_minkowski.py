"""
Test script for checking that Minkowski coordinates compile.
"""

# Modules
import scripts.utils.athena as athena


# Prepare Athena++
def prepare(**kwargs):
    athena.configure('gtb',
                     prob='gr_shock_tube',
                     coord='minkowski',
                     flux='hlle', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    pass


# Analyze outputs
def analyze():
    return True
