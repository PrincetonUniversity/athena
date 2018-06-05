"""
Test script for checking that Kerr-Schild coordinates compile.
"""

# Modules
import scripts.utils.athena as athena


# Prepare Athena++
def prepare(**kwargs):
    athena.configure('gtb',
                     prob='gr_torus',
                     coord='kerr-schild',
                     flux='hlle', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    pass


# Analyze outputs
def analyze():
    return True
