"""
Test script for checking that Schwarzschild coordinates compile.
"""

# Modules
import logging
import scripts.utils.athena as athena


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure('gtb',
                     prob='gr_bondi',
                     coord='schwarzschild',
                     flux='hlle', **kwargs)
    athena.make()


# Run Athena++
def run(**kwargs):
    pass


# Analyze outputs
def analyze():
    return True
