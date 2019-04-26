"""
Test script for checking that Minkowski coordinates compile.
"""

# Modules
import logging
import scripts.utils.athena as athena
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
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
