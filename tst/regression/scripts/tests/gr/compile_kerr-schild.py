"""
Test script for checking that Kerr-Schild coordinates compile.
"""

# Modules
import logging
import scripts.utils.athena as athena
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
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
