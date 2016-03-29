"""
Test script for checking that Schwarzschild coordinates compile.
"""

# Modules
import scripts.utils.athena as athena

# Prepare Athena++
def prepare():
  athena.configure('gtb',
      prob='bh_spherical_accretion',
      coord='schwarzschild',
      flux='hlle')
  athena.make()

# Run Athena++
def run():
  pass

# Analyze outputs
def analyze():
  return True
