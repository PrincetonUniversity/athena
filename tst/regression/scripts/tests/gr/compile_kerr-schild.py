"""
Test script for checking that Kerr-Schild coordinates compile.
"""

# Modules
import scripts.utils.athena as athena

# Prepare Athena++
def prepare():
  athena.configure('gtb',
      prob='accretion_gr_fm_torus',
      coord='kerr-schild',
      flux='hlle')
  athena.make()

# Run Athena++
def run():
  pass

# Analyze outputs
def analyze():
  return True
