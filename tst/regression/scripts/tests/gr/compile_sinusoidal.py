"""
Test script for checking that sinusoidal coordinates compile.
"""

# Modules
import scripts.utils.athena as athena

# Prepare Athena++
def prepare():
  athena.configure('gtb',
      prob='gr_blast',
      coord='sinusoidal',
      flux='hlle')
  athena.make()

# Run Athena++
def run():
  pass

# Analyze outputs
def analyze():
  return True
