#!/usr/bin/env bash

# SCRIPT: cpplint_athena.sh
# AUTHOR: Kyle Gerard Felker - kfelker@princeton.edu
# DATE:   4/18/2018
# PURPOSE:  Wrapper script to ./cpplint.py application to check Athena++ src/ code
#           compliance with C++ style guildes. User's Python and/or Bash shell
#           implementation may not support recursive globbing of src/ subdirectories
#           and files, so this uses "find" cmd.
#
# USAGE: ./cpplint_athena.sh
#        Assumes this script is executed from ./tst/style/ with cpplint.py in
#        the same directory, and that CPPLINT.cfg is in root directory
# ========================================

find ../../src/ -type f \( -name "*.cpp" -o -name "*.hpp" \) -not -path "*/fft/plimpton/*" -print | xargs ./cpplint.py --counting=detailed
