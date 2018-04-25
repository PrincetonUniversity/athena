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

# src/plimpton/ should probably be removed from the src/ folder. Exclude from style checks for now.

# Apply Google C++ Style Linter to all source code files at once:
find ../../src/ -type f \( -name "*.cpp" -o -name "*.hpp" \) -not -path "*/fft/plimpton/*" -print | xargs ./cpplint.py --counting=detailed

# Ignoring inline comments, check that all sqrt() and cbrt() function calls reside in std::, not global namespace
set -e
echo "Starting std::sqrt(), std::cbrt() test"
find ../../src/ -type f \( -name "*.cpp" -o -name "*.hpp" \) -not -path "*/fft/plimpton/*" -print | while read -r file; do
    [ $(grep -ri "sqrt(" "$file" | grep -v "std::sqrt(" | grep -v "//" | wc -l) -eq 0 ]
    echo "Checking $file...."
    # echo $?  # silent return
    [ $(grep -ri "cbrt(" "$file" | grep -v "std::cbrt(" | grep -v "//" | wc -l) -eq 0 ]
    # ./cpplint.py --counting=detailed "$file" # for linting each src/ file separately
done

echo "End of std::sqrt(), std::cbrt() test"
set +e
