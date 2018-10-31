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
echo "Starting Google C++ Style cpplint.py test"
set -e
find ../../src/ -type f \( -name "*.cpp" -o -name "*.hpp" \) -not -path "*/fft/plimpton/*" -not -name "defs.hpp" -print | xargs ./cpplint.py --counting=detailed
set +e
echo "End of Google C++ Style cpplint.py test"

# Ignoring inline comments, check that all sqrt() and cbrt() function calls reside in std::, not global namespace
echo "Starting std::sqrt(), std::cbrt(), \t test"
while read -r file
do
    echo "Checking $file...."
    # sed -n '/\t/p' $file
    grep -n "$(printf '\t')" $file
    if [ $? -ne 1 ]; then echo "ERROR: Do not use \t tab characters"; exit 1; fi

    grep -ri "sqrt(" "$file" | grep -v "std::sqrt(" | grep -v "//"
    if [ $? -ne 1 ]; then echo "ERROR: Use std::sqrt(), not sqrt()"; exit 1; fi

    grep -ri "cbrt(" "$file" | grep -v "std::cbrt(" | grep -v "//"
    if [ $? -ne 1 ]; then echo "ERROR: Use std::cbrt(), not cbrt()"; exit 1; fi
    # To lint each src/ file separately, use:
    # ./cpplint.py --counting=detailed "$file"
done < <(find ../../src/ -type f \( -name "*.cpp" -o -name "*.hpp" \) -not -path "*/fft/plimpton/*" -print)

echo "End of std::sqrt(), std::cbrt(), \t test"

# Search src/ C++ source code for trailing whitespace errors
# (Google C++ Style Linter does not check for this, but flake8 via pycodestyle warning W291 will check *.py)
echo "Checking for trailing whitespace in src/"
find ../../src/ -type f \( -name "*.cpp" -o -name "*.hpp*" \) -not -path "*/fft/plimpton/*" -exec grep -n -E " +$" {} +
if [ $? -ne 1 ]; then echo "ERROR: Found C++ files with trailing whitespace"; exit 1; fi
echo "End of trailing whitespace test"
