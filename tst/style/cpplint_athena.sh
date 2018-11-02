#!/usr/bin/env bash

# SCRIPT: cpplint_athena.sh
# AUTHOR: Kyle Gerard Felker - kfelker@princeton.edu
# DATE:   4/18/2018
# PURPOSE:  Wrapper script to ./cpplint.py application to check Athena++ src/ code
#           compliance with C++ style guildes. User's Python and/or Bash shell
#           implementation may not support recursive globbing of src/ subdirectories
#           and files, so this uses "find" cmd w/ non-POSIX Bash process substitution.
#
# USAGE: ./cpplint_athena.sh
#        Assumes this script is executed from ./tst/style/ with cpplint.py in
#        the same directory, and that CPPLINT.cfg is in root directory.
#        TODO: add explicit check of execution directory
# ========================================

# src/plimpton/ should probably be removed from the src/ folder. Exclude from style checks for now.

# no buffering issue with stdout of "git ls-files"; appears normally in Jenkins log
git ls-files ./
git --version
which git
# is this command's output buffered in Jenkins?-- Appears normally in Travis CI, regardless of linting outcome
git ls-tree -r HEAD ../../src 2>&1
git ls-tree -rz HEAD ../../src 2>&1
# from GNU coreutils 7.5 and later: "stdout -oL cmd" turns on line-buffering for cmd output, -o0 makes it unbuffered
# (no built-in counterpart available on macOS; 'script' is a possible alternative)
stdbuf -o0 git ls-tree -r HEAD ../../src | awk '{print substr($1,4,5), $4}'
stdbuf -o0 git ls-tree -r HEAD ../../src 2>&1 | awk '{print substr($1,4,5), $4}'
# gawk output is unbuffered by default for INTERACTIVE tty
# By default, grep output is line buffered when standard output is a terminal and block buffered otherwise.
# Use --line-buffered to force line-by-line buffering for non-terminal stdout

# Apply Google C++ Style Linter to all source code files at once:
echo "Starting Google C++ Style cpplint.py test"
set -e
# Use "python2 -u" to prevent buffering of sys.stdout,stderr.write() calls in cpplint.py and mix-up in Jenkins logs,
# and since the local version of cpplint.py is currently incompatible with Python 3. Monitor potential errors on macOS from O_NONBLOCK.
find ../../src/ -type f \( -name "*.cpp" -o -name "*.hpp" \) -not -path "*/fft/plimpton/*" -not -name "defs.hpp" -print | xargs python2 -u ./cpplint.py --counting=detailed
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
if [ $? -ne 1 ]; then echo "ERROR: Found C++ file(s) in src/ with trailing whitespace"; exit 1; fi
echo "End of trailing whitespace test"

# Check that all files in src/ have the correct, non-executable octal permission 644
# Git only tracks permission changes (when core.filemode=true) for the "user/owner" executable permissions bit,
# and ignores the user read/write and all "group" and "other", setting file modes to 100644 or 100755 (exec)
echo "Checking for correct file permissions in src/"

# Option 1: stat --- is not portable, e.g. BSD stat on macOS uses -f FORMAT flag
# - Directories have 755 global executable permissions to enable cd
# - No recursion into subdirectories
#stat -c '%a - %n' ../../src/*

# Option 2: find --- is portable, but tracks the local working direcotry NOT the Git working tree
# It will process local temp files which we don't care about.
#find ../../src/ -type f -printf '%m %p\n' | grep -v "644"

# Option 3: git ls-tree --- portable if the copy was cloned with git, only checks working tree
# (but not the staging area / index, unlike git ls-files)
# - First column of output is 6 octal digit UNIX file mode: 2x file type, 1x sticky bits, 3x permissions
# - As Git tree objects, directories have 040000 file mode--- permissions are ignored
# - Recurse into sub-trees to get only blob (file) entries:
git ls-tree -r HEAD ../../src
git ls-tree -r HEAD ../../src | awk '{print substr($1,4,5), $4}'
git ls-tree -r HEAD ../../src | awk '{print substr($1,4,5), $4}' | grep -v "644"
grep_code=$?
echo $grep_code
# | sort -r | tee >(head -n1) | tail -n1
if [ $grep_code -ne 1 ]; then echo "ERROR: Found C++ file(s) in src/ with executable permission"; exit 1; fi
echo "End of file permissions test"
