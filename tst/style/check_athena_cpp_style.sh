#!/usr/bin/env bash

# SCRIPT: check_athena_cpp_style.sh
# AUTHOR: Kyle Gerard Felker - kfelker@princeton.edu
# DATE:   4/18/2018
# PURPOSE:  Wrapper script to ./cpplint.py application to check Athena++ src/ code
#           compliance with C++ style guildes. User's Python and/or Bash shell
#           implementation may not support recursive globbing of src/ subdirectories
#           and files, so this uses "find" cmd w/ non-POSIX Bash process substitution.
#
# USAGE: ./check_athena_cpp_style.sh
#        Assumes this script is executed from ./tst/style/ with cpplint.py in
#        the same directory, and that CPPLINT.cfg is in root directory.
#        TODO: add explicit check of execution directory
# ========================================

# src/plimpton/ should probably be removed from the src/ folder. Exclude from style checks for now.

# Apply Google C++ Style Linter to all source code files at once:
echo "Starting Google C++ Style cpplint.py test"
set -e
# Use "python[23] -u" to prevent buffering of sys.stdout,stderr.write() calls in cpplint.py and mix-up in Jenkins logs,
find ../../src/ -type f \( -name "*.cpp" -o -name "*.hpp" \) -not -path "*/fft/plimpton/*" -not -name "defs.hpp" -print | xargs python -u ./cpplint.py --counting=detailed
set +e
echo "End of Google C++ Style cpplint.py test"

# Begin custom Athena++ style rules and checks:
echo "Starting std::sqrt(), std::cbrt(), \t test"
while read -r file
do
    echo "Checking $file...."
    # sed -n '/\t/p' $file

    # TYPE 1: may cause bugs, or introduces abhorrent style (e.g. mixing tabs and spaces).
    # --------------------------
    grep -n "$(printf '\t')" $file
    if [ $? -ne 1 ]; then echo "ERROR: Do not use \t tab characters"; exit 1; fi

    # TYPE 2: strict ISO C++11 compilance and/or technical edge-cases.
    # Code would be fine for >95% of environments and libraries with these violations, but they may affect portability.
    # --------------------------

    # Ignoring inline comments, check that all sqrt() and cbrt() function calls reside in std::, not global namespace
    # Note, currently all such chained grep calls will miss violations if a comment is at the end of line, e.g.:
    #     }}  // this is a comment after a style error
    grep -nri "sqrt(" "$file" | grep -v "std::sqrt(" | grep -v "//"
    if [ $? -ne 1 ]; then echo "ERROR: Use std::sqrt(), not sqrt()"; exit 1; fi

    grep -nri "cbrt(" "$file" | grep -v "std::cbrt(" | grep -v "//"
    if [ $? -ne 1 ]; then echo "ERROR: Use std::cbrt(), not cbrt()"; exit 1; fi

    # TYPE 3: purely stylistic inconsistencies.
    # These errors would not cause any changes to code behavior if they were ignored, but they may affect readability.
    # --------------------------
    grep -nri "}}" "$file" | grep -v "//"
    if [ $? -ne 1 ]; then echo "ERROR: Use single closing brace '}}' per line"; exit 1; fi

    # GNU Grep Extended Regex (ERE) syntax:
    grep -nrEi '^\s+#pragma' "$file"
    if [ $? -ne 1 ]; then echo "ERROR: Left justify any #pragma statements"; exit 1; fi

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
# (but not the staging area/index, unlike git ls-files)
# - First column of output is 6 octal digit UNIX file mode: 2x file type, 1x sticky bits, 3x permissions
# - As Git tree objects, directories have 040000 file mode--- permissions are ignored
# - Recurse into sub-trees to get only blob (file) entries:

# As of 2018-11-02, Jenkins worker Perseus has system /usr/bin/git version 1.8.3.1 (2013-06-10), which does not
# have the fix from 1cf9952d (2014-11-30) first released in 2.3.0 that addressed duplicate path filtering in
# builtin/ls-tree.c. (on top of "pathspec" filtering in read_tree_recursive()), which broke the ability to
# specify "../" relative parent dirs. Use --full-tree to ignore current working directory when pattern matching

# Furthermore, even the latest version of "git ls-tree" will silently fail & return nothing if [<path>...]
# is in repository but does not match any tree contents.
git ls-tree -r --full-tree HEAD src/ | awk '{print substr($1,4,5), $4}' | grep -v "644"
if [ $? -ne 1 ]; then echo "ERROR: Found C++ file(s) in src/ with executable permission"; exit 1; fi
echo "End of file permissions test"
