#!/usr/bin/env bash

# SCRIPT: count_athena.sh
# AUTHOR: Kyle Gerard Felker - kfelker@princeton.edu
# DATE:   4/19/2018
#
# PURPOSE:  Count files in each major subdirectory of the Athena++ repo, and
#           sort each list by the number of lines. Excludes data files and
#           external library source code that is packaged with Athena++
#
# USAGE: ./count_athena.sh
#        Assumes this script is executed from ./tst/style/
# ========================================

cd ../../src
echo "Counting src/ files...."
git ls-files | grep -v "fft/plimpton/" | xargs wc -l | sort -k1 -r

cd ../vis/
echo "Counting vis/ files...."
# Exclude vis/visit/ .xml files from count
git ls-files | grep -vE ".xml" | xargs wc -l | sort -k1 -r

cd ../tst/
echo "Counting tst/ files...."
# Exclude tst/regression/data/ and Google C++ Style linter cpplint.py from count
git ls-files | grep -vE "\.vtk|\.hst|\.tab|cpplint\.py" | xargs wc -l | sort -k1 -r
