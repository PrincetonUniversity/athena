#!/usr/bin/env bash

# Define list of compiler warning flags for Jenkins and Travis CI regression tests
set_warning_cflag () {
    # Flags common to all compilers
    warn_flags=("-Wall"
		"-Wextra"
		"-Werror"
	        "-pedantic"
	        "-pedantic-errors"
	       )
    # Suppress or add warnings based on C++ compiler
    if [ "$1" == "g++" ] || [ "$1" == "gcc" ]; then
	warn_flags+=("-Wno-unused-variable"
		     "-Wno-unused-but-set-variable"
		     "-Wno-unused-parameter"
		     "-Wno-unused-but-set-parameter"
		     "-Wno-unknown-pragmas"
		     "-Wno-unused-function")
    elif [ "$1" == "clang++" ] || [ "$1" == "clang" ]; then
	warn_flags+=("-Wno-unused-private-field"  # Added to Clang ~v3.2 in 2012
		     "-Wno-unused-variable"
		     "-Wno-unused-parameter"
		     "-Wno-unknown-pragmas"
		     "-Wno-unused-function"
		     # Add even more warnings:
		     "-Wshorten-64-to-32")
    elif [ "$1" == "icpc" ] || [ "$1" == "icc" ]; then
	# TODO(felker): Intel compiler remarks, warnings, and errors are not well-documented.
	# Intel's -Wall preset is minimalistic compared to GCC's. Elevate default -w1 to -w2 or -w3
	# and add individual warnings (or just try "-w3 -diag-disable:remark"

	# Does ICC have -pedantic?? will it error out?
		     "-Wno-unused-variable"
		     "-Wno-unused-function")
	# There are still illegal narrowing warnings with icc that need to be addressed
		     #"-Wshorten-64-to-32")
    else
	echo "Unknown CXX=$1 compiler"
	return 1
    fi
    echo ${warn_flags[*]}
    return 0
}
set_warning_cflag $1
