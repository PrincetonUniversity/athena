#!/usr/bin/env bash

# Define list of compiler warning flags for Jenkins and Travis CI regression tests
set_warning_cflag () {
    # Flags common to all compilers
    warn_flags=("-Wall"
		"-Wextra"
		"-Werror")
    # Suppress or add warnings based on C++ compiler
    if [ "$1" = "g++" ]; then
	warn_flags+=("-Wno-unused-private-field"
		     "-Wno-maybe-uninitialized"
		     "-Wno-address"
		     "-Wno-unused-but-set-variable"
		     "-Wno-comment"
		     "-Wno-unused-variable"
		     "-Wno-unused-parameter"
		     "-Wno-unknown-pragmas"
		     "-Wno-unused-function")
    elif [ "$1" == "clang++" ]; then
	warn_flags+=("-Wno-unused-private-field"
		     "-Wno-address"
		     "-Wno-unused-variable"
		     "-Wno-unused-parameter"
		     "-Wno-unknown-pragmas"
		     "-Wno-unused-function"
		     "-Wshorten-64-to-32")
    elif [ "$1" == "icc" ]; then
	warn_flags+=("-diag-disable=175"
		     "-Wno-unused-variable"
		     "-Wno-unused-function")
	# There are still illegal narrowing warnings with icc to be fixed
		     #"-Wshorten-64-to-32")
    else
	echo "Unknown CXX=$1 compiler"
	return 1
    fi
    echo ${warn_flags[*]}
    return 0
}
set_warning_cflag $1
