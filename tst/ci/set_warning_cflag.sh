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
		     "-Wno-unused-but-set-variable"  # separate case vs. previous warning
		     "-Wno-unused-parameter"
		     "-Wno-unused-but-set-parameter"
		     "-Wno-unknown-pragmas"
		     "-Wno-unused-function"
		     # Add even more warnings:
		     #"-Wconversion"      # also controls sign-,float-conversion
		     # TODO(felker): all floating-point literals need toggleable "real"=float/double
		     # suffix to avoid the -Wfloat-conversion warnings
		    )
    elif [ "$1" == "clang++" ] || [ "$1" == "clang" ]; then
	warn_flags+=("-Wno-unused-private-field"  # Added to Clang ~v3.2 in 2012
		     "-Wno-unused-variable"
		     "-Wno-unused-parameter"
		     "-Wno-unknown-pragmas"
		     "-Wno-unused-function"
		     # Add even more warnings:
		     # "-Wconversion"       # also controls sign-,float-conversion
		     "-Wshorten-64-to-32" # (controlled by above flag; not avail in GCC)
		    )
    elif [ "$1" == "icpc" ] || [ "$1" == "icc" ]; then
	# TODO(felker): Intel Classic C/++ compiler remarks, warnings, and errors are not
	# well-documented.  Intel's -Wall preset is minimalistic compared to
	# GCC's. Elevate default -w1 to -w2 or -w3 and add individual warnings (or just
	# try "-w3 -diag-disable:remark"

	# Does ICC have -pedantic?? will it error out?
	warn_flags+=(#"-w3"
	             "-diag-disable=175" # "error #175: subscript out of range" for IBY during Hydro, e.g.
		     "-Wno-unused-variable"
		     "-Wno-unused-function"
		     # Add even more warnings:
		     # "-Wconversion"
		     "-Wshorten-64-to-32" # illegal narrowing warnings (#2259)
	            )
    elif [ "$1" == "icpx" ] || [ "$1" == "icx" ]; then
	# oneAPI compiler drivers accept ALL Clang compiler options (interpreted directly)
	# plus most Intel C++ Compiler Classic options (translated to Clang equiv)
	# Undocumented options from Intel classic compilers are not implemented
	warn_flags+=("-Wno-unused-private-field"
		     "-Wno-unused-variable"
		     "-Wno-unused-parameter"
		     "-Wunknown-pragmas"
		     #"-Wno-unknown-pragmas"  # only a subset of ICC pragmas supported in ICX
		     "-Wno-unused-function"
		     # Add even more warnings:
		     # "-Wconversion"
		     "-Wshorten-64-to-32"
		     "-Wno-array-bounds"  # false positives for isothermal problems with wli[4] in hlle.cpp
		     "-Wno-pass-failed"  # #pragma unroll (in mesh_refinement.cpp in ICX 2021.1.2), #pragma distribute_point (HLLC outside loop only)
		     "-fp-model=precise"  # Or use -Wno-tautological-constant-compare
	             # std::isnan() calls in default -fp-model=fast throw errors with -Wall,-Wextra,-Werror
		    )
    else
	echo "Unknown CXX=$1 compiler"
	return 1
    fi
    echo ${warn_flags[*]}
    return 0
}
set_warning_cflag $1
