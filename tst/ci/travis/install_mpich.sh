#!/bin/bash

# MPICH 3.2.1 = released 2017-11-10
version_str=3.2.1
# for macOS builds, install MPICH from Homebrew
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    export HOMEBREW_NO_AUTO_UPDATE=1
    cd mpich
    brew unlink open-mpi || true
    # MPICH and dependencies
    #brew update > /dev/null  # Travis CI's Homebrew was out of date such that MPICH gcc vs. gfortran deps are broken
    # Depending on macOS image, may need to install Xcode CLI (workaround for "missing stdio.h" header):
    #sudo softwareupdate -i "Command Line Tools (macOS High Sierra version 10.13) for Xcode-9.4"
    # Always install dependencies from pre-compiled bottles before attempting to build-from-source:
    brew install gcc  # Homebrew dependency for MPICH due to gfrotran (would take 1 hr to build from source)
    rm '/usr/local/include/c++' || true # recommended by Homebrew
    brew link --overwrite gcc || true
    brew link mpich || true
    # check to see if MPICH executable is cached from previous build
    if [ -f "bin/mpiexec" ]; then
	echo "Using cached MPICH"
    else
        echo "Installing MPICH with Homebrew"
	HOMEBREW_TEMP=$TRAVIS_BUILD_DIR/mpich
	# There is no libtoolize in macOS system. Homebrew uses "g" prefix
	ln -s `which glibtoolize`  /usr/local/opt/libtool/bin/libtoolize
	brew install mpich # --HEAD
        # brew install --build-from-source --cc=clang mpich # --cc=gcc-7
	# Unlink MPICH to allow for simultaneous installation with OpenMPI
	brew unlink mpich
	# /usr/local/opt symlinks to Cellar are preserved, use these:
	# (note, /usr/local/bin/bin -> /usr/local/bin = circular symlink)
	ln -s /usr/local/opt/mpich/bin bin
	ln -s /usr/local/opt/mpich/include include
	ln -s /usr/local/opt/mpich/lib lib
	ln -s /usr/local/opt/mpich/share share
    fi
else
    # for Ubuntu builds, install MPICH from source
    # check to see if MPICH is cached from previous build
    if [ -f mpich/lib/libmpich.so ]; then
	echo "libmpich.so found -- nothing to build."
    else
	# download, configure, and compile MPICH
	echo "Downloading mpich source."
	wget http://www.mpich.org/static/downloads/${version_str}/mpich-${version_str}.tar.gz
	tar xfz mpich-${version_str}.tar.gz
	rm mpich-${version_str}.tar.gz
	cd mpich-${version_str}
	echo "configuring and building mpich."
	# Disabled Fortran bindings to shorten MPICH install time when building from source
	# Need to enable romio for MPI-I/O
	./configure \
            --prefix=`pwd`/../mpich \
            --enable-static=false \
            --enable-alloca=true \
            --disable-long-double \
            --enable-threads=multiple \
	    --enable-fortran=no \
	    --enable-romio=yes \
	    --enable-fast=all \
            --enable-g=none \
	    --enable-timing=none
	make -j4
	make install
	cd -
	# (Optional) Delete MPICH build directory
	#	rm -rf mpich-3.2.1
    fi
fi
