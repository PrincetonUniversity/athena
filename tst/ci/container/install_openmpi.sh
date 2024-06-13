#!/bin/bash

# OpenMPI 3.0.2 = released 2018-06-01
# will need to manually adjust below URL after major version 3.0
version_str=3.0.2
# for macOS builds, install OpenMPI from Homebrew:
if [ "$OS_NAME" == "osx" ]; then
    export HOMEBREW_NO_AUTO_UPDATE=1
    cd openmpi
    brew unlink mpich || true
    # OpenMPI and dependencies:
    # Depending on macOS image, may need to install Xcode CLI (workaround for "missing stdio.h" header):
    #sudo softwareupdate -i "Command Line Tools (macOS High Sierra version 10.13) for Xcode-9.4"
    # Always install dependencies from pre-compiled bottles before attempting to build-from-source:
    brew install gcc
    brew link libevent || true
    brew install libevent || true
    brew link hwloc || true
    brew install hwloc || true
    brew link open-mpi || true
    # check to see if OpenMPI is cached from previous build
    if [ -f "bin/mpirun" ]; then
	echo "Using cached OpenMPI"
    else
        echo "Installing OpenMPI with Homebrew"
	HOMEBREW_TEMP=$BUILD_DIR/openmpi
	# brew update
        brew install open-mpi #--cc=gcc-7 --build-from-source --without-fortran
	brew unlink open-mpi
	# /usr/local/opt symlinks to Cellar are preserved, use these:
	# (note, /usr/local/bin/bin -> /usr/local/bin = circular symlink)
	ln -s /usr/local/opt/open-mpi/bin bin
	ln -s /usr/local/opt/open-mpi/include include
	ln -s /usr/local/opt/open-mpi/lib lib
	ln -s /usr/local/opt/open-mpi/share share
    fi
else
    # for Ubuntu builds, install OpenMPI from source
    # check to see if OpenMPI is cached from previous build
    if [ -f "openmpi/bin/mpirun" ] && [ -f "openmpi-${version_str}/config.log" ]; then
	echo "Using cached OpenMPI"
    else
        # download, configure, and compile OpenMPI
	echo "Downloading OpenMPI Source"
	# Using 2.1.1 on Homebrew for osx Travis builds
	wget https://download.open-mpi.org/release/open-mpi/v3.0/openmpi-${version_str}.tar.gz
	tar zxf openmpi-${version_str}.tar.gz
	cd openmpi-${version_str}
	echo "Configuring and building OpenMPI"
	# The configure output is not printed to the Travis CI log due to the redirect
	./configure --prefix=$BUILD_DIR/openmpi &> openmpi.configure # CC=$C_COMPILER CXX=$CXX_COMPILER --without-fortran
	make -j4 &> openmpi.make
	make install &> openmpi.install
	cd ..
    fi
fi
