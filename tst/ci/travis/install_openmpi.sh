#!/bin/bash
# for macOS builds, install OpenMPI from Homebrew:
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    export HOMEBREW_NO_AUTO_UPDATE=1
    cd openmpi
    brew unlink mpich || true
    # OpenMPI and dependencies:
    # always install dependencies from pre-compiled bottles before attempting to build-from-source
    brew install gcc
    brew link libevent || true
    brew install libevent || true
    brew link open-mpi || true
    # check to see if OpenMPI is cached from previous build
    if [ -f "bin/mpirun" ]; then
	echo "Using cached OpenMPI"
    else
        echo "Installing OpenMPI with Homebrew"
	HOMEBREW_TEMP=$TRAVIS_BUILD_DIR/openmpi
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
    if [ -f "openmpi/bin/mpirun" ] && [ -f "openmpi-3.0.2/config.log" ]; then
	echo "Using cached OpenMPI"
    else
        # install OpenMPI from source
	echo "Downloading OpenMPI Source"
	# Using 2.1.1 on Homebrew for osx Travis builds
	wget https://download.open-mpi.org/release/open-mpi/v3.0/openmpi-3.0.2.tar.gz
	tar zxf openmpi-3.0.2.tar.gz
	echo "Configuring and building OpenMPI"
	cd openmpi-3.0.2
	# The configure output is not printed to the Travis CI log due to the redirect
	./configure --prefix=$TRAVIS_BUILD_DIR/openmpi &> openmpi.configure # CC=$C_COMPILER CXX=$CXX_COMPILER --without-fortran
	make -j4 &> openmpi.make
	make install &> openmpi.install
	cd ..
    fi
fi
