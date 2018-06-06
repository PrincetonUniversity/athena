#!/bin/bash
# for macOS builds, install MPICH from Homebrew
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    export HOMEBREW_NO_AUTO_UPDATE=1
    cd mpich
    brew unlink open-mpi || true
    # MPICH and dependencies
    brew install gcc # Homebrew dependency due to gfrotran
    rm '/usr/local/include/c++' || true # recommended by Homebrew
    brew link --overwrite gcc || true
    brew link mpich || true
    # check to see if MPICH executable is cached from previous build
    if [ -f "bin/mpiexec" ]; then
	echo "Using cached MPICH"
    else
        echo "Installing MPICH with Homebrew"
	HOMEBREW_TEMP=$TRAVIS_BUILD_DIR/mpich
	# brew update
	# There is no libtoolize in macOS system. Homebrew uses "g" prefix
	ln -s `which glibtoolize`  /usr/local/opt/libtool/bin/libtoolize
	brew install mpich # --HEAD
        # brew install --build-from-source --cc=clang mpich # --without-fortran # --cc=gcc-7
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
	echo "Downloading mpich source."
	wget http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz
	tar xfz mpich-3.2.1.tar.gz
	rm mpich-3.2.1.tar.gz
	echo "configuring and building mpich."
	cd mpich-3.2.1
	# Disabled fortran to shorten MPICH install time when building from source
	# Need to enable romio for MPI-I/O
	./configure \
            --prefix=`pwd`/../mpich \
            --enable-static=false \
            --enable-alloca=true \
            --disable-long-double \
            --enable-threads=single \
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
    # Recommended by Travis CI documentation to unset these for MPI builds
    # (put in .travis.yml before_install section)
    # test -n $CC && unset CC
    # test -n $CXX && unset CXX
fi
