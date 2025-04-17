#!/bin/bash

# HDF5 1.10.4 = released 2018-10-09
# HDF5 1.10.3 = released 2018-08-22
# HDF5 1.10.2 = released 2018-03-30
# HDF5 1.10.1 = released 2018-08-22
version_str=1.10.4
# for macOS builds, install HDF5 from Homebrew
if [ "$OS_NAME" == "osx" ]; then
    export HOMEBREW_NO_AUTO_UPDATE=1
    cd hdf5
    # HDF5 dependencies
    brew link szip || true
    brew install szip || true
    brew link hdf5 || true
    # check to see if HDF5 is cached from previous build
    if [ -f "bin/h5stat" ]; then
	echo "Using cached HDF5"
    else
        echo "Installing HDF5 with Homebrew"
	HOMEBREW_TEMP=$BUILD_DIR/hdf5
	# brew update
        brew install hdf5 # --with-mpi
	brew unlink hdf5
	# /usr/local/opt symlinks to Cellar are preserved, use these:
	# (note, /usr/local/bin/bin -> /usr/local/bin = circular symlink)
	ln -s /usr/local/opt/hdf5/bin bin || true
	ln -s /usr/local/opt/hdf5/lib lib || true
	ln -s /usr/local/opt/hdf5/include include || true
	ln -s /usr/local/opt/hdf5/share share || true
    fi
else
    # for Ubuntu builds, install HDF5 from source
    # check to see if HDF5 is cached from previous build
    if [ -f "hdf5/bin/h5stat" ]; then
	echo "Using cached HDF5"
    else
        # download, configure, and compile HDF5
	echo "Downloading HDF5 Source"
	wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-${version_str}/src/hdf5-${version_str}.tar.gz
	tar zxf hdf5-${version_str}.tar.gz
	cd hdf5-${version_str}
	echo "Configuring and building HDF5"
	./configure --prefix=$BUILD_DIR/hdf5 &> hdf5.configure
	# For parallel MPI build of HDF5 library:
	# CC=/Users/kfelker/mpich-install/bin/mpicc ./configure --enable-parallel
	make -j4 &> hdf5.make
	# make check &> hdf5.makecheck # (expensive tests)
	make install &> hdf5.install
	cd ..
    fi
fi
