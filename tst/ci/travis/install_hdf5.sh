#!/bin/bash
# for macOS builds, install HDF5 from Homebrew
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
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
	HOMEBREW_TEMP=$TRAVIS_BUILD_DIR/hdf5
	brew update
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
        # install HDF5 from source
	echo "Downloading HDF5 Source"
	# Using 1.10.1 on Homebrew for osx Travis builds
	wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.1.tar.gz
	tar zxf hdf5-1.10.1.tar.gz
	echo "Configuring and building HDF5"
	cd hdf5-1.10.1
	./configure --prefix=$TRAVIS_BUILD_DIR/hdf5 CC=$C_COMPILER CXX=$CXX_COMPILER &> hdf5.configure
	# CC=/Users/kfelker/mpich-install/bin/mpicc ./configure --enable-parallel
	make -j4 &> hdf5.make
	#make check &> hdf5.makecheck
	make install &> hdf5.install
	cd ..
    fi
fi
