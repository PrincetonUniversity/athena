#!/bin/bash

# FFTW 3.3.8 = released 2018-05-28
# FFTW 3.3.7 = released 2017-10-29
version_str=3.3.8
# for macOS builds, install FFTW from Homebrew
if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    export HOMEBREW_NO_AUTO_UPDATE=1
    cd fftw
    brew link fftw || true
    # check to see if FFTW is cached from previous build
    if [ -f "bin/fft-wisdom" ]; then
	echo "Using cached FFTW"
    else
        echo "Installing FFTW with Homebrew"
	HOMEBREW_TEMP=$TRAVIS_BUILD_DIR/fftw
	# brew update
        brew install fftw
	brew unlink fftw
	# /usr/local/opt symlinks to Cellar are preserved, use these:
	# (note, /usr/local/bin/bin -> /usr/local/bin = circular symlink)
	ln -s /usr/local/opt/fftw/bin bin || true
	ln -s /usr/local/opt/fftw/lib lib || true
	ln -s /usr/local/opt/fftw/include include || true
	ln -s /usr/local/opt/fftw/share share || true
    fi
else
    # for Ubuntu builds, install FFTW from source
    # check to see if FFTW is cached from previous build
    if [ -f "fftw/bin/fft-wisdom" ]; then # && [ -f "fftw-${version_str}/config.log" ]
	echo "Using cached FFTW"
    else
	# download, configure, and compile FFTW3
	echo "Downloading FFTW Source"
	wget http://www.fftw.org/fftw-${version_str}.tar.gz
	tar zxf fftw-${version_str}.tar.gz
	cd fftw-${version_str}
	echo "Configuring and building FFTW"
	# --enable-mpi: not needed, since Plimpton's FFT library is used in Athena++ for MPI+FFT grav
	./configure --prefix=$TRAVIS_BUILD_DIR/fftw &> fftw.configure
	# make -j4 &> fftw.make
	make &> fftw.make
	make install &> fftw.install
	cd ..
    fi
fi
