### FFT wrapper
The Athena FFT wrapper is designed to use [FFTW3](http://www.fftw.org) for serial FFT and [Plimpton's library](http://www.sandia.gov/~sjplimp/docs/fft/README.html) for parallel FFT (again utilizing FFTW3 for the underlying serial FFT).
* Currently, only works on uniform grids without refinement
* Can have multiple `MeshBlocks` per processor
* Supports any domain decomposition, but the best performance is expected for pencil decomposition
* [FFTW3 library](http://www.fftw.org) is required whenever using the `-fft` configuration flag. The library does not need to be built with `--with-mpi`, since only the serial functionality is used in Athena++

To use FFT, configure the code with the `-fft` flag, and specify the exact `--fftw_path [path_to_FFTW_library]` if the library is not in the default linker search path. For example,
```
    > ./configure.py --prob=pgen_name -fft --fftw_path=path_to_FFTW_library
```

### FFT Driver and Block
The Athena FFT wrapper uses Athena's `Mesh` and `MeshBlock` information. In particular, the `FFTDriver` class is initialized with a pointer to the `Mesh` object. Conceptually, an `FFTBlock` is similar to a `MeshBlock`, but for the parallel FFT solver **each MPI rank can have only one `FFTBlock`** that spans every `MeshBlock` owned by the MPI rank. Also, the shape of this `FFTBlock` should be a cuboid (or **the number of `MeshBlocks` per processor should be a power of 2** for the Z-ordering). During the `FFTDriver` object construction, the code automatically checks that each MPI rank possess a single cuboid `FFTBlock`.

`LoadSource()/RetrieveResult()` member functions can automatically populate an input/output array to/from `FFTBlock` with `LogicalLocation` and `RegionSize` information of `MeshBlock` (see below).

### Example Calling Sequence (MeshBlock/processor=1)
```c++
    FFTDriver *pfftd;
    pfftd = new FFTDriver(pMesh, pin);
    // initialize FFTBlocks. With true, it will set the `norm_factor_=1/gcnt_` that to be multiplied for the backward FFT.
    // can be done separtely `pfftd->pmy_fb->SetNormFactor(norm)`
    pfftd->InitializeFFTBlock(true);
    // automatic creation of forward/backward FFT plans
    pfftd->QuickCreatePlan();

    FFTBlock *pfft = pfftd->pmy_fb;

    LogicalLocation &loc = pblock->loc;
    RegionSize &block_size = pblock->block_size;

    // src array will be loaded to pfft->in_
    // 1 for real, 2 for complex (real, imaginary) src array
    pfft->LoadSource(src,1,NGHOST,loc,block_size);
    pfft->ExecuteForward();
    // specific kernel to be applied with option
    // default is simply swapping in_ and out_ arrays for preparing backward FFT that recovers source
    pfft->ApplyKernel(0);
    pfft->ExecuteBackward();
    // retrieve results in pfft->out_ to dst array
    // 1 retrieve only the real part, 2 for both real and imaginary parts
    pfft->RetrieveResult(dst,2,NGHOST,loc,block_size);
```
* See also `mesh/mesh.cpp` for `FFTDriver` constructor calls
* See below examples for more practical usages

### Use cases
* [[Self-Gravity with FFT]] `athena/src/gravity/fftgravity.cpp`
* [[Turbulence Driver]] `athena/src/fft/turblence.cpp`

The simplest example problem generator is located in `athena/src/pgen/fft.cpp`
```
    > ./configure.py --prob=fft -fft
    > ./athena -i ../inputs/hydro/athinput.fft
```
