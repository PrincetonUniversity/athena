### Self-Gravity with FFT
A new `Gravity` class is constructed as a member of `MeshBlock`. The `Gravity` class contains the gravitational potential array `phi`. For the FFT gravity solver, the boundary condition is also separately set for `phi` (which is unnecessary for the multigrid gravity solver; see [[Self-Gravity with Multigrid]]).

When self-gravity is enabled, `phi` is automatically included in any [[Outputs]] that include either the `prim` or `cons` variables. Alternatively, it can be explicitly included in the `variable` field as `phi`.

A new `FFTGravityDriver` class is constructed and solves Poisson's equation using [[FFT]] assuming periodic BCs in all directions. The `FFTGravityDriver::Solve()` function is called in the main loop at every time-integrator stage, and it performs the following steps:
1. load density from MeshBlocks to the FFT input array
2. execute forward FFT
3. multiply kernel (-k<sup>-2</sup>; `SetNormFactor()` takes the `four_pi_G` factor into account)
4. execute backward FFT
5. retrieve the real part of the FFT output, and store in the gravitational potential array `phi`
6. call gravity tasklist to send/receive boundary values

The momentum source term due to the gravitational force is included in the flux calculation by adding the gravitational tensor. See `athena/src/hydro/gravity_fluxes.cpp`.

The energy source term due to the gravitational energy flux is separately added as a source term. See `athena/src/hydro/srcterm/gravity.cpp`.

The gravitational constant (`four_pi_G`) should be set in `Mesh::InitMeshData` using `SetGravityConstant` or `SetFourPiG` function.

The background mean density (`grav_mean_rho`) should be set in `Mesh::InitMeshData` using `SetMeanDensity` functions. For problems with the Jeans' swindle (= fully periodic BCs), this should be the initial mean density (see `athena/src/pgen/jeans.cpp`; see also [#155](https://github.com/PrincetonUniversity/athena/issues/155)). For problems with other BCs (not implemented yet), this should be zero. 
 
### Regression Test (Jeans problem)
```
> cd tst/regression/
> ./run_tests.py grav
```

### Example
```
> ./configure.py --prob=poisson -fft --grav=fft
> ./athena -i ../inputs/hydro/athinput.poisson 
```