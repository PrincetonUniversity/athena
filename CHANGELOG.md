# Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) and this project (mostly) adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).
We will attempt to follow the guideline of incrementing only one of the `X.Y.Z` values for each tag and/or release:
- `X` (MAJOR) version when you make incompatible API changes,
- `Y` (MINOR) version when you add functionality in a backwards-compatible manner, and
- `Z` (PATCH) version when you make backwards-compatible bug fixes.

The `X` vs. `Y` division won’t be strictly followed for Athena++ versioning; certain backwards-compatible releases may be released as a new `X` value to signify major new physics capabilities, for example. As of `v1.1.0`, the Athena++ public API is only loosely documented in the GitHub Wiki, so the notion of backwards-compatibility is ambiguous. 

## [Unreleased](https://github.com/PrincetonUniversity/athena/tree/HEAD)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v1.1.0...HEAD)

#### Fixed bugs:

- SIMD vectorization disabled for Roe-type Riemann solvers [\#126](https://github.com/PrincetonUniversity/athena/issues/126)

## [v1.1.0](https://github.com/PrincetonUniversity/athena/tree/v1.1.0) (2018-05-23)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v1.1.0-dev...v1.1.0)

#### Merged pull requests:

- Removal of the Multigrid sover from the next public version [\#128](https://github.com/PrincetonUniversity/athena/pull/128) ([tomidakn](https://github.com/tomidakn))

## [v1.1.0-dev](https://github.com/PrincetonUniversity/athena/tree/v1.1.0-dev) (2018-05-23)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v1.0.1...v1.1.0-dev)

#### Implemented enhancements:

- HDF5 byte order [\#109](https://github.com/PrincetonUniversity/athena/issues/109)
- athena\_read.py delayed loading data [\#83](https://github.com/PrincetonUniversity/athena/issues/83)
- Double precision in HDF5 [\#28](https://github.com/PrincetonUniversity/athena/issues/28)
- improve vectorization for loops in ppm.cpp [\#121](https://github.com/PrincetonUniversity/athena/pull/121) ([beiwang2003](https://github.com/beiwang2003))
- Improve symmetry-preservation of HLLD floating point operations [\#120](https://github.com/PrincetonUniversity/athena/pull/120) ([felker](https://github.com/felker))
- Fix restarted simulation calculation of next\_time for writing outputs [\#116](https://github.com/PrincetonUniversity/athena/pull/116) ([felker](https://github.com/felker))
- Support double precision floating point HDF5 output [\#108](https://github.com/PrincetonUniversity/athena/pull/108) ([felker](https://github.com/felker))
- Add templates for Issue and PR; creating draft of CONTRIBUTING.md [\#93](https://github.com/PrincetonUniversity/athena/pull/93) ([felker](https://github.com/felker))
- Made athdf a class which will delay reading data until it's asked for. [\#84](https://github.com/PrincetonUniversity/athena/pull/84) ([msbc](https://github.com/msbc))

#### Fixed bugs:

- Regression tests are writing .pyc files [\#104](https://github.com/PrincetonUniversity/athena/issues/104)
- Compiler warnings [\#94](https://github.com/PrincetonUniversity/athena/issues/94)
- Segfault with PPM on spherical grid [\#76](https://github.com/PrincetonUniversity/athena/issues/76)
- Summed outputs have significant errors [\#73](https://github.com/PrincetonUniversity/athena/issues/73)
- restart/input conflicts for output files [\#62](https://github.com/PrincetonUniversity/athena/issues/62)
- FFTBlock.norm\_factor\_ not initialized properly [\#51](https://github.com/PrincetonUniversity/athena/issues/51)
- Regression test suite with MPI causes nonblocking IO write errors on macOS [\#47](https://github.com/PrincetonUniversity/athena/issues/47)
- UserWorkAfterLoop called before final outputs [\#45](https://github.com/PrincetonUniversity/athena/issues/45)
- Restart regression test fails due to FFT input parameter parsing [\#42](https://github.com/PrincetonUniversity/athena/issues/42)
- scopy\_ initialization [\#38](https://github.com/PrincetonUniversity/athena/issues/38)
- Better default behavior when setting start\_time [\#36](https://github.com/PrincetonUniversity/athena/issues/36)
- Fix SIMD vectorization and regression tests for diffusion module [\#127](https://github.com/PrincetonUniversity/athena/pull/127) ([jmshi](https://github.com/jmshi))
- Use simple average when computing cell-centered B  [\#122](https://github.com/PrincetonUniversity/athena/pull/122) ([felker](https://github.com/felker))
- Fix floating-point asymmetry of AMR cell-centered prolong/restrict operators [\#114](https://github.com/PrincetonUniversity/athena/pull/114) ([felker](https://github.com/felker))
- Fix asymmetry of Coordinates x1f for multiple MeshBlock grids [\#112](https://github.com/PrincetonUniversity/athena/pull/112) ([felker](https://github.com/felker))
- Fix general asymmetry of x1f, dx1v Coordinates terms; ensure exact float64 symmetry of PPM stencils [\#98](https://github.com/PrincetonUniversity/athena/pull/98) ([felker](https://github.com/felker))
- Reconstruction for energy shouldn't be done for the isothermal EOS. [\#37](https://github.com/PrincetonUniversity/athena/pull/37) ([changgoo](https://github.com/changgoo))

#### Closed issues:

- C++ source code does not conform to Athena++ Style Guide [\#95](https://github.com/PrincetonUniversity/athena/issues/95)
- Wiki page on OpenMP parallelization is out of date [\#92](https://github.com/PrincetonUniversity/athena/issues/92)
- Require code reviews before merging pull requests [\#88](https://github.com/PrincetonUniversity/athena/issues/88)
- accidentally merged pull request, please undo [\#87](https://github.com/PrincetonUniversity/athena/issues/87)
- Multigrid gravity NGHOST \> 2 causes artifacts [\#67](https://github.com/PrincetonUniversity/athena/issues/67)
- Slack workspace for Athena++ [\#61](https://github.com/PrincetonUniversity/athena/issues/61)
- Allow the regression test suite to accept platform-specific compiler and execution options [\#55](https://github.com/PrincetonUniversity/athena/issues/55)
- Regression test suite incompatible with Python 3 [\#46](https://github.com/PrincetonUniversity/athena/issues/46)
- Add optional runtime parameter for reducing/suppressing per-cycle summary info [\#39](https://github.com/PrincetonUniversity/athena/issues/39)
- Standardizing sqrt\(\) and cbrt\(\) calls [\#27](https://github.com/PrincetonUniversity/athena/issues/27)

#### Merged pull requests:

- pass meshblock in user defined diffusivity [\#130](https://github.com/PrincetonUniversity/athena/pull/130) ([jmshi](https://github.com/jmshi))
- Fix Bash error handling of C++ style checker [\#125](https://github.com/PrincetonUniversity/athena/pull/125) ([felker](https://github.com/felker))
- remove inverse h coeff [\#124](https://github.com/PrincetonUniversity/athena/pull/124) ([jmshi](https://github.com/jmshi))
- Pull request for diffusion physics [\#123](https://github.com/PrincetonUniversity/athena/pull/123) ([jmshi](https://github.com/jmshi))
- use default simdlen\(e.g.,4\) for a memory intensive loop in plm [\#119](https://github.com/PrincetonUniversity/athena/pull/119) ([beiwang2003](https://github.com/beiwang2003))
- use loop fusion and simdlen\(SIMD\_WIDTH\) for plm.cpp [\#118](https://github.com/PrincetonUniversity/athena/pull/118) ([beiwang2003](https://github.com/beiwang2003))
- Add default implementations of Laplacian operators to Coordinates class [\#117](https://github.com/PrincetonUniversity/athena/pull/117) ([felker](https://github.com/felker))
- C++11 standard compliance [\#113](https://github.com/PrincetonUniversity/athena/pull/113) ([felker](https://github.com/felker))
- improve vectorization for RiemannSolver \(mhd/hlld.cpp\) [\#110](https://github.com/PrincetonUniversity/athena/pull/110) ([beiwang2003](https://github.com/beiwang2003))
- Add Google C++ Style Guide linter and check with CI [\#107](https://github.com/PrincetonUniversity/athena/pull/107) ([felker](https://github.com/felker))
- Fix main GCC compiler warnings; strictly test code changes for new warnings [\#105](https://github.com/PrincetonUniversity/athena/pull/105) ([felker](https://github.com/felker))
- Use C++ style static\_cast operator instead of C style type casting [\#101](https://github.com/PrincetonUniversity/athena/pull/101) ([felker](https://github.com/felker))
- Revert "Revert "Fix general asymmetry of x1f, dx1v Coordinates terms; ensure exact float64 symmetry of PPM stencils"" [\#100](https://github.com/PrincetonUniversity/athena/pull/100) ([felker](https://github.com/felker))
- Revert "Fix general asymmetry of x1f, dx1v Coordinates terms; ensure exact float64 symmetry of PPM stencils" [\#99](https://github.com/PrincetonUniversity/athena/pull/99) ([felker](https://github.com/felker))
- Add continuous integration with Travis CI and Jenkins [\#97](https://github.com/PrincetonUniversity/athena/pull/97) ([felker](https://github.com/felker))
- Add Bash wrapper to join\_vtk++ for many MeshBlocks, output steps [\#91](https://github.com/PrincetonUniversity/athena/pull/91) ([felker](https://github.com/felker))
- Add CODEOWNERS file to new .github/ root dir [\#90](https://github.com/PrincetonUniversity/athena/pull/90) ([felker](https://github.com/felker))
- adding shearingbox component   [\#89](https://github.com/PrincetonUniversity/athena/pull/89) ([jmshi](https://github.com/jmshi))
- Revert "Simple gravity source", pushed to wrong origin [\#86](https://github.com/PrincetonUniversity/athena/pull/86) ([xzackli](https://github.com/xzackli))
- Simple gravity source [\#85](https://github.com/PrincetonUniversity/athena/pull/85) ([xzackli](https://github.com/xzackli))
- Improve floating point precision consistency throughout code [\#82](https://github.com/PrincetonUniversity/athena/pull/82) ([felker](https://github.com/felker))
- Add primitive variable flooring functions to the EquationOfState class  [\#81](https://github.com/PrincetonUniversity/athena/pull/81) ([felker](https://github.com/felker))
- Improve TimeIntegratorTaskList performance [\#80](https://github.com/PrincetonUniversity/athena/pull/80) ([felker](https://github.com/felker))
- Various PPM fixes [\#79](https://github.com/PrincetonUniversity/athena/pull/79) ([felker](https://github.com/felker))
- Use AVX-512 registers for vectorized loops on Skylake; enable OpenMP 4.0 SIMD with newer GCC versions [\#77](https://github.com/PrincetonUniversity/athena/pull/77) ([felker](https://github.com/felker))
- Add MG and FFT convergence regression tests based on 3D linear Jeans instability [\#75](https://github.com/PrincetonUniversity/athena/pull/75) ([alwinm](https://github.com/alwinm))
- Fix to self-gravity FFT and MG Solve [\#72](https://github.com/PrincetonUniversity/athena/pull/72) ([pdmullen](https://github.com/pdmullen))
- Require negative xrat for user-defined mesh generation [\#71](https://github.com/PrincetonUniversity/athena/pull/71) ([felker](https://github.com/felker))
- Add missing space in hdf5\_path string to fix bug [\#70](https://github.com/PrincetonUniversity/athena/pull/70) ([msbc](https://github.com/msbc))
- Add cylindrical and spherical polar coordinates limiter terms to PPM [\#69](https://github.com/PrincetonUniversity/athena/pull/69) ([felker](https://github.com/felker))
- Replace "\#pragma simd" with "\#pragma omp simd" for icc 18 [\#66](https://github.com/PrincetonUniversity/athena/pull/66) ([felker](https://github.com/felker))
- Fix to self gravity energy source term [\#63](https://github.com/PrincetonUniversity/athena/pull/63) ([pdmullen](https://github.com/pdmullen))
- Read athena [\#60](https://github.com/PrincetonUniversity/athena/pull/60) ([msbc](https://github.com/msbc))
- Read athinput [\#59](https://github.com/PrincetonUniversity/athena/pull/59) ([msbc](https://github.com/msbc))
- Init FFTBlock::norm\_factor\_ to 1 [\#57](https://github.com/PrincetonUniversity/athena/pull/57) ([msbc](https://github.com/msbc))
- Redesign the TimeIntegratorTask list to enable time integration methods with order, nsub\_steps, nregisters \> 2  [\#49](https://github.com/PrincetonUniversity/athena/pull/49) ([felker](https://github.com/felker))
- Python2-3 compatible regression test suite [\#48](https://github.com/PrincetonUniversity/athena/pull/48) ([felker](https://github.com/felker))
- Intel AVX2 vectorization improvements for hydrodynamics [\#44](https://github.com/PrincetonUniversity/athena/pull/44) ([felker](https://github.com/felker))
- Nonuniform reconstruction with Mignone PPM limiter [\#43](https://github.com/PrincetonUniversity/athena/pull/43) ([felker](https://github.com/felker))
- Add ncycle\_out optional input parameter [\#41](https://github.com/PrincetonUniversity/athena/pull/41) ([felker](https://github.com/felker))
- fft gravity solver update [\#40](https://github.com/PrincetonUniversity/athena/pull/40) ([changgoo](https://github.com/changgoo))
- FFT gravity [\#35](https://github.com/PrincetonUniversity/athena/pull/35) ([changgoo](https://github.com/changgoo))

## [v1.0.1](https://github.com/PrincetonUniversity/athena/tree/v1.0.1) (2017-08-22)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v1.0.0-beta...v1.0.1)

#### Implemented enhancements:

- Code structure diagrams now exist [\#34](https://github.com/PrincetonUniversity/athena/issues/34)

#### Fixed bugs:

- Alfven wave test failure [\#29](https://github.com/PrincetonUniversity/athena/issues/29)

#### Closed issues:

- Cannot compile properly on tiger [\#8](https://github.com/PrincetonUniversity/athena/issues/8)

#### Merged pull requests:

- Do not access de-allocated memory in SliceOutputData [\#32](https://github.com/PrincetonUniversity/athena/pull/32) ([dradice](https://github.com/dradice))

## [v1.0.0-beta](https://github.com/PrincetonUniversity/athena/tree/v1.0.0-beta) (2016-11-28)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v0.3.0...v1.0.0-beta)

#### Fixed bugs:

- blast problem missing header needed to compile [\#25](https://github.com/PrincetonUniversity/athena/issues/25)
- Segfault with mpi/hdf5/ghost zones [\#23](https://github.com/PrincetonUniversity/athena/issues/23)
- 1D Newtonian MHD broken with periodic boundary [\#16](https://github.com/PrincetonUniversity/athena/issues/16)
- Problem with fields at physical boundaries in 1D [\#15](https://github.com/PrincetonUniversity/athena/issues/15)
- Compilation error -- accessing invalid member of RegionSize [\#11](https://github.com/PrincetonUniversity/athena/issues/11)
- Need OMP-ready scratch array for relativistic MHD Riemann solver [\#6](https://github.com/PrincetonUniversity/athena/issues/6)

#### Closed issues:

- Compilation failure when using HDF5 without MPI [\#22](https://github.com/PrincetonUniversity/athena/issues/22)
- ParameterInput should be generally accessible [\#20](https://github.com/PrincetonUniversity/athena/issues/20)
- User physical src term cannot be enrolled when restarting [\#18](https://github.com/PrincetonUniversity/athena/issues/18)
- Primitive boundary functions [\#17](https://github.com/PrincetonUniversity/athena/issues/17)
- Line missing in src/pgen/jet.cpp [\#12](https://github.com/PrincetonUniversity/athena/issues/12)
- Flux correction for source terms [\#9](https://github.com/PrincetonUniversity/athena/issues/9)

#### Merged pull requests:

- Small configure script enhancements [\#21](https://github.com/PrincetonUniversity/athena/pull/21) ([jzuhone](https://github.com/jzuhone))

## [v0.3.0](https://github.com/PrincetonUniversity/athena/tree/v0.3.0) (2015-08-31)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v0.2.0...v0.3.0)

#### Fixed bugs:

- Uninitialized variable root\_level is used [\#10](https://github.com/PrincetonUniversity/athena/issues/10)

## [v0.2.0](https://github.com/PrincetonUniversity/athena/tree/v0.2.0) (2015-04-10)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v0.1.0...v0.2.0)

#### Closed issues:

- Restart output problem [\#7](https://github.com/PrincetonUniversity/athena/issues/7)
- MHD block boundary is needed [\#4](https://github.com/PrincetonUniversity/athena/issues/4)



\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*
