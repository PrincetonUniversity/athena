# Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) and this project (mostly) adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).
We will attempt to follow the guideline of incrementing only one of the `X.Y.Z` values for each tag and/or release:
- `X` (MAJOR) version when you make incompatible API changes,
- `Y` (MINOR) version when you add functionality in a backwards-compatible manner, and
- `Z` (PATCH) version when you make backwards-compatible bug fixes.

The `X` vs. `Y` division won’t be strictly followed for Athena++ releases; for example, certain backwards-compatible versions may be released as a new `X` value to signify major new physics capabilities. As of `v1.1.0`, the Athena++ public API is only loosely documented in the GitHub Wiki, so the notion of backwards-compatibility is ambiguous. Nevertheless, versions with major changes to existing Athena++ core classes and  functions will generally be released under a new `X` value.

All major changes to the Athena++ private repository between each version/tag are summarized in this `CHANGLEOG.md` document. Each version has an **Issues and Pull Requests** section, whose subsections are automatically populated from the issue/PR labels. The list entries contain links to the private repository issue tracker `#N` id.

Additionally, the changes are **manually** summarized using the following categories:
- **Added:** for brand new features or extended capabilities
- **Fixed/Changed:** for bug fixes or modified behavior
- **Removed:** for removed functionality

The automatically-generated content should be used for reference when writing these sections. At this time, both the private and public [GitHub Release Notes](https://help.github.com/articles/creating-releases/) are started by copy/pasting from these sections.

<!-- "Implemented enhancements" (enhancement) vs. "Merged pull requests" (feature request, etc.) division doesn't make a ton of sense-->
<!-- Eventually, need to add label for "backwards-incompatible" and announce "BREAKING CHANGES" -->

## [Unreleased](https://github.com/PrincetonUniversity/athena/tree/HEAD)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v1.1.1-dev...HEAD)

### Added
Feature branches to merge to `master`:
- Chemistry (`chemistry`)
- Fourth-order solvers
  - Hydrodynamics (`hydro4`)
  - MHD (`mhd4`, `mhd4_3D`)

## [v1.1.1](https://github.com/PrincetonUniversity/athena/tree/v1.1.1) (2018-07-31)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v1.1.1-dev...v1.1.1)

### Removed
- Multigrid solver and related regression tests (removed directly from `v1.1.1-dev`)

## [v1.1.1-dev](https://github.com/PrincetonUniversity/athena/tree/v1.1.1-dev) (2018-07-31)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v1.1.0-dev...v1.1.1-dev)

### Added
- Performance improvements for SR calculations
- Pragmas and clauses from OpenMP 4.0 and 4.5 for SIMD-enabled functions inside vectorized loops (instead of depending on the compiler to inline the function calls)
- LLF Riemann solver for Newtonian hydrodynamics
- Rules for checking Python style of Athena++ scripts with `flake8`
- Interactive plotting for spherical polar coordinates results
- Several new input files and problem generators

### Fixed/Changed
- Re-enabled SIMD vectorization for Roe-type Riemann solvers and rewrote eigenmatrix calculations to improve performance
- Improved readability of `TimeIntegratorTaskList`
- Fixed spherical coordinates terms for non-ideal MHD
- Fixed reflective symmetry preservation for hydrodynamic viscosity calculations
- Eliminated small floating-point errors when analyzing uniform grid results using included Python HDF5 reader
- Changed turbulence driving switches to avoid possible bug during initial cycle
- Plugged MPI resource leaks in Multigrid

### Issues and Pull Requests:

#### Fixed bugs:

- Gravity boundary functions not set when enrolled [\#148](https://github.com/PrincetonUniversity/athena/issues/148)
- Viscosity terms break Rayleigh-Taylor symmetry [\#144](https://github.com/PrincetonUniversity/athena/issues/144)
- SIMD vectorization disabled for Roe-type Riemann solvers [\#126](https://github.com/PrincetonUniversity/athena/issues/126)
- Memory leak in jeans\_3d.py test MPI run [\#115](https://github.com/PrincetonUniversity/athena/issues/115)
- athena\_read.py athdf\(\) results in small floating-point errors [\#111](https://github.com/PrincetonUniversity/athena/issues/111)
- Remove unused GravityBoundaryFunction\_\[\] array from Mesh class [\#149](https://github.com/PrincetonUniversity/athena/pull/149) ([felker](https://github.com/felker))
- Cleanup minor issues before v1.1.1 release; fix viscosity asymmetry  [\#147](https://github.com/PrincetonUniversity/athena/pull/147) ([felker](https://github.com/felker))
- Return exact floating-point values when reading HDF5 coordinates [\#145](https://github.com/PrincetonUniversity/athena/pull/145) ([c-white](https://github.com/c-white))
- Fixed spherical coordinates for non-ideal MHD [\#142](https://github.com/PrincetonUniversity/athena/pull/142) ([tomidakn](https://github.com/tomidakn))
- Added interactive spherical plotting: [\#139](https://github.com/PrincetonUniversity/athena/pull/139) ([c-white](https://github.com/c-white))
- Use \#pragma omp declare simd for functions called in SIMD loops [\#138](https://github.com/PrincetonUniversity/athena/pull/138) ([felker](https://github.com/felker))
- Plugged MPI resource leaks in Multigrid. [\#137](https://github.com/PrincetonUniversity/athena/pull/137) ([tomidakn](https://github.com/tomidakn))
- AMR Fix [\#133](https://github.com/PrincetonUniversity/athena/pull/133) ([tomidakn](https://github.com/tomidakn))
- Bug fix for initial turbulence driving  [\#132](https://github.com/PrincetonUniversity/athena/pull/132) ([changgoo](https://github.com/changgoo))

#### Closed issues:

- Python scripts don't comply with PEP 8 style; vis scripts need documentation [\#96](https://github.com/PrincetonUniversity/athena/issues/96)

#### Merged pull requests:

- Improve performance for all SR problems [\#140](https://github.com/PrincetonUniversity/athena/pull/140) ([beiwang2003](https://github.com/beiwang2003))
- Redefine and initialize grav\_mean\_rho directly in Mesh class [\#136](https://github.com/PrincetonUniversity/athena/pull/136) ([changgoo](https://github.com/changgoo))
- Upgrade MPICH and OpenMPI in Travis CI build environments [\#135](https://github.com/PrincetonUniversity/athena/pull/135) ([felker](https://github.com/felker))
- Add Python style checker and fix existing violations [\#134](https://github.com/PrincetonUniversity/athena/pull/134) ([felker](https://github.com/felker))
- Improve CHANGELOG, release process, and CI testing [\#131](https://github.com/PrincetonUniversity/athena/pull/131) ([felker](https://github.com/felker))

## [v1.1.0](https://github.com/PrincetonUniversity/athena/tree/v1.1.0) (2018-05-23)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v1.1.0-dev...v1.1.0)

<!-- ### Added
### Fixed/Changed -->

### Removed
- Multigrid solver and related regression tests (removed directly from `v1.1.0-dev`)

## [v1.1.0-dev](https://github.com/PrincetonUniversity/athena/tree/v1.1.0-dev) (2018-05-23)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v1.0.1...v1.1.0-dev)

### Added
- Self-gravity
  - FFT
  - Multigrid
- Shearing box
- Viscosity, resistivity, and conduction
  - Anisotropic viscosity is not implemented, but it can be added as user-defined viscosity function. Can copy function's source code from Athena 4.2.
- Piecewise parabolic method (PPM) with `time/xorder=3` runtime option
  - Robust flooring of reconstructed states
  - Curvilinear and nonuniform mesh terms
- Characteristic variable reconstruction (PLM and PPM)  with `time/xorder=2c` or `3c` runtime option
- Redesign of time-integrator to support high-order schemes
- Turbulence driving
- Double precision floating-point HDF5 output
- Software development tools (mostly) exclusive to private repository:
  - Continuous integration (Jenkins and Travis CI)
  - C++ style checker
  - Expanded regression test suite flexibility and code coverage
  - Code Reviews, Issue/PR categories, `probot` automated closing of stale Issues/PRs, protected branch status for `master`
  - Slack workspace, `CONTRIBUTING.md` guide, Issue and PR templates.

### Fixed/Changed
- Changed OpenMP setup to coarse threading over `MeshBlock`; thread-safe MPI now in use.
- Performance optimizations: improved vectorization, decreased memory traffic, etc.
- Fixed reflective symmetry preservation (exact to double precision for hydrodynamics)
- Fixed all compiler warnings and ensured C++11 compliance

<!-- ### Removed -->

### Issues and Pull Requests:

#### Implemented enhancements:

- Pass MeshBlock to user defined diffusivity function [\#130](https://github.com/PrincetonUniversity/athena/pull/130) ([jmshi](https://github.com/jmshi))
- Use native HDF5 byte order [\#109](https://github.com/PrincetonUniversity/athena/issues/109)
- Conflicting restart and input parameters for writing output [\#62](https://github.com/PrincetonUniversity/athena/issues/62)
- Add double precision option for HDF5 output [\#28](https://github.com/PrincetonUniversity/athena/issues/28)
- Improve vectorization for loops in PPM [\#121](https://github.com/PrincetonUniversity/athena/pull/121) ([beiwang2003](https://github.com/beiwang2003))
- Improve symmetry-preservation of HLLD floating-point operations [\#120](https://github.com/PrincetonUniversity/athena/pull/120) ([felker](https://github.com/felker))
- Fix restarted simulation calculation of next\_time for writing outputs [\#116](https://github.com/PrincetonUniversity/athena/pull/116) ([felker](https://github.com/felker))
- Support double precision floating-point HDF5 output [\#108](https://github.com/PrincetonUniversity/athena/pull/108) ([felker](https://github.com/felker))
- Add templates for Issue and PR; creating draft of CONTRIBUTING.md [\#93](https://github.com/PrincetonUniversity/athena/pull/93) ([felker](https://github.com/felker))
- Make athdf a class which will delay reading until data is requested [\#84](https://github.com/PrincetonUniversity/athena/pull/84) ([msbc](https://github.com/msbc))
- Improve TimeIntegratorTaskList performance [\#80](https://github.com/PrincetonUniversity/athena/pull/80) ([felker](https://github.com/felker))
- Use AVX-512 registers for vectorized loops on Skylake; enable OpenMP 4.0 SIMD with newer GCC versions [\#77](https://github.com/PrincetonUniversity/athena/pull/77) ([felker](https://github.com/felker))
- Require negative xrat for user-defined mesh generation [\#71](https://github.com/PrincetonUniversity/athena/pull/71) ([felker](https://github.com/felker))
- Replace "\#pragma simd" with "\#pragma omp simd" for icc 18 [\#66](https://github.com/PrincetonUniversity/athena/pull/66) ([felker](https://github.com/felker))
- Make athinput reader Python 3 friendly in athena\_read.py [\#60](https://github.com/PrincetonUniversity/athena/pull/60) ([msbc](https://github.com/msbc))
- Add athinput reader to athena\_read.py [\#59](https://github.com/PrincetonUniversity/athena/pull/59) ([msbc](https://github.com/msbc))
- Ensure Python 2 and 3 compatibility of regression test suite [\#48](https://github.com/PrincetonUniversity/athena/pull/48) ([felker](https://github.com/felker))
- Improve Intel AVX2 vectorization of hydrodynamics routines [\#44](https://github.com/PrincetonUniversity/athena/pull/44) ([felker](https://github.com/felker))

#### Fixed bugs:

- Regression tests are writing .pyc files [\#104](https://github.com/PrincetonUniversity/athena/issues/104)
- Need to address g++, clang++, icc compiler warnings [\#94](https://github.com/PrincetonUniversity/athena/issues/94)
- PPM on spherical grid causes segmentation fault [\#76](https://github.com/PrincetonUniversity/athena/issues/76)
- Summed outputs have significant errors [\#73](https://github.com/PrincetonUniversity/athena/issues/73)
- Multigrid gravity NGHOST \> 2 causes artifacts [\#67](https://github.com/PrincetonUniversity/athena/issues/67)
- FFTBlock.norm\_factor\_ may be uninitialized [\#51](https://github.com/PrincetonUniversity/athena/issues/51)
- Regression test suite with MPI causes nonblocking IO write errors on macOS [\#47](https://github.com/PrincetonUniversity/athena/issues/47)
- UserWorkAfterLoop is called before final outputs [\#45](https://github.com/PrincetonUniversity/athena/issues/45)
- Restart regression test fails due to FFT input parameter parsing [\#42](https://github.com/PrincetonUniversity/athena/issues/42)
- Incorrect initialization of scopy\_ used for AthenaArray deep copy [\#38](https://github.com/PrincetonUniversity/athena/issues/38)
- Need better default behavior when setting start\_time [\#36](https://github.com/PrincetonUniversity/athena/issues/36)
- Fix SIMD vectorization and regression tests for diffusion module [\#127](https://github.com/PrincetonUniversity/athena/pull/127) ([jmshi](https://github.com/jmshi))
- Use simple average when computing cell-centered B  [\#122](https://github.com/PrincetonUniversity/athena/pull/122) ([felker](https://github.com/felker))
- Fix floating-point asymmetry of AMR cell-centered prolong/restrict operators [\#114](https://github.com/PrincetonUniversity/athena/pull/114) ([felker](https://github.com/felker))
- Fix asymmetry of Coordinates x1f for multiple MeshBlock grids [\#112](https://github.com/PrincetonUniversity/athena/pull/112) ([felker](https://github.com/felker))
- Fix general asymmetry of x1f, dx1v Coordinates terms; ensure exact float64 symmetry of PPM stencils [\#98](https://github.com/PrincetonUniversity/athena/pull/98) ([felker](https://github.com/felker))
- Implement various PPM fixes [\#79](https://github.com/PrincetonUniversity/athena/pull/79) ([felker](https://github.com/felker))
- Fix to self-gravity FFT and MG Solve [\#72](https://github.com/PrincetonUniversity/athena/pull/72) ([pdmullen](https://github.com/pdmullen))
- Add missing space in hdf5\_path string to fix bug [\#70](https://github.com/PrincetonUniversity/athena/pull/70) ([msbc](https://github.com/msbc))
- Fix self-gravity energy source term [\#63](https://github.com/PrincetonUniversity/athena/pull/63) ([pdmullen](https://github.com/pdmullen))
- Initialize FFTBlock::norm\_factor\_ to 1 [\#57](https://github.com/PrincetonUniversity/athena/pull/57) ([msbc](https://github.com/msbc))
- Improve stability of FFT gravity solver [\#40](https://github.com/PrincetonUniversity/athena/pull/40) ([changgoo](https://github.com/changgoo))
- Reconstruction for energy shouldn't be done for the isothermal EOS. [\#37](https://github.com/PrincetonUniversity/athena/pull/37) ([changgoo](https://github.com/changgoo))

#### Closed issues:

- C++ source code does not conform to Athena++ Style Guide [\#95](https://github.com/PrincetonUniversity/athena/issues/95)
- Wiki page on OpenMP parallelization is out of date [\#92](https://github.com/PrincetonUniversity/athena/issues/92)
- Require code reviews before merging pull requests [\#88](https://github.com/PrincetonUniversity/athena/issues/88)
- Create Slack workspace for Athena++ [\#61](https://github.com/PrincetonUniversity/athena/issues/61)
- Allow the regression test suite to accept platform-specific compiler and execution options [\#55](https://github.com/PrincetonUniversity/athena/issues/55)
- Regression test suite incompatible with Python 3 [\#46](https://github.com/PrincetonUniversity/athena/issues/46)
- Add optional runtime parameter for reducing/suppressing per-cycle summary info [\#39](https://github.com/PrincetonUniversity/athena/issues/39)
- Standardize sqrt\(\) and cbrt\(\) calls for C++11 compliance and portability [\#27](https://github.com/PrincetonUniversity/athena/issues/27)

#### Merged pull requests:

- Fix Bash error handling of C++ style checker [\#125](https://github.com/PrincetonUniversity/athena/pull/125) ([felker](https://github.com/felker))
- Remove inverse h\_coeff from Coordinates class [\#124](https://github.com/PrincetonUniversity/athena/pull/124) ([jmshi](https://github.com/jmshi))
- Add diffusion physics capabilities  [\#123](https://github.com/PrincetonUniversity/athena/pull/123) ([jmshi](https://github.com/jmshi))
- Use default simdlen for a memory intensive loop in PLM [\#119](https://github.com/PrincetonUniversity/athena/pull/119) ([beiwang2003](https://github.com/beiwang2003))
- Use loop fusion and define simdlen\(SIMD\_WIDTH\) for PLM [\#118](https://github.com/PrincetonUniversity/athena/pull/118) ([beiwang2003](https://github.com/beiwang2003))
- Add default implementations of Laplacian operators to Coordinates class [\#117](https://github.com/PrincetonUniversity/athena/pull/117) ([felker](https://github.com/felker))
- Codify and fix C++11 standard compliance [\#113](https://github.com/PrincetonUniversity/athena/pull/113) ([felker](https://github.com/felker))
- Improve HLLD vectorization [\#110](https://github.com/PrincetonUniversity/athena/pull/110) ([beiwang2003](https://github.com/beiwang2003))
- Add Google C++ Style Guide linter and check with CI [\#107](https://github.com/PrincetonUniversity/athena/pull/107) ([felker](https://github.com/felker))
- Fix main GCC compiler warnings; strictly test code changes for new warnings [\#105](https://github.com/PrincetonUniversity/athena/pull/105) ([felker](https://github.com/felker))
- Use C++ style static\_cast operator instead of C style type casting [\#101](https://github.com/PrincetonUniversity/athena/pull/101) ([felker](https://github.com/felker))
- Add continuous integration with Travis CI and Jenkins [\#97](https://github.com/PrincetonUniversity/athena/pull/97) ([felker](https://github.com/felker))
- Add Bash wrapper to join\_vtk++ for many MeshBlocks, output steps [\#91](https://github.com/PrincetonUniversity/athena/pull/91) ([felker](https://github.com/felker))
- Add CODEOWNERS file to new .github/ root dir [\#90](https://github.com/PrincetonUniversity/athena/pull/90) ([felker](https://github.com/felker))
- Add ShearingBox component   [\#89](https://github.com/PrincetonUniversity/athena/pull/89) ([jmshi](https://github.com/jmshi))
- Improve floating-point precision consistency throughout code [\#82](https://github.com/PrincetonUniversity/athena/pull/82) ([felker](https://github.com/felker))
- Add primitive variable flooring functions to the EquationOfState class  [\#81](https://github.com/PrincetonUniversity/athena/pull/81) ([felker](https://github.com/felker))
- Add MG and FFT convergence regression tests based on 3D linear Jeans instability [\#75](https://github.com/PrincetonUniversity/athena/pull/75) ([alwinm](https://github.com/alwinm))
- Add cylindrical and spherical polar coordinates limiter terms to PPM [\#69](https://github.com/PrincetonUniversity/athena/pull/69) ([felker](https://github.com/felker))
- Redesign the TimeIntegratorTask list to enable time integration methods with order, nsub\_steps, nregisters \> 2  [\#49](https://github.com/PrincetonUniversity/athena/pull/49) ([felker](https://github.com/felker))
- Enable nonuniform reconstruction with Mignone PPM limiter [\#43](https://github.com/PrincetonUniversity/athena/pull/43) ([felker](https://github.com/felker))
- Add ncycle\_out optional input parameter [\#41](https://github.com/PrincetonUniversity/athena/pull/41) ([felker](https://github.com/felker))
- Add FFT gravity solver [\#35](https://github.com/PrincetonUniversity/athena/pull/35) ([changgoo](https://github.com/changgoo))

## [v1.0.1](https://github.com/PrincetonUniversity/athena/tree/v1.0.1) (2017-08-22)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v1.0.0-beta...v1.0.1)

### Added
- Scripts and documentation for releasing from private to public repository

### Fixed/Changed
- Upgraded left and right reconstructed primitive state `AthenaArray` in `Hydro` class to 3D (from `x1`-sliced 1D arrays)
- Riemann solvers now directly fill electric field arrays

<!-- ### Removed -->

### Issues and Pull Requests:

#### Implemented enhancements:

- Add code structure diagrams to Wiki [\#34](https://github.com/PrincetonUniversity/athena/issues/34)

#### Fixed bugs:

- Alfven wave test fails due to incorrect error thresholds [\#29](https://github.com/PrincetonUniversity/athena/issues/29)
- Cannot compile properly on tiger [\#8](https://github.com/PrincetonUniversity/athena/issues/8)
- Do not access de-allocated memory in SliceOutputData [\#32](https://github.com/PrincetonUniversity/athena/pull/32) ([dradice](https://github.com/dradice))

## [v1.0.0-beta](https://github.com/PrincetonUniversity/athena/tree/v1.0.0-beta) (2016-11-28)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v0.3.0...v1.0.0-beta)

### Added
- BSD license

### Fixed/Changed
- Redesigned `Coordinates` class
<!-- ### Removed -->

### Issues and Pull Requests:

#### Implemented enhancements:

- ParameterInput should be generally accessible [\#20](https://github.com/PrincetonUniversity/athena/issues/20)
- Need OMP-ready scratch array for relativistic MHD Riemann solver [\#6](https://github.com/PrincetonUniversity/athena/issues/6)
- Small configure script enhancements [\#21](https://github.com/PrincetonUniversity/athena/pull/21) ([jzuhone](https://github.com/jzuhone))

#### Fixed bugs:

- Blast problem missing header needed to compile [\#25](https://github.com/PrincetonUniversity/athena/issues/25)
- MPI + HDF5 writer with ghost zone output causes segmentation fault [\#23](https://github.com/PrincetonUniversity/athena/issues/23)
- Compilation failure when using HDF5 without MPI [\#22](https://github.com/PrincetonUniversity/athena/issues/22)
- 1D Newtonian MHD broken with periodic boundary [\#16](https://github.com/PrincetonUniversity/athena/issues/16)
- Problem with fields at physical boundaries in 1D [\#15](https://github.com/PrincetonUniversity/athena/issues/15)
- Line missing in src/pgen/jet.cpp [\#12](https://github.com/PrincetonUniversity/athena/issues/12)
- Compilation error -- accessing invalid member of RegionSize [\#11](https://github.com/PrincetonUniversity/athena/issues/11)

#### Closed issues:

- Add  boundary functions based on primitive variables [\#17](https://github.com/PrincetonUniversity/athena/issues/17)
- Need flux correction for source terms [\#9](https://github.com/PrincetonUniversity/athena/issues/9)

## [v0.3.0](https://github.com/PrincetonUniversity/athena/tree/v0.3.0) (2015-08-31)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v0.2.0...v0.3.0)

### Added

### Fixed/Changed

<!-- ### Removed -->

### Issues and Pull Requests:

#### Fixed bugs:

- Uninitialized variable root\_level is used [\#10](https://github.com/PrincetonUniversity/athena/issues/10)

## [v0.2.0](https://github.com/PrincetonUniversity/athena/tree/v0.2.0) (2015-04-10)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v0.1.0...v0.2.0)

### Added

### Fixed/Changed

<!-- ### Removed -->

### Issues and Pull Requests:

#### Fixed bugs:

- Restart output problem due to next\_time=0 [\#7](https://github.com/PrincetonUniversity/athena/issues/7)

#### Closed issues:

- Add MeshBlock boundary condition for MHD [\#4](https://github.com/PrincetonUniversity/athena/issues/4)

## [v0.1.0](https://github.com/PrincetonUniversity/athena/tree/v0.1.0) (2015-02-02)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/efaebe9aa67738beb74879d3cd8b7cdaecab441c...v0.1.0)
### Added

### Fixed/Changed

<!-- ### Removed -->

### Issues and Pull Requests:

#### Fixed bugs:

- Problem with boundary functions in 1D hydro shock tube [\#5](https://github.com/PrincetonUniversity/athena/issues/5)

#### Closed issues:

- Add C++ style document to wiki [\#1](https://github.com/PrincetonUniversity/athena/issues/1)

\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*
