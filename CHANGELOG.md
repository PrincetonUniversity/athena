# Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/) and this project (mostly) adheres to Calendar Versioning.

All major changes to the Athena++ private repository between each version/tag are summarized in this `CHANGLEOG.md` document. Each version has an **Issues and Pull Requests** section, whose subsections are automatically populated from the issue/PR labels. The list entries contain links to the private repository issue tracker `#N` id.

Additionally, the changes are **manually** summarized using the following categories:
- **Added:** for brand new features or extended capabilities
- **Fixed/Changed:** for bug fixes or modified behavior
- **Removed:** for removed functionality

The automatically-generated content should be used for reference when writing these sections. At this time, both the private and public [GitHub Release Notes](https://help.github.com/articles/creating-releases/) are started by copy/pasting from these sections.

<!-- "Implemented enhancements" (enhancement) vs. "Merged pull requests" (feature request, etc.) division doesn't make a ton of sense-->
<!-- Eventually, need to add label for "backwards-incompatible" and announce "BREAKING CHANGES" -->

## [Unreleased](https://github.com/PrincetonUniversity/athena/tree/HEAD)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v21.0...HEAD)

## [v21.0](https://github.com/PrincetonUniversity/athena/tree/v21.0) (2021-01-06)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v21.0-dev...v21.0)

### Removed
- Multigrid solver and related regression tests (removed directly from `v21.0-dev`)

## [v21.0-dev](https://github.com/PrincetonUniversity/athena/tree/v21.0-dev) (2021-01-06)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v19.0...v21.0-dev)

#### Implemented enhancements:

- Need more flexible task control in STS [\#270](https://github.com/PrincetonUniversity/athena/issues/270)
- Enable AMR/SMR and source terms with -sts [\#250](https://github.com/PrincetonUniversity/athena/issues/250)
- Robustness improvements to relativity [\#316](https://github.com/PrincetonUniversity/athena/pull/316) ([c-white](https://github.com/c-white))
- Change primitive velocities in SR [\#303](https://github.com/PrincetonUniversity/athena/pull/303) ([c-white](https://github.com/c-white))
- Update passive scalars to have user source terms and work with relativity [\#297](https://github.com/PrincetonUniversity/athena/pull/297) ([c-white](https://github.com/c-white))

#### Fixed bugs:

- Spurious waves emerging from corners using the intel 19.1 compiler and MPI? [\#340](https://github.com/PrincetonUniversity/athena/issues/340)
- EOS Riemann regression test fails [\#335](https://github.com/PrincetonUniversity/athena/issues/335)
- Unused and uninitialized variables with orbital advection [\#332](https://github.com/PrincetonUniversity/athena/issues/332)
- regression tests for diffusion do not work & failed [\#328](https://github.com/PrincetonUniversity/athena/issues/328)
- Reconstruction of passive scalars uses wrong array size [\#306](https://github.com/PrincetonUniversity/athena/issues/306)
- malloc error with gcc compiler, spherical polar coordinate and xorder=3 [\#304](https://github.com/PrincetonUniversity/athena/issues/304)
- RK3 integrator breaks conservation [\#300](https://github.com/PrincetonUniversity/athena/issues/300)
- Runtime error "pure virtual method called" with intel compiler [\#289](https://github.com/PrincetonUniversity/athena/issues/289)
- Test GR/SR compatibility with passive scalars [\#281](https://github.com/PrincetonUniversity/athena/issues/281)
- New Limitation of MeshBlock Size [\#343](https://github.com/PrincetonUniversity/athena/pull/343) ([tomo-ono](https://github.com/tomo-ono))
- Move Poisson solve to after ptlist-\>DoTaskListOneStage\(\) [\#309](https://github.com/PrincetonUniversity/athena/pull/309) ([pdmullen](https://github.com/pdmullen))
- Add support for NSCALARS \> NWAVE [\#307](https://github.com/PrincetonUniversity/athena/pull/307) ([msbc](https://github.com/msbc))
- Fix logic in bvars\_sts designation [\#305](https://github.com/PrincetonUniversity/athena/pull/305) ([pdmullen](https://github.com/pdmullen))

#### Closed issues:

- bugs in using std::min or std::max for single precision [\#326](https://github.com/PrincetonUniversity/athena/issues/326)
- Mismatch of the function argument in WeightedAve for FaceField [\#318](https://github.com/PrincetonUniversity/athena/issues/318)

#### Merged pull requests:

- Implemented the source term formula for gravitational acceleration [\#344](https://github.com/PrincetonUniversity/athena/pull/344) ([tomidakn](https://github.com/tomidakn))
- Enable STS in Shearing Box with Orbital Advection [\#342](https://github.com/PrincetonUniversity/athena/pull/342) ([tomo-ono](https://github.com/tomo-ono))
- Logical location \(long\)is used to set array index \(int\) -- Jenkins passed! [\#341](https://github.com/PrincetonUniversity/athena/pull/341) ([changgoo](https://github.com/changgoo))
- Fix lack of deallocating dynamic memory [\#338](https://github.com/PrincetonUniversity/athena/pull/338) ([tomo-ono](https://github.com/tomo-ono))
- fix dependency in tasklist [\#337](https://github.com/PrincetonUniversity/athena/pull/337) ([tomo-ono](https://github.com/tomo-ono))
- Increased thresholds for eos\_riemann test [\#336](https://github.com/PrincetonUniversity/athena/pull/336) ([msbc](https://github.com/msbc))
- fix EMF flux correction in shearing box [\#334](https://github.com/PrincetonUniversity/athena/pull/334) ([tomo-ono](https://github.com/tomo-ono))
- fix unused and uninitialized variables [\#333](https://github.com/PrincetonUniversity/athena/pull/333) ([tomo-ono](https://github.com/tomo-ono))
- Fix conditionals in TimeIntegratorTaskList fns for compatibility with STS [\#331](https://github.com/PrincetonUniversity/athena/pull/331) ([pdmullen](https://github.com/pdmullen))
- Docstring Style Correction [\#330](https://github.com/PrincetonUniversity/athena/pull/330) ([changgoo](https://github.com/changgoo))
- Correct the input block name in the turbulence regression test and make driving isotropic regardless of the box shape [\#329](https://github.com/PrincetonUniversity/athena/pull/329) ([changgoo](https://github.com/changgoo))
- put cast in the min or max function [\#327](https://github.com/PrincetonUniversity/athena/pull/327) ([tomo-ono](https://github.com/tomo-ono))
- Implement a New Discretization of the Gravitational Stress Tensor [\#325](https://github.com/PrincetonUniversity/athena/pull/325) ([pdmullen](https://github.com/pdmullen))
- Change in turbulence input style [\#324](https://github.com/PrincetonUniversity/athena/pull/324) ([changgoo](https://github.com/changgoo))
- Fix Bug with Mesh Structure Output [\#322](https://github.com/PrincetonUniversity/athena/pull/322) ([tomo-ono](https://github.com/tomo-ono))
- Orbital Advection & New Shearing Box [\#321](https://github.com/PrincetonUniversity/athena/pull/321) ([tomo-ono](https://github.com/tomo-ono))
- Fix typo: bounday\_flag.cpp to boundary\_flag.cpp [\#320](https://github.com/PrincetonUniversity/athena/pull/320) ([changgoo](https://github.com/changgoo))
- Fix mismatch between WeightedAve prototype and definition [\#319](https://github.com/PrincetonUniversity/athena/pull/319) ([pdmullen](https://github.com/pdmullen))
- Improve performance of turbulence driver [\#317](https://github.com/PrincetonUniversity/athena/pull/317) ([changgoo](https://github.com/changgoo))
- Update to bug fix on fft pgen [\#315](https://github.com/PrincetonUniversity/athena/pull/315) ([changgoo](https://github.com/changgoo))
- Fix typo in HLLC and bug in Hydrogen EOS [\#314](https://github.com/PrincetonUniversity/athena/pull/314) ([msbc](https://github.com/msbc))
- Bug fixes for FFT; New regression test [\#313](https://github.com/PrincetonUniversity/athena/pull/313) ([changgoo](https://github.com/changgoo))
- Modifying the KH problem and a input file for MHD KH [\#312](https://github.com/PrincetonUniversity/athena/pull/312) ([tomo-ono](https://github.com/tomo-ono))

## [v19.0](https://github.com/PrincetonUniversity/athena/tree/v19.0) (2019-08-06)

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v19.0-dev...v19.0)

### Removed
- Multigrid solver and related regression tests (removed directly from `v19.0-dev`)

## [v19.0-dev](https://github.com/PrincetonUniversity/athena/tree/v19.0-dev) (2019-08-06)

**Switched to using CalVer (Calendar Versioning) scheme from SemVer (Semantic Versioning)**. `v1.1.1-dev` was the previous tag and version.

[Full Changelog](https://github.com/PrincetonUniversity/athena/compare/v1.1.1...v19.0-dev)

#### Implemented enhancements:

- Extend user-defined history outputs beyond parallel summed quantities  [\#237](https://github.com/PrincetonUniversity/athena/issues/237)
- Add test coverage for Roe Riemann solvers [\#231](https://github.com/PrincetonUniversity/athena/issues/231)
- Add "logging" module to all regression tests to replace calls to print\(\) [\#203](https://github.com/PrincetonUniversity/athena/issues/203)
- Warn users about proper polar vs. polar\_wedge boundary condition flag   [\#196](https://github.com/PrincetonUniversity/athena/issues/196)
- Loading SMR data at native resolution [\#173](https://github.com/PrincetonUniversity/athena/issues/173)
- Create cleaner failure mode for incorrect input parameter file text formatting [\#168](https://github.com/PrincetonUniversity/athena/issues/168)
- Extend vis/python/uniform.py to 2D HDF5 files [\#164](https://github.com/PrincetonUniversity/athena/issues/164)
- \(Non-Zero\) Mean Density Must Be Set for FFT Self-Gravity with Periodic BC's [\#155](https://github.com/PrincetonUniversity/athena/issues/155)
- include input parameters in hdf5 output? [\#143](https://github.com/PrincetonUniversity/athena/issues/143)
- Extend TaskList to support more than 64 task IDs [\#286](https://github.com/PrincetonUniversity/athena/pull/286) ([tomidakn](https://github.com/tomidakn))
- Add timing to regression tests [\#284](https://github.com/PrincetonUniversity/athena/pull/284) ([msbc](https://github.com/msbc))
- Throw error when reading HDF5 files with ghost zones if their number hasn't been specified. [\#280](https://github.com/PrincetonUniversity/athena/pull/280) ([msbc](https://github.com/msbc))
- Add user ability to specify names of libraries when linking [\#277](https://github.com/PrincetonUniversity/athena/pull/277) ([c-white](https://github.com/c-white))
- Remove MeshBlock dependency of FFTDriver [\#276](https://github.com/PrincetonUniversity/athena/pull/276) ([changgoo](https://github.com/changgoo))
- Move MGGravity from GravityDriver to MeshBlock [\#275](https://github.com/PrincetonUniversity/athena/pull/275) ([tomidakn](https://github.com/tomidakn))
- Add input option to output more timestep diagnostics [\#273](https://github.com/PrincetonUniversity/athena/pull/273) ([felker](https://github.com/felker))
- Refactor MeshBlockTree to save memory footprint [\#262](https://github.com/PrincetonUniversity/athena/pull/262) ([tomidakn](https://github.com/tomidakn))
- Replace AthenaArray\<T\>::InitWithShallowCopy with C++ references [\#261](https://github.com/PrincetonUniversity/athena/pull/261) ([felker](https://github.com/felker))
- Rewrite shearing box capabilities  [\#260](https://github.com/PrincetonUniversity/athena/pull/260) ([felker](https://github.com/felker))
- Remove do-nothing DeleteAthenaArray calls [\#258](https://github.com/PrincetonUniversity/athena/pull/258) ([felker](https://github.com/felker))
- Refactor bvals/ and Mesh::AdaptiveMeshRefinement\(\) [\#257](https://github.com/PrincetonUniversity/athena/pull/257) ([felker](https://github.com/felker))
- Replace AthenaFFTComplex with std::complex\<Real\> [\#256](https://github.com/PrincetonUniversity/athena/pull/256) ([felker](https://github.com/felker))
- Deduplicate 4x STS regression tests [\#252](https://github.com/PrincetonUniversity/athena/pull/252) ([pdmullen](https://github.com/pdmullen))
- Adds logging module to all regression tests to replace calls to print\(\) [\#251](https://github.com/PrincetonUniversity/athena/pull/251) ([msbc](https://github.com/msbc))
- Global initialization of turbulence realization [\#241](https://github.com/PrincetonUniversity/athena/pull/241) ([changgoo](https://github.com/changgoo))
- Enable FFT for an entirely refined mesh [\#240](https://github.com/PrincetonUniversity/athena/pull/240) ([changgoo](https://github.com/changgoo))
- Optimize performance of SR/GR calculations [\#235](https://github.com/PrincetonUniversity/athena/pull/235) ([beiwang2003](https://github.com/beiwang2003))
- Improve all enumerated types [\#234](https://github.com/PrincetonUniversity/athena/pull/234) ([felker](https://github.com/felker))
- Add test coverage for MHD and Hydro Roe Riemann solvers  [\#233](https://github.com/PrincetonUniversity/athena/pull/233) ([msbc](https://github.com/msbc))
- Replace holdovers from C with C++11 alternatives [\#219](https://github.com/PrincetonUniversity/athena/pull/219) ([felker](https://github.com/felker))
- Update ShowConfig for general EOS [\#218](https://github.com/PrincetonUniversity/athena/pull/218) ([msbc](https://github.com/msbc))
- Extend line plotting to .athdf files [\#211](https://github.com/PrincetonUniversity/athena/pull/211) ([c-white](https://github.com/c-white))
- Make C++ exception handling optional; standardize more C++11 headers and namespace usage [\#190](https://github.com/PrincetonUniversity/athena/pull/190) ([felker](https://github.com/felker))
- Extend Kelvin-Helmholtz problem generator  [\#187](https://github.com/PrincetonUniversity/athena/pull/187) ([felker](https://github.com/felker))
- Use C++ compiler front ends for main/documented "configure.py --cxx" options but support C front ends; add built-in support for Apple+LLVM with OpenMP workaround [\#186](https://github.com/PrincetonUniversity/athena/pull/186) ([felker](https://github.com/felker))

#### Fixed bugs:

- athena\_read.py silently fails when ghost zones are saved to file [\#279](https://github.com/PrincetonUniversity/athena/issues/279)
- MPI + continuous turbulence driving with OU process is broken  [\#267](https://github.com/PrincetonUniversity/athena/issues/267)
- Invalid read in AMR [\#229](https://github.com/PrincetonUniversity/athena/issues/229)
- NaN'ing outputs may not be caught by some regression tests [\#223](https://github.com/PrincetonUniversity/athena/issues/223)
- Bug with permuting fft directions in athena\_fft.cpp [\#217](https://github.com/PrincetonUniversity/athena/issues/217)
- 1D AMR not refining third component of magnetic field correctly [\#212](https://github.com/PrincetonUniversity/athena/issues/212)
- MPI/AMR segfaulting [\#199](https://github.com/PrincetonUniversity/athena/issues/199)
- RECV\_\* boundary buffer receive tasks are missing dependencies on matching INT\_\* integration tasks [\#198](https://github.com/PrincetonUniversity/athena/issues/198)
- Possible incorrect memory access in spherical\_polar coordinates: private & public version [\#192](https://github.com/PrincetonUniversity/athena/issues/192)
- athena\_read.athdf\(\) handles ghost zones incorrectly [\#184](https://github.com/PrincetonUniversity/athena/issues/184)
- Possible race condition in OpenMP  [\#183](https://github.com/PrincetonUniversity/athena/issues/183)
- SMR causes issues with cell sizes difference warning in spherical coordinates with user mesh generator [\#180](https://github.com/PrincetonUniversity/athena/issues/180)
- uniform.py now writes all 0s when overwriting a file [\#179](https://github.com/PrincetonUniversity/athena/issues/179)
- OpenMP is broken; add regression tests to prevent this in the future [\#175](https://github.com/PrincetonUniversity/athena/issues/175)
- Some regression tests never fail [\#174](https://github.com/PrincetonUniversity/athena/issues/174)
- Python reader fails for slices taken across multiple blocks [\#159](https://github.com/PrincetonUniversity/athena/issues/159)
- Python HDF5 reader no longer returns refinement levels [\#156](https://github.com/PrincetonUniversity/athena/issues/156)
- Is ApplyPrimitiveFloors\(\) being called too often? [\#153](https://github.com/PrincetonUniversity/athena/issues/153)
- Fix bug when saving zero scalars to restart file [\#285](https://github.com/PrincetonUniversity/athena/pull/285) ([msbc](https://github.com/msbc))
- Fix bug when saving multiple scalars to restart file [\#283](https://github.com/PrincetonUniversity/athena/pull/283) ([msbc](https://github.com/msbc))
- Fix minor issues with Python HDF5 reader [\#278](https://github.com/PrincetonUniversity/athena/pull/278) ([c-white](https://github.com/c-white))
- Fix curvilinear reconstruction, add regression test scripts \(passive scalars, MPI shearing box, etc.\), and other QoL improvements [\#274](https://github.com/PrincetonUniversity/athena/pull/274) ([felker](https://github.com/felker))
- Revert GR optimizations [\#272](https://github.com/PrincetonUniversity/athena/pull/272) ([c-white](https://github.com/c-white))
- Fix a NaN bug \(\#267\) caused by the incorrect out-of-place FFT [\#268](https://github.com/PrincetonUniversity/athena/pull/268) ([changgoo](https://github.com/changgoo))
- Fix bug in spherical plotting script [\#265](https://github.com/PrincetonUniversity/athena/pull/265) ([c-white](https://github.com/c-white))
- Fix shearing box + MPI [\#264](https://github.com/PrincetonUniversity/athena/pull/264) ([felker](https://github.com/felker))
- Fix FFT permutation bug [\#238](https://github.com/PrincetonUniversity/athena/pull/238) ([changgoo](https://github.com/changgoo))
- Improve NaN handling in Python scripts [\#224](https://github.com/PrincetonUniversity/athena/pull/224) ([felker](https://github.com/felker))
- Fix fail criterion for diffusion/\*\_diffusion\* tests [\#222](https://github.com/PrincetonUniversity/athena/pull/222) ([pdmullen](https://github.com/pdmullen))
- Fix typo in 1D prolongation of magnetic field [\#216](https://github.com/PrincetonUniversity/athena/pull/216) ([c-white](https://github.com/c-white))
- Fix general EOS bugs and add eos\_test pgen [\#214](https://github.com/PrincetonUniversity/athena/pull/214) ([msbc](https://github.com/msbc))
- 1D MHD AMR fix [\#213](https://github.com/PrincetonUniversity/athena/pull/213) ([c-white](https://github.com/c-white))
- Correct energy bug in PrimitiveToConserved in general\_hydro.cpp [\#208](https://github.com/PrincetonUniversity/athena/pull/208) ([msbc](https://github.com/msbc))
- Apply minor fixes to Python scripts [\#201](https://github.com/PrincetonUniversity/athena/pull/201) ([c-white](https://github.com/c-white))
- Fix memory access issue with spherical diffusion [\#195](https://github.com/PrincetonUniversity/athena/pull/195) ([c-white](https://github.com/c-white))
- Improve post-processing and plotting scripts [\#185](https://github.com/PrincetonUniversity/athena/pull/185) ([c-white](https://github.com/c-white))
- Fix minor OpenMP error and duplicate MPI regression test for OpenMP, hybrid MPI+OpenMP [\#178](https://github.com/PrincetonUniversity/athena/pull/178) ([felker](https://github.com/felker))
- Fix and extend Python scripts for analysis of HDF5 files [\#167](https://github.com/PrincetonUniversity/athena/pull/167) ([c-white](https://github.com/c-white))
- Fix Python HDF5 reader's ability to return refinement levels [\#158](https://github.com/PrincetonUniversity/athena/pull/158) ([msbc](https://github.com/msbc))
- Eliminate redundant calls to ApplyPrimitiveFloors\(\) [\#154](https://github.com/PrincetonUniversity/athena/pull/154) ([felker](https://github.com/felker))

#### Closed issues:

- Add optional time/sts\_nstage\_out input parameter [\#249](https://github.com/PrincetonUniversity/athena/issues/249)
- Deduplicate 4x STS regression test scripts in diffusion/ [\#242](https://github.com/PrincetonUniversity/athena/issues/242)
- "register" keyword used in src/fft/plimpton/ files is deprecated and was removed in C++17 [\#230](https://github.com/PrincetonUniversity/athena/issues/230)
- Make `athinput` a file extension as opposed to a file prefix. [\#226](https://github.com/PrincetonUniversity/athena/issues/226)
- 1D MHD SMR fields wrong at refinement boundaries with more than 2 ghost zones [\#209](https://github.com/PrincetonUniversity/athena/issues/209)
- Add interface for operator splitting and super-time-stepping [\#172](https://github.com/PrincetonUniversity/athena/issues/172)

#### Merged pull requests:

- Add user-defined passive scalar floor option in input file [\#291](https://github.com/PrincetonUniversity/athena/pull/291) ([munan](https://github.com/munan))
- Add MHD capabilities to the general EOS module [\#282](https://github.com/PrincetonUniversity/athena/pull/282) ([msbc](https://github.com/msbc))
- Add passive scalars [\#263](https://github.com/PrincetonUniversity/athena/pull/263) ([felker](https://github.com/felker))
- Extend automatic and manual load balancing beyond AMR [\#259](https://github.com/PrincetonUniversity/athena/pull/259) ([tomidakn](https://github.com/tomidakn))
- Fix edge-case bug in typecast function. [\#255](https://github.com/PrincetonUniversity/athena/pull/255) ([msbc](https://github.com/msbc))
- Implement OU Process for driven turbulence and Helmholtz decomposition [\#244](https://github.com/PrincetonUniversity/athena/pull/244) ([changgoo](https://github.com/changgoo))
- Remove "register" keyword from 2x .cpp files in plimpton/ [\#232](https://github.com/PrincetonUniversity/athena/pull/232) ([felker](https://github.com/felker))
- Fix general EOS hydro function names and prepare for general EOS MHD extension [\#228](https://github.com/PrincetonUniversity/athena/pull/228) ([msbc](https://github.com/msbc))
- Add compatibility checks for polar, polar\_wedge boundary conditions to class constructor [\#200](https://github.com/PrincetonUniversity/athena/pull/200) ([felker](https://github.com/felker))
- Generalize the hydrodynamics equation of state capabilities for non-ideal gases [\#197](https://github.com/PrincetonUniversity/athena/pull/197) ([msbc](https://github.com/msbc))
- Revert to using 1D pencil arrays for performance [\#193](https://github.com/PrincetonUniversity/athena/pull/193) ([tomidakn](https://github.com/tomidakn))
- Update input params for diffusion regression problems [\#188](https://github.com/PrincetonUniversity/athena/pull/188) ([pdmullen](https://github.com/pdmullen))
- Integrate code coverage analysis into regression test framework and continuous integration [\#181](https://github.com/PrincetonUniversity/athena/pull/181) ([felker](https://github.com/felker))
- Implement RKL1 Super-Time-Stepping [\#176](https://github.com/PrincetonUniversity/athena/pull/176) ([pdmullen](https://github.com/pdmullen))
- Add fourth-order hydrodynamics solver [\#157](https://github.com/PrincetonUniversity/athena/pull/157) ([felker](https://github.com/felker))
- Add new script for plotting 2D slices [\#151](https://github.com/PrincetonUniversity/athena/pull/151) ([c-white](https://github.com/c-white))
- Add HDF5 reader and pgen for initializing from HDF5 [\#146](https://github.com/PrincetonUniversity/athena/pull/146) ([c-white](https://github.com/c-white))

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
