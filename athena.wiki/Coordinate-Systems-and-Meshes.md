### Select a Coordinate System
Currently the following coordinate systems are available:
* Cartesian (default; <var>x</var>, <var>y</var>, <var>z</var>)
* Cylindrical coordinates (<var>r</var>, <var>&#981;</var>, <var>z</var>)
* Spherical-polar coordinates (<var>r</var>, <var>&theta;</var>, <var>&#981;</var>)
* Various general-relativistic coordinates:
  * Minkowski
  * Schwarzschild
  * Kerr-Schild
  * User-defined metric

To select a coordinate system, configure the code with the `--coord [coord]` option.
```
> python ./configure.py --prob sphtorus --coord spherical_polar
```
The names of coordinate systems are corresponding to the file names in `src/coordinates`, and are listed in the help message (use the `-h` option).

#### 2D Polar Coordinates
Note, if a complete circular mesh in 2D is desired, the appropriate choice for the coordinate system is `--coord=cylindrical`, **not** `spherical_polar`. The `x2` polar angle coordinate is strictly limited to the [0, π] range in both 2D and 3D spherical-polar coordinates. The `x2` azimuthal angle coordinate can take on any value in [0, 2π) for the cylindrical coordinate system in 2D and 3D.

In Athena++, 2D spherical-like coordinates are used to simulate a very thin "orange slice" of the domain.

### Mesh Spacing
By default, the cell size is simply set evenly: `dx = L / Nx`. It is sometimes convenient to use geometric spacing, i.e. the cell size increases (or decreases) by a constant factor <var>r</var>. This is especially useful in spherical-polar coordinates so that it can maintain the cell aspect ratio. To enable this, set `x1rat, x2rat`, or `x3rat` in the `<mesh>` block in the input file.
```
    <mesh>
    ...
    x1rat = 1.01
```
While Athena++ uses formulae consistent with nonuniform mesh spacing, too large a mesh ratio should not be used; a non-fatal warning is raised if the value is not between 0.9 and 1.1 everywhere. The reconstruction algorithm reduces to first-order if the grid spacing changes too much from cell-to-cell, and moreover, the nonuniform truncation error can produce anomalous reflections of waves.

Note, the geometric spacing is equivalent to logarithmic spacing (uniform `dx'` in logspace given by the coordinate transformation <var>x'</var> = <span>log<sub>b</sub></span>(<var>x<sub>1</sub></var>) for some <span><var>b>0 </var>, <var>b</var> &ne; 1</span>) only for a particular choice of (`x1rat`, `x1min`, and `x1max`).

Note that this ratio should be adjusted when resolution is manually changed. For example, when the number of cells in the direction is doubled, the ratio should be adjusted to square root of the original one to maintain consistency. This adjustment is done automatically whenever mesh refinement is invoked.

### User-Defined Mesh Generator
More generally, any kind of mesh spacing can be used through the user-defined mesh generator feature. Similar to the user-defined boundary conditions, define the MeshGenerator function first:
```c++
    Real MyMeshSpacingX2(Real x, RegionSize rs);
```
And then enroll it in `Mesh::InitUserMeshData`:
```c++
    EnrollUserMeshGenerator(X2DIR, MyMeshSpacingX2);
```
Finally, the input variable `mesh/x2rat=-1.0` must be set in [[The Input File]] or on the command line in order to run a simulation with a user-defined MeshGenerator function in that direction.

The MeshGenerator function in x1, for example, must be a monotonic, 1:1 mapping from [0:1] to `[rs.x1min:rs.x1max]`, which returns a coordinate location between `rs.x1min:rs.x1max` for given `x = i/Nx` (the logical location of the cell).

For example, let us consider a mesh generator which produces high resolution near the midplane in the <var>&theta;</var> (x2) direction in spherical-polar coordinates.
```c++
    Real CompressedX2(Real x, RegionSize rs)
    {
      Real t=2.0*x-1.0;
      Real w=0.25*(t*(t*t+1.0))+0.5;
      return w*rs.x2max+(1.0-w) * rs.x2min;
    }
```

This function gives mesh spacing proportional to (3<var>&psi;</var><sup>2</sup> + 1) where <var>&psi;</var> is the angle from the midplane, and the mesh spacing near the pole is about 4 times larger than that near the midplane. The code will give a warning if the resulting mesh contains a steep change in the mesh size.

<!-- TODO: consider splitting the remainder of this page into separate wiki page-->
### Symmetry preservation in Athena++
Hydrodynamics and MHD test problems with physical symmetries are discriminating tests of computational astrophysics codes. In particular, linear instabilities may amplify any small floating-point differences (as small as 1 [unit in the last place (ULP)](https://en.wikipedia.org/wiki/Unit_in_the_last_place)), resulting in visible asymmetries in late-time solutions. The issues typically become more pronounced when using schemes with high-order accuracy; we have observed solutions with symmetry errors many orders of magnitude larger when using `time/xorder=3` PPM than when using `time/xorder=2` PLM reconstruction, for example.

While the underlying integration algorithm of Athena++ is dimensionally unsplit and is well-suited to preserving symmetries *in theory*, there are several (related) challenges when attempting to preserve symmetry with floating-point operations *in practice*:
* Both absolute and relative round-off errors are asymmetric about any value != 0.0
* Floating-point arithmetic is non-associative
* Compilers may enable value-unsafe optimizations of floating-point operations

Although symmetry-preservation is a desirable solver behavior, it may not be worth the developer effort (or even feasible) to maintain certain complex symmetry relationships for all possible test problems and physics. The design choices can be difficult, especially if the required changes would degrade the clarity, generality, or performance of the overall solver (see [[Style Guide]]).

<!-- Language note: reflection in/along a coordinate direction; reflection about an axis or plane -->
Athena++ has been designed to preserve certain test problem symmetries **exactly** to double precision accuracy. This capability is automatic; no additional options when [[Configuring]] or [[Running the Code]] are required to preserve reflective symmetry in a coordinate direction (e.g. transformation in `x1` about `x2-x3` plane), as long as the following constraints are satisfied:
* Uniform mesh using Cartesian coordinates
  * Confirmed working with AMR/SMR without MPI
* Symmetric initial condition
* Hydrodynamics problems only
  * Constrained transport (CT) magnetic field updates currently break symmetry due to non-associativity of operations in `field/calculate_corner_e.cpp`
* Plane of symmetry located at `x1=0.0`; this restriction is sufficient for any number of `MeshBlock`s and `nx1` (even or odd valued ---> reflection either about cell faces `x1f=0.0` or volume-centers `x1v=0.0`, respectively)
* *Known bug:* Configured without MPI. Solvers using unrefined grids with MPI might preserve symmetry; AMR with MPI likely won't.
* Compiled with a value-safe floating-point arithmetic mode

This has been tested up to 64-bit floating-point precision for a variety of resolutions, additional physics, Riemann solvers, and reconstruction methods. 

<!-- Certain optional capabilities are known/suspected to introduce asymmetric errors when enabled:
* N/A -->

*Miscellaneous notes for users*:
- When performing symmetry tests and using HDF5 output, the code should be configured with `-h5double` for maximum precision.
- The [[Analysis Tools]] provided in the repository includes the `athdf()` class in `athena_read.py` for HDF5 data. This reader was designed to analyze both AMR/SMR and unrefined results. When reading AMR/SMR results into Python, the script interpolates quantities to a common grid and may introduce asymmetries. However, when provided with unrefined, uniform grid HDF5 data, the reader is guaranteed to preserve any symmetries present in the raw floating-point data. Alternatively, the `h5py` module can be used to directly load the `.athdf` files.

<!--
Several examples of 2D and 3D problem generators in Athena++ with physical symmetries are discussed.

#### Grid-aligned reflection symmetry
If the 1D axis of symmetry (2D domain) or 2D plane of symmetry (3D domain) of a problem is aligned along one of the coordinate directions, we refer to the symmetry relationship as 1D reflection symmetry

Solutions to the 2D and 3D (single-mode `iprob=1`) Rayleigh-Taylor problems ought to be symmetric about the midplane at `x1=0` for all time.

See also: Kelvin-Hemholtz (x2)

#### Diagonal reflection symmetry
The only oblique reflection symmetry that Athena++ can preserve is when the axis/plane of symmetry is along a grid diagonal.

The 2D Liska-Wendroff hydrodynamics implosion test has an axis of symmetry along the grid diagonal `x1=x2`. 2D reflective symmetry across the coordinate diagonal is a special type of reflective symmetry that has different concerns than 1D reflection symmetry.

#### Discrete rotational symmetry of order 2 / point reflection
In two dimensions, a point reflection is the same as a rotation of 180 degrees. I

The Orszag-Tang multidimensional vortex problem is an example of a problem whose density solution should always remain identical when rotated by 180 degrees. The solution to the standard MHD formulation of the problem will not exhibit perfect symmetry in Athena++ due to the current  

However, by modifying the initial condition to specify `B=0`, we can test the rotational symmetry preservation of the hydrodynamics solver, which is exactly preserved.  

#### Axisymmetry -- impossible on rectilinear grid
Spherical hydro blast wave

#### Translation symmetry
TODO: Discuss translation symmetry and Galilean invariance, See https://princetonuniversity.github.io/Athena-Cversion/AthenaDocsGalileo
-->

#### Compiler-dependent floating-point math modes
The Intel compiler is more aggressive with math optimizations than either the Clang or GNU Compiler Collection (GCC). Both `clang++` and `g++` [disable `-ffast-math`](https://gcc.gnu.org/onlinedocs/gcc-4.5.2/gcc/Optimize-Options.html) at all optimization levels, by default. In contrast, the default `icc` math mode is `-fp-model fast=1`, which is **value-unsafe**. Studies of the above symmetric problems in Athena++ have shown that compiling with `-fp-model strict` is required to preserve symmetry with the Intel compiler. However, this flag prevents the use of OpenMP 4.0 SIMD vectorization extensions; thus, the performance of the code drops dramatically. The less restrictive `-fp-model precise` option is value-safe, but is insufficient to prevent the creation of asymmetries in solutions to the 2D Rayleigh-Taylor problem, for example.

According to [Intel fp-model documentation](https://software.intel.com/en-us/cpp-compiler-18.0-developer-guide-and-reference-fp-model-fp), the only differences between `strict` and `precise` is that the former 1) also enables floating-point exception semantics and 2) disables floating-point contractions, including fused multiply-add (FMA) instructions. Both modes ensure reproducibility by prohibiting any reordering/change of associativity of floating-point operations:
> disables optimizations that can change the result of floating-point calculations, which is required for strict ANSI conformance.

Future work will ensure symmetry-preservation with the Intel compiler while enabling SIMD vectorization. The ability of OpenMP threading to reorder floating-point operations of inside SIMD loops can be toggled via the environment variable `KMP_DETERMINISTIC_REDUCTION`; setting it to `0` may be sufficient to maintain problem symmetry with `-fp-model precise`:
> Enables (1) or disables (0) the use of a specific ordering of the reduction operations for implementing the reduction clause for an OpenMP* parallel region. This has the effect that, for a given number of threads, in a given parallel region, for a given data set and reduction operation, a floating point reduction done for an OpenMP reduction clause will have a consistent floating point result from run to run, since round-off errors will be identical.

#### Tips for Athena++ developers
When adding new physics or changing existing Athena++ source code, it is important to keep the following concerns for preserving symmetry in mind:
- User-specified initial conditions, source terms, refinement conditions, etc. that depend on a function of position should refer to the cell-interface positions, if possible. For example, `f(x1f)` is preferred to `f(x1v)`.
- Any operation that is expressed in separate per-direction calculations, such as functions `OperatorX1()`, `OperatorX2()`, and `OperatorX3()`, should use identical sequences of floating-point operations. For example, earlier versions of the `PiecewiseLinearX*()` functions directly divided a quantity for reconstructions along `x1`, whereas an inverse was precomputed and stored for reconstructions along `x2`:
```c++
// PLMx1()
Real& dx_im2 = pco->dx1v(i-2);
dql = (q(nin,k,j,i-1) - q(nin,k,j,i-2))/dx_im2;
// PLMx2()
Real dx2jm2i = 1.0/pco->dx2v(j-2);
dql = (q(nin,k,j-1,i) - q(nin,k,j-2,i))*dx2jm2i;
```
The latter sequence incurs an intermediate rounding error that does not occur in the first formulation. This difference created visible artifacts in the Liska-Wendroff implosion problem at high-resolutions.
- The [[Regression Testing]] suite should be periodically extended to check symmetry maintenance of new problems. Unlike issues with error convergence, symmetry violation bugs can be introduced by very minor changes and may require significant developer effort to debug.
- Developers should be sure that stencils/expressions that reference mesh quantities indexed across a range of cells in `x1`, `x2`, or `x3` greater than +/- 1 are written **without biased floating-point associativity**. For example, in `reconstruct/ppm.cpp`, the original calculation of approximations to the second-derivative of input `AthenaArray q` along the `x1` direction (indexed by `i` in variable names) was expressed as:
```c++
d2qc_im1(i) = q_im2(n,i) - 2.0*q_im1(n,i) + q    (n,i);
d2qc    (i) = q_im1(n,i) - 2.0*q    (n,i) + q_ip1(n,i); //(CD eq 85a) (no 1/2)
d2qc_ip1(i) = q    (n,i) - 2.0*q_ip1(n,i) + q_ip2(n,i);
```
Even though these stencils are "centered" (`d2qc` expression equally weights `im1` and `ip1` quantities) in infinite precision arithmetic, they are **not centered** in finite-precision arithmetic. The solver always adds/subtracts the lower indexed quantities first due to implicit left-to-right associativity in floating-point arithmetic. Hence, the round-off error from subtracting the first two terms (lower `i`) "dominates" the round-off error in the final result and introduces asymmetries in the solution. Simply replacing the above lines with:
```c++
d2qc_im1(i) = q_im2(n,i) + q    (n,i) - 2.0*q_im1(n,i) ;
d2qc    (i) = q_im1(n,i) + q_ip1(n,i) - 2.0*q(n,i); //(CD eq 85a) (no 1/2)
d2qc_ip1(i) = q    (n,i) + q_ip2(n,i) - 2.0*q_ip1(n,i) ;
```
guarantees that the `ip1` and `im1` quantities can be swapped in the middle expression and return the identical floating-point value, for example.
<!-- they can commute w.r.t. round-off error in ternary finite-precision operator f(x,y,z)-->

**Implementation in `Coordinates` class**

The class is currently designed such that the calculation of cell face positions in the `Coordinates()` class constructor uses a special inlined function `UniformMeshGeneratorX1()` in `mesh.hpp`. This function parameterizes the real cell faces in the calling `MeshBlock` by [-0.5, 0.5] instead of the [0,1] range used by the `DefaultMeshGeneratorX1()` function for user-defined mesh generators. The same function is used to compute the boundary positions of each `MeshBlock`. It was designed to preserve exact floating point symmetry of `x1f` if `x1min=-x1max` and work for `nx1` even (central face at `x1f=0`) or odd (central cell has `x1v=0`). Furthermore, `dx1f=dx1v` are guaranteed to be identical and constant when `x1rat=1.0` for a Cartesian mesh.
<!-- recommended: use Nx1 as a power of 2.-->
