### Input File
The boundary condition flags must be specified in the input file.
```
    <mesh>
    ix1_bc     = user          # inner-X1 boundary flag
    ox1_bc     = user          # outer-X1 boundary flag
    ix2_bc     = outflow       # inner-X2 boundary flag
    ox2_bc     = outflow       # outer-X2 boundary flag
    ix3_bc     = periodic      # inner-X3 boundary flag
    ox3_bc     = periodic      # outer-X3 boundary flag
```
The following flags are available:
* `outflow` : outflowing (zero-gradient) boundary condition
* `reflecting` : reflecting boundary condition
* `periodic` : periodic boundary condition
* `polar` : geometric boundary condition at the pole, available only for the x2-direction in 3D spherical-polar and similar coordinates (see below)
* `polar_wedge` : geometric boundary condition away from the pole, available only for the x2-direction in 2D and 3D (partial wedge in `x3` only) spherical-polar and similar coordinates (see below)
* `user` : user-defined boundary condition

For cylindrical- and spherical-like coordinates one often has the azimuthal angle run from 0 to 2π with periodic boundaries. In addition, spherical-like coordinates often have the polar angle run from 0 to π with polar boundaries. See the next section for restrictions on polar boundary conditions.

### Polar Boundary Conditions
While Athena++ supports having a reflecting wall or even outflow conditions along the polar axes in spherical-like coordinates, it also supports physically connecting the domain across the axes, allowing fluid to flow from one side to the other with no artificial boundary. This is done by specifying either `polar` or `polar_wedge` as the boundary condition for the polar (`x2`) angle. The following restrictions affect these modes:
* `polar`:
  * The domain must be 3D.
  * The azimuthal angle `x3` must extend from 0 to 2π (6.283185307179586), and it must have both inner and outer boundaries be periodic.
* `polar_wedge`:
  * The domain can either be 2D or 3D.
  * If the domain is 3D, `polar_wedge` should only be used if and only if the Mesh's azimuthal angle `x3` does NOT extend from exactly 0 to 2π (6.283185307179586). Otherwise, use the `polar` flag for either or both domain coordinate limits boundaries `ix3_bc`, `ox3_bc`.
  * Furthermore, in 3D each MeshBlock should span an azimuthal length of `dx3`= π/n, for some integer `n`. Therefore, uniform [Mesh Spacing](Coordinate-Systems-and-Meshes#mesh-spacing) along `x3` is required. 
* The polar angle `x2` must extend to exactly 0 if connecting across the inner (North) pole, and to exactly π (3.141592653589793) if connecting across the outer (South) pole.
* There must be an even number of cells in the azimuthal direction.
* The number of MeshBlocks in the azimuthal direction must be either 1 or an even number.
* All MeshBlocks surrounding a given pole at the same radius must be at the same level of refinement.

Polar boundaries work with [[Static Mesh Refinement]], but the user should double-check that the refinement scheme does not violate the same-level requirement. Currently, [[Adaptive Mesh Refinement]] does not protect against the grid evolving to an invalid refinement state.

### User-Defined Boundary Conditions
In order to use the user-defined boundary conditions, they must be defined and then enrolled in `Mesh::InitUserMeshData()` in the problem generator file.

First, the boundary condition function must be defined and match the following function signature:
```c++
void MyBoundary_ix1(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh);
```
The custom function can then be enrolled by calling:
```c++
EnrollUserBoundaryFunction(BoundaryFace::inner_x1, MyBoundary_ix1);
```
In Athena++, the boundary conditions are applied to both the cell-centered primitive variables and the face-centered magnetic fields. The conservative variables and cell-centered magnetic fields are automatically calculated from those updated quantities.

The following tables list the loop limits for the boundary cells that should be filled. Each table specifies a physical variable, and each row corresponds to a boundary direction. For the loop limits (`il`, `iu`, etc.) and the number of the ghost cells (`ngh`), the local values provided as the function parameters must be used instead of `pmb->is`, `pmb->ie`, etc. or `NGHOST`, because these values are adjusted accordingly for mesh refinement.

| Primitives | i-lower | i-upper | j-lower | j-upper | k-lower | k-upper |
| ----:|:----:|:----:|:----:|:----:|:----:|:----:|
| ix1  | il-ngh | il-1 | jl   | ju   | kl   | ku   |
| ox1  | iu+1 | iu+ngh | jl   | ju   | kl   | ku   |
| ix2  | il   | iu   | jl-ngh | jl-1 | kl   | ku   |
| ox2  | il   | iu   | ju+1 | ju+ngh | kl   | ku   |
| ix3  | il   | iu   | jl   | ju   | kl-ngh | kl-1 |
| ox3  | il   | iu   | jl   | ju   | ku+1 | ku+ngh |


| x1 B-fields | i-lower | i-upper | j-lower | j-upper | k-lower | k-upper |
| ----:|:----:|:----:|:----:|:----:|:----:|:----:|
| ix1  | il-ngh | il-1 | jl   | ju   | kl   | ku   |
| ox1  | iu+2 | iu+ngh+1 | jl   | ju   | kl   | ku   |
| ix2  | il   | iu+1 | jl-ngh | jl-1 | kl   | ku   |
| ox2  | il   | iu+1 | ju+1 | ju+ngh | kl   | ku   |
| ix3  | il   | iu+1 | jl   | ju   | kl-ngh | kl-1 |
| ox3  | il   | iu+1 | jl   | ju   | ku+1 | ku+ngh |


| x2 B-fields | i-lower | i-upper | j-lower | j-upper | k-lower | k-upper |
| ----:|:----:|:----:|:----:|:----:|:----:|:----:|
| ix1  | il-ngh | il-1 | jl   | ju+1 | kl   | ku   |
| ox1  | iu+1 | iu+ngh | jl   | ju+1 | kl   | ku   |
| ix2  | il   | iu   | jl-ngh | jl-1 | kl   | ku   |
| ox2  | il   | iu   | ju+2 | ju+ngh+1 | kl   | ku   |
| ix3  | il   | iu   | jl   | ju+1 | kl-ngh | kl-1 |
| ox3  | il   | iu   | jl   | ju+1 | ku+1 | ku+ngh |


| x3 B-fields | i-lower | i-upper | j-lower | j-upper | k-lower | k-upper |
| ----:|:----:|:----:|:----:|:----:|:----:|:----:|
| ix1  | il-ngh | il-1 | jl   | ju   | kl   | ku+1 |
| ox1  | iu+1 | iu+ngh | jl   | ju   | kl   | ku+1 |
| ix2  | il   | iu   | jl-ngh | jl   | kl   | ku+1 |
| ox2  | il   | iu   | ju+1 | ju+ngh | kl   | ku+1 |
| ix3  | il   | iu   | jl   | ju   | kl-ngh | kl-1 |
| ox3  | il   | iu   | jl   | ju   | ku+2 | ku+ngh+1 |


For details, see `src/bvals/outflow.cpp` (or `reflect.cpp`) and the Programmer Guide. Note that the `time` and `dt` parameters are intended for calculating time-dependent boundary conditions.
