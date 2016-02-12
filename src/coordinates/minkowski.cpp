// Minkowski spacetime, Minkowski (Cartesian) coordinates
// Notes:
//   coordinates: t, x, y, z
//   metric: ds^2 = -dt^2 + dx^2 + dy^2 + dz^2

// Primary header
#include "coordinates.hpp"

// C++ headers
#include <cmath>  // sqrt()

// Athena headers
#include "../athena.hpp"           // enums, macros, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../mesh.hpp"             // MeshBlock

//--------------------------------------------------------------------------------------

// Constructor
// Inputs:
//   pmb: pointer to block containing this grid
//   pin: pointer to runtime inputs
Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin, int flag)
{
  // Set pointer to host MeshBlock and note active zone boundaries
  pmy_block = pmb;
  cflag=flag;
  int is, ie, js, je, ks, ke, ng;
  if(cflag==0) {
    is = pmb->is; js = pmb->js; ks = pmb->ks;
    ie = pmb->ie; je = pmb->je; ke = pmb->ke;
    ng=NGHOST;
  }
  else {
    is = pmb->cis; js = pmb->cjs; ks = pmb->cks;
    ie = pmb->cie; je = pmb->cje; ke = pmb->cke;
    ng=pmb->cnghost;
  }

  // Set face-centered positions and distances
  AllocateAndSetBasicCoordinates();

  // Initialize volume-averaged positions and spacings: x-direction
  for (int i = is-NGHOST; i <= ie+NGHOST; ++i)
    x1v(i) = 0.5 * (x1f(i) + x1f(i+1));
  for (int i = is-NGHOST; i <= ie+NGHOST-1; ++i)
    dx1v(i) = x1v(i+1) - x1v(i);

  // Initialize volume-averaged positions and spacings: y-direction
  if (pmb->block_size.nx2 == 1)  // no extent
  {
    x2v(js) = 0.5 * (x2f(js) + x2f(js+1));
    dx2v(js) = dx2f(js);
  }
  else  // extended
  {
    for (int j = js-NGHOST; j <= je+NGHOST; ++j)
      x2v(j) = 0.5 * (x2f(j) + x2f(j+1));
    for (int j = js-NGHOST; j <= je+NGHOST-1; ++j)
      dx2v(j) = x2v(j+1) - x2v(j);
  }

  // Initialize volume-averaged positions and spacings: z-direction
  if (pmb->block_size.nx3 == 1)  // no extent
  {
    x3v(ks) = 0.5 * (x3f(ks) + x3f(ks+1));
    dx3v(ks) = dx3f(ks);
  }
  else  // extended
  {
    for (int k = ks-NGHOST; k <= ke+NGHOST; ++k)
      x3v(k) = 0.5 * (x3f(k) + x3f(k+1));
    for (int k = ks-NGHOST; k <= ke+NGHOST-1; ++k)
      dx3v(k) = x3v(k+1) - x3v(k);
  }

  // Prepare for MHD mesh refinement
  if (pmb->pmy_mesh->multilevel == true && MAGNETIC_FIELDS_ENABLED)
  {
    for (int i = is-NGHOST; i <= ie+NGHOST; ++i)
      x1s2(i) = x1s3(i) = x1v(i);
    if (pmb->block_size.nx2 == 1)
      x2s1(js) = x2s3(js) = x2v(js);
    else
      for (int j = js-NGHOST; j <= je+NGHOST; ++j)
        x2s1(j) = x2s3(j) = x2v(j);
    if (pmb->block_size.nx3 == 1)
      x3s1(ks) = x3s2(ks) = x3v(ks);
    else
      for (int k = ks-NGHOST; k <= ke+NGHOST; ++k)
        x3s1(k) = x3s2(k) = x3v(k);
  }
}

//--------------------------------------------------------------------------------------

// Destructor
Coordinates::~Coordinates()
{
  DeleteBasicCoordinates();
}

//--------------------------------------------------------------------------------------

// Function for computing cell volumes
// Inputs:
//   k: z-index
//   j: y-index
//   il,iu: x-index bounds
// Outputs:
//   volumes: 1D array of cell volumes
// Notes:
//   \Delta V = \Delta x * \Delta y * \Delta z
//   cf. GetCellVolume()
void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &volumes)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    volumes(i) = GetCellVolume(k, j, i);
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing single cell volume
// Inputs:
//   k,j,i: z-, y-, and x-indices
// Outputs:
//   returned value: cell volume
// Notes:
//   \Delta V = \Delta x * \Delta y * \Delta z
//   cf. CellVolume()
Real Coordinates::GetCellVolume(const int k, const int j, const int i)
{
  return dx1f(i) * dx2f(j) * dx3f(k);
}

//--------------------------------------------------------------------------------------

// Function for computing areas orthogonal to x
// Inputs:
//   k: z-index
//   j: y-index
//   il,iu: x-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to x
// Notes:
//   \Delta A = \Delta y * \Delta z
//   cf. GetFace1Area()
void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    areas(i) = GetFace1Area(k, j, i);
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing single area orthogonal to x
// Inputs:
//   k,j,i: z-, y-, and x-indices
// Outputs:
//   returned value: interface area orthogonal to x
// Notes:
//   \Delta A = \Delta y * \Delta z
//   cf. Face1Area()
Real Coordinates::GetFace1Area(const int k, const int j, const int i)
{
  return dx2f(j) * dx3f(k);
}

//--------------------------------------------------------------------------------------

// Function for computing areas orthogonal to y
// Inputs:
//   k: z-index
//   j: y-index
//   il,iu: x-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to y
// Notes:
//   \Delta A = \Delta x * \Delta z
void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    areas(i) = dx1f(i) * dx3f(k);
  return;
}


//--------------------------------------------------------------------------------------

// Function for computing areas orthogonal to z
// Inputs:
//   k: z-index
//   j: y-index
//   il,iu: x-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to z
// Notes:
//   \Delta A = \Delta x * \Delta y
void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    areas(i) = dx1f(i) * dx2f(j);
  return;
}


//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the x-direction
// Inputs:
//   k: z-index (unused)
//   j: y-index (unused)
//   il,iu: x-index bounds
// Outputs:
//   lengths: 1D array of edge lengths along x
// Notes:
//   \Delta L = \Delta x
void Coordinates::Edge1Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &lengths)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    lengths(i) = dx1f(i);
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the y-direction
// Inputs:
//   k: z-index (unused)
//   j: y-index
//   il,iu: x-index bounds
// Outputs:
//   lengths: 1D array of edge lengths along y
// Notes:
//   \Delta L = \Delta y
//   cf. GetEdge2Length()
void Coordinates::Edge2Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &lengths)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    lengths(i) = GetEdge2Length(k, j, i);
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing single length of edge in the y-direction
// Inputs:
//   k,i: z- and x-indices (unused)
//   j: y-index
// Outputs:
//   returned value: length of edge along y
// Notes:
//   \Delta L = \Delta y
//   cf. Edge2Length()
Real Coordinates::GetEdge2Length(const int k, const int j, const int i)
{
  return dx2f(j);
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the z-direction
// Inputs:
//   k: z-index
//   j: y-index (unused)
//   il,iu: x-index bounds
// Outputs:
//   lengths: 1D array of edge lengths along z
// Notes:
//   \Delta L = \Delta z
//   cf. GetEdge3Length()
void Coordinates::Edge3Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &lengths)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    lengths(i) = GetEdge3Length(k, j, i);
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing single length of edge in the z-direction
// Inputs:
//   k: z-index
//   j,i: y- and x-indices (unused)
// Outputs:
//   returned value: length of edge along z
// Notes:
//   \Delta L = \Delta z
//   cf. Edge3Length()
Real Coordinates::GetEdge3Length(const int k, const int j, const int i)
{
  return dx3f(k);
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the x-direction
// Inputs:
//   k: z-index (unused)
//   j: y-index (unused)
//   i: x-index
// Outputs:
//   returned value: x-width of cell (i,j,k)
// Notes:
//   \Delta W = \Delta x
Real Coordinates::CenterWidth1(const int k, const int j, const int i)
{
  return dx1f(i);
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the y-direction
// Inputs:
//   k: z-index (unused)
//   j: y-index
//   i: x-index (unused)
// Outputs:
//   returned value: y-width of cell (i,j,k)
// Notes:
//   \Delta W = \Delta y
Real Coordinates::CenterWidth2(const int k, const int j, const int i)
{
  return dx2f(j);
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the z-direction
// Inputs:
//   k: z-index
//   j: y-index (unused)
//   i: x-index (unused)
// Outputs:
//   returned value: z-width of cell (i,j,k)
// Notes:
//   \Delta W = \Delta z
Real Coordinates::CenterWidth3(const int k, const int j, const int i)
{
  return dx3f(k);
}

//--------------------------------------------------------------------------------------

// Function for computing source terms using fluxes
// Inputs:
//   dt: size of timestep
//   flux: 1D array of x-fluxes
//   prim: 3D array of primitive values at beginning of half timestep
//   bb_cc: 3D array of cell-centered magnetic fields
// Outputs:
//   cons: source terms added to 3D array of conserved variables
// Notes:
//   source terms all vanish identically
void Coordinates::CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bb_cc,
    AthenaArray<Real> &cons)
{
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing cell-centered metric coefficients
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::CellMetric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    g(I00,i) = -1.0;
    g(I11,i) = 1.0;
    g(I22,i) = 1.0;
    g(I33,i) = 1.0;
    g_inv(I00,i) = -1.0;
    g_inv(I11,i) = 1.0;
    g_inv(I22,i) = 1.0;
    g_inv(I33,i) = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: x-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face1Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    g(I00,i) = -1.0;
    g(I11,i) = 1.0;
    g(I22,i) = 1.0;
    g(I33,i) = 1.0;
    g_inv(I00,i) = -1.0;
    g_inv(I11,i) = 1.0;
    g_inv(I22,i) = 1.0;
    g_inv(I33,i) = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: y-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face2Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    g(I00,i) = -1.0;
    g(I11,i) = 1.0;
    g(I22,i) = 1.0;
    g(I33,i) = 1.0;
    g_inv(I00,i) = -1.0;
    g_inv(I11,i) = 1.0;
    g_inv(I22,i) = 1.0;
    g_inv(I33,i) = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: z-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face3Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    g(I00,i) = -1.0;
    g(I11,i) = 1.0;
    g(I22,i) = 1.0;
    g(I33,i) = 1.0;
    g_inv(I00,i) = -1.0;
    g_inv(I11,i) = 1.0;
    g_inv(I22,i) = 1.0;
    g_inv(I33,i) = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: x-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   bb1: 3D array of normal components B^1 of magnetic field, in global coordinates
//   prim_l: 1D array of left primitives, using global coordinates
//   prim_r: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
//   bbx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   transformation is trivial
void Coordinates::PrimToLocal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb1, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &bbx)
{
  if (MAGNETIC_FIELDS_ENABLED)
  {
    #pragma simd
    for (int i = il; i <= iu; ++i)
      bbx(i) = bb1(k,j,i);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: y-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   bb2: 3D array of normal components B^2 of magnetic field, in global coordinates
//   prim_l: 1D array of left primitives, using global coordinates
//   prim_r: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
//   bbx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   transformation is trivial
void Coordinates::PrimToLocal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb2, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &bbx)
{
  if (MAGNETIC_FIELDS_ENABLED)
  {
    #pragma simd
    for (int i = il; i <= iu; ++i)
      bbx(i) = bb2(k,j,i);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: z-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   bb3: 3D array of normal components B^3 of magnetic field, in global coordinates
//   prim_l: 1D array of left primitives, using global coordinates
//   prim_r: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
//   bbx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   transformation is trivial
void Coordinates::PrimToLocal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb3, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &bbx)
{
  if (MAGNETIC_FIELDS_ENABLED)
  {
    #pragma simd
    for (int i = il; i <= iu; ++i)
      bbx(i) = bb3(k,j,i);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: x-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates (unused)
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   transformation is trivial except for sign change from lowering time index
void Coordinates::FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &txt = flux(IEN,i);
    Real &t10 = flux(IEN,i);
    t10 = -txt;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: y-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates (unused)
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   transformation is trivial except for sign change from lowering time index
void Coordinates::FluxToGlobal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &tyt = flux(IEN,i);
    Real &t20 = flux(IEN,i);
    t20 = -tyt;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: z-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates (unused)
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   transformation is trivial except for sign change from lowering time index
void Coordinates::FluxToGlobal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &tzt = flux(IEN,i);
    Real &t30 = flux(IEN,i);
    t30 = -tzt;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for calculating distance between two points
// Inputs:
//   a1,a2,a3: global coordinates of first point
//   bx,by,bz: Minkowski coordinates of second point
// Outputs:
//   returned value: Euclidean distance between a and b
Real Coordinates::DistanceBetweenPoints(Real a1, Real a2, Real a3, Real bx, Real by,
    Real bz)
{
  Real ax = a1;
  Real ay = a2;
  Real az = a3;
  return std::sqrt(SQR(ax-bx) + SQR(ay-by) + SQR(az-bz));
}

//--------------------------------------------------------------------------------------

// Function for calculating Minkowski coordinates of cell
// Inputs:
//   x0,x1,x2,x3: Minkowski coordinates
// Outputs:
//   pt,px,py,pz: Minkowski coordinate values set
// Notes:
//   transformation is trivial
void Coordinates::MinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3, Real *pt,
    Real *px, Real *py, Real *pz)
{
  *pt = x0;
  *px = x1;
  *py = x2;
  *pz = x3;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: cell-centered
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j,i: z-, y-, and x-indices (unused)
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial
void Coordinates::TransformVectorCell(Real at, Real ax, Real ay, Real az, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = at;
  *pa1 = ax;
  *pa2 = ay;
  *pa3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: x-interface
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j,i: z-, y-, and x-indices (unused)
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial
void Coordinates::TransformVectorFace1(Real at, Real ax, Real ay, Real az, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = at;
  *pa1 = ax;
  *pa2 = ay;
  *pa3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: y-interface
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j,i: z-, y-, and x-indices (unused)
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial
void Coordinates::TransformVectorFace2(Real at, Real ax, Real ay, Real az, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = at;
  *pa1 = ax;
  *pa2 = ay;
  *pa3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: z-interface
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j,i: z-, y-, and x-indices (unused)
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial
void Coordinates::TransformVectorFace3(Real at, Real ax, Real ay, Real az, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = at;
  *pa1 = ax;
  *pa2 = ay;
  *pa3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for raising covariant components of a vector
// Inputs:
//   a_0,a_1,a_2,a_3: covariant components of vector
//   k,j,i: indices of cell in which transformation is desired
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to contravariant 4-vector components
void Coordinates::RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = -a_0;
  *pa1 = a_1;
  *pa2 = a_2;
  *pa3 = a_3;
  return;
}

//--------------------------------------------------------------------------------------

// Function for lowering contravariant components of a vector
// Inputs:
//   a0,a1,a2,a3: contravariant components of vector
//   k,j,i: indices of cell in which transformation is desired
// Outputs:
//   pa_0,pa_1,pa_2,pa_3: pointers to covariant 4-vector components
void Coordinates::LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j,
    int i, Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3)
{
  *pa_0 = -a0;
  *pa_1 = a1;
  *pa_2 = a2;
  *pa_3 = a3;
  return;
}
