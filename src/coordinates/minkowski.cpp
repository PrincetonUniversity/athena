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
Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin)
{
  // Set pointer to host MeshBlock
  pmy_block = pmb;

  // Set face centered positions and distances
  AllocateAndSetBasicCoordinates();

  // Initialize volume-averaged positions and spacings: x-direction
  for (int i = pmb->is-NGHOST; i <= pmb->ie+NGHOST; ++i)
    x1v(i) = 0.5 * (x1f(i) + x1f(i+1));
  for (int i = pmb->is-NGHOST; i <= pmb->ie+NGHOST-1; ++i)
    dx1v(i) = x1v(i+1) - x1v(i);

  // Initialize volume-averaged positions and spacings: y-direction
  if (pmb->block_size.nx2 == 1)  // no extent
  {
    x2v(pmb->js) = 0.5 * (x2f(pmb->js) + x2f(pmb->js+1));
    dx2v(pmb->js) = dx2f(pmb->js);
  }
  else  // extended
  {
    for (int j = pmb->js-NGHOST; j <= pmb->je+NGHOST; ++j)
      x2v(j) = 0.5 * (x2f(j) + x2f(j+1));
    for (int j = pmb->js-NGHOST; j <= pmb->je+NGHOST-1; ++j)
      dx2v(j) = x2v(j+1) - x2v(j);
  }

  // Initialize volume-averaged positions and spacings: z-direction
  if (pmb->block_size.nx3 == 1)  // no extent
  {
    x3v(pmb->ks) = 0.5 * (x3f(pmb->ks) + x3f(pmb->ks+1));
    dx3v(pmb->ks) = dx3f(pmb->ks);
  }
  else  // extended
  {
    for (int k = pmb->ks-NGHOST; k <= pmb->ke+NGHOST; ++k)
      x3v(k) = 0.5 * (x3f(k) + x3f(k+1));
    for (int k = pmb->ks-NGHOST; k <= pmb->ke+NGHOST-1; ++k)
      dx3v(k) = x3v(k+1) - x3v(k);
  }

  // Calculate coarse coordinates
  if (pmb->pmy_mesh->multilevel == true)
  {
    int cis = pmb->cis;
    int cie = pmb->cie;
    int cjs = pmb->cjs;
    int cje = pmb->cje;
    int cks = pmb->cks;
    int cke = pmb->cke;
    for (int i = cis-(pmb->cnghost); i <= cie+(pmb->cnghost); ++i)
      coarse_x1v(i) = 0.5 * (coarse_x1f(i+1) + coarse_x1f(i));
    for (int i = cis-(pmb->cnghost); i <= cie+(pmb->cnghost)-1; ++i)
      coarse_dx1v(i) = coarse_x1v(i+1) - coarse_x1v(i);
    if (pmb->block_size.nx2 == 1)
    {
      coarse_x2v(cjs) = 0.5 * (coarse_x2f(cjs+1) + coarse_x2f(cjs));
      coarse_dx2v(cjs) = coarse_dx2f(cjs);
    }
    else
    {
      for (int j = cjs-(pmb->cnghost); j <= cje+(pmb->cnghost); ++j)
        coarse_x2v(j) = 0.5 * (coarse_x2f(j+1) + coarse_x2f(j));
      for (int j = cjs-(pmb->cnghost); j <= cje+(pmb->cnghost)-1; ++j)
        coarse_dx2v(j) = coarse_x2v(j+1) - coarse_x2v(j);
    }
    if (pmb->block_size.nx3 == 1)
    {
      coarse_x3v(cks) = 0.5 * (coarse_x3f(cks+1) + coarse_x3f(cks));
      coarse_dx3v(cks) = coarse_dx3f(cks);
    }
    else
    {
      for (int k = cks-(pmb->cnghost); k <= cke+(pmb->cnghost); ++k)
        coarse_x3v(k) = 0.5 * (coarse_x3f(k+1) + coarse_x3f(k));
      for (int k = cks-(pmb->cnghost); k <= cke+(pmb->cnghost)-1; ++k)
        coarse_dx3v(k) = coarse_x3v(k+1) - coarse_x3v(k);
    }
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
    volumes(i) = dx1f(i) * dx2f(j) * dx3f(k);
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
    areas(i) = dx2f(j) * dx3f(k);
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
void Coordinates::Edge2Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &lengths)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    lengths(i) = dx2f(j);
  return;
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
void Coordinates::Edge3Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &lengths)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    lengths(i) = dx3f(k);
  return;
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

// Function for computing source terms using x-fluxes
// Inputs:
//   k,j: z- and y-indices
//   dt: size of timestep
//   flux: 1D array of x-fluxes
//   prim: 3D array of primitive values at beginning of half timestep
//   bb_cc: 3D array of cell-centered magnetic fields
// Outputs:
//   cons: source terms added to k,j-slice of 3D array of conserved variables
// Notes:
//   source terms all vanish identically
void Coordinates::CoordSrcTermsX1(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flux, const AthenaArray<Real> &prim,
  const AthenaArray<Real> &bb_cc, AthenaArray<Real> &cons)
{
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing source terms using y-fluxes
// Inputs:
//   k,j: z- and y-indices
//   dt: size of timestep
//   flux_j: 1D array of y-fluxes left of cells j
//   flux_jp1: 1D array of y-fluxes right of cells j
//   prim: 3D array of primitive values at beginning of half timestep
//   bcc: 3D array of cell-centered magnetic fields
// Outputs:
//   cons: source terms added to k,j-slice of 3D array of conserved variables
// Notes:
//   source terms all vanish identically
void Coordinates::CoordSrcTermsX2(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flux_j, const AthenaArray<Real> &flux_jp1,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing source terms using z-fluxes
// Inputs:
//   k,j: z- and y-indices
//   dt: size of timestep
//   flux_k: 2D array of z-fluxes left of cells k
//   flux_kp1: 2D array of z-fluxes right of cells k
//   prim: 3D array of primitive values at beginning of half timestep
//   bcc: 3D array of cell-centered magnetic fields
// Outputs:
//   cons: source terms added to k,j-slice of 3D array of conserved variables
// Notes:
//   source terms all vanish identically
void Coordinates::CoordSrcTermsX3(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flux_k, const AthenaArray<Real> &flux_kp1,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
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
void Coordinates::MinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3,
    Real *pt, Real *px, Real *py, Real *pz)
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
void Coordinates::TransformVectorCell(
    Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
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
void Coordinates::TransformVectorFace1(
    Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
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
void Coordinates::TransformVectorFace2(
    Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
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
void Coordinates::TransformVectorFace3(
    Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = at;
  *pa1 = ax;
  *pa2 = ay;
  *pa3 = az;
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
