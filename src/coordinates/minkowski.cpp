// Minkowski spacetime, Minkowski (Cartesian) coordinates
// Notes:
//   coordinates: t, x, y, z
//   metric: ds^2 = -dt^2 + dx^2 + dy^2 + dz^2

// Primary header
#include "coordinates.hpp"

// C++ headers
#include <cmath>  // sqrt()

// Athena headers
#include "../athena.hpp"         // enums, macros, Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../mesh.hpp"           // MeshBlock

//--------------------------------------------------------------------------------------

// Constructor
// Inputs:
//   pb: pointer to block containing this grid
//   pin: pointer to runtime inputs
Coordinates::Coordinates(MeshBlock *pb, ParameterInput *pin)
{
  // Set pointer to host MeshBlock
  pmy_block = pb;

  // Initialize volume-averaged positions and spacings: x-direction
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST; ++i)
    pb->x1v(i) = 0.5 * (pb->x1f(i) + pb->x1f(i+1));
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST-1; ++i)
    pb->dx1v(i) = pb->x1v(i+1) - pb->x1v(i);

  // Initialize volume-averaged positions and spacings: y-direction
  if (pb->block_size.nx2 == 1)  // no extent
  {
    pb->x2v(pb->js) = 0.5 * (pb->x2f(pb->js) + pb->x2f(pb->js+1));
    pb->dx2v(pb->js) = pb->dx2f(pb->js);
  }
  else  // extended
  {
    for (int j = pb->js-NGHOST; j <= pb->je+NGHOST; ++j)
      pb->x2v(j) = 0.5 * (pb->x2f(j) + pb->x2f(j+1));
    for (int j = pb->js-NGHOST; j <= pb->je+NGHOST-1; ++j)
      pb->dx2v(j) = pb->x2v(j+1) - pb->x2v(j);
  }

  // Initialize volume-averaged positions and spacings: z-direction
  if (pb->block_size.nx3 == 1)  // no extent
  {
    pb->x3v(pb->ks) = 0.5 * (pb->x3f(pb->ks) + pb->x3f(pb->ks+1));
    pb->dx3v(pb->ks) = pb->dx3f(pb->ks);
  }
  else  // extended
  {
    for (int k = pb->ks-NGHOST; k <= pb->ke+NGHOST; ++k)
      pb->x3v(k) = 0.5 * (pb->x3f(k) + pb->x3f(k+1));
    for (int k = pb->ks-NGHOST; k <= pb->ke+NGHOST-1; ++k)
      pb->dx3v(k) = pb->x3v(k+1) - pb->x3v(k);
  }
}

//--------------------------------------------------------------------------------------

// Destructor
Coordinates::~Coordinates()
{
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
void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &volumes)
{
  const Real &delta_y = pmy_block->dx2f(j);
  const Real &delta_z = pmy_block->dx3f(k);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &delta_x = pmy_block->dx1f(i);
    Real &volume = volumes(i);
    volume = delta_x * delta_y * delta_z;
  }
  return;
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
void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  const Real &delta_y = pmy_block->dx2f(j);
  const Real &delta_z = pmy_block->dx3f(k);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    Real &area = areas(i);
    area = delta_y * delta_z;
  }
  return;
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
  const Real &delta_z = pmy_block->dx3f(k);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &delta_x = pmy_block->dx1f(i);
    Real &area = areas(i);
    area = delta_x * delta_z;
  }
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
  const Real &delta_y = pmy_block->dx2f(j);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &delta_x = pmy_block->dx1f(i);
    Real &area = areas(i);
    area = delta_x * delta_y;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the x-direction
// Inputs:
//   k: z-index (unused)
//   j: y-index (unused)
//   il,iu: x-index bounds
// Outputs:
//   len: 1D array of edge lengths along x
// Notes:
//   \Delta L = \Delta x
void Coordinates::Edge1Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    len(i) = pmy_block->dx1f(i);
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the y-direction
// Inputs:
//   k: z-index (unused)
//   j: y-index
//   il,iu: x-index bounds
// Outputs:
//   len: 1D array of edge lengths along y
// Notes:
//   \Delta L = \Delta y
void Coordinates::Edge2Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len)
{
  const Real &length = pmy_block->dx2f(j);
  #pragma simd
  for (int i = il; i <= iu; ++i)
    len(i) = length;
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the z-direction
// Inputs:
//   k: z-index
//   j: y-index (unused)
//   il,iu: x-index bounds
// Outputs:
//   len: 1D array of edge lengths along z
// Notes:
//   \Delta L = \Delta z
void Coordinates::Edge3Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &len)
{
  const Real &length = pmy_block->dx3f(k);
  #pragma simd
  for (int i = il; i <= iu; ++i)
    len(i) = length;
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the x-direction
// Inputs:
//   k: z-index (unused)
//   j: y-index (unused)
//   i: x-index
// Outputs:
//   returned value: width of cell (i,j,k)
// Notes:
//   \Delta W = \Delta x
Real Coordinates::CenterWidth1(const int k, const int j, const int i)
{
  return pmy_block->dx1f(i);
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the y-direction
// Inputs:
//   k: z-index (unused)
//   j: y-index
//   i: x-index (unused)
// Outputs:
//   returned value: width of cell (i,j,k)
// Notes:
//   \Delta W = \Delta y
Real Coordinates::CenterWidth2(const int k, const int j, const int i)
{
  return pmy_block->dx2f(j);
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the z-direction
// Inputs:
//   k: z-index
//   j: y-index (unused)
//   i: x-index (unused)
// Outputs:
//   returned value: width of cell (i,j,k)
// Notes:
//   \Delta W = \Delta z
Real Coordinates::CenterWidth3(const int k, const int j, const int i)
{
  return pmy_block->dx3f(k);
}

//--------------------------------------------------------------------------------------

// Function for computing source terms
// Inputs:
//   dt: size of timestep
//   prim: full grid of primitive values at beginning of half timestep
//   cons: full grid of conserved variables at end of half timestep
// Outputs:
//   cons: source terms added
// Notes:
//   source terms all vanish identically
void Coordinates::CoordinateSourceTerms(Real dt, const AthenaArray<Real> &prim,
    AthenaArray<Real> &cons)
{
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing cell-centered metric coefficients
// Inputs:
//   k: z-index
//   j: y-index
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::CellMetric(const int k, const int j, AthenaArray<Real> &g,
    AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = pmy_block->is-NGHOST; i <= pmy_block->ie+NGHOST; ++i)
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

// Function for computing face-centered metric coefficients: r-interface
// Inputs:
//   k: phi-index
//   j: theta-index
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face1Metric(const int k, const int j, AthenaArray<Real> &g,
    AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie+1; i++)
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

// Function for computing face-centered metric coefficients: theta-interface
// Inputs:
//   k: phi-index
//   j: theta-index
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face2Metric(const int k, const int j, AthenaArray<Real> &g,
    AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
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

// Function for computing face-centered metric coefficients: phi-interface
// Inputs:
//   k: phi-index
//   j: theta-index
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face3Metric(const int k, const int j, AthenaArray<Real> &g,
    AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
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
//   k: z-index
//   j: y-index
//   b1_vals: 3D array of normal components B^1 of magnetic field, in global coordinates
//   prim_left: 1D array of left primitives, using global coordinates
//   prim_right: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_left: values overwritten in local coordinates
//   prim_right: values overwritten in local coordinates
//   bx: 1D array of longitudinal magnetic fields, in local coordinates
// Notes:
//   transformation is trivial
void Coordinates::PrimToLocal1(const int k, const int j,
    const AthenaArray<Real> &b1_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx)
{
  if (MAGNETIC_FIELDS_ENABLED)
  {
    #pragma simd
    for (int i = pmy_block->is; i <= pmy_block->ie+1; i++)
      bx(i) = b1_vals(k,j,i);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: y-interface
// Inputs:
//   k: z-index
//   j: y-index
//   b2_vals: 3D array of normal components B^2 of magnetic field, in global coordinates
//   prim_left: 1D array of left primitives, using global coordinates
//   prim_right: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_left: values overwritten in local coordinates
//   prim_right: values overwritten in local coordinates
//   by: 1D array of longitudinal magnetic fields, in local coordinates
// Notes:
//   transformation is trivial
void Coordinates::PrimToLocal2(const int k, const int j,
    const AthenaArray<Real> &b2_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &by)
{
  if (MAGNETIC_FIELDS_ENABLED)
  {
    #pragma simd
    for (int i = pmy_block->is; i <= pmy_block->ie; i++)
      by(i) = b2_vals(k,j,i);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: z-interface
// Inputs:
//   k: z-index
//   j: y-index
//   b3_vals: 3D array of normal components B^3 of magnetic field, in global coordinates
//   prim_left: 1D array of left primitives, using global coordinates
//   prim_right: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_left: values overwritten in local coordinates
//   prim_right: values overwritten in local coordinates
//   bz: 1D array of longitudinal magnetic fields, in local coordinates
// Notes:
//   transformation is trivial
void Coordinates::PrimToLocal3(const int k, const int j,
    const AthenaArray<Real> &b3_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bz)
{
  if (MAGNETIC_FIELDS_ENABLED)
  {
    #pragma simd
    for (int i = pmy_block->is; i <= pmy_block->ie; i++)
      bz(i) = b3_vals(k,j,i);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: x-interface
// Inputs:
//   k: z-index
//   j: y-index
//   pflux: pointer to array of fluxes in 1D, using local coordinates
// Outputs:
//   pflux: pointer to values overwritten in global coordinates
// Notes:
//   transformation is trivial except for sign change from lowering time index
void Coordinates::FluxToGlobal1(const int k, const int j, AthenaArray<Real> &flux)
{
  #pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie+1; ++i)
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
//   k: z-index
//   j: y-index
//   pflux: pointer to array of fluxes in 1D, using local coordinates
// Outputs:
//   pflux: pointer to values overwritten in global coordinates
// Notes:
//   transformation is trivial except for sign change from lowering time index
void Coordinates::FluxToGlobal2(const int k, const int j, AthenaArray<Real> &flux)
{
  #pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; ++i)
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
//   k: z-index
//   j: y-index
//   pflux: pointer to array of fluxes in 1D, using local coordinates
// Outputs:
//   pflux: pointer to values overwritten in global coordinates
// Notes:
//   transformation is trivial except for sign change from lowering time index
void Coordinates::FluxToGlobal3(const int k, const int j, AthenaArray<Real> &flux)
{
  #pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; ++i)
  {
    const Real &tzt = flux(IEN,i);
    Real &t30 = flux(IEN,i);
    t30 = -tzt;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for converting all primitives to conserved variables
// Inputs:
//   prim: 3D array of primitives
//   b: 3D array of cell-centered magnetic fields
//   gamma_adi_red: \Gamma/(\Gamma-1) for ratio of specific heats \Gamma
// Outputs:
//   cons: 3D array of conserved variables
void Coordinates::PrimToCons(
    const AthenaArray<Real> &prim, const AthenaArray<Real> &b, Real gamma_adi_red,
    AthenaArray<Real> &cons)
{
  // Prepare index bounds
  int il = pmy_block->is - NGHOST;
  int iu = pmy_block->ie + NGHOST;
  int jl = pmy_block->js;
  int ju = pmy_block->je;
  if (pmy_block->block_size.nx2 > 1)
  {
    jl -= (NGHOST);
    ju += (NGHOST);
  }
  int kl = pmy_block->ks;
  int ku = pmy_block->ke;
  if (pmy_block->block_size.nx3 > 1)
  {
    kl -= (NGHOST);
    ku += (NGHOST);
  }

  // Go through all cells
  for (int k = kl; k <= ku; k++)
    for (int j = jl; j <= ju; j++)
    {
      #pragma simd
      for (int i = il; i <= iu; i++)
      {
        // Extract geometric quantities
        const Real g00 = -1.0;
        const Real g11 = 1.0;
        const Real g22 = 1.0;
        const Real g33 = 1.0;

        // Extract primitives
        const Real &rho = prim(IDN,k,j,i);
        const Real &pgas = prim(IEN,k,j,i);
        const Real &v1 = prim(IVX,k,j,i);
        const Real &v2 = prim(IVY,k,j,i);
        const Real &v3 = prim(IVZ,k,j,i);

        // Extract magnetic fields
        Real b1 = 0.0, b2 = 0.0, b3 = 0.0;
        if (MAGNETIC_FIELDS_ENABLED)
        {
          b1 = b(IB1,k,j,i);
          b2 = b(IB2,k,j,i);
          b3 = b(IB3,k,j,i);
        }

        // Calculate 4-velocity
        Real u0 = std::sqrt(-1.0 / (g00 + g11*v1*v1 + g22*v2*v2 + g33*v3*v3));
        Real u1 = u0 * v1;
        Real u2 = u0 * v2;
        Real u3 = u0 * v3;
        Real u_0 = g00*u0;
        Real u_1 = g11*u1;
        Real u_2 = g22*u2;
        Real u_3 = g33*u3;

        // Calculate 4-magnetic field
        Real bcon0 = g11*b1*u1 + g22*b2*u2 + g33*b3*u3;
        Real bcon1 = 1.0/u0 * (b1 + bcon0 * u1);
        Real bcon2 = 1.0/u0 * (b2 + bcon0 * u2);
        Real bcon3 = 1.0/u0 * (b3 + bcon0 * u3);
        Real bcov0 = g00*bcon0;
        Real bcov1 = g11*bcon1;
        Real bcov2 = g22*bcon2;
        Real bcov3 = g33*bcon3;
        Real b_sq = bcon0*bcov0 + bcon1*bcov1 + bcon2*bcov2 + bcon3*bcov3;

        // Extract conserved quantities
        Real &rho_u0 = cons(IDN,k,j,i);
        Real &t0_0 = cons(IEN,k,j,i);
        Real &t0_1 = cons(IM1,k,j,i);
        Real &t0_2 = cons(IM2,k,j,i);
        Real &t0_3 = cons(IM3,k,j,i);

        // Set conserved quantities
        Real w = rho + gamma_adi_red * pgas + b_sq;
        Real ptot = pgas + 0.5*b_sq;
        rho_u0 = rho * u0;
        t0_0 = w * u0 * u_0 - bcon0 * bcov0 + ptot;
        t0_1 = w * u0 * u_1 - bcon0 * bcov1;
        t0_2 = w * u0 * u_2 - bcon0 * bcov2;
        t0_3 = w * u0 * u_3 - bcon0 * bcov3;
      }
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

// Function for transforming 4-vector from Minkowski to global: cell-centered
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j: z- and y-indices (unused)
//   i: x-index
// Outputs:
//   a0,a1,a2,a3: upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial
void Coordinates::TransformVectorCell(Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *a0, Real *a1, Real *a2, Real *a3)
{
  *a0 = at;
  *a1 = ax;
  *a2 = ay;
  *a3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: x-interface
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j: z- and y-indices (unused)
//   i: x-index
// Outputs:
//   a0,a1,a2,a3: upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial
void Coordinates::TransformVectorFace1(Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *a0, Real *a1, Real *a2, Real *a3)
{
  *a0 = at;
  *a1 = ax;
  *a2 = ay;
  *a3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: y-interface
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j: z- and y-indices (unused)
//   i: x-index
// Outputs:
//   a0,a1,a2,a3: upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial
void Coordinates::TransformVectorFace2(Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *a0, Real *a1, Real *a2, Real *a3)
{
  *a0 = at;
  *a1 = ax;
  *a2 = ay;
  *a3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: z-interface
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j: z- and y-indices (unused)
//   i: x-index
// Outputs:
//   a0,a1,a2,a3: upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial
void Coordinates::TransformVectorFace3(Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *a0, Real *a1, Real *a2, Real *a3)
{
  *a0 = at;
  *a1 = ax;
  *a2 = ay;
  *a3 = az;
  return;
}
