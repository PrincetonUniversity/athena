// Minkowski spacetime, tilted coordinates
// Notes:
//   coordinates: t', x', y, z
//   parameters: a (velocity)
//   metric:
//     ds^2 = -(1-a^2)/(1+a^2) dt'^2 + 2a/(1+a^2) (dt'dx' + dx'dt')
//        + (1-a^2)/(1+a^2) dx'^2 + dy^2 + dz^2
//   relation to Minkowski (t, x, y, z):
//     t' = (t + a x) / \sqrt{1 + a^2}
//     x' = (x - a t) / \sqrt{1 + a^2}

// Primary header
#include "coordinates.hpp"

// C++ headers
#include <cmath>  // NAN, sqrt()

// Athena headers
#include "../athena.hpp"           // enums, macros, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../fluid/eos/eos.hpp"    // FluidEqnOfState
#include "../fluid/fluid.hpp"      // Fluid
#include "../mesh.hpp"             // MeshBlock
#include "../parameter_input.hpp"  // ParameterInput

// Global variables
static Real alpha;  // \sqrt{1+a^2}
static Real beta;   // \sqrt{1-a^2}

//--------------------------------------------------------------------------------------

// Constructor
// Inputs:
//   pb: pointer to block containing this grid
//   pin: pointer to runtime inputs
Coordinates::Coordinates(MeshBlock *pb, ParameterInput *pin)
{
  // Set face centered positions and distances
  AllocateAndSetBasicCoordinates();

  // Set parameters
  tilted_a_ = pin->GetReal("coord", "a");
  const Real &a = tilted_a_;
  alpha = std::sqrt(1.0 + SQR(a));
  beta = std::sqrt(1.0 - SQR(a));

  // Set pointer to host MeshBlock
  pmy_block = pb;

  // Initialize volume-averaged positions and spacings: x'-direction
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST; ++i)
    x1v(i) = 0.5 * (x1f(i) + x1f(i+1));
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST-1; ++i)
    dx1v(i) = x1v(i+1) - x1v(i);

  // Initialize volume-averaged positions and spacings: y-direction
  if (pb->block_size.nx2 == 1)  // no extent
  {
    x2v(pb->js) = 0.5 * (x2f(pb->js) + x2f(pb->js+1));
    dx2v(pb->js) = dx2f(pb->js);
  }
  else  // extended
  {
    for (int j = pb->js-NGHOST; j <= pb->je+NGHOST; ++j)
      x2v(j) = 0.5 * (x2f(j) + x2f(j+1));
    for (int j = pb->js-NGHOST; j <= pb->je+NGHOST-1; ++j)
      dx2v(j) = x2v(j+1) - x2v(j);
  }

  // Initialize volume-averaged positions and spacings: z-direction
  if (pb->block_size.nx3 == 1)  // no extent
  {
    x3v(pb->ks) = 0.5 * (x3f(pb->ks) + x3f(pb->ks+1));
    dx3v(pb->ks) = dx3f(pb->ks);
  }
  else  // extended
  {
    for (int k = pb->ks-NGHOST; k <= pb->ke+NGHOST; ++k)
      x3v(k) = 0.5 * (x3f(k) + x3f(k+1));
    for (int k = pb->ks-NGHOST; k <= pb->ke+NGHOST-1; ++k)
      dx3v(k) = x3v(k+1) - x3v(k);
  }

  if(pmb->pmy_mesh->multilevel==true) { // calc coarse coodinates
    int cis = pmb->cis; int cjs = pmb->cjs; int cks = pmb->cks;
    int cie = pmb->cie; int cje = pmb->cje; int cke = pmb->cke;
    for (int i=cis-(pmb->cnghost); i<=cie+(pmb->cnghost); ++i) {
      coarse_x1v(i) = 0.5*(coarse_x1f(i+1) + coarse_x1f(i));
    }
    if (pmb->block_size.nx2 == 1) {
      coarse_x2v(cjs) = 0.5*(coarse_x2f(cjs+1) + coarse_x2f(cjs));
    } else {
      for (int j=cjs-(pmb->cnghost); j<=cje+(pmb->cnghost); ++j) {
        coarse_x2v(j) = 0.5*(coarse_x2f(j+1) + coarse_x2f(j));
      }
    }
    if (pmb->block_size.nx3 == 1) {
      coarse_x3v(cks) = 0.5*(coarse_x3f(cks+1) + coarse_x3f(cks));
    } else {
      for (int k=cks-(pmb->cnghost); k<=cke+(pmb->cnghost); ++k) {
        coarse_x3v(k) = 0.5*(coarse_x3f(k+1) + coarse_x3f(k));
      }
    }

    if (MAGNETIC_FIELDS_ENABLED) {
      for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i)
        x1s2(i) = x1s3(i) = x1v(i);
      for (int i=cis-(pmb->cnghost); i<=cie+(pmb->cnghost); ++i)
        coarse_x1s2(i) = coarse_x1s3(i) = coarse_x1v(i);
      if (pmb->block_size.nx2 == 1) {
        x2s1(js) = x2s3(js) = x2v(js);
        coarse_x2s1(js) = coarse_x2s3(js) = coarse_x2v(js);
      }
      else {
        for (int j=js-(NGHOST); j<=je+(NGHOST); ++j)
          x2s1(j) = x2s3(j) = x2v(j);
        for (int j=cjs-(pmb->cnghost); j<=cje+(pmb->cnghost); ++j)
          coarse_x2s1(j) = coarse_x2s3(j) = coarse_x2v(j);
      }
      if (pmb->block_size.nx3 == 1) {
        x3s1(ks) = x3s2(ks) = x3v(ks);
        coarse_x3s1(ks) = coarse_x3s2(ks) = coarse_x3v(ks);
      }
      else {
        for (int k=ks-(NGHOST); k<=ke+(NGHOST); ++k)
          x3s1(k) = x3s2(k) = x3v(k);
        for (int k=cks-(pmb->cnghost); k<=cke+(pmb->cnghost); ++k)
          coarse_x3s1(k) = coarse_x3s2(k) = coarse_x3v(k);
      }
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
//   il,iu: x'-index bounds
// Outputs:
//   volumes: 1D array of cell volumes
// Notes:
//   \Delta V = \Delta x' * \Delta y * \Delta z
void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &volumes)
{
  const Real &delta_y = dx2f(j);
  const Real &delta_z = dx3f(k);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &delta_x = dx1f(i);
    Real &volume = volumes(i);
    volume = delta_x * delta_y * delta_z;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing areas orthogonal to x'
// Inputs:
//   k: z-index
//   j: y-index
//   il,iu: x'-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to x'
// Notes:
//   \Delta A = \Delta y * \Delta z
void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  const Real &delta_y = dx2f(j);
  const Real &delta_z = dx3f(k);
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
//   il,iu: x'-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to y
// Notes:
//   \Delta A = \Delta x' * \Delta z
void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  const Real &delta_z = dx3f(k);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &delta_x = dx1f(i);
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
//   il,iu: x'-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to z
// Notes:
//   \Delta A = \Delta x' * \Delta y
void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  const Real &delta_y = dx2f(j);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &delta_x = dx1f(i);
    Real &area = areas(i);
    area = delta_x * delta_y;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the x'-direction
// Inputs:
//   k: z-index (unused)
//   j: y-index (unused)
//   il,iu: x'-index bounds
// Outputs:
//   lengths: 1D array of edge lengths along x'
// Notes:
//   \Delta L = \Delta x'
void Coordinates::Edge1Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &lengths)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &delta_x = dx1f(i);
    Real &length = lengths(i);
    length = delta_x;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the y-direction
// Inputs:
//   k: z-index (unused)
//   j: y-index
//   il,iu: x'-index bounds
// Outputs:
//   lengths: 1D array of edge lengths along y
// Notes:
//   \Delta L = \Delta y
void Coordinates::Edge2Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &lengths)
{
  const Real &delta_y = dx2f(j);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    Real &length = lengths(i);
    length = delta_y;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the z-direction
// Inputs:
//   k: z-index
//   j: y-index (unused)
//   il,iu: x'-index bounds
// Outputs:
//   lengths: 1D array of edge lengths along z
// Notes:
//   \Delta L = \Delta z
void Coordinates::Edge3Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &lengths)
{
  const Real &delta_z = dx3f(k);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    Real &length = lengths(i);
    length = delta_z;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the x'-direction
// Inputs:
//   k: z-index (unused)
//   j: y-index (unused)
//   i: x'-index
// Outputs:
//   returned value: x'-width of cell (i,j,k)
// Notes:
//   \Delta W = (\beta/\alpha) \Delta x'
Real Coordinates::CenterWidth1(const int k, const int j, const int i)
{
  return (beta/alpha) * dx1f(i);
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the y-direction
// Inputs:
//   k: z-index (unused)
//   j: y-index
//   i: x'-index (unused)
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
//   i: x'-index (unused)
// Outputs:
//   returned value: z-width of cell (i,j,k)
// Notes:
//   \Delta W = \Delta z
Real Coordinates::CenterWidth3(const int k, const int j, const int i)
{
  return dx3f(k);
}

//--------------------------------------------------------------------------------------

// Function for computing source terms using x'-fluxes
// Inputs:
//   k,j: z- and y-indices
//   dt: size of timestep
//   flux: 1D array of x'-fluxes
//   prim: 3D array of primitive values at beginning of half timestep
//   bcc: 3D array of cell-centered magnetic fields
// Outputs:
//   cons: source terms added to k,j-slice of 3D array of conserved variables
// Notes:
//   source terms all vanish identically
void Coordinates::CoordSrcTermsX1(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flux, const AthenaArray<Real> &prim,
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
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
    g(I00,i) = g_inv(I00,i) = -SQR(beta) / SQR(alpha);
    g(I01,i) = g_inv(I01,i) = 2.0*tilted_a_ / SQR(alpha);
    g(I11,i) = g_inv(I11,i) = SQR(beta) / SQR(alpha);
    g(I22,i) = g_inv(I22,i) = 1.0;
    g(I33,i) = g_inv(I33,i) = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: x'-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x'-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face1Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    g(I00,i) = g_inv(I00,i) = -SQR(beta) / SQR(alpha);
    g(I01,i) = g_inv(I01,i) = 2.0*tilted_a_ / SQR(alpha);
    g(I11,i) = g_inv(I11,i) = SQR(beta) / SQR(alpha);
    g(I22,i) = g_inv(I22,i) = 1.0;
    g(I33,i) = g_inv(I33,i) = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: y-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x'-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face2Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    g(I00,i) = g_inv(I00,i) = -SQR(beta) / SQR(alpha);
    g(I01,i) = g_inv(I01,i) = 2.0*tilted_a_ / SQR(alpha);
    g(I11,i) = g_inv(I11,i) = SQR(beta) / SQR(alpha);
    g(I22,i) = g_inv(I22,i) = 1.0;
    g(I33,i) = g_inv(I33,i) = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: z-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x'-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face3Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    g(I00,i) = g_inv(I00,i) = -SQR(beta) / SQR(alpha);
    g(I01,i) = g_inv(I01,i) = 2.0*tilted_a_ / SQR(alpha);
    g(I11,i) = g_inv(I11,i) = SQR(beta) / SQR(alpha);
    g(I22,i) = g_inv(I22,i) = 1.0;
    g(I33,i) = g_inv(I33,i) = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: x'-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x'-index bounds
//   b1_vals: 3D array of normal components B^1 of magnetic field, in global coordinates
//   prim_left: 1D array of left primitives, using global coordinates
//   prim_right: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_left: values overwritten in local coordinates
//   prim_right: values overwritten in local coordinates
//   bx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   expects v1/v2/v3 in IVX/IVY/IVZ slots
//   expects B1 in b1_vals
//   expects B2/B3 in IBY/IBZ slots
//   puts vx/vy/vz in IVX/IVY/IVZ
//   puts Bx in bx
//   puts By/Bz in IBY/IBZ slots
void Coordinates::PrimToLocal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b1_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx)
{
  // Prepare geometric quantities
  const Real g_00 = -SQR(beta) / SQR(alpha);
  const Real g_01 = 2.0*tilted_a_ / SQR(alpha);
  const Real g_10 = g_01;
  const Real g_11 = SQR(beta) / SQR(alpha);
  const Real g_22 = 1.0;
  const Real g_33 = 1.0;
  const Real mt_0 = alpha/beta;
  const Real mx_0 = 2.0*tilted_a_ / (alpha*beta);
  const Real mx_1 = beta/alpha;
  const Real my_2 = 1.0;
  const Real mz_3 = 1.0;

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract global 3-velocities
    Real &v1l = prim_left(IVX,i);
    Real &v2l = prim_left(IVY,i);
    Real &v3l = prim_left(IVZ,i);
    Real &v1r = prim_right(IVX,i);
    Real &v2r = prim_right(IVY,i);
    Real &v3r = prim_right(IVZ,i);

    // Construct global 4-velocities
    Real u0l = std::sqrt(-1.0 /
        (g_00 + 2.0*g_01*v1l + g_11*v1l*v1l + g_22*v2l*v2l + g_33*v3l*v3l));
    Real u1l = u0l * v1l;
    Real u2l = u0l * v2l;
    Real u3l = u0l * v3l;
    Real u0r = std::sqrt(-1.0 /
        (g_00 + 2.0*g_01*v1r + g_11*v1r*v1r + g_22*v2r*v2r + g_33*v3r*v3r));
    Real u1r = u0r * v1r;
    Real u2r = u0r * v2r;
    Real u3r = u0r * v3r;

    // Transform 4-velocities
    Real utl = mt_0*u0l;
    Real uxl = mx_0*u0l + mx_1*u1l;
    Real uyl = my_2*u2l;
    Real uzl = mz_3*u3l;
    Real utr = mt_0*u0r;
    Real uxr = mx_0*u0r + mx_1*u1r;
    Real uyr = my_2*u2r;
    Real uzr = mz_3*u3r;

    // Set local 3-velocities
    v1l = uxl / utl;
    v2l = uyl / utl;
    v3l = uzl / utl;
    v1r = uxr / utr;
    v2r = uyr / utr;
    v3r = uzr / utr;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      // Extract global magnetic fields
      const Real &bb1 = b1_vals(k,j,i);
      Real &bb2l = prim_left(IBY,i);
      Real &bb3l = prim_left(IBZ,i);
      Real &bb2r = prim_right(IBY,i);
      Real &bb3r = prim_right(IBZ,i);

      // Construct global contravariant magnetic fields
      Real b0l = g_10*bb1*u0l + g_11*bb1*u1l + g_22*bb2l*u2l + g_33*bb3l*u3l;
      Real b1l = (bb1 + b0l * u1l) / u0l;
      Real b2l = (bb2l + b0l * u2l) / u0l;
      Real b3l = (bb3l + b0l * u3l) / u0l;
      Real b0r = g_10*bb1*u0r + g_11*bb1*u1r + g_22*bb2r*u2r + g_33*bb3r*u3r;
      Real b1r = (bb1 + b0r * u1r) / u0r;
      Real b2r = (bb2r + b0r * u2r) / u0r;
      Real b3r = (bb3r + b0r * u3r) / u0r;

      // Transform contravariant magnetic fields
      Real btl = mt_0*b0l;
      Real bxl = mx_0*b0l + mx_1*b1l;
      Real byl = my_2*b2l;
      Real bzl = mz_3*b3l;
      Real btr = mt_0*b0r;
      Real bxr = mx_0*b0r + mx_1*b1r;
      Real byr = my_2*b2r;
      Real bzr = mz_3*b3r;

      // Set local magnetic fields
      Real bbxl = utl * bxl - uxl * btl;
      Real bbxr = utr * bxr - uxr * btr;
      bx(i) = 0.5 * (bbxl + bbxr);
      bb2l = utl * byl - uyl * btl;
      bb3l = utl * bzl - uzl * btl;
      bb2r = utr * byr - uyr * btr;
      bb3r = utr * bzr - uzr * btr;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: y-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x'-index bounds
//   b2_vals: 3D array of normal components B^2 of magnetic field, in global coordinates
//   prim_left: 1D array of left primitives, using global coordinates
//   prim_right: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_left: values overwritten in local coordinates
//   prim_right: values overwritten in local coordinates
//   by: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   expects v1/v2/v3 in IVX/IVY/IVZ slots
//   expects B2 in b2_vals
//   expects B3/B1 in IBY/IBZ slots
//   puts vx/vy/vz in IVY/IVZ/IVX slots
//   puts Bx in bx
//   puts By/Bz in IBY/IBZ slots
void Coordinates::PrimToLocal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b2_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx)
{
  // Prepare geometric quantities
  const Real g_00 = -SQR(beta) / SQR(alpha);
  const Real g_01 = 2.0*tilted_a_ / SQR(alpha);
  const Real g_10 = g_01;
  const Real g_11 = SQR(beta) / SQR(alpha);
  const Real g_22 = 1.0;
  const Real g_33 = 1.0;
  const Real mt_0 = alpha/beta;
  const Real mx_2 = 1.0;
  const Real my_3 = 1.0;
  const Real mz_0 = 2.0*tilted_a_ / (alpha*beta);
  const Real mz_1 = beta/alpha;

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract global 3-velocities
    Real &v1l = prim_left(IVX,i);
    Real &v2l = prim_left(IVY,i);
    Real &v3l = prim_left(IVZ,i);
    Real &v1r = prim_right(IVX,i);
    Real &v2r = prim_right(IVY,i);
    Real &v3r = prim_right(IVZ,i);

    // Construct global 4-velocities
    Real u0l = std::sqrt(-1.0 /
        (g_00 + 2.0*g_01*v1l + g_11*v1l*v1l + g_22*v2l*v2l + g_33*v3l*v3l));
    Real u1l = u0l * v1l;
    Real u2l = u0l * v2l;
    Real u3l = u0l * v3l;
    Real u0r = std::sqrt(-1.0 /
        (g_00 + 2.0*g_01*v1r + g_11*v1r*v1r + g_22*v2r*v2r + g_33*v3r*v3r));
    Real u1r = u0r * v1r;
    Real u2r = u0r * v2r;
    Real u3r = u0r * v3r;

    // Transform 4-velocities
    Real utl = mt_0*u0l;
    Real uxl = mx_2*u2l;
    Real uyl = my_3*u3l;
    Real uzl = mz_0*u0l + mz_1*u1l;
    Real utr = mt_0*u0r;
    Real uxr = mx_2*u2r;
    Real uyr = my_3*u3r;
    Real uzr = mz_0*u0r + mz_1*u1r;

    // Set local 3-velocities
    v2l = uxl / utl;
    v3l = uyl / utl;
    v1l = uzl / utl;
    v2r = uxr / utr;
    v3r = uyr / utr;
    v1r = uzr / utr;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      // Extract global magnetic fields
      const Real &bb2 = b2_vals(k,j,i);
      Real &bb3l = prim_left(IBY,i);
      Real &bb1l = prim_left(IBZ,i);
      Real &bb3r = prim_right(IBY,i);
      Real &bb1r = prim_right(IBZ,i);

      // Construct global contravariant magnetic fields
      Real b0l = g_10*bb1l*u0l + g_11*bb1l*u1l + g_22*bb2*u2l + g_33*bb3l*u3l;
      Real b1l = (bb1l + b0l * u1l) / u0l;
      Real b2l = (bb2 + b0l * u2l) / u0l;
      Real b3l = (bb3l + b0l * u3l) / u0l;
      Real b0r = g_10*bb1r*u0r + g_11*bb1r*u1r + g_22*bb2*u2r + g_33*bb3r*u3r;
      Real b1r = (bb1r + b0r * u1r) / u0r;
      Real b2r = (bb2 + b0r * u2r) / u0r;
      Real b3r = (bb3r + b0r * u3r) / u0r;

      // Transform contravariant magnetic fields
      Real btl = mt_0*b0l;
      Real bxl = mx_2*b2l;
      Real byl = my_3*b3l;
      Real bzl = mz_0*b0l + mz_1*b1l;
      Real btr = mt_0*b0r;
      Real bxr = mx_2*b2r;
      Real byr = my_3*b3r;
      Real bzr = mz_0*b0r + mz_1*b1r;

      // Set local magnetic fields
      Real bbxl = utl * bxl - uxl * btl;
      Real bbxr = utr * bxr - uxr * btr;
      bx(i) = 0.5 * (bbxl + bbxr);
      bb3l = utl * byl - uyl * btl;
      bb1l = utl * bzl - uzl * btl;
      bb3r = utr * byr - uyr * btr;
      bb1r = utr * bzr - uzr * btr;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: z-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x'-index bounds
//   b3_vals: 3D array of normal components B^3 of magnetic field, in global coordinates
//   prim_left: 1D array of left primitives, using global coordinates
//   prim_right: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_left: values overwritten in local coordinates
//   prim_right: values overwritten in local coordinates
//   bz: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   expects v1/v2/v3 in IVX/IVY/IVZ slots
//   expects B3 in b3_vals
//   expects B1/B2 in IBY/IBZ slots
//   puts vx/vy/vz in IVZ/IVX/IVY slots
//   puts Bx in bx
//   puts By/Bz in IBY/IBZ slots
void Coordinates::PrimToLocal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &b3_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx)
{
  // Prepare geometric quantities
  const Real g_00 = -SQR(beta) / SQR(alpha);
  const Real g_01 = 2.0*tilted_a_ / SQR(alpha);
  const Real g_10 = g_01;
  const Real g_11 = SQR(beta) / SQR(alpha);
  const Real g_22 = 1.0;
  const Real g_33 = 1.0;
  const Real mt_0 = alpha/beta;
  const Real mx_3 = 1.0;
  const Real my_0 = 2.0*tilted_a_ / (alpha*beta);
  const Real my_1 = beta/alpha;
  const Real mz_2 = 1.0;

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract global 3-velocities
    Real &v1l = prim_left(IVX,i);
    Real &v2l = prim_left(IVY,i);
    Real &v3l = prim_left(IVZ,i);
    Real &v1r = prim_right(IVX,i);
    Real &v2r = prim_right(IVY,i);
    Real &v3r = prim_right(IVZ,i);

    // Construct global 4-velocities
    Real u0l = std::sqrt(-1.0 /
        (g_00 + 2.0*g_01*v1l + g_11*v1l*v1l + g_22*v2l*v2l + g_33*v3l*v3l));
    Real u1l = u0l * v1l;
    Real u2l = u0l * v2l;
    Real u3l = u0l * v3l;
    Real u0r = std::sqrt(-1.0 /
        (g_00 + 2.0*g_01*v1r + g_11*v1r*v1r + g_22*v2r*v2r + g_33*v3r*v3r));
    Real u1r = u0r * v1r;
    Real u2r = u0r * v2r;
    Real u3r = u0r * v3r;

    // Transform 4-velocities
    Real utl = mt_0*u0l;
    Real uxl = mx_3*u3l;
    Real uyl = my_0*u0l + my_1*u1l;
    Real uzl = mz_2*u2l;
    Real utr = mt_0*u0r;
    Real uxr = mx_3*u3r;
    Real uyr = my_0*u0r + my_1*u1r;
    Real uzr = mz_2*u2r;

    // Set local 3-velocities
    v3l = uxl / utl;
    v1l = uyl / utl;
    v2l = uzl / utl;
    v3r = uxr / utr;
    v1r = uyr / utr;
    v2r = uzr / utr;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      // Extract global magnetic fields
      const Real &bb3 = b3_vals(k,j,i);
      Real &bb1l = prim_left(IBY,i);
      Real &bb2l = prim_left(IBZ,i);
      Real &bb1r = prim_right(IBY,i);
      Real &bb2r = prim_right(IBZ,i);

      // Construct global contravariant magnetic fields
      Real b0l = g_10*bb1l*u0l + g_11*bb1l*u1l + g_22*bb2l*u2l + g_33*bb3*u3l;
      Real b1l = (bb1l + b0l * u1l) / u0l;
      Real b2l = (bb2l + b0l * u2l) / u0l;
      Real b3l = (bb3 + b0l * u3l) / u0l;
      Real b0r = g_10*bb1r*u0r + g_11*bb1r*u1r + g_22*bb2r*u2r + g_33*bb3*u3r;
      Real b1r = (bb1r + b0r * u1r) / u0r;
      Real b2r = (bb2r + b0r * u2r) / u0r;
      Real b3r = (bb3 + b0r * u3r) / u0r;

      // Transform contravariant magnetic fields
      Real btl = mt_0*b0l;
      Real bxl = mx_3*b3l;
      Real byl = my_0*b0l + my_1*b1l;
      Real bzl = mz_2*b2l;
      Real btr = mt_0*b0r;
      Real bxr = mx_3*b3r;
      Real byr = my_0*b0r + my_1*b1r;
      Real bzr = mz_2*b2r;

      // Set local magnetic fields
      Real bbxl = utl * bxl - uxl * btl;
      Real bbxr = utr * bxr - uxr * btr;
      bx(i) = 0.5 * (bbxl + bbxr);
      bb1l = utl * byl - uyl * btl;
      bb2l = utl * bzl - uzl * btl;
      bb1r = utr * byr - uyr * btr;
      bb2r = utr * bzr - uzr * btr;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: x'-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x'-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates
//   bx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM1/IM2/IM3 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts x1-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts x1-fluxes of B2/B3 in IBY/IBZ slots
void Coordinates::FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx, AthenaArray<Real> &flux)
{
  // Prepare geometric quantities
  const Real g_00 = -SQR(beta) / SQR(alpha);
  const Real g_01 = 2.0*tilted_a_ / SQR(alpha);
  const Real g_10 = g_01;
  const Real g_11 = SQR(beta) / SQR(alpha);
  const Real g_22 = 1.0;
  const Real g_33 = 1.0;
  const Real m0_t = beta/alpha;
  const Real m1_t = -2.0*tilted_a_ / (alpha*beta);
  const Real m1_x = alpha/beta;
  const Real m2_y = 1.0;
  const Real m3_z = 1.0;

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract local conserved quantities and fluxes
    const Real &dt = cons(IDN,i);
    const Real &ttt = cons(IEN,i);
    const Real &ttx = cons(IM1,i);
    const Real &tty = cons(IM2,i);
    const Real &ttz = cons(IM3,i);
    Real dx = flux(IDN,i);
    Real txt = flux(IEN,i);
    Real txx = flux(IM1,i);
    Real txy = flux(IM2,i);
    Real txz = flux(IM3,i);

    // Transform stress-energy tensor
    Real t10 = m1_t*m0_t*ttt + m1_x*m0_t*txt;
    Real t11 = m1_t*m1_t*ttt + m1_t*m1_x*ttx + m1_x*m1_t*txt + m1_x*m1_x*txx;
    Real t12 = m1_t*m2_y*tty + m1_x*m2_y*txy;
    Real t13 = m1_t*m3_z*ttz + m1_x*m3_z*txz;

    // Extract global fluxes
    Real &d1 = flux(IDN,i);
    Real &t1_0 = flux(IEN,i);
    Real &t1_1 = flux(IM1,i);
    Real &t1_2 = flux(IM2,i);
    Real &t1_3 = flux(IM3,i);

    // Set fluxes
    d1 = m1_t*dt + m1_x*dx;
    t1_0 = g_00*t10 + g_01*t11;
    t1_1 = g_10*t10 + g_11*t11;
    t1_2 = g_22*t12;
    t1_3 = g_33*t13;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      const Real &fyt = cons(IBY,i);
      const Real &fzt = cons(IBZ,i);
      Real fyx = flux(IBY,i);
      Real fzx = flux(IBZ,i);
      Real &f21 = flux(IBY,i);
      Real &f31 = flux(IBZ,i);
      f21 = m2_y*m1_t*fyt + m2_y*m1_x*fyx;
      f31 = m3_z*m1_t*fzt + m3_z*m1_x*fzx;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: y-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x'-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates
//   bx: 1D array of longitudinal magnetic fields, in local coordinates
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM2/IM3/IM1 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts x2-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts x2-fluxes of B3/B1 in IBY/IBZ slots
void Coordinates::FluxToGlobal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx, AthenaArray<Real> &flux)
{
  // Prepare geometric quantities
  const Real g_00 = -SQR(beta) / SQR(alpha);
  const Real g_01 = 2.0*tilted_a_ / SQR(alpha);
  const Real g_10 = g_01;
  const Real g_11 = SQR(beta) / SQR(alpha);
  const Real g_22 = 1.0;
  const Real g_33 = 1.0;
  const Real m0_t = beta/alpha;
  const Real m1_t = -2.0*tilted_a_ / (alpha*beta);
  const Real m1_z = alpha/beta;
  const Real m2_x = 1.0;
  const Real m3_y = 1.0;

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract local conserved quantities and fluxes
    const Real &dt = cons(IDN,i);
    const Real &ttt = cons(IEN,i);
    const Real &ttx = cons(IM2,i);
    const Real &tty = cons(IM3,i);
    const Real &ttz = cons(IM1,i);
    Real dx = flux(IDN,i);
    Real txt = flux(IEN,i);
    Real txx = flux(IM2,i);
    Real txy = flux(IM3,i);
    Real txz = flux(IM1,i);

    // Transform stress-energy tensor
    Real t20 = m2_x*m0_t*txt;
    Real t21 = m2_x*m1_t*txt + m2_x*m1_z*txz;
    Real t22 = m2_x*m2_x*txx;
    Real t23 = m2_x*m3_y*txy;

    // Extract global fluxes
    Real &d2 = flux(IDN,i);
    Real &t2_0 = flux(IEN,i);
    Real &t2_1 = flux(IM1,i);
    Real &t2_2 = flux(IM2,i);
    Real &t2_3 = flux(IM3,i);

    // Set fluxes
    d2 = m2_x*dx;
    t2_0 = g_00*t20 + g_01*t21;
    t2_1 = g_10*t20 + g_11*t21;
    t2_2 = g_22*t22;
    t2_3 = g_33*t23;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      Real ftx = -bx(i);
      Real fyx = flux(IBY,i);
      Real fzx = flux(IBZ,i);
      Real &f32 = flux(IBY,i);
      Real &f12 = flux(IBZ,i);
      f32 = m3_y*m2_x*fyx;
      f12 = m1_t*m2_x*ftx + m1_z*m2_x*fzx;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: z-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x'-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates
//   bx: 1D array of longitudinal magnetic fields, in local coordinates
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM3/IM1/IM2 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts x3-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts x3-fluxes of B1/B2 in IBY/IBZ slots
void Coordinates::FluxToGlobal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx, AthenaArray<Real> &flux)
{
  // Prepare geometric quantities
  const Real g_00 = -SQR(beta) / SQR(alpha);
  const Real g_01 = 2.0*tilted_a_ / SQR(alpha);
  const Real g_10 = g_01;
  const Real g_11 = SQR(beta) / SQR(alpha);
  const Real g_22 = 1.0;
  const Real g_33 = 1.0;
  const Real m0_t = beta/alpha;
  const Real m1_t = -2.0*tilted_a_ / (alpha*beta);
  const Real m1_y = alpha/beta;
  const Real m2_z = 1.0;
  const Real m3_x = 1.0;

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract local conserved quantities and fluxes
    const Real &dt = cons(IDN,i);
    const Real &ttt = cons(IEN,i);
    const Real &ttx = cons(IM3,i);
    const Real &tty = cons(IM1,i);
    const Real &ttz = cons(IM2,i);
    Real dx = flux(IDN,i);
    Real txt = flux(IEN,i);
    Real txx = flux(IM3,i);
    Real txy = flux(IM1,i);
    Real txz = flux(IM2,i);

    // Transform stress-energy tensor
    Real t30 = m3_x*m0_t*txt;
    Real t31 = m3_x*m1_t*txt + m3_x*m1_y*txy;
    Real t32 = m3_x*m2_z*txz;
    Real t33 = m3_x*m3_x*txx;

    // Extract global fluxes
    Real &d3 = flux(IDN,i);
    Real &t3_0 = flux(IEN,i);
    Real &t3_1 = flux(IM1,i);
    Real &t3_2 = flux(IM2,i);
    Real &t3_3 = flux(IM3,i);

    // Set fluxes
    d3 = m3_x*dx;
    t3_0 = g_00*t30 + g_01*t31;
    t3_1 = g_10*t30 + g_11*t31;
    t3_2 = g_22*t32;
    t3_3 = g_33*t33;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      Real ftx = -bx(i);
      Real fyx = flux(IBY,i);
      Real fzx = flux(IBZ,i);
      Real &f13 = flux(IBY,i);
      Real &f23 = flux(IBZ,i);
      f13 = m1_t*m3_x*ftx + m1_y*m3_x*fyx;
      f23 = m2_z*m3_x*fzx;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for calculating distance between two points
// Inputs:
//   a1,a2,a3: global coordinates of first point (unused)
//   bx,by,bz: Minkowski coordinates of second point (unused)
// Outputs:
//   returned value: NAN
// Notes:
//   should not be called, as constant-t' surfaces are not constant-t surfaces
Real Coordinates::DistanceBetweenPoints(Real a1, Real a2, Real a3, Real bx, Real by,
    Real bz)
{
  return NAN;
}

//--------------------------------------------------------------------------------------

// Function for calculating Minkowski coordinates of cell
// Inputs:
//   x0,x1,x2,x3: tilted coordinates
// Outputs:
//   pt,px,py,pz: Minkowski coordinate values set
// Notes:
//   transformation given by:
//     t = (t' - a x') / \alpha
//     x = (x' + a t') / \alpha
void Coordinates::MinkowskiCoordinates(Real x0, Real x1, Real x2, Real x3,
    Real *pt, Real *px, Real *py, Real *pz)
{
  *pt = (x0 - tilted_a_ * x1) / alpha;
  *px = (x1 + tilted_a_ * x0) / alpha;
  *py = x2;
  *pz = x3;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: cell-centered
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j,i: z-, y-, and x'-indices (unused)
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial
void Coordinates::TransformVectorCell(
    Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = 1.0/alpha * at + tilted_a_/alpha * ax;
  *pa1 = -tilted_a_/alpha * at + 1.0/alpha * ax;
  *pa2 = ay;
  *pa3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: x-interface
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j,i: z-, y-, and x'-indices (unused)
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial
void Coordinates::TransformVectorFace1(
    Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = 1.0/alpha * at + tilted_a_/alpha * ax;
  *pa1 = -tilted_a_/alpha * at + 1.0/alpha * ax;
  *pa2 = ay;
  *pa3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: y-interface
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j,i: z-, y-, and x'-indices (unused)
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial
void Coordinates::TransformVectorFace2(
    Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = 1.0/alpha * at + tilted_a_/alpha * ax;
  *pa1 = -tilted_a_/alpha * at + 1.0/alpha * ax;
  *pa2 = ay;
  *pa3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: z-interface
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j,i: z-, y-, and x'-indices (unused)
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
// Notes:
//   transformation is trivial
void Coordinates::TransformVectorFace3(
    Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = 1.0/alpha * at + tilted_a_/alpha * ax;
  *pa1 = -tilted_a_/alpha * at + 1.0/alpha * ax;
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
void Coordinates::LowerVectorCell(
    Real a0, Real a1, Real a2, Real a3, int k, int j, int i,
    Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3)
{
  // Extract geometric quantities
  Real g_00 = -SQR(beta) / SQR(alpha);
  Real g_01 = 2.0 * tilted_a_ / SQR(alpha);
  Real g_10 = g_01;
  Real g_11 = SQR(beta) / SQR(alpha);
  Real g_22 = 1.0;
  Real g_33 = 1.0;

  // Set lowered components
  *pa_0 = g_00*a0 + g_01*a1;
  *pa_1 = g_10*a0 + g_11*a1;
  *pa_2 = g_22*a2;
  *pa_3 = g_33*a3;
  return;
}
