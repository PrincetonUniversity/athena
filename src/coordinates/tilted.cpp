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
#include "../mesh.hpp"             // MeshBlock
#include "../parameter_input.hpp"  // ParameterInput
#include "../fluid/fluid.hpp"      // Fluid
#include "../fluid/eos/eos.hpp"    // FluidEqnOfState

// Global variables
static Real alpha;  // \sqrt{1+a^2}
static Real beta;   // \sqrt{1-a^2}

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

  // Set parameters
  tilted_a_ = pin->GetReal("coord", "a");
  const Real &a = tilted_a_;
  alpha = std::sqrt(1.0 + SQR(a));
  beta = std::sqrt(1.0 - SQR(a));

  // Initialize volume-averaged positions and spacings: x'-direction
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

  // Allocate arrays for intermediate geometric quantities: x'-direction
  int n_cells_1 = pmb->block_size.nx1 + 2*NGHOST;
  g_.NewAthenaArray(NMETRIC, n_cells_1);
  gi_.NewAthenaArray(NMETRIC, n_cells_1);
}

//--------------------------------------------------------------------------------------

// Destructor
Coordinates::~Coordinates()
{
  DeleteBasicCoordinates();
  g_.DeleteAthenaArray();
  gi_.DeleteAthenaArray();
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
//   cf. GetCellVolume()
void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &volumes)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    volumes(i) = dx1f(i) * dx2f(j) * dx3f(k);
  return;
}

// GetCellVolume returns only one CellVolume at i
Real Coordinates::GetCellVolume(const int k, const int j, const int i)
{
  return dx1f(i)*dx2f(j)*dx3f(k);
}

//--------------------------------------------------------------------------------------

// Function for computing single cell volume
// Inputs:
//   k,j,i: z-, y-, and x'-indices
// Outputs:
//   returned value: cell volume
// Notes:
//   \Delta V = \Delta x' * \Delta y * \Delta z
//   cf. CellVolume()
Real Coordinates::GetCellVolume(const int k, const int j, const int i)
{
    return dx1f(i) * dx2f(j) * dx3f(k);
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

// Function for computing single area orthogonal to x'
// Inputs:
//   k,j,i: z-, y-, and x'-indices
// Outputs:
//   returned value: interface area orthogonal to x'
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
//   il,iu: x'-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to y
// Notes:
//   \Delta A = \Delta x' * \Delta z
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
//   il,iu: x'-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to z
// Notes:
//   \Delta A = \Delta x' * \Delta y
void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    areas(i) = dx1f(i) * dx2f(j);
  return;
}


// GetFace1Area returns only one Face1Area at i
Real Coordinates::GetFace1Area(const int k, const int j, const int i)
{
  return dx2f(j)*dx3f(k);
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
    lengths(i) = dx1f(i);
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
//   il,iu: x'-index bounds
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

// GetEdge?Length functions: return one edge length at i
Real Coordinates::GetEdge2Length(const int k, const int j, const int i)
{
  return dx2f(j);
}

Real Coordinates::GetEdge3Length(const int k, const int j, const int i)
{
  return dx3f(k);
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
//   bb1: 3D array of normal components B^1 of magnetic field, in global coordinates
//   prim_l: 1D array of left primitives, using global coordinates
//   prim_r: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
//   bbx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   expects \tilde{u}^1/\tilde{u}^2/\tilde{u}^3 in IVX/IVY/IVZ slots
//   expects B^1 in bb1
//   expects B^2/B^3 in IBY/IBZ slots
//   puts \tilde{u}^x/\tilde{u}^y/\tilde{u}^z in IVX/IVY/IVZ
//   puts B^x in bbx
//   puts B^y/B^z in IBY/IBZ slots
//   \tilde{u}^\hat{i} = u^\hat{i}
void Coordinates::PrimToLocal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb1, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &bbx)
{
  // Calculate metric coefficients
  if (MAGNETIC_FIELDS_ENABLED)
    Face1Metric(k, j, il, iu, g_, gi_);

  // Extract transformation coefficients
  const Real mt_0 = alpha/beta;
  const Real mx_0 = 2.0*tilted_a_ / (alpha*beta);
  const Real mx_1 = beta/alpha;
  const Real my_2 = 1.0;
  const Real mz_3 = 1.0;

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract global projected 4-velocities
    Real uu0_l = 0.0;
    Real uu1_l = prim_l(IVX,i);
    Real uu2_l = prim_l(IVY,i);
    Real uu3_l = prim_l(IVZ,i);
    Real uu0_r = 0.0;
    Real uu1_r = prim_r(IVX,i);
    Real uu2_r = prim_r(IVY,i);
    Real uu3_r = prim_r(IVZ,i);

    // Transform projected 4-velocities
    Real uux_l = mx_0*uu0_l + mx_1*uu1_l;
    Real uuy_l = my_2*uu2_l;
    Real uuz_l = mz_3*uu3_l;
    Real uux_r = mx_0*uu0_r + mx_1*uu1_r;
    Real uuy_r = my_2*uu2_r;
    Real uuz_r = mz_3*uu3_r;

    // Set local projected 4-velocities
    prim_l(IVX,i) = uux_l;
    prim_l(IVY,i) = uuy_l;
    prim_l(IVZ,i) = uuz_l;
    prim_r(IVX,i) = uux_r;
    prim_r(IVY,i) = uuy_r;
    prim_r(IVZ,i) = uuz_r;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      // Extract metric coefficients
      const Real &g_00 = g_(I00,i);
      const Real &g_01 = g_(I01,i);
      const Real &g_02 = 0.0;
      const Real &g_03 = 0.0;
      const Real &g_10 = g_(I01,i);
      const Real &g_11 = g_(I11,i);
      const Real &g_12 = 0.0;
      const Real &g_13 = 0.0;
      const Real &g_20 = 0.0;
      const Real &g_21 = 0.0;
      const Real &g_22 = g_(I22,i);
      const Real &g_23 = 0.0;
      const Real &g_30 = 0.0;
      const Real &g_31 = 0.0;
      const Real &g_32 = 0.0;
      const Real &g_33 = g_(I33,i);
      const Real &g01 = gi_(I01,i);
      const Real &g02 = 0.0;
      const Real &g03 = 0.0;
      Real alpha = std::sqrt(-1.0/gi_(I00,i));

      // Calculate global 4-velocities
      Real tmp = g_11*uu1_l*uu1_l + 2.0*g_12*uu1_l*uu2_l + 2.0*g_13*uu1_l*uu3_l
               + g_22*uu2_l*uu2_l + 2.0*g_23*uu2_l*uu3_l
               + g_33*uu3_l*uu3_l;
      Real gamma_l = std::sqrt(1.0 + tmp);
      Real u0_l = gamma_l / alpha;
      Real u1_l = uu1_l - alpha * gamma_l * g01;
      Real u2_l = uu2_l - alpha * gamma_l * g02;
      Real u3_l = uu3_l - alpha * gamma_l * g03;
      tmp = g_11*uu1_r*uu1_r + 2.0*g_12*uu1_r*uu2_r + 2.0*g_13*uu1_r*uu3_r
          + g_22*uu2_r*uu2_r + 2.0*g_23*uu2_r*uu3_r
          + g_33*uu3_r*uu3_r;
      Real gamma_r = std::sqrt(1.0 + tmp);
      Real u0_r = gamma_r / alpha;
      Real u1_r = uu1_r - alpha * gamma_r * g01;
      Real u2_r = uu2_r - alpha * gamma_r * g02;
      Real u3_r = uu3_r - alpha * gamma_r * g03;

      // Extract global magnetic fields
      const Real &bb1_l = bb1(k,j,i);
      const Real &bb1_r = bb1(k,j,i);
      Real &bb2_l = prim_l(IBY,i);
      Real &bb3_l = prim_l(IBZ,i);
      Real &bb2_r = prim_r(IBY,i);
      Real &bb3_r = prim_r(IBZ,i);

      // Calculate global 4-magnetic fields
      Real b0_l = g_10*bb1_l*u0_l + g_11*bb1_l*u1_l + g_12*bb1_l*u2_l + g_13*bb1_l*u3_l
                + g_20*bb2_l*u0_l + g_21*bb2_l*u1_l + g_22*bb2_l*u2_l + g_23*bb2_l*u3_l
                + g_30*bb3_l*u0_l + g_31*bb3_l*u1_l + g_32*bb3_l*u2_l + g_33*bb3_l*u3_l;
      Real b1_l = (bb1_l + b0_l * u1_l) / u0_l;
      Real b2_l = (bb2_l + b0_l * u2_l) / u0_l;
      Real b3_l = (bb3_l + b0_l * u3_l) / u0_l;
      Real b0_r = g_10*bb1_r*u0_r + g_11*bb1_r*u1_r + g_12*bb1_r*u2_r + g_13*bb1_r*u3_r
                + g_20*bb2_r*u0_r + g_21*bb2_r*u1_r + g_22*bb2_r*u2_r + g_23*bb2_r*u3_r
                + g_30*bb3_r*u0_r + g_31*bb3_r*u1_r + g_32*bb3_r*u2_r + g_33*bb3_r*u3_r;
      Real b1_r = (bb1_r + b0_r * u1_r) / u0_r;
      Real b2_r = (bb2_r + b0_r * u2_r) / u0_r;
      Real b3_r = (bb3_r + b0_r * u3_r) / u0_r;

      // Transform 4-velocities
      Real ut_l = mt_0*u0_l;
      Real ux_l = mx_0*u0_l + mx_1*u1_l;
      Real uy_l = my_2*u2_l;
      Real uz_l = mz_3*u3_l;
      Real ut_r = mt_0*u0_r;
      Real ux_r = mx_0*u0_r + mx_1*u1_r;
      Real uy_r = my_2*u2_r;
      Real uz_r = mz_3*u3_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_0*b0_l + mx_1*b1_l;
      Real by_l = my_2*b2_l;
      Real bz_l = mz_3*b3_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_0*b0_r + mx_1*b1_r;
      Real by_r = my_2*b2_r;
      Real bz_r = mz_3*b3_r;

      // Set local magnetic fields
      Real bbx_l = ut_l * bx_l - ux_l * bt_l;
      Real bbx_r = ut_r * bx_r - ux_r * bt_r;
      bbx(i) = 0.5 * (bbx_l + bbx_r);
      bb2_l = ut_l * by_l - uy_l * bt_l;
      bb3_l = ut_l * bz_l - uz_l * bt_l;
      bb2_r = ut_r * by_r - uy_r * bt_r;
      bb3_r = ut_r * bz_r - uz_r * bt_r;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: y-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x'-index bounds
//   bb2: 3D array of normal components B^2 of magnetic field, in global coordinates
//   prim_l: 1D array of left primitives, using global coordinates
//   prim_r: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
//   bbx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   expects \tilde{u}^1/\tilde{u}^2/\tilde{u}^3 in IVX/IVY/IVZ slots
//   expects B^2 in bb2
//   expects B^3/B^1 in IBY/IBZ slots
//   puts \tilde{u}^x/\tilde{u}^y/\tilde{u}^z in IVY/IVZ/IVX slots
//   puts B^x in bbx
//   puts B^y/B^z in IBY/IBZ slots
//   \tilde{u}^\hat{i} = u^\hat{i}
void Coordinates::PrimToLocal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb2, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &bbx)
{
  // Calculate metric coefficients
  if (MAGNETIC_FIELDS_ENABLED)
    Face2Metric(k, j, il, iu, g_, gi_);

  // Extract transformation coefficients
  const Real mt_0 = alpha/beta;
  const Real mx_2 = 1.0;
  const Real my_3 = 1.0;
  const Real mz_0 = 2.0*tilted_a_ / (alpha*beta);
  const Real mz_1 = beta/alpha;

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract global projected 4-velocities
    Real uu0_l = 0.0;
    Real uu1_l = prim_l(IVX,i);
    Real uu2_l = prim_l(IVY,i);
    Real uu3_l = prim_l(IVZ,i);
    Real uu0_r = 0.0;
    Real uu1_r = prim_r(IVX,i);
    Real uu2_r = prim_r(IVY,i);
    Real uu3_r = prim_r(IVZ,i);

    // Transform projected 4-velocities
    Real uux_l = mx_2*uu2_l;
    Real uuy_l = my_3*uu3_l;
    Real uuz_l = mz_0*uu0_l + mz_1*uu1_l;
    Real uux_r = mx_2*uu2_r;
    Real uuy_r = my_3*uu3_r;
    Real uuz_r = mz_0*uu0_r + mz_1*uu1_r;

    // Set local projected 4-velocities
    prim_l(IVY,i) = uux_l;
    prim_l(IVZ,i) = uuy_l;
    prim_l(IVX,i) = uuz_l;
    prim_r(IVY,i) = uux_r;
    prim_r(IVZ,i) = uuy_r;
    prim_r(IVX,i) = uuz_r;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      // Extract metric coefficients
      const Real &g_00 = g_(I00,i);
      const Real &g_01 = g_(I01,i);
      const Real &g_02 = 0.0;
      const Real &g_03 = 0.0;
      const Real &g_10 = g_(I01,i);
      const Real &g_11 = g_(I11,i);
      const Real &g_12 = 0.0;
      const Real &g_13 = 0.0;
      const Real &g_20 = 0.0;
      const Real &g_21 = 0.0;
      const Real &g_22 = g_(I22,i);
      const Real &g_23 = 0.0;
      const Real &g_30 = 0.0;
      const Real &g_31 = 0.0;
      const Real &g_32 = 0.0;
      const Real &g_33 = g_(I33,i);
      const Real &g01 = gi_(I01,i);
      const Real &g02 = 0.0;
      const Real &g03 = 0.0;
      Real alpha = std::sqrt(-1.0/gi_(I00,i));

      // Calculate global 4-velocities
      Real tmp = g_11*uu1_l*uu1_l + 2.0*g_12*uu1_l*uu2_l + 2.0*g_13*uu1_l*uu3_l
               + g_22*uu2_l*uu2_l + 2.0*g_23*uu2_l*uu3_l
               + g_33*uu3_l*uu3_l;
      Real gamma_l = std::sqrt(1.0 + tmp);
      Real u0_l = gamma_l / alpha;
      Real u1_l = uu1_l - alpha * gamma_l * g01;
      Real u2_l = uu2_l - alpha * gamma_l * g02;
      Real u3_l = uu3_l - alpha * gamma_l * g03;
      tmp = g_11*uu1_r*uu1_r + 2.0*g_12*uu1_r*uu2_r + 2.0*g_13*uu1_r*uu3_r
          + g_22*uu2_r*uu2_r + 2.0*g_23*uu2_r*uu3_r
          + g_33*uu3_r*uu3_r;
      Real gamma_r = std::sqrt(1.0 + tmp);
      Real u0_r = gamma_r / alpha;
      Real u1_r = uu1_r - alpha * gamma_r * g01;
      Real u2_r = uu2_r - alpha * gamma_r * g02;
      Real u3_r = uu3_r - alpha * gamma_r * g03;

      // Extract global magnetic fields
      const Real &bb2_l = bb2(k,j,i);
      const Real &bb2_r = bb2(k,j,i);
      Real &bb3_l = prim_l(IBY,i);
      Real &bb1_l = prim_l(IBZ,i);
      Real &bb3_r = prim_r(IBY,i);
      Real &bb1_r = prim_r(IBZ,i);

      // Calculate global 4-magnetic fields
      Real b0_l = g_10*bb1_l*u0_l + g_11*bb1_l*u1_l + g_12*bb1_l*u2_l + g_13*bb1_l*u3_l
                + g_20*bb2_l*u0_l + g_21*bb2_l*u1_l + g_22*bb2_l*u2_l + g_23*bb2_l*u3_l
                + g_30*bb3_l*u0_l + g_31*bb3_l*u1_l + g_32*bb3_l*u2_l + g_33*bb3_l*u3_l;
      Real b1_l = (bb1_l + b0_l * u1_l) / u0_l;
      Real b2_l = (bb2_l + b0_l * u2_l) / u0_l;
      Real b3_l = (bb3_l + b0_l * u3_l) / u0_l;
      Real b0_r = g_10*bb1_r*u0_r + g_11*bb1_r*u1_r + g_12*bb1_r*u2_r + g_13*bb1_r*u3_r
                + g_20*bb2_r*u0_r + g_21*bb2_r*u1_r + g_22*bb2_r*u2_r + g_23*bb2_r*u3_r
                + g_30*bb3_r*u0_r + g_31*bb3_r*u1_r + g_32*bb3_r*u2_r + g_33*bb3_r*u3_r;
      Real b1_r = (bb1_r + b0_r * u1_r) / u0_r;
      Real b2_r = (bb2_r + b0_r * u2_r) / u0_r;
      Real b3_r = (bb3_r + b0_r * u3_r) / u0_r;

      // Transform 4-velocities
      Real ut_l = mt_0*u0_l;
      Real ux_l = mx_2*u2_l;
      Real uy_l = my_3*u3_l;
      Real uz_l = mz_0*u0_l + mz_1*u1_l;
      Real ut_r = mt_0*u0_r;
      Real ux_r = mx_2*u2_r;
      Real uy_r = my_3*u3_r;
      Real uz_r = mz_0*u0_r + mz_1*u1_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_2*b2_l;
      Real by_l = my_3*b3_l;
      Real bz_l = mz_0*b0_l + mz_1*b1_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_2*b2_r;
      Real by_r = my_3*b3_r;
      Real bz_r = mz_0*b0_r + mz_1*b1_r;

      // Set local magnetic fields
      Real bbx_l = ut_l * bx_l - ux_l * bt_l;
      Real bbx_r = ut_r * bx_r - ux_r * bt_r;
      bbx(i) = 0.5 * (bbx_l + bbx_r);
      bb3_l = ut_l * by_l - uy_l * bt_l;
      bb1_l = ut_l * bz_l - uz_l * bt_l;
      bb3_r = ut_r * by_r - uy_r * bt_r;
      bb1_r = ut_r * bz_r - uz_r * bt_r;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: z-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x'-index bounds
//   bb3: 3D array of normal components B^3 of magnetic field, in global coordinates
//   prim_l: 1D array of left primitives, using global coordinates
//   prim_r: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
//   bbx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   expects \tilde{u}^1/\tilde{u}^2/\tilde{u}^3 in IVX/IVY/IVZ slots
//   expects B^3 in bb3
//   expects B^1/B^2 in IBY/IBZ slots
//   puts \tilde{u}^x/\tilde{u}^y/\tilde{u}^z in IVZ/IVX/IVY slots
//   puts B^x in bbx
//   puts B^y/B^z in IBY/IBZ slots
//   \tilde{u}^\hat{i} = u^\hat{i}
void Coordinates::PrimToLocal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb3, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &bbx)
{
  // Calculate metric coefficients
  if (MAGNETIC_FIELDS_ENABLED)
    Face3Metric(k, j, il, iu, g_, gi_);

  // Extract transformation coefficients
  const Real mt_0 = alpha/beta;
  const Real mx_3 = 1.0;
  const Real my_0 = 2.0*tilted_a_ / (alpha*beta);
  const Real my_1 = beta/alpha;
  const Real mz_2 = 1.0;

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract global projected 4-velocities
    Real uu0_l = 0.0;
    Real uu1_l = prim_l(IVX,i);
    Real uu2_l = prim_l(IVY,i);
    Real uu3_l = prim_l(IVZ,i);
    Real uu0_r = 0.0;
    Real uu1_r = prim_r(IVX,i);
    Real uu2_r = prim_r(IVY,i);
    Real uu3_r = prim_r(IVZ,i);

    // Transform projected 4-velocities
    Real uux_l = mx_3*uu3_l;
    Real uuy_l = my_0*uu0_l + my_1*uu1_l;
    Real uuz_l = mz_2*uu2_l;
    Real uux_r = mx_3*uu3_r;
    Real uuy_r = my_0*uu0_r + my_1*uu1_r;
    Real uuz_r = mz_2*uu2_r;

    // Set local projected 4-velocities
    prim_l(IVZ,i) = uux_l;
    prim_l(IVX,i) = uuy_l;
    prim_l(IVY,i) = uuz_l;
    prim_r(IVZ,i) = uux_r;
    prim_r(IVX,i) = uuy_r;
    prim_r(IVY,i) = uuz_r;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      // Extract metric coefficients
      const Real &g_00 = g_(I00,i);
      const Real &g_01 = g_(I01,i);
      const Real &g_02 = 0.0;
      const Real &g_03 = 0.0;
      const Real &g_10 = g_(I01,i);
      const Real &g_11 = g_(I11,i);
      const Real &g_12 = 0.0;
      const Real &g_13 = 0.0;
      const Real &g_20 = 0.0;
      const Real &g_21 = 0.0;
      const Real &g_22 = g_(I22,i);
      const Real &g_23 = 0.0;
      const Real &g_30 = 0.0;
      const Real &g_31 = 0.0;
      const Real &g_32 = 0.0;
      const Real &g_33 = g_(I33,i);
      const Real &g01 = gi_(I01,i);
      const Real &g02 = 0.0;
      const Real &g03 = 0.0;
      Real alpha = std::sqrt(-1.0/gi_(I00,i));

      // Calculate global 4-velocities
      Real tmp = g_11*uu1_l*uu1_l + 2.0*g_12*uu1_l*uu2_l + 2.0*g_13*uu1_l*uu3_l
               + g_22*uu2_l*uu2_l + 2.0*g_23*uu2_l*uu3_l
               + g_33*uu3_l*uu3_l;
      Real gamma_l = std::sqrt(1.0 + tmp);
      Real u0_l = gamma_l / alpha;
      Real u1_l = uu1_l - alpha * gamma_l * g01;
      Real u2_l = uu2_l - alpha * gamma_l * g02;
      Real u3_l = uu3_l - alpha * gamma_l * g03;
      tmp = g_11*uu1_r*uu1_r + 2.0*g_12*uu1_r*uu2_r + 2.0*g_13*uu1_r*uu3_r
          + g_22*uu2_r*uu2_r + 2.0*g_23*uu2_r*uu3_r
          + g_33*uu3_r*uu3_r;
      Real gamma_r = std::sqrt(1.0 + tmp);
      Real u0_r = gamma_r / alpha;
      Real u1_r = uu1_r - alpha * gamma_r * g01;
      Real u2_r = uu2_r - alpha * gamma_r * g02;
      Real u3_r = uu3_r - alpha * gamma_r * g03;

      // Extract global magnetic fields
      const Real &bb3_l = bb3(k,j,i);
      const Real &bb3_r = bb3(k,j,i);
      Real &bb1_l = prim_l(IBY,i);
      Real &bb2_l = prim_l(IBZ,i);
      Real &bb1_r = prim_r(IBY,i);
      Real &bb2_r = prim_r(IBZ,i);

      // Calculate global 4-magnetic fields
      Real b0_l = g_10*bb1_l*u0_l + g_11*bb1_l*u1_l + g_12*bb1_l*u2_l + g_13*bb1_l*u3_l
                + g_20*bb2_l*u0_l + g_21*bb2_l*u1_l + g_22*bb2_l*u2_l + g_23*bb2_l*u3_l
                + g_30*bb3_l*u0_l + g_31*bb3_l*u1_l + g_32*bb3_l*u2_l + g_33*bb3_l*u3_l;
      Real b1_l = (bb1_l + b0_l * u1_l) / u0_l;
      Real b2_l = (bb2_l + b0_l * u2_l) / u0_l;
      Real b3_l = (bb3_l + b0_l * u3_l) / u0_l;
      Real b0_r = g_10*bb1_r*u0_r + g_11*bb1_r*u1_r + g_12*bb1_r*u2_r + g_13*bb1_r*u3_r
                + g_20*bb2_r*u0_r + g_21*bb2_r*u1_r + g_22*bb2_r*u2_r + g_23*bb2_r*u3_r
                + g_30*bb3_r*u0_r + g_31*bb3_r*u1_r + g_32*bb3_r*u2_r + g_33*bb3_r*u3_r;
      Real b1_r = (bb1_r + b0_r * u1_r) / u0_r;
      Real b2_r = (bb2_r + b0_r * u2_r) / u0_r;
      Real b3_r = (bb3_r + b0_r * u3_r) / u0_r;

      // Transform 4-velocities
      Real ut_l = mt_0*u0_l;
      Real ux_l = mx_3*u3_l;
      Real uy_l = my_0*u0_l + my_1*u1_l;
      Real uz_l = mz_2*u2_l;
      Real ut_r = mt_0*u0_r;
      Real ux_r = mx_3*u3_r;
      Real uy_r = my_0*u0_r + my_1*u1_r;
      Real uz_r = mz_2*u2_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_3*b3_l;
      Real by_l = my_0*b0_l + my_1*b1_l;
      Real bz_l = mz_2*b2_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_3*b3_r;
      Real by_r = my_0*b0_r + my_1*b1_r;
      Real bz_r = mz_2*b2_r;

      // Set local magnetic fields
      Real bbx_l = ut_l * bx_l - ux_l * bt_l;
      Real bbx_r = ut_r * bx_r - ux_r * bt_r;
      bbx(i) = 0.5 * (bbx_l + bbx_r);
      bb1_l = ut_l * by_l - uy_l * bt_l;
      bb2_l = ut_l * bz_l - uz_l * bt_l;
      bb1_r = ut_r * by_r - uy_r * bt_r;
      bb2_r = ut_r * bz_r - uz_r * bt_r;
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
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM1/IM2/IM3 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts x1-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts x1-fluxes of B2/B3 in IBY/IBZ slots
void Coordinates::FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux)
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
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM2/IM3/IM1 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts x2-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts x2-fluxes of B3/B1 in IBY/IBZ slots
void Coordinates::FluxToGlobal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux)
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
      Real ftx = -bbx(i);
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
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM3/IM1/IM2 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts x3-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts x3-fluxes of B1/B2 in IBY/IBZ slots
void Coordinates::FluxToGlobal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux)
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
      Real ftx = -bbx(i);
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
void Coordinates::LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j,
    int i, Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3)
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
