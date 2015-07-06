// Minkowski spacetime, sinusoidal ("snake") coordinates
// Notes:
//   coordinates: t, x, y, z
//   parameters: a, k
//   metric:
//     ds^2 = -dt^2 + \alpha^2 dx^2 - 2 \beta dx dy + dy^2 + dz^2
//     alpha = \sqrt(1 + a^2 k^2 \cos^2(k x))
//     \beta = a k \cos(k x)
//   relation to Minkowski (Cartesian) coordinates:
//     t = t_m
//     x = x_m
//     y = y_m + a \sin(k x_m)
//     z = z_m

// Primary header
#include "coordinates.hpp"

// C++ headers
#include <cmath>  // cos(), sin(), sqrt()

// Athena headers
#include "../athena.hpp"           // enums, macros, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../fluid/eos/eos.hpp"    // FluidEqnofState
#include "../fluid/fluid.hpp"      // Fluid
#include "../mesh.hpp"             // MeshBlock
#include "../parameter_input.hpp"  // ParameterInput

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
  sinu_amplitude_ = pin->GetReal("coord", "a");
  sinu_wavenumber_ = pin->GetReal("coord", "k");
  const Real &a = sinu_amplitude_;
  const Real &k = sinu_wavenumber_;

  // Set pointer to host MeshBlock
  pmy_block = pb;

  // Initialize volume-averaged positions and spacings: x-direction
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

  if(pb->pmy_mesh->multilevel==true) { // calc coarse coodinates
    int cis = pb->cis; int cjs = pb->cjs; int cks = pb->cks;
    int cie = pb->cie; int cje = pb->cje; int cke = pb->cke;
    for (int i=cis-(pb->cnghost); i<=cie+(pb->cnghost); ++i) {
      coarse_x1v(i) = 0.5*(coarse_x1f(i+1) + coarse_x1f(i));
    }
    for (int i=cis-(pb->cnghost); i<=cie+(pb->cnghost)-1; ++i) {
      coarse_dx1v(i) = coarse_x1v(i+1) - coarse_x1v(i);
    }
    if (pb->block_size.nx2 == 1) {
      coarse_x2v(cjs) = 0.5*(coarse_x2f(cjs+1) + coarse_x2f(cjs));
      coarse_dx2v(cjs) = coarse_dx2f(cjs);
    } else {
      for (int j=cjs-(pb->cnghost); j<=cje+(pb->cnghost); ++j) {
        coarse_x2v(j) = 0.5*(coarse_x2f(j+1) + coarse_x2f(j));
      }
      for (int j=cjs-(pb->cnghost); j<=cje+(pb->cnghost)-1; ++j) {
        coarse_dx2v(j) = coarse_x2v(j+1) - coarse_x2v(j);
      }
    }
    if (pb->block_size.nx3 == 1) {
      coarse_x3v(cks) = 0.5*(coarse_x3f(cks+1) + coarse_x3f(cks));
      coarse_dx3v(cks) = coarse_dx3f(cks);
    } else {
      for (int k=cks-(pb->cnghost); k<=cke+(pb->cnghost); ++k) {
        coarse_x3v(k) = 0.5*(coarse_x3f(k+1) + coarse_x3f(k));
      }
      for (int k=cks-(pb->cnghost); k<=cke+(pb->cnghost)-1; ++k) {
        coarse_dx3v(k) = coarse_x3v(k+1) - coarse_x3v(k);
      }
    }
  }

  // Allocate arrays for intermediate geometric quantities: x-direction
  int n_cells_1 = pb->block_size.nx1 + 2*NGHOST;
  coord_src_i1_.NewAthenaArray(n_cells_1);
  metric_cell_i1_.NewAthenaArray(n_cells_1);
  metric_cell_i2_.NewAthenaArray(n_cells_1);
  metric_face1_i1_.NewAthenaArray(n_cells_1);
  metric_face1_i2_.NewAthenaArray(n_cells_1);
  metric_face2_i1_.NewAthenaArray(n_cells_1);
  metric_face2_i2_.NewAthenaArray(n_cells_1);
  metric_face3_i1_.NewAthenaArray(n_cells_1);
  metric_face3_i2_.NewAthenaArray(n_cells_1);
  trans_face1_i2_.NewAthenaArray(n_cells_1);
  trans_face2_i1_.NewAthenaArray(n_cells_1);
  trans_face2_i2_.NewAthenaArray(n_cells_1);
  trans_face3_i2_.NewAthenaArray(n_cells_1);

  // Calculate intermediate geometric quantities: x-direction
  #pragma simd
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST; ++i)
  {
    // Useful quantities
    Real r_c = x1v(i);
    Real r_m = x1f(i);
    Real r_p = x1f(i+1);
    Real sin_2m = std::sin(2.0*k*r_m);
    Real sin_2p = std::sin(2.0*k*r_p);
    Real cos_c = std::cos(k*r_c);
    Real cos_m = std::cos(k*r_m);
    Real cos_p = std::cos(k*r_p);
    Real alpha_sq_c = 1.0 + SQR(a)*SQR(k) * SQR(cos_c);
    Real alpha_sq_m = 1.0 + SQR(a)*SQR(k) * SQR(cos_m);
    Real alpha_c = std::sqrt(alpha_sq_c);
    Real beta_c = a*k * cos_c;
    Real beta_m = a*k * cos_m;
    Real beta_p = a*k * cos_p;

    // Source terms
    coord_src_i1_(i) = (beta_m - beta_p) / dx1f(i);

    // Cell-centered metric
    metric_cell_i1_(i) = alpha_sq_c;
    metric_cell_i2_(i) = beta_c;

    // Face-centered metric
    metric_face1_i1_(i) = alpha_sq_m;
    metric_face1_i2_(i) = beta_m;
    metric_face2_i1_(i) = alpha_sq_c;
    metric_face2_i2_(i) = beta_c;
    metric_face3_i1_(i) = alpha_sq_c;
    metric_face3_i2_(i) = beta_c;

    // Coordinate transformations
    trans_face1_i2_(i) = beta_m;
    trans_face2_i1_(i) = alpha_c;
    trans_face2_i2_(i) = beta_c;
    trans_face3_i2_(i) = beta_m;
  }
}

//--------------------------------------------------------------------------------------

// Destructor
Coordinates::~Coordinates()
{
  DeleteBasicCoordinates();

  coord_src_i1_.DeleteAthenaArray();
  metric_cell_i1_.DeleteAthenaArray();
  metric_cell_i2_.DeleteAthenaArray();
  metric_face1_i1_.DeleteAthenaArray();
  metric_face1_i2_.DeleteAthenaArray();
  metric_face2_i1_.DeleteAthenaArray();
  metric_face2_i2_.DeleteAthenaArray();
  metric_face3_i1_.DeleteAthenaArray();
  metric_face3_i2_.DeleteAthenaArray();
  trans_face1_i2_.DeleteAthenaArray();
  trans_face2_i1_.DeleteAthenaArray();
  trans_face2_i2_.DeleteAthenaArray();
  trans_face3_i2_.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------

// Function for computing cell volumes
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
// Outputs:
//   volumes: 1D array of cell volumes
// Notes:
//   \Delta V = \Delta x * \Delta y * \Delta z
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

Real Coordinates::GetCellVolume(const int k, const int j, const int i)
{
  return dx1f(i)*dx2f(j)*dx3f(k);
}


//--------------------------------------------------------------------------------------

// Function for computing areas orthogonal to x
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to x
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
//   k,j: z- and y-indices
//   il,iu: x-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to y
// Notes:
//   \Delta A = \Delta x * \Delta z
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
//   k,j: z- and y-indices
//   il,iu: x-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to z
// Notes:
//   \Delta A = \Delta x * \Delta y
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


Real Coordinates::GetFace1Area(const int k, const int j, const int i)
{
  return dx2f(j)*dx3f(k);
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the x-direction
// Inputs:
//   k,j: z- and y-indices (unused)
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
//   il,iu: x-index bounds
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
//   il,iu: x-index bounds
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

// Function for computing widths of cells in the x-direction
// Inputs:
//   k,j: z- and y-indices (unused)
//   i: x-index
// Outputs:
//   returned value: x-width of cell (i,j,k)
// Notes:
//   \Delta W >= \Delta x
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
//   j,i: y- and x-indices (unused)
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
//   bcc: 3D array of cell-centered magnetic fields
// Outputs:
//   cons: source terms added to k,j-slice of 3D array of conserved variables
// Notes:
//   all source terms computed in this function
void Coordinates::CoordSrcTermsX1(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flux, const AthenaArray<Real> &prim,
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  // Extract ratio of specific heats
  const Real gamma_adi = pmy_block->pfluid->pf_eos->GetGamma();
  const Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Go through cells
  #pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; ++i)
  {
    // Extract geometric quantities
    const Real g00 = -1.0;
    const Real &g11 = metric_cell_i1_(i);
    const Real g12 = -metric_cell_i2_(i);
    const Real g21 = -metric_cell_i2_(i);
    const Real g22 = 1.0;
    const Real g33 = 1.0;
    const Real &gamma2_11 = coord_src_i1_(i);

    // Extract primitives
    const Real &rho = prim(IDN,k,j,i);
    const Real &pgas = prim(IEN,k,j,i);
    const Real &v1 = prim(IVX,k,j,i);
    const Real &v2 = prim(IVY,k,j,i);
    const Real &v3 = prim(IVZ,k,j,i);

    // Calculate 4-velocity
    Real u0 = std::sqrt(-1.0 /
        (g00 + g11*v1*v1 + g12*v1*v2 + g21*v2*v1 + g22*v2*v2 + g33*v3*v3));
    Real u1 = u0 * v1;
    Real u2 = u0 * v2;
    Real u3 = u0 * v3;
    Real u_2 = g21*u1 + g22*u2;

    // Extract and calculate magnetic field
    Real bcon1 = 0.0, bcov2 = 0.0, b_sq = 0.0;
    if (MAGNETIC_FIELDS_ENABLED)
    {
      const Real &b1 = bcc(IB1,k,j,i);
      const Real &b2 = bcc(IB2,k,j,i);
      const Real &b3 = bcc(IB3,k,j,i);
      Real bcon0 = g11*b1*u1 + g12*b1*u2 + g21*b2*u1 + g22*b2*u2 + g33*b3*u3;
      bcon1 = (b1 + bcon0 * u1) / u0;
      Real bcon2 = (b2 + bcon0 * u2) / u0;
      Real bcon3 = (b3 + bcon0 * u3) / u0;
      Real bcov0 = g00*bcon0;
      Real bcov1 = g11*bcon1 + g12*bcon2;
      bcov2 = g21*bcon1 + g22*bcon2;
      Real bcov3 = g33*bcon3;
      b_sq = bcov0*bcon0 + bcov1*bcon1 + bcov2*bcon2 + bcov3*bcon3;
    }

    // Calculate stress-energy tensor
    Real w = rho + gamma_adi_red * pgas + b_sq;
    Real t1_2 = w*u1*u_2 - bcon1*bcov2;

    // Calculate source terms
    Real s1 = gamma2_11 * t1_2;

    // Extract conserved quantities
    Real &m1 = cons(IM1,k,j,i);

    // Add source terms to conserved quantities
    m1 += dt * s1;
  }
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
//   not using this function
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
//   not using this function
void Coordinates::CoordSrcTermsX3(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flux_k, const AthenaArray<Real> &flux_kp1,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing cell-centered metric coefficients
// Inputs:
//   k,j: z- and y-indices (unused)
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
    // Extract geometric quantities
    const Real &alpha_sq = metric_cell_i1_(i);
    const Real &beta = metric_cell_i2_(i);

    // Extract metric terms
    Real &g00 = g(I00,i);
    Real &g11 = g(I11,i);
    Real &g12 = g(I12,i);
    Real &g22 = g(I22,i);
    Real &g33 = g(I33,i);
    Real &gi00 = g_inv(I00,i);
    Real &gi11 = g_inv(I11,i);
    Real &gi12 = g_inv(I12,i);
    Real &gi22 = g_inv(I22,i);
    Real &gi33 = g_inv(I33,i);

    // Set metric terms
    g00 = -1.0;
    g11 = alpha_sq;
    g12 = -beta;
    g22 = 1.0;
    g33 = 1.0;
    gi00 = -1.0;
    gi11 = 1.0;
    gi12 = beta;
    gi22 = alpha_sq;
    gi33 = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: x-interface
// Inputs:
//   k,j: z- and y-indices (unused)
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
    // Extract geometric quantities
    const Real &alpha_sq = metric_face1_i1_(i);
    const Real &beta = metric_face1_i2_(i);

    // Extract metric terms
    Real &g00 = g(I00,i);
    Real &g11 = g(I11,i);
    Real &g12 = g(I12,i);
    Real &g22 = g(I22,i);
    Real &g33 = g(I33,i);
    Real &gi00 = g_inv(I00,i);
    Real &gi11 = g_inv(I11,i);
    Real &gi12 = g_inv(I12,i);
    Real &gi22 = g_inv(I22,i);
    Real &gi33 = g_inv(I33,i);

    // Set metric terms
    g00 = -1.0;
    g11 = alpha_sq;
    g12 = -beta;
    g22 = 1.0;
    g33 = 1.0;
    gi00 = -1.0;
    gi11 = 1.0;
    gi12 = beta;
    gi22 = alpha_sq;
    gi33 = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: y-interface
// Inputs:
//   k,j: z- and y-indices (unused)
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
    // Extract geometric quantities
    const Real &alpha_sq = metric_face2_i1_(i);
    const Real &beta = metric_face2_i2_(i);

    // Extract metric terms
    Real &g00 = g(I00,i);
    Real &g11 = g(I11,i);
    Real &g12 = g(I12,i);
    Real &g22 = g(I22,i);
    Real &g33 = g(I33,i);
    Real &gi00 = g_inv(I00,i);
    Real &gi11 = g_inv(I11,i);
    Real &gi12 = g_inv(I12,i);
    Real &gi22 = g_inv(I22,i);
    Real &gi33 = g_inv(I33,i);

    // Set metric terms
    g00 = -1.0;
    g11 = alpha_sq;
    g12 = -beta;
    g22 = 1.0;
    g33 = 1.0;
    gi00 = -1.0;
    gi11 = 1.0;
    gi12 = beta;
    gi22 = alpha_sq;
    gi33 = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: z-interface
// Inputs:
//   k,j: z- and y-indices (unused)
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
    // Extract geometric quantities
    const Real &alpha_sq = metric_face3_i1_(i);
    const Real &beta = metric_face3_i2_(i);

    // Extract metric terms
    Real &g00 = g(I00,i);
    Real &g11 = g(I11,i);
    Real &g12 = g(I12,i);
    Real &g22 = g(I22,i);
    Real &g33 = g(I33,i);
    Real &gi00 = g_inv(I00,i);
    Real &gi11 = g_inv(I11,i);
    Real &gi12 = g_inv(I12,i);
    Real &gi22 = g_inv(I22,i);
    Real &gi33 = g_inv(I33,i);

    // Set metric terms
    g00 = -1.0;
    g11 = alpha_sq;
    g12 = -beta;
    g22 = 1.0;
    g33 = 1.0;
    gi00 = -1.0;
    gi11 = 1.0;
    gi12 = beta;
    gi22 = alpha_sq;
    gi33 = 1.0;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: x-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
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
  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract geometric quantities
    const Real g00 = -1.0;
    const Real &g11 = metric_face1_i1_(i);
    const Real g12 = -metric_face1_i2_(i);
    const Real g22 = 1.0;
    const Real g33 = 1.0;
    const Real mt0 = 1.0;
    const Real mx1 = 1.0;
    const Real my1 = -trans_face1_i2_(i);
    const Real my2 = 1.0;
    const Real mz3 = 1.0;

    // Extract global 3-velocities
    Real &v1l = prim_left(IVX,i);
    Real &v2l = prim_left(IVY,i);
    Real &v3l = prim_left(IVZ,i);
    Real &v1r = prim_right(IVX,i);
    Real &v2r = prim_right(IVY,i);
    Real &v3r = prim_right(IVZ,i);

    // Construct global 4-velocities
    Real u0l = std::sqrt(-1.0 /
        (g00 + g11*v1l*v1l + 2.0*g12*v1l*v2l + g22*v2l*v2l + g33*v3l*v3l));
    Real u1l = u0l * v1l;
    Real u2l = u0l * v2l;
    Real u3l = u0l * v3l;
    Real u0r = std::sqrt(-1.0 /
        (g00 + g11*v1r*v1r + 2.0*g12*v1r*v2r + g22*v2r*v2r + g33*v3r*v3r));
    Real u1r = u0r * v1r;
    Real u2r = u0r * v2r;
    Real u3r = u0r * v3r;

    // Transform 4-velocities
    Real utl = mt0*u0l;
    Real uxl = mx1*u1l;
    Real uyl = my1*u1l + my2*u2l;
    Real uzl = mz3*u3l;
    Real utr = mt0*u0r;
    Real uxr = mx1*u1r;
    Real uyr = my1*u1r + my2*u2r;
    Real uzr = mz3*u3r;

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
      const Real &b1 = b1_vals(k,j,i);
      Real &b2l = prim_left(IBY,i);
      Real &b3l = prim_left(IBZ,i);
      Real &b2r = prim_right(IBY,i);
      Real &b3r = prim_right(IBZ,i);

      // Construct global contravariant magnetic fields
      Real bcon0l = g11*b1*u1l + g12*(b1*u2l+b2l*u1l) + g22*b2l*u2l + g33*b3l*u3l;
      Real bcon1l = (b1 + bcon0l * u1l) / u0l;
      Real bcon2l = (b2l + bcon0l * u2l) / u0l;
      Real bcon3l = (b3l + bcon0l * u3l) / u0l;
      Real bcon0r = g11*b1*u1r + g12*(b1*u2r+b2r*u1r) + g22*b2r*u2r + g33*b3r*u3r;
      Real bcon1r = (b1 + bcon0r * u1r) / u0r;
      Real bcon2r = (b2r + bcon0r * u2r) / u0r;
      Real bcon3r = (b3r + bcon0r * u3r) / u0r;

      // Transform contravariant magnetic fields
      Real bcontl = mt0*bcon0l;
      Real bconxl = mx1*bcon1l;
      Real bconyl = my1*bcon1l + my2*bcon2l;
      Real bconzl = mz3*bcon3l;
      Real bcontr = mt0*bcon0r;
      Real bconxr = mx1*bcon1r;
      Real bconyr = my1*bcon1r + my2*bcon2r;
      Real bconzr = mz3*bcon3r;

      // Set local magnetic fields
      Real bxl = utl * bconxl - uxl * bcontl;
      Real bxr = utr * bconxr - uxr * bcontr;
      bx(i) = 0.5 * (bxl + bxr);
      b2l = utl * bconyl - uyl * bcontl;
      b3l = utl * bconzl - uzl * bcontl;
      b2r = utr * bconyr - uyr * bcontr;
      b3r = utr * bconzr - uzr * bcontr;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: y-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
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
  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract geometric quantities
    const Real g00 = -1.0;
    const Real &g11 = metric_face2_i1_(i);
    const Real g12 = -metric_face2_i2_(i);
    const Real g22 = 1.0;
    const Real g33 = 1.0;
    const Real mt0 = 1.0;
    const Real mx2 = 1.0 / trans_face2_i1_(i);
    const Real my3 = 1.0;
    const Real &mz1 = trans_face2_i1_(i);
    const Real mz2 = -trans_face2_i2_(i) / trans_face2_i1_(i);

    // Extract global 3-velocities
    Real &v1l = prim_left(IVX,i);
    Real &v2l = prim_left(IVY,i);
    Real &v3l = prim_left(IVZ,i);
    Real &v1r = prim_right(IVX,i);
    Real &v2r = prim_right(IVY,i);
    Real &v3r = prim_right(IVZ,i);

    // Construct global 4-velocities
    Real u0l = std::sqrt(-1.0 /
        (g00 + g11*v1l*v1l + 2.0*g12*v1l*v2l + g22*v2l*v2l + g33*v3l*v3l));
    Real u1l = u0l * v1l;
    Real u2l = u0l * v2l;
    Real u3l = u0l * v3l;
    Real u0r = std::sqrt(-1.0 /
        (g00 + g11*v1r*v1r + 2.0*g12*v1r*v2r + g22*v2r*v2r + g33*v3r*v3r));
    Real u1r = u0r * v1r;
    Real u2r = u0r * v2r;
    Real u3r = u0r * v3r;

    // Transform 4-velocities
    Real utl = mt0*u0l;
    Real uxl = mx2*u2l;
    Real uyl = my3*u3l;
    Real uzl = mz1*u1l + mz2*u2l;
    Real utr = mt0*u0r;
    Real uxr = mx2*u2r;
    Real uyr = my3*u3r;
    Real uzr = mz1*u1r + mz2*u2r;

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
      const Real &b2 = b2_vals(k,j,i);
      Real &b3l = prim_left(IBY,i);
      Real &b1l = prim_left(IBZ,i);
      Real &b3r = prim_right(IBY,i);
      Real &b1r = prim_right(IBZ,i);

      // Construct global contravariant magnetic fields
      Real bcon0l = g11*b1l*u1l + g12*(b1l*u2l+b2*u1l) + g22*b2*u2l + g33*b3l*u3l;
      Real bcon1l = (b1l + bcon0l * u1l) / u0l;
      Real bcon2l = (b2 + bcon0l * u2l) / u0l;
      Real bcon3l = (b3l + bcon0l * u3l) / u0l;
      Real bcon0r = g11*b1r*u1r + g12*(b1r*u2r+b2*u1r) + g22*b2*u2r + g33*b3r*u3r;
      Real bcon1r = (b1r + bcon0r * u1r) / u0r;
      Real bcon2r = (b2 + bcon0r * u2r) / u0r;
      Real bcon3r = (b3r + bcon0r * u3r) / u0r;

      // Transform contravariant magnetic fields
      Real bcontl = mt0*bcon0l;
      Real bconxl = mx2*bcon2l;
      Real bconyl = my3*bcon3l;
      Real bconzl = mz1*bcon1l + mz2*bcon2l;
      Real bcontr = mt0*bcon0r;
      Real bconxr = mx2*bcon2r;
      Real bconyr = my3*bcon3r;
      Real bconzr = mz1*bcon1r + mz2*bcon2r;

      // Set local magnetic fields
      Real bxl = utl * bconxl - uxl * bcontl;
      Real bxr = utr * bconxr - uxr * bcontr;
      bx(i) = 0.5 * (bxl + bxr);
      b3l = utl * bconyl - uyl * bcontl;
      b1l = utl * bconzl - uzl * bcontl;
      b3r = utr * bconyr - uyr * bcontr;
      b1r = utr * bconzr - uzr * bcontr;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: z-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
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
  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract geometric quantities
    const Real g00 = -1.0;
    const Real &g11 = metric_face3_i1_(i);
    const Real g12 = -metric_face3_i2_(i);
    const Real g22 = 1.0;
    const Real g33 = 1.0;
    const Real mt0 = 1.0;
    const Real mx3 = 1.0;
    const Real my1 = 1.0;
    const Real mz1 = -trans_face3_i2_(i);
    const Real mz2 = 1.0;

    // Extract global 3-velocities
    Real &v1l = prim_left(IVX,i);
    Real &v2l = prim_left(IVY,i);
    Real &v3l = prim_left(IVZ,i);
    Real &v1r = prim_right(IVX,i);
    Real &v2r = prim_right(IVY,i);
    Real &v3r = prim_right(IVZ,i);

    // Construct global 4-velocities
    Real u0l = std::sqrt(-1.0 /
        (g00 + g11*v1l*v1l + 2.0*g12*v1l*v2l + g22*v2l*v2l + g33*v3l*v3l));
    Real u1l = u0l * v1l;
    Real u2l = u0l * v2l;
    Real u3l = u0l * v3l;
    Real u0r = std::sqrt(-1.0 /
        (g00 + g11*v1r*v1r + 2.0*g12*v1r*v2r + g22*v2r*v2r + g33*v3r*v3r));
    Real u1r = u0r * v1r;
    Real u2r = u0r * v2r;
    Real u3r = u0r * v3r;

    // Transform 4-velocities
    Real utl = mt0*u0l;
    Real uxl = mx3*u3l;
    Real uyl = my1*u1l;
    Real uzl = mz1*u1l + mz2*u2l;
    Real utr = mt0*u0r;
    Real uxr = mx3*u3r;
    Real uyr = my1*u1r;
    Real uzr = mz1*u1r + mz2*u2r;

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
      const Real &b3 = b3_vals(k,j,i);
      Real &b1l = prim_left(IBY,i);
      Real &b2l = prim_left(IBZ,i);
      Real &b1r = prim_right(IBY,i);
      Real &b2r = prim_right(IBZ,i);

      // Construct global contravariant magnetic fields
      Real bcon0l = g11*b1l*u1l + g12*(b1l*u2l+b2l*u1l) + g22*b2l*u2l + g33*b3*u3l;
      Real bcon1l = (b1l + bcon0l * u1l) / u0l;
      Real bcon2l = (b2l + bcon0l * u2l) / u0l;
      Real bcon3l = (b3 + bcon0l * u3l) / u0l;
      Real bcon0r = g11*b1r*u1r + g12*(b1r*u2r+b2r*u1r) + g22*b2r*u2r + g33*b3*u3r;
      Real bcon1r = (b1r + bcon0r * u1r) / u0r;
      Real bcon2r = (b2r + bcon0r * u2r) / u0r;
      Real bcon3r = (b3 + bcon0r * u3r) / u0r;

      // Transform contravariant magnetic fields
      Real bcontl = mt0*bcon0l;
      Real bconxl = mx3*bcon3l;
      Real bconyl = my1*bcon1l;
      Real bconzl = mz1*bcon1l + mz2*bcon2l;
      Real bcontr = mt0*bcon0r;
      Real bconxr = mx3*bcon3r;
      Real bconyr = my1*bcon1r;
      Real bconzr = mz1*bcon1r + mz2*bcon2r;

      // Set local magnetic fields
      Real bxl = utl * bconxl - uxl * bcontl;
      Real bxr = utr * bconxr - uxr * bcontr;
      bx(i) = 0.5 * (bxl + bxr);
      b1l = utl * bconyl - uyl * bcontl;
      b2l = utl * bconzl - uzl * bcontl;
      b1r = utr * bconyr - uyr * bcontr;
      b2r = utr * bconzr - uzr * bcontr;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: x-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates (unused)
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
  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract geometric quantities
    const Real g00 = -1.0;
    const Real &g11 = metric_face1_i1_(i);
    const Real g12 = -metric_face1_i2_(i);
    const Real g22 = 1.0;
    const Real g33 = 1.0;
    const Real m0t = 1.0;
    const Real m1x = 1.0;
    const Real &m2x = trans_face1_i2_(i);
    const Real m2y = 1.0;
    const Real m3z = 1.0;

    // Extract local conserved quantities and fluxes
    const Real dx = flux(IDN,i);
    const Real txt = flux(IEN,i);
    const Real txx = flux(IM1,i);
    const Real txy = flux(IM2,i);
    const Real txz = flux(IM3,i);

    // Transform stress-energy tensor
    const Real tcon_10 = m1x * m0t * txt;
    const Real tcon_11 = m1x * m1x * txx;
    const Real tcon_12 = m1x * (m2x * txx + m2y * txy);
    const Real tcon_13 = m1x * m3z * txz;

    // Extract global fluxes
    Real &d1 = flux(IDN,i);
    Real &t10 = flux(IEN,i);
    Real &t11 = flux(IM1,i);
    Real &t12 = flux(IM2,i);
    Real &t13 = flux(IM3,i);

    // Set fluxes
    d1 = m1x*dx;
    t10 = g00*tcon_10;
    t11 = g11*tcon_11 + g12*tcon_12;
    t12 = g12*tcon_11 + g22*tcon_12;
    t13 = g33*tcon_13;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      const Real fyx = flux(IBY,i);
      const Real fzx = flux(IBZ,i);
      Real &f21 = flux(IBY,i);
      Real &f31 = flux(IBZ,i);
      f21 = m1x * m2y * fyx;
      f31 = m1x * m3z * fzx;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: y-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates (unused)
//   bx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
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
  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract geometric quantities
    const Real g00 = -1.0;
    const Real &g11 = metric_face2_i1_(i);
    const Real g12 = -metric_face2_i2_(i);
    const Real g22 = 1.0;
    const Real g33 = 1.0;
    const Real m0t = 1.0;
    const Real m1x = trans_face2_i2_(i) / trans_face2_i1_(i);
    const Real m1z = 1.0 / trans_face2_i1_(i);
    const Real &m2x = trans_face2_i1_(i);
    const Real m3y = 1.0;

    // Extract local conserved quantities and fluxes
    const Real dx = flux(IDN,i);
    const Real txt = flux(IEN,i);
    const Real txx = flux(IM2,i);
    const Real txy = flux(IM3,i);
    const Real txz = flux(IM1,i);

    // Transform stress-energy tensor
    const Real tcon_20 = m2x * m0t * txt;
    const Real tcon_21 = m2x * (m1x * txx + m1z * txz);
    const Real tcon_22 = m2x * m2x * txx;
    const Real tcon_23 = m2x * m3y * txy;

    // Extract global fluxes
    Real &d2 = flux(IDN,i);
    Real &t20 = flux(IEN,i);
    Real &t21 = flux(IM1,i);
    Real &t22 = flux(IM2,i);
    Real &t23 = flux(IM3,i);

    // Set fluxes
    d2 = m2x*dx;
    t20 = g00*tcon_20;
    t21 = g11*tcon_21 + g12*tcon_22;
    t22 = g12*tcon_21 + g22*tcon_22;
    t23 = g33*tcon_23;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      const Real fyx = flux(IBY,i);
      const Real fzx = flux(IBZ,i);
      Real &f32 = flux(IBY,i);
      Real &f12 = flux(IBZ,i);
      f32 = m3y * m2x * fyx;
      f12 = m2x * m1z * fzx;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: z-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates (unused)
//   bx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
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
  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract geometric quantities
    const Real g00 = -1.0;
    const Real &g11 = metric_face3_i1_(i);
    const Real g12 = -metric_face3_i2_(i);
    const Real g22 = 1.0;
    const Real g33 = 1.0;
    const Real m0t = 1.0;
    const Real m1y = 1.0;
    const Real &m2y = trans_face3_i2_(i);
    const Real m2z = 1.0;
    const Real m3x = 1.0;

    // Extract local conserved quantities and fluxes
    const Real dx = flux(IDN,i);
    const Real txt = flux(IEN,i);
    const Real txx = flux(IM3,i);
    const Real txy = flux(IM1,i);
    const Real txz = flux(IM2,i);

    // Transform stress-energy tensor
    const Real tcon_30 = m3x * m0t * txt;
    const Real tcon_31 = m3x * m1y * txy;
    const Real tcon_32 = m3x * (m2y * txy + m2z * txz);
    const Real tcon_33 = m3x * m3x * txx;

    // Extract global fluxes
    Real &d3 = flux(IDN,i);
    Real &t30 = flux(IEN,i);
    Real &t31 = flux(IM1,i);
    Real &t32 = flux(IM2,i);
    Real &t33 = flux(IM3,i);

    // Set fluxes
    d3 = m3x*dx;
    t30 = g00*tcon_30;
    t31 = g11*tcon_31 + g12*tcon_32;
    t32 = g12*tcon_31 + g22*tcon_32;
    t33 = g33*tcon_33;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      const Real fyx = flux(IBY,i);
      const Real fzx = flux(IBZ,i);
      Real &f13 = flux(IBY,i);
      Real &f23 = flux(IBZ,i);
      f13 = m1y * m3x * fyx;
      f23 = m3x * (m2y * fyx + m2z * fzx);
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
  const Real &a = sinu_amplitude_;
  const Real &k = sinu_wavenumber_;
  Real ax = a1;
  Real ay = a2 - a * std::sin(k * a1);
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
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
void Coordinates::TransformVectorCell(
    Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  const Real &beta = metric_cell_i2_(i);
  *pa0 = at;
  *pa1 = ax;
  *pa2 = beta * ax + ay;
  *pa3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: x-interface
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j: z- and y-indices (unused)
//   i: x-index
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
void Coordinates::TransformVectorFace1(
    Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  const Real &beta = metric_face1_i2_(i);
  *pa0 = at;
  *pa1 = ax;
  *pa2 = beta * ax + ay;
  *pa3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: y-interface
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j: z- and y-indices (unused)
//   i: x-index
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
void Coordinates::TransformVectorFace2(
    Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  const Real &beta = metric_face2_i2_(i);
  *pa0 = at;
  *pa1 = ax;
  *pa2 = beta * ax + ay;
  *pa3 = az;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Minkowski to global: z-interface
// Inputs:
//   at,ax,ay,az: upper 4-vector components in Minkowski coordinates
//   k,j: z- and y-indices (unused)
//   i: x-index
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in global coordinates
void Coordinates::TransformVectorFace3(
    Real at, Real ax, Real ay, Real az, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  const Real &beta = metric_face3_i2_(i);
  *pa0 = at;
  *pa1 = ax;
  *pa2 = beta * ax + ay;
  *pa3 = az;
  return;
}
