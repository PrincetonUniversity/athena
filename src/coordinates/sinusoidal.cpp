// Minkowski spacetime, sinusoidal ("snake") coordinates
// Notes:
//   coordinates: t, x, y, z
//   parameters: a (aa in code), k (kk in code)
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
#include "../mesh.hpp"             // MeshBlock
#include "../parameter_input.hpp"  // ParameterInput
#include "../hydro/hydro.hpp"
#include "../hydro/eos/eos.hpp"

//--------------------------------------------------------------------------------------

// Constructor
// Inputs:
//   pmb: pointer to block containing this grid
//   pin: pointer to runtime inputs
Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin)
{
  // Set pointer to host MeshBlock and note active zone boundaries
  pmy_block = pmb;
  int is = pmb->is;
  int ie = pmb->ie;
  int js = pmb->js;
  int je = pmb->je;
  int ks = pmb->ks;
  int ke = pmb->ke;

  // Set face centered positions and distances
  AllocateAndSetBasicCoordinates();

  // Set parameters
  sinu_amplitude_ = pin->GetReal("coord", "a");
  sinu_wavenumber_ = pin->GetReal("coord", "k");
  const Real &aa = sinu_amplitude_;
  const Real &kk = sinu_wavenumber_;

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

  // Allocate arrays for intermediate geometric quantities: x-direction
  int n_cells_1 = pmb->block_size.nx1 + 2*NGHOST;
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
  g_.NewAthenaArray(NMETRIC, n_cells_1);
  gi_.NewAthenaArray(NMETRIC, n_cells_1);

  // Calculate intermediate geometric quantities: x-direction
  #pragma simd
  for (int i = is-NGHOST; i <= ie+NGHOST; ++i)
  {
    // Useful quantities
    Real r_c = x1v(i);
    Real r_m = x1f(i);
    Real r_p = x1f(i+1);
    Real sin_2m = std::sin(2.0*kk*r_m);
    Real sin_2p = std::sin(2.0*kk*r_p);
    Real cos_c = std::cos(kk*r_c);
    Real cos_m = std::cos(kk*r_m);
    Real cos_p = std::cos(kk*r_p);
    Real alpha_sq_c = 1.0 + SQR(aa)*SQR(kk) * SQR(cos_c);
    Real alpha_sq_m = 1.0 + SQR(aa)*SQR(kk) * SQR(cos_m);
    Real alpha_c = std::sqrt(alpha_sq_c);
    Real beta_c = aa*kk * cos_c;
    Real beta_m = aa*kk * cos_m;
    Real beta_p = aa*kk * cos_p;

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
  g_.DeleteAthenaArray();
  gi_.DeleteAthenaArray();
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
//   k,j: z- and y-indices
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
//   k,j: z- and y-indices
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
//   k,j: z- and y-indices
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
    lengths(i) = dx2f(j);
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
    lengths(i) = dx3f(k);
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
//   bb_cc: 3D array of cell-centered magnetic fields
// Outputs:
//   cons: source terms added to k,j-slice of 3D array of conserved variables
// Notes:
//   all source terms computed in this function
void Coordinates::CoordSrcTermsX1(const int k, const int j, const Real dt,
  const AthenaArray<Real> &flux, const AthenaArray<Real> &prim,
  const AthenaArray<Real> &bb_cc, AthenaArray<Real> &cons)
{
  // Extract ratio of specific heats
  const Real gamma_adi = pmy_block->phydro->pf_eos->GetGamma();

  // Calculate metric coefficients
  CellMetric(k, j, pmy_block->is, pmy_block->ie, g_, gi_);

  // Go through cells
  #pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; ++i)
  {
    // Extract metric coefficients
    const Real &g_00 = g_(I00,i);
    const Real &g_01 = 0.0;
    const Real &g_02 = 0.0;
    const Real &g_03 = 0.0;
    const Real &g_10 = 0.0;
    const Real &g_11 = g_(I11,i);
    const Real &g_12 = g_(I12,i);
    const Real &g_13 = 0.0;
    const Real &g_20 = 0.0;
    const Real &g_21 = g_(I12,i);
    const Real &g_22 = g_(I22,i);
    const Real &g_23 = 0.0;
    const Real &g_30 = 0.0;
    const Real &g_31 = 0.0;
    const Real &g_32 = 0.0;
    const Real &g_33 = g_(I33,i);
    const Real &g01 = 0.0;
    const Real &g02 = 0.0;
    const Real &g03 = 0.0;
    Real alpha = std::sqrt(-1.0/gi_(I00,i));

    // Extract primitives
    const Real &rho = prim(IDN,k,j,i);
    const Real &pgas = prim(IEN,k,j,i);
    const Real &uu1 = prim(IVX,k,j,i);
    const Real &uu2 = prim(IVY,k,j,i);
    const Real &uu3 = prim(IVZ,k,j,i);

    // Calculate 4-velocity
    Real tmp = g_11*uu1*uu1 + 2.0*g_12*uu1*uu2 + 2.0*g_13*uu1*uu3
             + g_22*uu2*uu2 + 2.0*g_23*uu2*uu3
             + g_33*uu3*uu3;
    Real gamma = std::sqrt(1.0 + tmp);
    Real u0 = gamma / alpha;
    Real u1 = uu1 - alpha * gamma * g01;
    Real u2 = uu2 - alpha * gamma * g02;
    Real u3 = uu3 - alpha * gamma * g03;
    Real u_0 = g_00*u0 + g_01*u1 + g_02*u2 + g_03*u3;
    Real u_1 = g_10*u0 + g_11*u1 + g_12*u2 + g_13*u3;
    Real u_2 = g_20*u0 + g_21*u1 + g_22*u2 + g_23*u3;
    Real u_3 = g_30*u0 + g_31*u1 + g_32*u2 + g_33*u3;

    // Extract and calculate magnetic field
    Real b0 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0;
    Real b_0 = 0.0, b_1 = 0.0, b_2 = 0.0, b_3 = 0.0;
    Real b_sq = 0.0;
    if (MAGNETIC_FIELDS_ENABLED)
    {
      const Real &bb1 = bb_cc(IB1,k,j,i);
      const Real &bb2 = bb_cc(IB2,k,j,i);
      const Real &bb3 = bb_cc(IB3,k,j,i);
      b0 =
            g_10*bb1*u0 + g_11*bb1*u1 + g_12*bb1*u2 + g_13*bb1*u3
          + g_20*bb2*u0 + g_21*bb2*u1 + g_22*bb2*u2 + g_23*bb2*u3
          + g_30*bb3*u0 + g_31*bb3*u1 + g_32*bb3*u2 + g_33*bb3*u3;
      b1 = (bb1 + b0 * u1) / u0;
      b2 = (bb2 + b0 * u2) / u0;
      b3 = (bb3 + b0 * u3) / u0;
      b_0 = g_00*b0 + g_01*b1 + g_02*b2 + g_03*b3;
      b_1 = g_10*b0 + g_11*b1 + g_12*b2 + g_13*b3;
      b_2 = g_20*b0 + g_21*b1 + g_22*b2 + g_23*b3;
      b_3 = g_30*b0 + g_31*b1 + g_32*b2 + g_33*b3;
      b_sq = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3;
    }

    // Calculate stress-energy tensor
    Real wtot = rho + gamma_adi/(gamma_adi-1.0) * pgas + b_sq;
    Real ptot = pgas + 0.5*b_sq;
    Real t1_2 = wtot*u1*u_2 - b1*b_2;

    // Calculate connection coefficients
    const Real &gamma2_11 = coord_src_i1_(i);

    // Calculate source terms
    Real s_1 = gamma2_11*t1_2;

    // Extract conserved quantities
    Real &m_1 = cons(IM1,k,j,i);

    // Add source terms to conserved quantities
    m_1 += dt * s_1;
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

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract transformation coefficients
    const Real mt_0 = 1.0;
    const Real mx_1 = 1.0;
    const Real my_1 = -trans_face1_i2_(i);
    const Real my_2 = 1.0;
    const Real mz_3 = 1.0;

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
    Real uux_l = mx_1*uu1_l;
    Real uuy_l = my_1*uu1_l + my_2*uu2_l;
    Real uuz_l = mz_3*uu3_l;
    Real uux_r = mx_1*uu1_r;
    Real uuy_r = my_1*uu1_r + my_2*uu2_r;
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
      const Real &g_01 = 0.0;
      const Real &g_02 = 0.0;
      const Real &g_03 = 0.0;
      const Real &g_10 = 0.0;
      const Real &g_11 = g_(I11,i);
      const Real &g_12 = g_(I12,i);
      const Real &g_13 = 0.0;
      const Real &g_20 = 0.0;
      const Real &g_21 = g_(I12,i);
      const Real &g_22 = g_(I22,i);
      const Real &g_23 = 0.0;
      const Real &g_30 = 0.0;
      const Real &g_31 = 0.0;
      const Real &g_32 = 0.0;
      const Real &g_33 = g_(I33,i);
      const Real &g01 = 0.0;
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
      Real ux_l = mx_1*u1_l;
      Real uy_l = my_1*u1_l + my_2*u2_l;
      Real uz_l = mz_3*u3_l;
      Real ut_r = mt_0*u0_r;
      Real ux_r = mx_1*u1_r;
      Real uy_r = my_1*u1_r + my_2*u2_r;
      Real uz_r = mz_3*u3_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_1*b1_l;
      Real by_l = my_1*b1_l + my_2*b2_l;
      Real bz_l = mz_3*b3_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_1*b1_r;
      Real by_r = my_1*b1_r + my_2*b2_r;
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
//   il,iu: x-index bounds
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

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract transformation coefficients
    const Real mt_0 = 1.0;
    const Real mx_2 = 1.0 / trans_face2_i1_(i);
    const Real my_3 = 1.0;
    const Real &mz_1 = trans_face2_i1_(i);
    const Real mz_2 = -trans_face2_i2_(i) / trans_face2_i1_(i);

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
    Real uuz_l = mz_1*uu1_l + mz_2*uu2_l;
    Real uux_r = mx_2*uu2_r;
    Real uuy_r = my_3*uu3_r;
    Real uuz_r = mz_1*uu1_r + mz_2*uu2_r;

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
      const Real &g_01 = 0.0;
      const Real &g_02 = 0.0;
      const Real &g_03 = 0.0;
      const Real &g_10 = 0.0;
      const Real &g_11 = g_(I11,i);
      const Real &g_12 = g_(I12,i);
      const Real &g_13 = 0.0;
      const Real &g_20 = 0.0;
      const Real &g_21 = g_(I12,i);
      const Real &g_22 = g_(I22,i);
      const Real &g_23 = 0.0;
      const Real &g_30 = 0.0;
      const Real &g_31 = 0.0;
      const Real &g_32 = 0.0;
      const Real &g_33 = g_(I33,i);
      const Real &g01 = 0.0;
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
      Real uz_l = mz_1*u1_l + mz_2*u2_l;
      Real ut_r = mt_0*u0_r;
      Real ux_r = mx_2*u2_r;
      Real uy_r = my_3*u3_r;
      Real uz_r = mz_1*u1_r + mz_2*u2_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_2*b2_l;
      Real by_l = my_3*b3_l;
      Real bz_l = mz_1*b1_l + mz_2*b2_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_2*b2_r;
      Real by_r = my_3*b3_r;
      Real bz_r = mz_1*b1_r + mz_2*b2_r;

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
//   il,iu: x-index bounds
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

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract transformation coefficients
    const Real mt_0 = 1.0;
    const Real mx_3 = 1.0;
    const Real my_1 = 1.0;
    const Real mz_1 = -trans_face3_i2_(i);
    const Real mz_2 = 1.0;

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
    Real uuy_l = my_1*uu1_l;
    Real uuz_l = mz_1*uu1_l + mz_2*uu2_l;
    Real uux_r = mx_3*uu3_r;
    Real uuy_r = my_1*uu1_r;
    Real uuz_r = mz_1*uu1_r + mz_2*uu2_r;

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
      const Real &g_01 = 0.0;
      const Real &g_02 = 0.0;
      const Real &g_03 = 0.0;
      const Real &g_10 = 0.0;
      const Real &g_11 = g_(I11,i);
      const Real &g_12 = g_(I12,i);
      const Real &g_13 = 0.0;
      const Real &g_20 = 0.0;
      const Real &g_21 = g_(I12,i);
      const Real &g_22 = g_(I22,i);
      const Real &g_23 = 0.0;
      const Real &g_30 = 0.0;
      const Real &g_31 = 0.0;
      const Real &g_32 = 0.0;
      const Real &g_33 = g_(I33,i);
      const Real &g01 = 0.0;
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
      Real uy_l = my_1*u1_l;
      Real uz_l = mz_1*u1_l + mz_2*u2_l;
      Real ut_r = mt_0*u0_r;
      Real ux_r = mx_3*u3_r;
      Real uy_r = my_1*u1_r;
      Real uz_r = mz_1*u1_r + mz_2*u2_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_3*b3_l;
      Real by_l = my_1*b1_l;
      Real bz_l = mz_1*b1_l + mz_2*b2_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_3*b3_r;
      Real by_r = my_1*b1_r;
      Real bz_r = mz_1*b1_r + mz_2*b2_r;

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
//   expects values and x-fluxes of Mx/My/Mz in IM1/IM2/IM3 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts x1-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts x1-fluxes of B2/B3 in IBY/IBZ slots
void Coordinates::FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux)
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
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
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
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
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
  Real g_00 = -1.0;
  const Real &g_11 = metric_cell_i1_(i);
  Real g_12 = -metric_cell_i2_(i);
  const Real &g_21 = g_12;
  Real g_22 = 1.0;
  Real g_33 = 1.0;
  *pa_0 = g_00*a0;
  *pa_1 = g_11*a1 + g_12*a2;
  *pa_2 = g_21*a1 + g_22*a2;
  *pa_3 = g_33*a3;
  return;
}
