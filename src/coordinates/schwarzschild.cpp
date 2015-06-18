// Schwarzschild spacetime, spherical coordinates
// Notes:
//   coordinates: t, r, theta, phi
//   parameters: M (mass)
//   metric:
//     ds^2 = -\alpha^2 dt^2 + 1/\alpha^2 * dr^2
//            + r^2 (d\theta^2 + \sin^2\theta d\phi^2)
//     \alpha = \sqrt(1 - 2M/r)

// Primary header
#include "coordinates.hpp"

// C++ headers
#include <cmath>  // acos(), cos(), log(), sin(), sqrt()

// Athena headers
#include "../athena.hpp"           // enums, macros, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../fluid/eos/eos.hpp"    // FluidEqnOfState
#include "../fluid/fluid.hpp"      // Fluid
#include "../mesh.hpp"             // MeshBlock
#include "../parameter_input.hpp"  // ParameterInput

// Constructor
// Inputs:
//   pb: pointer to MeshBlock containing this grid
//   pin: pointer to runtime inputs
Coordinates::Coordinates(MeshBlock *pb, ParameterInput *pin)
{
  // Set pointer to host MeshBlock
  pmy_block = pb;

  // Set parameters
  bh_mass_ = pin->GetReal("coord", "m");
  bh_spin_ = 0.0;
  const Real &m = bh_mass_;

  // Initialize volume-averaged positions and spacings: r-direction
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST; ++i)
  {
    Real r_m = pb->x1f(i);
    Real r_p = pb->x1f(i+1);
    pb->x1v(i) = std::pow(0.5 * (r_m*r_m*r_m + r_p*r_p*r_p), 1.0/3.0);
  }
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST-1; ++i)
    pb->dx1v(i) = pb->x1v(i+1) - pb->x1v(i);

  // Initialize volume-averaged positions and spacings: theta-direction
  if (pb->block_size.nx2 == 1)  // no extent
  {
    Real theta_m = pb->x2f(pb->js);
    Real theta_p = pb->x2f(pb->js+1);
    pb->x2v(pb->js) = std::acos(0.5 * (std::cos(theta_m) + std::cos(theta_p)));
    pb->dx2v(pb->js) = pb->dx2f(pb->js);
  }
  else  // extended
  {
    for (int j = pb->js-NGHOST; j <= pb->je+NGHOST; j++)
    {
      Real theta_m = pb->x2f(j);
      Real theta_p = pb->x2f(j+1);
      pb->x2v(j) = std::acos(0.5 * (std::cos(theta_m) + std::cos(theta_p)));
    }
    for (int j = pb->js-NGHOST; j <= pb->je+NGHOST-1; j++)
      pb->dx2v(j) = pb->x2v(j+1) - pb->x2v(j);
  }

  // Initialize volume-averaged positions and spacings: phi-direction
  if (pb->block_size.nx3 == 1)  // no extent
  {
    Real phi_m = pb->x3f(pb->ks);
    Real phi_p = pb->x3f(pb->ks+1);
    pb->x3v(pb->ks) = 0.5 * (phi_m + phi_p);
    pb->dx3v(pb->ks) = pb->dx3f(pb->ks);
  }
  else  // extended
  {
    for (int k = pb->ks-NGHOST; k <= pb->ke+NGHOST; k++)
    {
      Real phi_m = pb->x3f(k);
      Real phi_p = pb->x3f(k+1);
      pb->x3v(k) = 0.5 * (phi_m + phi_p);
    }
    for (int k = pb->ks-NGHOST; k <= pb->ke+NGHOST-1; k++)
      pb->dx3v(k) = pb->x3v(k+1) - pb->x3v(k);
  }

  if(pb->pmy_mesh->multilevel==true) { // calc coarse coodinates
    int cis = pb->cis; int cjs = pb->cjs; int cks = pb->cks;
    int cie = pb->cie; int cje = pb->cje; int cke = pb->cke;
    for (int i=cis-(pb->cnghost); i<=cie+(pb->cnghost); ++i) {
      Real r_m = pb->coarse_x1f(i);
      Real r_p = pb->coarse_x1f(i+1);
      pb->coarse_x1v(i) = std::pow(0.5 * (r_m*r_m*r_m + r_p*r_p*r_p), 1.0/3.0);
    }
    for (int i=cis-(pb->cnghost); i<=cie+(pb->cnghost)-1; ++i) {
      pb->coarse_dx1v(i) = pb->coarse_x1v(i+1) - pb->coarse_x1v(i);
    }
    if (pb->block_size.nx2 == 1) {
      pb->coarse_x2v(cjs) = 0.5*(pb->coarse_x2f(cjs+1) + pb->coarse_x2f(cjs));
      pb->coarse_dx2v(cjs) = pb->coarse_dx2f(cjs);
    } else {
      for (int j=cjs-(pb->cnghost); j<=cje+(pb->cnghost); ++j) {
        Real theta_m = pb->coarse_x2f(j);
        Real theta_p = pb->coarse_x2f(j+1);
        pb->coarse_x2v(j) = std::acos(0.5 * (std::cos(theta_m) + std::cos(theta_p)));
      }
      for (int j=cjs-(pb->cnghost); j<=cje+(pb->cnghost)-1; ++j) {
        pb->coarse_dx2v(j) = pb->coarse_x2v(j+1) - pb->coarse_x2v(j);
      }
    }
    if (pb->block_size.nx3 == 1) {
      pb->coarse_x3v(cks) = 0.5*(pb->coarse_x3f(cks+1) + pb->coarse_x3f(cks));
      pb->coarse_dx3v(cks) = pb->coarse_dx3f(cks);
    } else {
      for (int k=cks-(pb->cnghost); k<=cke+(pb->cnghost); ++k) {
        pb->coarse_x3v(k) = 0.5*(pb->coarse_x3f(k+1) + pb->coarse_x3f(k));
      }
      for (int k=cks-(pb->cnghost); k<=cke+(pb->cnghost)-1; ++k) {
        pb->coarse_dx3v(k) = pb->coarse_x3v(k+1) - pb->coarse_x3v(k);
      }
    }
  }

  // Allocate arrays for intermediate geometric quantities: r-direction
  int n_cells_1 = pb->block_size.nx1 + 2*NGHOST;
  coord_vol_i1_.NewAthenaArray(n_cells_1);
  coord_area1_i1_.NewAthenaArray(n_cells_1);
  coord_area2_i1_.NewAthenaArray(n_cells_1);
  coord_area3_i1_.NewAthenaArray(n_cells_1);
  coord_len1_i1_.NewAthenaArray(n_cells_1);
  coord_len2_i1_.NewAthenaArray(n_cells_1);
  coord_len3_i1_.NewAthenaArray(n_cells_1);
  coord_width1_i1_.NewAthenaArray(n_cells_1);
  coord_src_i1_.NewAthenaArray(n_cells_1);
  coord_src_i2_.NewAthenaArray(n_cells_1);
  coord_src_i3_.NewAthenaArray(n_cells_1);
  coord_src_i4_.NewAthenaArray(n_cells_1);
  metric_cell_i1_.NewAthenaArray(n_cells_1);
  metric_face1_i1_.NewAthenaArray(n_cells_1);
  metric_face2_i1_.NewAthenaArray(n_cells_1);
  metric_face3_i1_.NewAthenaArray(n_cells_1);
  trans_face1_i1_.NewAthenaArray(n_cells_1);
  trans_face2_i1_.NewAthenaArray(n_cells_1);
  trans_face3_i1_.NewAthenaArray(n_cells_1);

  // Allocate arrays for intermediate geometric quantities: theta-direction
  int n_cells_2 = (pb->block_size.nx2 > 1) ? pb->block_size.nx2 + 2*NGHOST : 1;
  coord_vol_j1_.NewAthenaArray(n_cells_2);
  coord_area1_j1_.NewAthenaArray(n_cells_2);
  coord_area2_j1_.NewAthenaArray(n_cells_2);
  coord_area3_j1_.NewAthenaArray(n_cells_2);
  coord_len1_j1_.NewAthenaArray(n_cells_2);
  coord_len2_j1_.NewAthenaArray(n_cells_2);
  coord_len3_j1_.NewAthenaArray(n_cells_2);
  coord_width3_j1_.NewAthenaArray(n_cells_2);
  coord_src_j1_.NewAthenaArray(n_cells_2);
  coord_src_j2_.NewAthenaArray(n_cells_2);
  coord_src_j3_.NewAthenaArray(n_cells_2);
  metric_cell_j1_.NewAthenaArray(n_cells_2);
  metric_face1_j1_.NewAthenaArray(n_cells_2);
  metric_face2_j1_.NewAthenaArray(n_cells_2);
  metric_face3_j1_.NewAthenaArray(n_cells_2);
  trans_face1_j1_.NewAthenaArray(n_cells_2);
  trans_face2_j1_.NewAthenaArray(n_cells_2);
  trans_face3_j1_.NewAthenaArray(n_cells_2);

  // Calculate intermediate geometric quantities: r-direction
  #pragma simd
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST; ++i)
  {
    // Useful quantities
    Real r_c = pb->x1v(i);
    Real r_m = pb->x1f(i);
    Real r_p = pb->x1f(i+1);
    Real alpha_c = std::sqrt(1.0 - 2.0*m/r_c);
    Real alpha_m = std::sqrt(1.0 - 2.0*m/r_m);
    Real alpha_p = std::sqrt(1.0 - 2.0*m/r_p);
    Real r_p_cu = r_p*r_p*r_p;
    Real r_m_cu = r_m*r_m*r_m;

    // Volumes, areas, lengths, and widths
    coord_vol_i1_(i) = 1.0/3.0 * (r_p_cu - r_m_cu);
    coord_area1_i1_(i) = SQR(r_m);
    coord_area2_i1_(i) = coord_vol_i1_(i);
    coord_area3_i1_(i) = coord_vol_i1_(i);
    coord_len1_i1_(i) = coord_vol_i1_(i);
    coord_len2_i1_(i) = coord_area1_i1_(i);
    coord_len3_i1_(i) = coord_area1_i1_(i);
    coord_width1_i1_(i) = r_p*alpha_p - r_m*alpha_m
        + m * std::log((r_p*(1.0+alpha_p)-m) / (r_m*(1.0+alpha_m)-m));

    // Source terms
    coord_src_i1_(i) = 3.0*m / (r_p_cu - r_m_cu)
        * (r_p - r_m + 2.0*m * std::log((r_p-2.0*m) / (r_m-2.0*m)));
    coord_src_i2_(i) = 3.0*m / (r_p_cu - r_m_cu)
        * (r_p - r_m - 2.0*m * std::log(r_p/r_m));
    coord_src_i3_(i) = 2.0*m - 3.0/4.0 * (r_m + r_p) * (SQR(r_m) + SQR(r_p))
        / (SQR(r_m) + r_m * r_p + SQR(r_p));
    coord_src_i4_(i) = 3.0/2.0 * (r_m + r_p) / (SQR(r_m) + r_m * r_p + SQR(r_p));

    // Metric coefficients
    metric_cell_i1_(i) = SQR(alpha_c);
    metric_face1_i1_(i) = SQR(alpha_m);
    metric_face2_i1_(i) = SQR(alpha_c);
    metric_face3_i1_(i) = SQR(alpha_c);

    // Coordinate transformations
    trans_face1_i1_(i) = alpha_m;
    trans_face2_i1_(i) = alpha_c;
    trans_face3_i1_(i) = alpha_c;
  }

  // Calculate intermediate geometric quantities: theta-direction
  if (n_cells_2 > 1)  // extended
  {
#pragma simd
    for (int j = pb->js-NGHOST; j <= pb->je+NGHOST; j++)
    {
      // Useful quantities
      Real theta_c = pb->x2v(j);
      Real theta_m = pb->x2f(j);
      Real theta_p = pb->x2f(j+1);
      Real sin_c = std::sin(theta_c);
      Real sin_m = std::sin(theta_m);
      Real sin_p = std::sin(theta_p);
      Real cos_m = std::cos(theta_m);
      Real cos_p = std::cos(theta_p);
      Real sin_c_sq = SQR(sin_c);
      Real sin_m_sq = SQR(sin_m);
      Real sin_m_cu = sin_m_sq*sin_m;
      Real sin_p_cu = SQR(sin_p)*sin_p;

      // Volumes, areas, lengths, and widths
      coord_vol_j1_(j) = cos_m - cos_p;
      coord_area1_j1_(j) = coord_vol_j1_(j);
      coord_area2_j1_(j) = sin_m;
      coord_area3_j1_(j) = coord_vol_j1_(j);
      coord_len1_j1_(j) = coord_area2_j1_(j);
      coord_len2_j1_(j) = coord_vol_j1_(j);
      coord_len3_j1_(j) = coord_area2_j1_(j);
      coord_width3_j1_(j) = sin_c;

      // Source terms
      coord_src_j1_(j) = 1.0/6.0
          * (4.0 - std::cos(2.0*theta_m)
          - 2.0 * cos_m * cos_p - std::cos(2.0*theta_p));
      coord_src_j2_(j) = 1.0/3.0 * (sin_m_cu - sin_p_cu)
          / (cos_m - cos_p);
      coord_src_j3_(j) = (sin_p-sin_m) / (cos_m-cos_p);    // cot((theta_p+theta_m)/2)

      // Metric coefficients
      metric_cell_j1_(j) = sin_c_sq;
      metric_face1_j1_(j) = sin_c_sq;
      metric_face2_j1_(j) = sin_m_sq;
      metric_face3_j1_(j) = sin_c_sq;

      // Coordinate transformations
      trans_face1_j1_(j) = sin_c;
      trans_face2_j1_(j) = sin_m;
      trans_face3_j1_(j) = sin_c;
    }
  }
  else  // no extent
  {
    // Useful quantities
    Real theta_c = pb->x2v(pb->js);
    Real theta_m = pb->x2f(pb->js);
    Real theta_p = pb->x2f(pb->js+1);
    Real sin_c = std::sin(theta_c);
    Real sin_m = std::sin(theta_m);
    Real sin_p = std::sin(theta_p);
    Real cos_m = std::cos(theta_m);
    Real cos_p = std::cos(theta_p);
    Real sin_c_sq = SQR(sin_c);
    Real sin_m_sq = SQR(sin_m);
    Real sin_m_cu = sin_m_sq*sin_m;
    Real sin_p_cu = SQR(sin_p)*sin_p;

    // Volumes and areas
    coord_vol_j1_(pb->js) = cos_m - cos_p;
    coord_area1_j1_(pb->js) = coord_vol_j1_(pb->js);
    coord_area2_j1_(pb->js) = sin_m;
    coord_area3_j1_(pb->js) = coord_vol_j1_(pb->js);
    coord_len1_j1_(pb->js) = coord_area2_j1_(pb->js);
    coord_len2_j1_(pb->js) = coord_vol_j1_(pb->js);
    coord_len3_j1_(pb->js) = coord_area2_j1_(pb->js);

    // Source terms
    coord_src_j1_(pb->js) = 1.0/6.0
        * (4.0 - std::cos(2.0*theta_m)
        - 2.0 * cos_m * cos_p - std::cos(2.0*theta_p));
    coord_src_j2_(pb->js) = 1.0/3.0 * (sin_m_cu - sin_p_cu)
        / (cos_m - cos_p);
    coord_src_j3_(pb->js) = (sin_p-sin_m) / (cos_m-cos_p);   // cot((theta_p+theta_m)/2)

    // Metric coefficients
    metric_cell_j1_(pb->js) = sin_c_sq;
    metric_face1_j1_(pb->js) = sin_c_sq;
    metric_face2_j1_(pb->js) = sin_m_sq;
    metric_face3_j1_(pb->js) = sin_c_sq;

    // Coordinate transformations
    trans_face1_j1_(pb->js) = sin_c;
    trans_face2_j1_(pb->js) = sin_m;
    trans_face3_j1_(pb->js) = sin_c;
  }
}

// Destructor
Coordinates::~Coordinates()
{
  coord_vol_i1_.DeleteAthenaArray();
  coord_area1_i1_.DeleteAthenaArray();
  coord_area2_i1_.DeleteAthenaArray();
  coord_area3_i1_.DeleteAthenaArray();
  coord_len1_i1_.DeleteAthenaArray();
  coord_len2_i1_.DeleteAthenaArray();
  coord_len3_i1_.DeleteAthenaArray();
  coord_width1_i1_.DeleteAthenaArray();
  coord_src_i1_.DeleteAthenaArray();
  coord_src_i2_.DeleteAthenaArray();
  coord_src_i3_.DeleteAthenaArray();
  coord_src_i4_.DeleteAthenaArray();
  coord_vol_j1_.DeleteAthenaArray();
  coord_area1_j1_.DeleteAthenaArray();
  coord_area2_j1_.DeleteAthenaArray();
  coord_area3_j1_.DeleteAthenaArray();
  coord_len1_j1_.DeleteAthenaArray();
  coord_len2_j1_.DeleteAthenaArray();
  coord_len3_j1_.DeleteAthenaArray();
  coord_width3_j1_.DeleteAthenaArray();
  coord_src_j1_.DeleteAthenaArray();
  coord_src_j2_.DeleteAthenaArray();
  coord_src_j3_.DeleteAthenaArray();
  metric_cell_i1_.DeleteAthenaArray();
  metric_cell_j1_.DeleteAthenaArray();
  metric_face1_i1_.DeleteAthenaArray();
  metric_face1_j1_.DeleteAthenaArray();
  metric_face2_i1_.DeleteAthenaArray();
  metric_face2_j1_.DeleteAthenaArray();
  metric_face3_i1_.DeleteAthenaArray();
  metric_face3_j1_.DeleteAthenaArray();
  trans_face1_i1_.DeleteAthenaArray();
  trans_face1_j1_.DeleteAthenaArray();
  trans_face2_i1_.DeleteAthenaArray();
  trans_face2_j1_.DeleteAthenaArray();
  trans_face3_i1_.DeleteAthenaArray();
  trans_face3_j1_.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------

// Function for computing cell volumes
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   volumes: 1D array of cell volumes
// Notes:
//   \Delta V = 1/3 * \Delta(r^3) (-\Delta\cos\theta) \Delta\phi
void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &volumes)
{
  const Real &neg_delta_cos_theta = coord_vol_j1_(j);
  const Real &delta_phi = pmy_block->dx3f(k);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &third_delta_r_cb = coord_vol_i1_(i);
    Real &volume = volumes(i);
    volume = third_delta_r_cb * neg_delta_cos_theta * delta_phi;
  }
  return;
}

Real Coordinates::GetCellVolume(const int k, const int j, const int i)
{
  return coord_vol_i_(i)*coord_vol_j_(j)*pmy_block->dx3f(k);
}

//--------------------------------------------------------------------------------------

// Function for computing areas orthogonal to r
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to r
// Notes:
//   \Delta A = r^2 (-\Delta\cos\theta) \Delta\phi
void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  const Real &neg_delta_cos_theta = coord_area1_j1_(j);
  const Real &delta_phi = pmy_block->dx3f(k);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &r_sq = coord_area1_i1_(i);
    Real &area = areas(i);
    area = r_sq * neg_delta_cos_theta * delta_phi;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing areas orthogonal to theta
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to theta
// Notes:
//   \Delta A = 1/3 \Delta(r^3) \sin\theta \Delta\phi
void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  const Real &sin_theta = coord_area2_j1_(j);
  const Real &delta_phi = pmy_block->dx3f(k);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &third_delta_r_cb = coord_area2_i1_(i);
    Real &area = areas(i);
    area = third_delta_r_cb * sin_theta * delta_phi;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing areas orthogonal to phi
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to phi
// Notes:
//   \Delta A = 1/3 \Delta(r^3) (-\Delta\cos\theta)
void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  const Real &neg_delta_cos_theta = coord_area3_j1_(j);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &third_delta_r_cb = coord_area3_i1_(i);
    Real &area = areas(i);
    area = third_delta_r_cb * neg_delta_cos_theta;
  }
  return;
}


Real Coordinates::GetFace1Area(const int k, const int j, const int i)
{
  return coord_area1_i_(i)*coord_area1_j_(j)*pmy_block->dx3f(k);
}


//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the r-direction
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   lengths: 1D array of edge lengths along r
// Notes:
//   \Delta L = 1/3 \Delta(r^3) \sin\theta
void Coordinates::Edge1Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &lengths)
{
  const Real &sin_theta = coord_len1_j1_(j);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &third_delta_r_cb = coord_len1_i1_(i);
    Real &length = lengths(i);
    length = third_delta_r_cb * sin_theta;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the theta-direction
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   lengths: 1D array of edge lengths along r
// Notes:
//   \Delta L = r^2 (-\Delta\cos\theta)
void Coordinates::Edge2Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &lengths)
{
  const Real &neg_delta_cos_theta = coord_len2_j1_(j);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &r_sq = coord_len2_i1_(i);
    Real &length = lengths(i);
    length = r_sq * neg_delta_cos_theta;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the phi-direction
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   lengths: 1D array of edge lengths along r
// Notes:
//   \Delta L = r^2 \sin\theta \Delta\phi
void Coordinates::Edge3Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &lengths)
{
  const Real &sin_theta = coord_len3_j1_(j);
  const Real &delta_phi = pmy_block->dx3f(k);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    const Real &r_sq = coord_len3_i1_(i);
    Real &length = lengths(i);
    length = r_sq * sin_theta * delta_phi;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the r-direction
// Inputs:
//   k,j: phi- and theta-indices (unused)
//   i: r-index
// Outputs:
//   returned value: r-width of cell (i,j,k)
// Notes:
//   \Delta W = \Delta(r \alpha) + M \Delta\log(r(1+\alpha)-M)
Real Coordinates::CenterWidth1(const int k, const int j, const int i)
{
  return coord_width1_i1_(i);
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the theta-direction
// Inputs:
//   k: phi-index (unused)
//   j,i: theta- and r-indices
// Outputs:
//   returned value: theta-width of cell (i,j,k)
// Notes:
//   \Delta W = r \Delta\theta
Real Coordinates::CenterWidth2(const int k, const int j, const int i)
{
  const Real &r = pmy_block->x1v(i);
  const Real &delta_theta = pmy_block->dx1f(j);
  return r * delta_theta;
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the phi-direction
// Inputs:
//   k,j,i: phi-, theta-, and r-indices
// Outputs:
//   returned value: phi-width of cell (i,j,k)
// Notes:
//   \Delta W = r \sin\theta \Delta\phi
Real Coordinates::CenterWidth3(const int k, const int j, const int i)
{
  const Real &r = pmy_block->x1v(i);
  const Real &sin_theta = coord_width3_j1_(j);
  const Real &delta_phi = pmy_block->dx1f(k);
  return r * sin_theta * delta_phi;
}

//--------------------------------------------------------------------------------------

// Function for computing source terms using r-fluxes
// Inputs:
//   k,j: phi- and theta-indices
//   dt: size of timestep
//   flux: 1D array of r-fluxes
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

  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_cell_j1_(j);
  const Real &gamma2_33 = coord_src_j2_(j);
  const Real &gamma3_23 = coord_src_j3_(j);
  const Real &gamma3_32 = gamma3_23;

  // Go through cells
  #pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; ++i)
  {
    // Extract remaining geometric quantities
    const Real &alpha_sq = metric_cell_i1_(i);
    const Real &r = pmy_block->x1v(i);
    const Real r_sq = SQR(r);
    const Real g00 = -alpha_sq;
    const Real &g11 = 1.0/alpha_sq;
    const Real &g22 = r_sq;
    const Real g33 = r_sq * sin_sq_theta;
    const Real &gamma0_01 = coord_src_i1_(i);
    const Real &gamma0_10 = gamma0_01;
    const Real &gamma1_00 = coord_src_i2_(i);
    const Real gamma1_11 = -gamma0_01;
    const Real &gamma1_22 = coord_src_i3_(i);
    const Real gamma1_33 = gamma1_22 * coord_src_j1_(j);
    const Real &gamma2_12 = coord_src_i4_(i);
    const Real &gamma2_21 = gamma2_12;
    const Real &gamma3_13 = gamma2_12;
    const Real &gamma3_31 = gamma3_13;

    // Extract primitives
    const Real &rho = prim(IDN,k,j,i);
    const Real &pgas = prim(IEN,k,j,i);
    const Real &v1 = prim(IVX,k,j,i);
    const Real &v2 = prim(IVY,k,j,i);
    const Real &v3 = prim(IVZ,k,j,i);

    // Calculate 4-velocity
    Real u0 = std::sqrt(-1.0 / (g00 + g11*v1*v1 + g22*v2*v2 + g33*v3*v3));
    Real u1 = u0 * v1;
    Real u2 = u0 * v2;
    Real u3 = u0 * v3;
    Real u_0 = g00 * u0;
    Real u_1 = g11 * u1;
    Real u_2 = g22 * u2;
    Real u_3 = g33 * u3;

    // Extract and calculate magnetic field
    Real bcon0 = 0.0, bcon1 = 0.0, bcon2 = 0.0, bcon3 = 0.0;
    Real bcov0 = 0.0, bcov1 = 0.0, bcov2 = 0.0, bcov3 = 0.0;
    Real b_sq = 0.0;
    if (MAGNETIC_FIELDS_ENABLED)
    {
      const Real &b1 = bcc(IB1,k,j,i);
      const Real &b2 = bcc(IB2,k,j,i);
      const Real &b3 = bcc(IB3,k,j,i);
      bcon0 = g11*b1*u1 + g22*b2*u2 + g33*b3*u3;
      bcon1 = (b1 + bcon0 * u1) / u0;
      bcon2 = (b2 + bcon0 * u2) / u0;
      bcon3 = (b3 + bcon0 * u3) / u0;
      bcov0 = g00*bcon0;
      bcov1 = g11*bcon1;
      bcov2 = g22*bcon2;
      bcov3 = g33*bcon3;
      b_sq = bcov0*bcon0 + bcov1*bcon1 + bcov2*bcon2 + bcov3*bcon3;
    }

    // Calculate stress-energy tensor
    Real wtot = rho + gamma_adi_red * pgas + b_sq;
    Real ptot = pgas + 0.5 * b_sq;
    Real t0_0 = wtot*u0*u_0 - bcon0*bcov0 + ptot;
    Real t0_1 = wtot*u0*u_1 - bcon0*bcov1;
    Real t1_0 = wtot*u1*u_0 - bcon1*bcov0;
    Real t1_1 = wtot*u1*u_1 - bcon1*bcov1 + ptot;
    Real t1_2 = wtot*u1*u_2 - bcon1*bcov2;
    Real t1_3 = wtot*u1*u_3 - bcon1*bcov3;
    Real t2_1 = wtot*u2*u_1 - bcon2*bcov1;
    Real t2_2 = wtot*u2*u_2 - bcon2*bcov2 + ptot;
    Real t2_3 = wtot*u2*u_3 - bcon2*bcov3;
    Real t3_1 = wtot*u3*u_1 - bcon3*bcov1;
    Real t3_2 = wtot*u3*u_2 - bcon3*bcov2;
    Real t3_3 = wtot*u3*u_3 - bcon3*bcov3 + ptot;

    // Calculate source terms
    Real s0 = gamma0_10 * t1_0 + gamma1_00 * t0_1;
    Real s1 = gamma0_01 * t0_0 + gamma1_11 * t1_1 + gamma2_21 * t2_2 + gamma3_31 * t3_3;
    Real s2 = gamma1_22 * t2_1 + gamma2_12 * t1_2 + gamma3_32 * t3_3;
    Real s3 = gamma1_33 * t3_1 + gamma2_33 * t3_2 + gamma3_13 * t1_3 + gamma3_23 * t2_3;

    // Extract conserved quantities
    Real &e = cons(IEN,k,j,i);
    Real &m1 = cons(IM1,k,j,i);
    Real &m2 = cons(IM2,k,j,i);
    Real &m3 = cons(IM3,k,j,i);

    // Add source terms to conserved quantities
    e += dt * s0;
    m1 += dt * s1;
    m2 += dt * s2;
    m3 += dt * s3;
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing source terms using theta-fluxes
// Inputs:
//   k,j: phi- and theta-indices
//   dt: size of timestep
//   flux_j: 1D array of theta-fluxes left of cells j
//   flux_jp1: 1D array of theta-fluxes right of cells j
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

// Function for computing source terms using phi-fluxes
// Inputs:
//   k,j: phi- and theta-indices
//   dt: size of timestep
//   flux_k: 2D array of phi-fluxes left of cells k
//   flux_kp1: 2D array of phi-fluxes right of cells k
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
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::CellMetric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_cell_j1_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract remaining geometric quantities
    const Real &alpha_sq = metric_cell_i1_(i);
    const Real &r = pmy_block->x1v(i);
    Real r_sq = SQR(r);

    // Extract metric terms
    Real &g00 = g(I00,i);
    Real &g11 = g(I11,i);
    Real &g22 = g(I22,i);
    Real &g33 = g(I33,i);
    Real &gi00 = g_inv(I00,i);
    Real &gi11 = g_inv(I11,i);
    Real &gi22 = g_inv(I22,i);
    Real &gi33 = g_inv(I33,i);

    // Set metric terms
    g00 = -alpha_sq;
    g11 = 1.0/alpha_sq;
    g22 = r_sq;
    g33 = r_sq * sin_sq_theta;
    gi00 = -1.0/alpha_sq;
    gi11 = alpha_sq;
    gi22 = 1.0/r_sq;
    gi33 = 1.0 / (r_sq * sin_sq_theta);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: r-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face1Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face1_j1_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract remaining geometric quantities
    const Real &alpha_sq = metric_face1_i1_(i);
    const Real &r = pmy_block->x1f(i);
    Real r_sq = SQR(r);

    // Extract metric terms
    Real &g00 = g(I00,i);
    Real &g11 = g(I11,i);
    Real &g22 = g(I22,i);
    Real &g33 = g(I33,i);
    Real &gi00 = g_inv(I00,i);
    Real &gi11 = g_inv(I11,i);
    Real &gi22 = g_inv(I22,i);
    Real &gi33 = g_inv(I33,i);

    // Set metric terms
    g00 = -alpha_sq;
    g11 = 1.0/alpha_sq;
    g22 = r_sq;
    g33 = r_sq * sin_sq_theta;
    gi00 = -1.0/alpha_sq;
    gi11 = alpha_sq;
    gi22 = 1.0/r_sq;
    gi33 = 1.0 / (r_sq * sin_sq_theta);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: theta-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face2Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face2_j1_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract remaining geometric quantities
    const Real &alpha_sq = metric_face2_i1_(i);
    const Real &r = pmy_block->x1v(i);
    Real r_sq = SQR(r);

    // Extract metric terms
    Real &g00 = g(I00,i);
    Real &g11 = g(I11,i);
    Real &g22 = g(I22,i);
    Real &g33 = g(I33,i);
    Real &gi00 = g_inv(I00,i);
    Real &gi11 = g_inv(I11,i);
    Real &gi22 = g_inv(I22,i);
    Real &gi33 = g_inv(I33,i);

    // Set metric terms
    g00 = -alpha_sq;
    g11 = 1.0/alpha_sq;
    g22 = r_sq;
    g33 = r_sq * sin_sq_theta;
    gi00 = -1.0/alpha_sq;
    gi11 = alpha_sq;
    gi22 = 1.0/r_sq;
    gi33 = 1.0 / (r_sq * sin_sq_theta);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: phi-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face3Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face3_j1_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract remaining geometric quantities
    const Real &alpha_sq = metric_face3_i1_(i);
    const Real &r = pmy_block->x1v(i);
    Real r_sq = SQR(r);

    // Extract metric terms
    Real &g00 = g(I00,i);
    Real &g11 = g(I11,i);
    Real &g22 = g(I22,i);
    Real &g33 = g(I33,i);
    Real &gi00 = g_inv(I00,i);
    Real &gi11 = g_inv(I11,i);
    Real &gi22 = g_inv(I22,i);
    Real &gi33 = g_inv(I33,i);

    // Set metric terms
    g00 = -alpha_sq;
    g11 = 1.0/alpha_sq;
    g22 = r_sq;
    g33 = r_sq * sin_sq_theta;
    gi00 = -1.0/alpha_sq;
    gi11 = alpha_sq;
    gi22 = 1.0/r_sq;
    gi33 = 1.0 / (r_sq * sin_sq_theta);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: r-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
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
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face1_j1_(j);
  const Real &sin_theta = trans_face1_j1_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract remaining geometric quantities
    const Real &alpha_sq = metric_face1_i1_(i);
    const Real &r = pmy_block->x1f(i);
    const Real r_sq = SQR(r);
    const Real &alpha = trans_face1_i1_(i);
    const Real g00 = -alpha_sq;
    const Real g11 = 1.0/alpha_sq;
    const Real g22 = r_sq;
    const Real g33 = r_sq * sin_sq_theta;
    const Real mt_0 = alpha;
    const Real mx_1 = 1.0/alpha;
    const Real my_2 = r;
    const Real mz_3 = r * sin_theta;

    // Extract global 3-velocities
    Real &v1l = prim_left(IVX,i);
    Real &v2l = prim_left(IVY,i);
    Real &v3l = prim_left(IVZ,i);
    Real &v1r = prim_right(IVX,i);
    Real &v2r = prim_right(IVY,i);
    Real &v3r = prim_right(IVZ,i);

    // Construct global 4-velocities
    Real u0l = std::sqrt(-1.0 / (g00 + g11*v1l*v1l + g22*v2l*v2l + g33*v3l*v3l));
    Real u1l = u0l * v1l;
    Real u2l = u0l * v2l;
    Real u3l = u0l * v3l;
    Real u0r = std::sqrt(-1.0 / (g00 + g11*v1r*v1r + g22*v2r*v2r + g33*v3r*v3r));
    Real u1r = u0r * v1r;
    Real u2r = u0r * v2r;
    Real u3r = u0r * v3r;

    // Transform 4-velocities
    Real utl = mt_0*u0l;
    Real uxl = mx_1*u1l;
    Real uyl = my_2*u2l;
    Real uzl = mz_3*u3l;
    Real utr = mt_0*u0r;
    Real uxr = mx_1*u1r;
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
      const Real &b1 = b1_vals(k,j,i);
      Real &b2l = prim_left(IBY,i);
      Real &b3l = prim_left(IBZ,i);
      Real &b2r = prim_right(IBY,i);
      Real &b3r = prim_right(IBZ,i);

      // Construct global contravariant magnetic fields
      Real bcon0l = g11*b1*u1l + g22*b2l*u2l + g33*b3l*u3l;
      Real bcon1l = (b1 + bcon0l * u1l) / u0l;
      Real bcon2l = (b2l + bcon0l * u2l) / u0l;
      Real bcon3l = (b3l + bcon0l * u3l) / u0l;
      Real bcon0r = g11*b1*u1r + g22*b2r*u2r + g33*b3r*u3r;
      Real bcon1r = (b1 + bcon0r * u1r) / u0r;
      Real bcon2r = (b2r + bcon0r * u2r) / u0r;
      Real bcon3r = (b3r + bcon0r * u3r) / u0r;

      // Transform contravariant magnetic fields
      Real bcontl = mt_0*bcon0l;
      Real bconxl = mx_1*bcon1l;
      Real bconyl = my_2*bcon2l;
      Real bconzl = mz_3*bcon3l;
      Real bcontr = mt_0*bcon0r;
      Real bconxr = mx_1*bcon1r;
      Real bconyr = my_2*bcon2r;
      Real bconzr = mz_3*bcon3r;

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

// Function for transforming primitives to locally flat frame: theta-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
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
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face2_j1_(j);
  const Real &sin_theta = trans_face2_j1_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract geometric quantities
    const Real &alpha_sq = metric_face2_i1_(i);
    const Real &r = pmy_block->x1v(i);
    const Real r_sq = SQR(r);
    const Real &alpha = trans_face2_i1_(i);
    const Real g00 = -alpha_sq;
    const Real g11 = 1.0/alpha_sq;
    const Real g22 = r_sq;
    const Real g33 = r_sq * sin_sq_theta;
    const Real mt_0 = alpha;
    const Real mx_2 = 1.0/r;
    const Real my_3 = r * sin_theta;
    const Real &mz_1 = 1.0/alpha;

    // Extract global 3-velocities
    Real &v1l = prim_left(IVX,i);
    Real &v2l = prim_left(IVY,i);
    Real &v3l = prim_left(IVZ,i);
    Real &v1r = prim_right(IVX,i);
    Real &v2r = prim_right(IVY,i);
    Real &v3r = prim_right(IVZ,i);

    // Construct global 4-velocities
    Real u0l = std::sqrt(-1.0 / (g00 + g11*v1l*v1l + g22*v2l*v2l + g33*v3l*v3l));
    Real u1l = u0l * v1l;
    Real u2l = u0l * v2l;
    Real u3l = u0l * v3l;
    Real u0r = std::sqrt(-1.0 / (g00 + g11*v1r*v1r + g22*v2r*v2r + g33*v3r*v3r));
    Real u1r = u0r * v1r;
    Real u2r = u0r * v2r;
    Real u3r = u0r * v3r;

    // Transform 4-velocities
    Real utl = mt_0*u0l;
    Real uxl = mx_2*u2l;
    Real uyl = my_3*u3l;
    Real uzl = mz_1*u1l;
    Real utr = mt_0*u0r;
    Real uxr = mx_2*u2r;
    Real uyr = my_3*u3r;
    Real uzr = mz_1*u1r;

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
      Real bcon0l = g11*b1l*u1l + g22*b2*u2l + g33*b3l*u3l;
      Real bcon1l = (b1l + bcon0l * u1l) / u0l;
      Real bcon2l = (b2 + bcon0l * u2l) / u0l;
      Real bcon3l = (b3l + bcon0l * u3l) / u0l;
      Real bcon0r = g11*b1r*u1r + g22*b2*u2r + g33*b3r*u3r;
      Real bcon1r = (b1r + bcon0r * u1r) / u0r;
      Real bcon2r = (b2 + bcon0r * u2r) / u0r;
      Real bcon3r = (b3r + bcon0r * u3r) / u0r;

      // Transform contravariant magnetic fields
      Real bcontl = mt_0*bcon0l;
      Real bconxl = mx_2*bcon2l;
      Real bconyl = my_3*bcon3l;
      Real bconzl = mz_1*bcon1l;
      Real bcontr = mt_0*bcon0r;
      Real bconxr = mx_2*bcon2r;
      Real bconyr = my_3*bcon3r;
      Real bconzr = mz_1*bcon1r;

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

// Function for transforming primitives to locally flat frame: phi-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
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
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face3_j1_(j);
  const Real &sin_theta = trans_face3_j1_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract geometric quantities
    const Real &alpha_sq = metric_face3_i1_(i);
    const Real &r = pmy_block->x1v(i);
    const Real r_sq = SQR(r);
    const Real &alpha = trans_face3_i1_(i);
    const Real g00 = -alpha_sq;
    const Real g11 = 1.0/alpha_sq;
    const Real g22 = r_sq;
    const Real g33 = r_sq * sin_sq_theta;
    const Real mt_0 = alpha;
    const Real mx_3 = r * sin_theta;
    const Real my_1 = 1.0/alpha;
    const Real mz_2 = r;

    // Extract global 3-velocities
    Real &v1l = prim_left(IVX,i);
    Real &v2l = prim_left(IVY,i);
    Real &v3l = prim_left(IVZ,i);
    Real &v1r = prim_right(IVX,i);
    Real &v2r = prim_right(IVY,i);
    Real &v3r = prim_right(IVZ,i);

    // Construct global 4-velocities
    Real u0l = std::sqrt(-1.0 / (g00 + g11*v1l*v1l + g22*v2l*v2l + g33*v3l*v3l));
    Real u1l = u0l * v1l;
    Real u2l = u0l * v2l;
    Real u3l = u0l * v3l;
    Real u0r = std::sqrt(-1.0 / (g00 + g11*v1r*v1r + g22*v2r*v2r + g33*v3r*v3r));
    Real u1r = u0r * v1r;
    Real u2r = u0r * v2r;
    Real u3r = u0r * v3r;

    // Transform 4-velocities
    Real utl = mt_0*u0l;
    Real uxl = mx_3*u3l;
    Real uyl = my_1*u1l;
    Real uzl = mz_2*u2l;
    Real utr = mt_0*u0r;
    Real uxr = mx_3*u3r;
    Real uyr = my_1*u1r;
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
      const Real &b3 = b3_vals(k,j,i);
      Real &b1l = prim_left(IBY,i);
      Real &b2l = prim_left(IBZ,i);
      Real &b1r = prim_right(IBY,i);
      Real &b2r = prim_right(IBZ,i);

      // Construct global contravariant magnetic fields
      Real bcon0l = g11*b1l*u1l + g22*b2l*u2l + g33*b3*u3l;
      Real bcon1l = (b1l + bcon0l * u1l) / u0l;
      Real bcon2l = (b2l + bcon0l * u2l) / u0l;
      Real bcon3l = (b3 + bcon0l * u3l) / u0l;
      Real bcon0r = g11*b1r*u1r + g22*b2r*u2r + g33*b3*u3r;
      Real bcon1r = (b1r + bcon0r * u1r) / u0r;
      Real bcon2r = (b2r + bcon0r * u2r) / u0r;
      Real bcon3r = (b3 + bcon0r * u3r) / u0r;

      // Transform contravariant magnetic fields
      Real bcontl = mt_0*bcon0l;
      Real bconxl = mx_3*bcon3l;
      Real bconyl = my_1*bcon1l;
      Real bconzl = mz_2*bcon2l;
      Real bcontr = mt_0*bcon0r;
      Real bconxr = mx_3*bcon3r;
      Real bconyr = my_1*bcon1r;
      Real bconzr = mz_2*bcon2r;

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

// Function for transforming fluxes to global frame: r-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates (unused)
//   bx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM1/IM2/IM3 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts r-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts r-fluxes of B2/B3 in IBY/IBZ slots
void Coordinates::FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx, AthenaArray<Real> &flux)
{
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face1_j1_(j);
  const Real &sin_theta = trans_face1_j1_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract geometric quantities
    const Real &alpha_sq = metric_face1_i1_(i);
    const Real &r = pmy_block->x1f(i);
    const Real r_sq = SQR(r);
    const Real &alpha = trans_face1_i1_(i);
    const Real g00 = -alpha_sq;
    const Real g11 = 1.0/alpha_sq;
    const Real g22 = r_sq;
    const Real g33 = r_sq * sin_sq_theta;
    const Real m0_t = 1.0/alpha;
    const Real m1_x = alpha;
    const Real m2_y = 1.0/r;
    const Real m3_z = 1.0 / (r * sin_theta);

    // Extract local conserved quantities and fluxes
    const Real dx = flux(IDN,i);
    const Real txt = flux(IEN,i);
    const Real txx = flux(IM1,i);
    const Real txy = flux(IM2,i);
    const Real txz = flux(IM3,i);

    // Transform stress-energy tensor
    Real t10 = m1_x*m0_t*txt;
    Real t11 = m1_x*m1_x*txx;
    Real t12 = m1_x*m2_y*txy;
    Real t13 = m1_x*m3_z*txz;

    // Extract global fluxes
    Real &d1 = flux(IDN,i);
    Real &t1_0 = flux(IEN,i);
    Real &t1_1 = flux(IM1,i);
    Real &t1_2 = flux(IM2,i);
    Real &t1_3 = flux(IM3,i);

    // Set fluxes
    d1 = m1_x*dx;
    t1_0 = g00*t10;
    t1_1 = g11*t11;
    t1_2 = g22*t12;
    t1_3 = g33*t13;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      const Real fyx = flux(IBY,i);
      const Real fzx = flux(IBZ,i);
      Real &f21 = flux(IBY,i);
      Real &f31 = flux(IBZ,i);
      f21 = m2_y*m1_x*fyx;
      f31 = m3_z*m1_x*fzx;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: theta-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates (unused)
//   bx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM2/IM3/IM1 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts theta-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts theta-fluxes of B3/B1 in IBY/IBZ slots
void Coordinates::FluxToGlobal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx, AthenaArray<Real> &flux)
{
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face2_j1_(j);
  const Real &sin_theta = trans_face2_j1_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract geometric quantities
    const Real &alpha_sq = metric_face2_i1_(i);
    const Real &r = pmy_block->x1v(i);
    const Real r_sq = SQR(r);
    const Real &alpha = trans_face2_i1_(i);
    const Real g00 = -alpha_sq;
    const Real g11 = 1.0/alpha_sq;
    const Real g22 = r_sq;
    const Real g33 = r_sq * sin_sq_theta;
    const Real m0_t = 1.0/alpha;
    const Real m1_z = alpha;
    const Real m2_x = 1.0/r;
    const Real m3_y = 1.0 / (r * sin_theta);

    // Extract local conserved quantities and fluxes
    const Real dx = flux(IDN,i);
    const Real txt = flux(IEN,i);
    const Real txx = flux(IM2,i);
    const Real txy = flux(IM3,i);
    const Real txz = flux(IM1,i);

    // Transform stress-energy tensor
    Real t20 = m2_x*m0_t*txt;
    Real t21 = m2_x*m1_z*txz;
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
    t2_0 = g00*t20;
    t2_1 = g11*t21;
    t2_2 = g22*t22;
    t2_3 = g33*t23;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      const Real fyx = flux(IBY,i);
      const Real fzx = flux(IBZ,i);
      Real &f32 = flux(IBY,i);
      Real &f12 = flux(IBZ,i);
      f32 = m3_y*m2_x*fyx;
      f12 = m1_z*m2_x*fzx;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: phi-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates (unused)
//   bx: 1D array of longitudinal magnetic fields, in local coordinates (unused)
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM3/IM1/IM2 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts phi-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts phi-fluxes of B1/B2 in IBY/IBZ slots
void Coordinates::FluxToGlobal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bx, AthenaArray<Real> &flux)
{
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face3_j1_(j);
  const Real &sin_theta = trans_face3_j1_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract geometric quantities
    const Real &alpha_sq = metric_face3_i1_(i);
    const Real &r = pmy_block->x1v(i);
    const Real r_sq = SQR(r);
    const Real &alpha = trans_face3_i1_(i);
    const Real g00 = -alpha_sq;
    const Real g11 = 1.0/alpha_sq;
    const Real g22 = r_sq;
    const Real g33 = r_sq * sin_sq_theta;
    const Real m0_t = 1.0/alpha;
    const Real m1_y = alpha;
    const Real m2_z = 1.0/r;
    const Real m3_x = 1.0 / (r * sin_theta);

    // Extract local conserved quantities and fluxes
    const Real dx = flux(IDN,i);
    const Real txt = flux(IEN,i);
    const Real txx = flux(IM3,i);
    const Real txy = flux(IM1,i);
    const Real txz = flux(IM2,i);

    // Transform stress-energy tensor
    Real t30 = m3_x*m0_t*txt;
    Real t31 = m3_x*m1_y*txy;
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
    t3_0 = g00*t30;
    t3_1 = g11*t31;
    t3_2 = g22*t32;
    t3_3 = g33*t33;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      const Real fyx = flux(IBY,i);
      const Real fzx = flux(IBZ,i);
      Real &f13 = flux(IBY,i);
      Real &f23 = flux(IBZ,i);
      f13 = m1_y*m3_x*fyx;
      f23 = m2_z*m3_x*fzx;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Boyer-Lindquist to Schwarzschild
// Inputs:
//   a0_bl,a1_bl,a2_bl,a3_bl: upper 4-vector components in Boyer-Lindquist coordinates
//   k,j,i: indices of cell in which transformation is desired
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in Schwarzschild coordinates
// Notes:
//   Schwarzschild coordinates match Boyer-Lindquist when a = 0
void Coordinates::TransformVectorCell(
    Real a0_bl, Real a1_bl, Real a2_bl, Real a3_bl, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = a0_bl;
  *pa1 = a1_bl;
  *pa2 = a2_bl;
  *pa3 = a3_bl;
  return;
}

// Function for transforming 4-vector from Boyer-Lindquist to Schwarzschild
// Inputs:
//   a0_bl,a1_bl,a2_bl,a3_bl: upper 4-vector components in Boyer-Lindquist coordinates
//   k,j,i: indices of x1-face in which transformation is desired
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in Schwarzschild coordinates
// Notes:
//   Schwarzschild coordinates match Boyer-Lindquist when a = 0
void Coordinates::TransformVectorFace1(
    Real a0_bl, Real a1_bl, Real a2_bl, Real a3_bl, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = a0_bl;
  *pa1 = a1_bl;
  *pa2 = a2_bl;
  *pa3 = a3_bl;
  return;
}

// Function for transforming 4-vector from Boyer-Lindquist to Schwarzschild
// Inputs:
//   a0_bl,a1_bl,a2_bl,a3_bl: upper 4-vector components in Boyer-Lindquist coordinates
//   k,j,i: indices of x2-face in which transformation is desired
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in Schwarzschild coordinates
// Notes:
//   Schwarzschild coordinates match Boyer-Lindquist when a = 0
void Coordinates::TransformVectorFace2(
    Real a0_bl, Real a1_bl, Real a2_bl, Real a3_bl, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = a0_bl;
  *pa1 = a1_bl;
  *pa2 = a2_bl;
  *pa3 = a3_bl;
  return;
}

// Function for transforming 4-vector from Boyer-Lindquist to Schwarzschild
// Inputs:
//   a0_bl,a1_bl,a2_bl,a3_bl: upper 4-vector components in Boyer-Lindquist coordinates
//   k,j,i: indices of x3-face in which transformation is desired
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in Schwarzschild coordinates
// Notes:
//   Schwarzschild coordinates match Boyer-Lindquist when a = 0
void Coordinates::TransformVectorFace3(
    Real a0_bl, Real a1_bl, Real a2_bl, Real a3_bl, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  *pa0 = a0_bl;
  *pa1 = a1_bl;
  *pa2 = a2_bl;
  *pa3 = a3_bl;
  return;
}

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
  const Real &sin_sq_theta = metric_cell_j1_(j);
  const Real &alpha_sq = metric_cell_i1_(i);
  const Real &r = pmy_block->x1v(i);
  Real r_sq = SQR(r);

  // Calculate metric terms
  Real g_00 = -alpha_sq;
  Real g_11 = 1.0/alpha_sq;
  Real g_22 = r_sq;
  Real g_33 = r_sq * sin_sq_theta;

  // Transform vector
  *pa_0 = g_00 * a0;
  *pa_1 = g_11 * a1;
  *pa_2 = g_22 * a2;
  *pa_3 = g_33 * a3;
  return;
}

// Function for returning Boyer-Lindquist coordinates of given cell
// Inputs:
//   x1,x2,x3: Schwarzschild coordinates to be converted
// Outputs:
//   pr: pointer to stored value of r
//   ptheta: pointer to stored value of theta
//   pphi: pointer to stored value of phi
// Notes:
//   Schwarzschild (x1,x2,x3) match Boyer-Lindquist (r,theta,phi) when a = 0
void Coordinates::GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3,
    Real *pr, Real *ptheta, Real *pphi)
{
  *pr = x1;
  *ptheta = x2;
  *pphi = x3;
  return;
}
