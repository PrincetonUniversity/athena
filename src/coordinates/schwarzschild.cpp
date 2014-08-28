// Schwarzschild spacetime, spherical coordinates
// Notes:
//   coordinates: t, r, theta, phi
//   parameters: M (mass)
//   metric: ds^2 = -(1-2M/r) dt^2 + 1/(1-2M/r) dr^2
//                  + r^2 (d\theta^2 + \sin^2\theta d\phi^2)

// Primary header
#include "coordinates.hpp"

// C++ headers
#include <cmath>  // acos(), cbrt(), cos(), cot, sin(), sqrt()

// Athena headers
#include "../athena.hpp"         // enums, macros, Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../mesh.hpp"           // MeshBlock

// TODO: find better input method
namespace globals
{
  const Real MASS = 1.0;
};
using namespace globals;

// Constructor
// Inputs:
//   pb: pointer to MeshBlock containing this grid
Coordinates::Coordinates(MeshBlock *pb)
{
  // Set pointer to host MeshBlock
  pmy_block = pb;

  // Initialize volume-averated positions and spacings: r-direction
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST; i++)
  {
    Real r_m = pb->x1f(i);
    Real r_p = pb->x1f(i+1);
    pb->x1v(i) = std::cbrt(0.5 * (r_m*r_m*r_m + r_p*r_p*r_p));
  }
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST-1; i++)
    pb->dx1v(i) = pb->x1v(i+1) - pb->x1v(i);

  // Initialize volume-averated positions and spacings: theta-direction
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

  // Initialize volume-averated positions and spacings: phi-direction
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

  // Allocate scratch arrays for integrator
  int n_cells_1 = pb->block_size.nx1 + 2*NGHOST;
  face_area.NewAthenaArray(n_cells_1);
  cell_volume.NewAthenaArray(n_cells_1);

  // Allocate arrays for intermediate geometric quantities: r-direction
  volume_i_.NewAthenaArray(n_cells_1);
  face1_area_i_.NewAthenaArray(n_cells_1);
  face2_area_i_.NewAthenaArray(n_cells_1);
  face3_area_i_.NewAthenaArray(n_cells_1);
  src_terms_i1_.NewAthenaArray(n_cells_1);
  src_terms_i2_.NewAthenaArray(n_cells_1);
  src_terms_i3_.NewAthenaArray(n_cells_1);
  src_terms_i4_.NewAthenaArray(n_cells_1);
  metric_cell_i1_.NewAthenaArray(n_cells_1);
  metric_cell_i2_.NewAthenaArray(n_cells_1);
  metric_cell_i3_.NewAthenaArray(n_cells_1);
  metric_cell_i4_.NewAthenaArray(n_cells_1);
  metric_cell_i5_.NewAthenaArray(n_cells_1);
  metric_cell_i6_.NewAthenaArray(n_cells_1);
  metric_face1_i1_.NewAthenaArray(n_cells_1);
  metric_face1_i2_.NewAthenaArray(n_cells_1);
  metric_face1_i3_.NewAthenaArray(n_cells_1);
  metric_face2_i1_.NewAthenaArray(n_cells_1);
  metric_face2_i2_.NewAthenaArray(n_cells_1);
  metric_face2_i3_.NewAthenaArray(n_cells_1);
  metric_face3_i1_.NewAthenaArray(n_cells_1);
  metric_face3_i2_.NewAthenaArray(n_cells_1);
  metric_face3_i3_.NewAthenaArray(n_cells_1);
  trans_face1_i1_.NewAthenaArray(n_cells_1);
  trans_face1_i2_.NewAthenaArray(n_cells_1);
  trans_face1_i3_.NewAthenaArray(n_cells_1);
  trans_face1_i4_.NewAthenaArray(n_cells_1);
  trans_face2_i1_.NewAthenaArray(n_cells_1);
  trans_face2_i2_.NewAthenaArray(n_cells_1);
  trans_face2_i3_.NewAthenaArray(n_cells_1);
  trans_face2_i4_.NewAthenaArray(n_cells_1);
  trans_face3_i1_.NewAthenaArray(n_cells_1);
  trans_face3_i2_.NewAthenaArray(n_cells_1);
  trans_face3_i3_.NewAthenaArray(n_cells_1);
  trans_face3_i4_.NewAthenaArray(n_cells_1);

  // Allocate arrays for intermediate geometric quantities: theta-direction
  int n_cells_2 = (pb->block_size.nx2 > 1) ? pb->block_size.nx2 + 2*NGHOST : 1;
  volume_j_.NewAthenaArray(n_cells_2);
  face1_area_j_.NewAthenaArray(n_cells_2);
  face2_area_j_.NewAthenaArray(n_cells_2);
  face3_area_j_.NewAthenaArray(n_cells_2);
  src_terms_j1_.NewAthenaArray(n_cells_2);
  src_terms_j2_.NewAthenaArray(n_cells_2);
  src_terms_j3_.NewAthenaArray(n_cells_2);
  metric_cell_j1_.NewAthenaArray(n_cells_2);
  metric_cell_j2_.NewAthenaArray(n_cells_2);
  metric_face1_j1_.NewAthenaArray(n_cells_2);
  metric_face2_j1_.NewAthenaArray(n_cells_2);
  metric_face3_j1_.NewAthenaArray(n_cells_2);
  trans_face1_j1_.NewAthenaArray(n_cells_2);
  trans_face1_j2_.NewAthenaArray(n_cells_2);
  trans_face2_j1_.NewAthenaArray(n_cells_2);
  trans_face2_j2_.NewAthenaArray(n_cells_2);
  trans_face3_j1_.NewAthenaArray(n_cells_2);
  trans_face3_j2_.NewAthenaArray(n_cells_2);

  // Calculate intermediate geometric quantities: r-direction
#pragma simd
  for (int i = is-NGHOST; i <= ie+NGHOST; i++)
  {
    // Useful quantities
    Real r_c = pb->x1v(i);
    Real r_m = pb->x1f(i);
    Real r_p = pb->x1f(i+1);

    // Volumes and areas
    volume_i_(i) = 1.0/3.0 * (r_p*r_p*r_p - r_m*r_m*r_m);
    face1_area_i_(i) = r_m*r_m;
    face2_area_i_(i) = volume_i_(i);
    face3_area_i_(i) = volume_i_(i);

    // Source terms
    src_terms_i1_(i) = 3.0*MASS / (r_p*r_p*r_p - r_m*r_m*r_m)
        * (r_p - r_m - 2.0*MASS * std::log((2.0*MASS - r_m)/(2.0*MASS - r_p)));
    src_terms_i2_(i) = 3.0*MASS / (r_p*r_p*r_p - r_m*r_m*r_m)
        * (r_p - r_m + 2.0*MASS * std::log(r_m/r_p));
    src_terms_i3_(i) = 2.0*MASS - 3.0/4.0 * (r_m + r_p) * (r_m*r_m + r_p*r_p)
        / (r_m*r_m + r_m * r_p + r_p*r_p);
    src_terms_i4_(i) = 3.0/2.0 * (r_m + r_p) / (r_m*r_m + r_m * r_p + r_p*r_p);

    // Cell-centered metric
    metric_cell_i1_(i) = -(1.0 - 2.0*MASS / r_c);
    metric_cell_i2_(i) = 1.0 / (1.0 - 2.0*MASS / r_c);
    metric_cell_i3_(i) = r_c*r_c;
    metric_cell_i4_(i) = 1.0 / metric_cell_i1_(i);
    metric_cell_i5_(i) = 1.0 / metric_cell_i2_(i);
    metric_cell_i6_(i) = 1.0 / metric_cell_i3_(i);

    // Face-centered metric
    metric_face1_i1_(i) = -(1.0 - 2.0*MASS / r_m);
    metric_face1_i2_(i) = 1.0 / (1.0 - 2.0*MASS / r_m);
    metric_face1_i3_(i) = r_m*r_m;
    metric_face2_i1_(i) = metric_cell_i1_(i);
    metric_face2_i2_(i) = metric_cell_i2_(i);
    metric_face2_i3_(i) = metric_cell_i3_(i);
    metric_face3_i1_(i) = metric_cell_i1_(i);
    metric_face3_i2_(i) = metric_cell_i2_(i);
    metric_face3_i3_(i) = metric_cell_i3_(i);

    // Coordinate transformations
    trans_face1_i1_(i) = std::sqrt(1.0 - 2.0*MASS/r_m);
    trans_face1_i2_(i) = 1.0 / trans_face1_i1_(i);
    trans_face1_i3_(i) = r_m;
    trans_face1_i4_(i) = 1.0 / r_m;
    trans_face2_i1_(i) = std::sqrt(1.0 - 2.0*MASS/r_c);
    trans_face2_i2_(i) = 1.0 / trans_face2_i1_(i);
    trans_face2_i3_(i) = r_c;
    trans_face2_i4_(i) = 1.0 / r_c;
    trans_face3_i1_(i) = trans_face2_i1_(i);
    trans_face3_i2_(i) = trans_face2_i2_(i);
    trans_face3_i3_(i) = trans_face2_i3_(i);
    trans_face3_i4_(i) = trans_face2_i4_(i);
  }

  // Calculate intermediate geometric quantities: theta-direction
  if (n_cells_2 > 1)  // extended
  {
#pragma simd
    for (int j = js-NGHOST; j <= js+NGHOST; j++)
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

      // Volumes and areas
      volume_j_(j) = cos_m - cos_p;
      face1_area_j_(j) = volume_j_(j);
      face2_area_j_(j) = sin_m;
      face3_area_j_(j) = volume_j_(j);

      // Source terms
      src_terms_j1_(j) = 1.0/6.0 * (4.0 - std::cos(2.0*theta_m)
          - 2.0 * cos_m * cos_p - std::cos(2.0*theta_p));
      src_terms_j2_(j) = 1.0/3.0 * (sin_m*sin_m*sin_m - sin_p*sin_p*sin_p)
          / (cos_m - cos_p);
      src_terms_j3_(j) = std::cot(0.5 * (theta_m + theta_p));

      // Cell-centered metric
      metric_cell_j1_(j) = sin_c*sin_c;
      metric_cell_j2_(j) = 1.0 / metric_cell_j1_(j);

      // Face-centered metric
      metric_face1_j1_(j) = metric_cell_j1_(j);
      metric_face2_j1_(j) = sin_m*sin_m;
      metric_face3_j1_(j) = metric_cell_j1_(j);

      // Coordinate transformations
      trans_face1_j1_(j) = sin_c;
      trans_face1_j2_(j) = 1.0 / sin_c;
      trans_face2_j1_(j) = sin_m;
      trans_face2_j2_(j) = 1.0 / sin_m;
      trans_face3_j1_(j) = sin_m;
      trans_face3_j2_(j) = trans_face2_j2_(j);
    }
  }
  else  // no extent
  {
    // Useful quantities
    Real theta_c = pb->x2v(js);
    Real theta_m = pb->x2f(js);
    Real theta_p = pb->x2f(js+1);
    Real sin_c = std::sin(theta_c);
    Real sin_m = std::sin(theta_m);
    Real sin_p = std::sin(theta_p);
    Real cos_m = std::cos(theta_m);
    Real cos_p = std::cos(theta_p);

    // Volumes and areas
    volume_j_(js) = cos_m - cos_p;
    face1_area_j_(js) = volume_j_(js);
    face2_area_j_(js) = sin_m;
    face3_area_j_(js) = volume_j_(js);

    // Source terms
    src_terms_j1_(js) = 1.0/6.0 * (4.0 - std::cos(2.0*theta_m)
        - 2.0 * cos_m * cos_p - std::cos(2.0*theta_p));
    src_terms_j2_(js) = 1.0/3.0 * (sin_m*sin_m*sin_m - sin_p*sin_p*sin_p)
        / (cos_m - cos_p);
    src_terms_j3_(js) = std::cot(0.5 * (theta_m + theta_p));

    // Cell-centered metric
    metric_cell_j1_(js) = sin_c*sin_c;
    metric_cell_j2_(js) = 1.0 / metric_cell_j1_(js);

    // Face-centered metric
    metric_face1_j1_(js) = metric_cell_j1_(js);
    metric_face2_j1_(js) = sin_m*sin_m;
    metric_face3_j1_(js) = metric_cell_j1_(js);

    // Coordinate transformations
    trans_face1_j1_(js) = sin_c;
    trans_face1_j2_(js) = 1.0 / sin_c;
    trans_face2_j1_(js) = sin_m;
    trans_face2_j2_(js) = 1.0 / sin_m;
    trans_face3_j1_(js) = sin_m;
    trans_face3_j2_(js) = trans_face2_j2_(js);
  }
}

// Destructor
Coordinates::~Coordinates()
{
  face_area.DeleteAthenaArray();
  cell_volume.DeleteAthenaArray();
  volume_i_.DeleteAthenaArray();
  face1_area_i_.DeleteAthenaArray();
  face2_area_i_.DeleteAthenaArray();
  face3_area_i_.DeleteAthenaArray();
  src_terms_i1_.DeleteAthenaArray();
  src_terms_i2_.DeleteAthenaArray();
  src_terms_i3_.DeleteAthenaArray();
  src_terms_i4_.DeleteAthenaArray();
  volume_j_.DeleteAthenaArray();
  face1_area_j_.DeleteAthenaArray();
  face2_area_j_.DeleteAthenaArray();
  face3_area_j_.DeleteAthenaArray();
  src_terms_j1_.DeleteAthenaArray();
  src_terms_j2_.DeleteAthenaArray();
  src_terms_j3_.DeleteAthenaArray();
  metric_cell_i1_.DeleteAthenaArray();
  metric_cell_i2_.DeleteAthenaArray();
  metric_cell_i3_.DeleteAthenaArray();
  metric_cell_i4_.DeleteAthenaArray();
  metric_cell_i5_.DeleteAthenaArray();
  metric_cell_i6_.DeleteAthenaArray();
  metric_cell_j1_.DeleteAthenaArray();
  metric_cell_j2_.DeleteAthenaArray();
  metric_face1_i1_.DeleteAthenaArray();
  metric_face1_i2_.DeleteAthenaArray();
  metric_face1_i3_.DeleteAthenaArray();
  metric_face1_j1_.DeleteAthenaArray();
  metric_face2_i1_.DeleteAthenaArray();
  metric_face2_i2_.DeleteAthenaArray();
  metric_face2_i3_.DeleteAthenaArray();
  metric_face2_j1_.DeleteAthenaArray();
  metric_face3_i1_.DeleteAthenaArray();
  metric_face3_i2_.DeleteAthenaArray();
  metric_face3_i3_.DeleteAthenaArray();
  metric_face3_j1_.DeleteAthenaArray();
  trans_face1_i1_.DeleteAthenaArray();
  trans_face1_i2_.DeleteAthenaArray();
  trans_face1_i3_.DeleteAthenaArray();
  trans_face1_i4_.DeleteAthenaArray();
  trans_face1_j1_.DeleteAthenaArray();
  trans_face1_j2_.DeleteAthenaArray();
  trans_face2_i1_.DeleteAthenaArray();
  trans_face2_i2_.DeleteAthenaArray();
  trans_face2_i3_.DeleteAthenaArray();
  trans_face2_i4_.DeleteAthenaArray();
  trans_face2_j1_.DeleteAthenaArray();
  trans_face2_j2_.DeleteAthenaArray();
  trans_face3_i1_.DeleteAthenaArray();
  trans_face3_i2_.DeleteAthenaArray();
  trans_face3_i3_.DeleteAthenaArray();
  trans_face3_i4_.DeleteAthenaArray();
  trans_face3_j1_.DeleteAthenaArray();
  trans_face3_j2_.DeleteAthenaArray();
}

// Function for computing areas orthogonal to r
// Inputs:
//   k: phi-index
//   j: theta-index
//   il, iu: r-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to r
// Notes:
//   \Delta A = r^2 (-\Delta\cos\theta) \Delta\phi
void Coordinates::Area1Face(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  Real &neg_delta_cos_theta = face1_area_j_(j);
  Real &delta_phi = pmy_block->dx3f(k);
#pragma simd
  for (int i = il; i <= iu; i++)
  {
    Real &area = areas(i);
    Real &r_sq = face1_area_i_(i);
    area = r_sq * neg_delta_cos_theta * delta_phi;
  }
  return;
}

// Function for computing areas orthogonal to theta
// Inputs:
//   k: phi-index
//   j: theta-index
//   il, iu: r-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to theta
// Notes:
//   \Delta A = 1/3 \Delta(r^3) \sin\theta \Delta\phi
void Coordinates::Area2Face(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  Real &sin_theta = face2_area_j_(j);
  Real &delta_phi = pmy_block->dx3f(k);
#pragma simd
  for (int i = il; i <= iu; i++)
  {
    Real &area = areas(i);
    Real &third_delta_r_cb = face2_area_i_(i);
    area = third_delta_r_cb * sin_theta * delta_phi;
  }
  return;
}

// Function for computing areas orthogonal to phi
// Inputs:
//   k: phi-index
//   j: theta-index
//   il, iu: r-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to phi
// Notes:
//   \Delta A = 1/3 \Delta(r^3) (-\Delta\cos\theta)
void Coordinates::Area3Face(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  Real &neg_delta_cos_theta = face3_area_j_(j);
#pragma simd
  for (int i = il; i <= iu; i++)
  {
    Real &area = areas(i);
    Real &third_delta_r_cb = face3_area_i_(i);
    area = third_delta_r_cb * neg_delta_cos_theta;
  }
  return;
}

// Function for computing cell volumes
// Inputs:
//   k: phi-index
//   j: theta-index
//   il, iu: r-index bounds
// Outputs:
//   volumes: 1D array of cell volumes
// Notes:
//   \Delta A = 1/3 \Delta(r^3) (-\Delta\cos\theta) \Delta\phi
void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &volumes)
{
  Real &neg_delta_cos_theta = volume_j_(j);
  Real &delta_phi = pmy_block->dx3f(k);
#pragma simd
  for (int i = il; i <= iu; i++)
  {
    Real &area = areas(i);
    Real &third_delta_r_cb = volume_i_(i);
    area = third_delta_r_cb * neg_delta_cos_theta * delta_phi;
  }
  return;
}

// Function for computing source terms
// Inputs:
//   k: z-index
//   j: y-index
//   prim: 1D array of primitive values in cells
// Outputs:
//   sources: array of source terms in 1D
// Notes:
//   source terms all vanish identically
//   sources assumed to be 0-initialized
void Coordinates::CoordinateSourceTerms(const int k, const int j,
    AthenaArray<Real> &prim, AthenaArray<Real> &sources)
{
  // Extract ratio of specific heats
  const Real gamma_adi = pmy_fluid->GetGamma();
  const Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Extract geometric quantities that do not depend on r
  Real &gamma_233 = src_terms_j2_(j);
  Real &gamma_323 = src_terms_j3_(j);
  Real &gamma_332 = gamma_323;

  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
  {
    // Extract remaining geometric quantities
    Real &g00 = metric_cell_i1_(i);
    Real &g11 = metric_cell_i2_(i);
    Real &g22 = metric_cell_i3_(i);
    Real g33 = metric_cell_i3_(i) * metric_cell_j1_(j);
    Real &gamma_001 = src_terms_i1_(i);
    Real &gamma_010 = gamma_001;
    Real &gamma_100 = src_terms_i2_(i);
    Real gamma_111 = -gamma_001;
    Real &gamma_122 = src_terms_i3_(i);
    Real gamma_133 = gamma_122 * src_terms_j1_(j);
    Real &gamma_212 = src_terms_i4_(i);
    Real &gamma_221 = gamma_212;
    Real &gamma_313 = gamma_212;
    Real &gamma_331 = gamma_313;

    // Extract primitives
    Real &rho = prim(IDN,k,j,i);
    Real &pgas = prim(IEN,k,j,i);
    Real &v1 = prim(IVX,k,j,i);
    Real &v2 = prim(IVY,k,j,i);
    Real &v3 = prim(IVZ,k,j,i);

    // Extract sources
    Real &s0 = sources(IEN,k,j,i);
    Real &s1 = sources(IM1,k,j,i);
    Real &s2 = sources(IM2,k,j,i);
    Real &s3 = sources(IM3,k,j,i);

    // Calculate 4-velocity
    Real u0 = std::sqrt(-1.0 / (g00 + g11*v1*v1 + g22*v2*v2 + g33*v3*v3));
    Real u1 = u0 * v1;
    Real u2 = u0 * v2;
    Real u3 = u0 * v3;
    Real u_lower_0 = g00 * u0;
    Real u_lower_1 = g11 * u1;
    Real u_lower_2 = g22 * u2;
    Real u_lower_3 = g33 * u3;

    // Calculate stress-energy tensor
    Real rho_h = rho + gamma_adi_red * pgas;
    Real t00 = rho_h * u0 * u_lower_0 + pgas;
    Real t01 = rho_h * u0 * u_lower_1;
    Real t10 = rho_h * u1 * u_lower_0;
    Real t11 = rho_h * u1 * u_lower_1 + pgas;
    Real t12 = rho_h * u1 * u_lower_2;
    Real t13 = rho_h * u1 * u_lower_3;
    Real t21 = rho_h * u2 * u_lower_1;
    Real t22 = rho_h * u2 * u_lower_2 + pgas;
    Real t23 = rho_h * u2 * u_lower_3;
    Real t31 = rho_h * u3 * u_lower_1;
    Real t32 = rho_h * u3 * u_lower_2;
    Real t33 = rho_h * u3 * u_lower_3 + pgas;

    // Set source terms
    s0 = gamma_010 * t10 + gamma_100 * t01;
    s1 = gamma_001 * t00 + gamma_111 * t11 + gamma_221 * t22 + gamma_331 * t33;
    s2 = gamma_122 * t21 + gamma_212 * t12 + gamma_332 * t33;
    s3 = gamma_133 * t31 + gamma_233 * t32 + gamma_313 * t13 + gamma_323 * t23;
  }
  return;
}

// Function for computing metric terms
// Inputs:
//   k: phi-index
//   j: theta-index
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::CellMetric(const int k, const int j, AthenaArray<Real> &g,
    AthenaArray<Real> &g_inv)
{
  // Extract geometric quantities that do not depend on r
  Real &sin_sq_theta = metric_cell_j1_(j);
  Real &csc_sq_theta = metric_cell_j2_(j);

  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is-NGHOST; i <= pmy_block->ie+NGHOST; i++)
  {
    // Extract remaining geometric quantities
    Real &t_factor = metric_cell_i1_(i);
    Real &r_factor = metric_cell_i2_(i);
    Real &r_sq = metric_cell_i3_(i);
    Real &t_inv_factor = metric_cell_i4_(i);
    Real &r_inv_factor = metric_cell_i5_(i);
    Real &r_inv_sq = metric_cell_i6_(i);

    // Extract metric terms
    Real &g00 = g(I00,i);
    Real &g11 = g(I11,i);
    Real &g22 = g(I22,i);
    Real &g33 = g(I33,i);
    Real &g_inv_00 = g_inv(I00,i);
    Real &g_inv_11 = g_inv(I11,i);
    Real &g_inv_22 = g_inv(I22,i);
    Real &g_inv_33 = g_inv(I33,i);

    // Set metric terms
    // TODO: should 0's be set explicitly?
    g00 = t_factor;
    g11 = r_factor;
    g22 = r_sq;
    g33 = r_sq * sin_sq_theta;
    g_inv_00 = t_inv_factor;
    g_inv_11 = r_inv_factor;
    g_inv_22 = r_inv_sq;
    g_inv_33 = r_inv_sq * csc_sq_theta;
  }
  return;
}

// Function for transforming primitives to locally flat frame: r-interface
// Inputs:
//   k: phi-index
//   j: theta-index
//   prim: array of primitives in 1D, using global coordinates
// Outputs:
//   prim: values overwritten in local coordinates
void Coordinates::PrimToLocal1(const int k, const int j, AthenaArray<Real> &prim)
{
  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie+1; i++)
  {
    // Extract primitives
    Real &v1 = prim(IVX,k,j,i);
    Real &v2 = prim(IVY,k,j,i);
    Real &v3 = prim(IVZ,k,j,i);

    // Extract geometric quantities
    Real &g00 = metric_face1_i1_(i);
    Real &g11 = metric_face1_i2_(i);
    Real &g22 = metric_face1_i3_(i);
    Real g33 = metric_face1_i3_(i) * metric_face1_j1_(j);
    Real &mt0 = trans_face1_i1_(i);
    Real &mx1 = trans_face1_i2_(i);
    Real &my2 = trans_face1_i3_(i);
    Real mz3 = trans_face1_i3_(i) * trans_face1_j1_(j);

    // Construct 4-velocity
    Real u0 = std::sqrt(-1.0 / (g00 + g11 * v1*v1 + g22 * v2*v2 + g33 * v3*v3));
    Real u1 = u0 * v1;
    Real u2 = u0 * v2;
    Real u3 = u0 * v3;

    // Transform 4-velocity
    Real u0_new = mt0 * u0;
    Real u1_new = mx1 * u1;
    Real u2_new = my2 * u2;
    Real u3_new = mz3 * u3;
    v1 = u1_new / u0_new;
    v2 = u2_new / u0_new;
    v2 = u3_new / u0_new;
  }
  return;
}

// Function for transforming primitives to locally flat frame: theta-interface
// Inputs:
//   k: phi-index
//   j: theta-index
//   prim: array of primitives in 1D, using global coordinates
// Outputs:
//   prim: values overwritten in local coordinates
void Coordinates::PrimToLocal2(const int k, const int j, AthenaArray<Real> &prim)
{
  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
  {
    // Extract primitives
    Real &v1 = prim(IVX,k,j,i);
    Real &v2 = prim(IVY,k,j,i);
    Real &v3 = prim(IVZ,k,j,i);

    // Extract geometric quantities
    Real &g00 = metric_face2_i1_(i);
    Real &g11 = metric_face2_i2_(i);
    Real &g22 = metric_face2_i3_(i);
    Real g33 = metric_face2_i3_(i) * metric_face2_j1_(j);
    Real &mt0 = trans_face2_i1_(i);
    Real &mx1 = trans_face2_i2_(i);
    Real &my2 = trans_face2_i3_(i);
    Real mz3 = trans_face2_i3_(i) * trans_face2_j1_(j);

    // Construct 4-velocity
    Real u0 = std::sqrt(-1.0 / (g00 + g11 * v1*v1 + g22 * v2*v2 + g33 * v3*v3));
    Real u1 = u0 * v1;
    Real u2 = u0 * v2;
    Real u3 = u0 * v3;

    // Transform 4-velocity
    Real u0_new = mt0 * u0;
    Real u1_new = mx1 * u1;
    Real u2_new = my2 * u2;
    Real u3_new = mz3 * u3;
    v1 = u1_new / u0_new;
    v2 = u2_new / u0_new;
    v2 = u3_new / u0_new;
  }
  return;
}

// Function for transforming primitives to locally flat frame: phi-interface
// Inputs:
//   k: phi-index
//   j: theta-index
//   prim: array of primitives in 1D, using global coordinates
// Outputs:
//   prim: values overwritten in local coordinates
void Coordinates::PrimToLocal3(const int k, const int j, AthenaArray<Real> &prim)
{
  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
  {
    // Extract primitives
    Real &v1 = prim(IVX,k,j,i);
    Real &v2 = prim(IVY,k,j,i);
    Real &v3 = prim(IVZ,k,j,i);

    // Extract geometric quantities
    Real &g00 = metric_face3_i1_(i);
    Real &g11 = metric_face3_i2_(i);
    Real &g22 = metric_face3_i3_(i);
    Real g33 = metric_face3_i3_(i) * metric_face3_j1_(j);
    Real &mt0 = trans_face3_i1_(i);
    Real &mx1 = trans_face3_i2_(i);
    Real &my2 = trans_face3_i3_(i);
    Real mz3 = trans_face3_i3_(i) * trans_face3_j1_(j);

    // Construct 4-velocity
    Real u0 = std::sqrt(-1.0 / (g00 + g11 * v1*v1 + g22 * v2*v2 + g33 * v3*v3));
    Real u1 = u0 * v1;
    Real u2 = u0 * v2;
    Real u3 = u0 * v3;

    // Transform 4-velocity
    Real u0_new = mt0 * u0;
    Real u1_new = mx1 * u1;
    Real u2_new = my2 * u2;
    Real u3_new = mz3 * u3;
    v1 = u1_new / u0_new;
    v2 = u2_new / u0_new;
    v2 = u3_new / u0_new;
  }
  return;
}

// Function for transforming fluxes to global frame: r-interface
// Inputs:
//   k: phi-index
//   j: theta-index
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
void Coordinates::FluxToGlobal1(const int k, const int j, AthenaArray<Real> &flux)
{
  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie+1; i++)
  {
    // Extract fluxes for reading
    Real &dx = flux(IDN,i);
    Real &txt = flux(IEN,i);
    Real &txx = flux(IM1,i);
    Real &txy = flux(IM2,i);
    Real &txz = flux(IM3,i);

    // Extract geometric quantities
    Real &g00 = metric_face1_i1_(i);
    Real &g11 = metric_face1_i2_(i);
    Real &g22 = metric_face1_i3_(i);
    Real g33 = metric_face1_i3_(i) * metric_face1_j1_(j);
    Real &m0t = trans_face1_i2_(i);
    Real &m1x = trans_face1_i1_(i);
    Real &m2y = trans_face1_i4_(i);
    Real m3z = trans_face1_i4_(i) * trans_face1_j2_(j);

    // Calculate new fluxes
    Real f_d = m1x * dx;
    Real f_m0 = m1x * g00 * m0t * txt;
    Real f_m1 = m1x * g11 * m1x * txx;
    Real f_m2 = m1x * g22 * m2y * txy;
    Real f_m3 = m1x * g33 * m3z * txz;

    // Extract fluxes for writing
    Real &d1 = flux(IDN,i);
    Real &t10 = flux(IEN,i);
    Real &t11 = flux(IM1,i);
    Real &t12 = flux(IM2,i);
    Real &t13 = flux(IM3,i);

    // Set fluxes
    d1 = f_d;
    t10 = f_m0;
    t11 = f_m1;
    t12 = f_m2;
    t13 = f_m3;
  }
  return;
}

// Function for transforming fluxes to global frame: theta-interface
// Inputs:
//   k: phi-index
//   j: theta-index
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
void Coordinates::FluxToGlobal2(const int k, const int j, AthenaArray<Real> &flux)
{
  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
  {
    // Extract fluxes for reading
    Real &dy = flux(IDN,i);
    Real &tyt = flux(IEN,i);
    Real &tyx = flux(IM1,i);
    Real &tyy = flux(IM2,i);
    Real &tyz = flux(IM3,i);

    // Extract geometric quantities
    Real &g00 = metric_face2_i1_(i);
    Real &g11 = metric_face2_i2_(i);
    Real &g22 = metric_face2_i3_(i);
    Real g33 = metric_face2_i3_(i) * metric_face2_j1_(j);
    Real &m0t = trans_face2_i2_(i);
    Real &m1x = trans_face2_i1_(i);
    Real &m2y = trans_face2_i4_(i);
    Real m3z = trans_face2_i4_(i) * trans_face2_j2_(j);

    // Calculate new fluxes
    Real f_d = m2y * dy;
    Real f_m0 = m2y * g00 * m0t * tyt;
    Real f_m1 = m2y * g11 * m1x * tyx;
    Real f_m2 = m2y * g22 * m2y * tyy;
    Real f_m3 = m2y * g33 * m3z * tyz;

    // Extract fluxes for writing
    Real &d2 = flux(IDN,i);
    Real &t20 = flux(IEN,i);
    Real &t21 = flux(IM1,i);
    Real &t22 = flux(IM2,i);
    Real &t23 = flux(IM3,i);

    // Set fluxes
    d2 = f_d;
    t20 = f_m0;
    t21 = f_m1;
    t22 = f_m2;
    t23 = f_m3;
  }
  return;
}

// Function for transforming fluxes to global frame: phi-interface
// Inputs:
//   k: phi-index
//   j: theta-index
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
void Coordinates::FluxToGlobal3(const int k, const int j, AthenaArray<Real> &flux)
{
  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
  {
    // Extract fluxes for reading
    Real &dz = flux(IDN,i);
    Real &tzt = flux(IEN,i);
    Real &tzx = flux(IM1,i);
    Real &tzy = flux(IM2,i);
    Real &tzz = flux(IM3,i);

    // Extract geometric quantities
    Real &g00 = metric_face3_i1_(i);
    Real &g11 = metric_face3_i2_(i);
    Real &g22 = metric_face3_i3_(i);
    Real g33 = metric_face3_i3_(i) * metric_face3_j1_(j);
    Real &m0t = trans_face3_i2_(i);
    Real &m1x = trans_face3_i1_(i);
    Real &m2y = trans_face3_i4_(i);
    Real m3z = trans_face3_i4_(i) * trans_face3_j2_(j);

    // Calculate new fluxes
    Real f_d = m3z * dz;
    Real f_m0 = m3z * g00 * m0t * tzt;
    Real f_m1 = m3z * g11 * m1x * tzx;
    Real f_m2 = m3z * g22 * m2y * tzy;
    Real f_m3 = m3z * g33 * m3z * tzz;

    // Extract fluxes for writing
    Real &d3 = flux(IDN,i);
    Real &t30 = flux(IEN,i);
    Real &t31 = flux(IM1,i);
    Real &t32 = flux(IM2,i);
    Real &t33 = flux(IM3,i);

    // Set fluxes
    d3 = f_d;
    t30 = f_m0;
    t31 = f_m1;
    t32 = f_m2;
    t33 = f_m3;
  }
  return;
}
