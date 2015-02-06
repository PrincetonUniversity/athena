// Schwarzschild spacetime, spherical coordinates
// Notes:
//   coordinates: t, r, theta, phi
//   parameters: M (mass)
//   metric: ds^2 = -(1-2M/r) dt^2 + 1/(1-2M/r) dr^2
//                  + r^2 (d\theta^2 + \sin^2\theta d\phi^2)

// Primary header
#include "coordinates.hpp"

// C++ headers
#include <cmath>  // acos(), cbrt(), cos(), sin(), sqrt()

// Athena headers
#include "../athena.hpp"         // enums, macros, Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../fluid/eos/eos.hpp"  // GetGamma()
#include "../fluid/fluid.hpp"    // Fluid
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
//   pin: pointer to runtime inputs
Coordinates::Coordinates(MeshBlock *pb, ParameterInput *pin)
{
  // Set pointer to host MeshBlock
  pmy_block = pb;

  // Initialize volume-averaged positions and spacings: r-direction
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST; i++)
  {
    Real r_m = pb->x1f(i);
    Real r_p = pb->x1f(i+1);
    pb->x1v(i) = cbrt(0.5 * (r_m*r_m*r_m + r_p*r_p*r_p));
  }
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST-1; i++)
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

  // Allocate arrays for intermediate geometric quantities: r-direction
  int n_cells_1 = pb->block_size.nx1 + 2*NGHOST;
  volume_i_.NewAthenaArray(n_cells_1);
  face1_area_i_.NewAthenaArray(n_cells_1);
  face2_area_i_.NewAthenaArray(n_cells_1);
  face3_area_i_.NewAthenaArray(n_cells_1);
  edge1_length_i_.NewAthenaArray(n_cells_1);
  edge2_length_i_.NewAthenaArray(n_cells_1);
  edge3_length_i_.NewAthenaArray(n_cells_1);
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
  metric_face1_i4_.NewAthenaArray(n_cells_1);
  metric_face1_i5_.NewAthenaArray(n_cells_1);
  metric_face1_i6_.NewAthenaArray(n_cells_1);
  metric_face2_i1_.NewAthenaArray(n_cells_1);
  metric_face2_i2_.NewAthenaArray(n_cells_1);
  metric_face2_i3_.NewAthenaArray(n_cells_1);
  metric_face2_i4_.NewAthenaArray(n_cells_1);
  metric_face2_i5_.NewAthenaArray(n_cells_1);
  metric_face2_i6_.NewAthenaArray(n_cells_1);
  metric_face3_i1_.NewAthenaArray(n_cells_1);
  metric_face3_i2_.NewAthenaArray(n_cells_1);
  metric_face3_i3_.NewAthenaArray(n_cells_1);
  metric_face3_i4_.NewAthenaArray(n_cells_1);
  metric_face3_i5_.NewAthenaArray(n_cells_1);
  metric_face3_i6_.NewAthenaArray(n_cells_1);
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
  edge1_length_j_.NewAthenaArray(n_cells_2);
  edge2_length_j_.NewAthenaArray(n_cells_2);
  edge3_length_j_.NewAthenaArray(n_cells_2);
  src_terms_j1_.NewAthenaArray(n_cells_2);
  src_terms_j2_.NewAthenaArray(n_cells_2);
  src_terms_j3_.NewAthenaArray(n_cells_2);
  metric_cell_j1_.NewAthenaArray(n_cells_2);
  metric_cell_j2_.NewAthenaArray(n_cells_2);
  metric_face1_j1_.NewAthenaArray(n_cells_2);
  metric_face1_j2_.NewAthenaArray(n_cells_2);
  metric_face2_j1_.NewAthenaArray(n_cells_2);
  metric_face2_j2_.NewAthenaArray(n_cells_2);
  metric_face3_j1_.NewAthenaArray(n_cells_2);
  metric_face3_j2_.NewAthenaArray(n_cells_2);
  trans_face1_j1_.NewAthenaArray(n_cells_2);
  trans_face1_j2_.NewAthenaArray(n_cells_2);
  trans_face2_j1_.NewAthenaArray(n_cells_2);
  trans_face2_j2_.NewAthenaArray(n_cells_2);
  trans_face3_j1_.NewAthenaArray(n_cells_2);
  trans_face3_j2_.NewAthenaArray(n_cells_2);

  // Calculate intermediate geometric quantities: r-direction
#pragma simd
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST; i++)
  {
    // Useful quantities
    Real r_c = pb->x1v(i);
    Real r_m = pb->x1f(i);
    Real r_p = pb->x1f(i+1);

    // Volumes, areas, and lengths
    volume_i_(i) = 1.0/3.0 * (r_p*r_p*r_p - r_m*r_m*r_m);
    face1_area_i_(i) = r_m*r_m;
    face2_area_i_(i) = volume_i_(i);
    face3_area_i_(i) = volume_i_(i);
    edge1_length_i_(i) = volume_i_(i);
    edge2_length_i_(i) = face1_area_i_(i);
    edge3_length_i_(i) = face1_area_i_(i);

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
    metric_face1_i4_(i) = 1.0 / metric_face1_i1_(i);
    metric_face1_i5_(i) = 1.0 / metric_face1_i2_(i);
    metric_face1_i6_(i) = 1.0 / metric_face1_i3_(i);
    metric_face2_i1_(i) = metric_cell_i1_(i);
    metric_face2_i2_(i) = metric_cell_i2_(i);
    metric_face2_i3_(i) = metric_cell_i3_(i);
    metric_face2_i4_(i) = metric_cell_i4_(i);
    metric_face2_i5_(i) = metric_cell_i5_(i);
    metric_face2_i6_(i) = metric_cell_i6_(i);
    metric_face3_i1_(i) = metric_cell_i1_(i);
    metric_face3_i2_(i) = metric_cell_i2_(i);
    metric_face3_i3_(i) = metric_cell_i3_(i);
    metric_face3_i4_(i) = metric_cell_i4_(i);
    metric_face3_i5_(i) = metric_cell_i5_(i);
    metric_face3_i6_(i) = metric_cell_i6_(i);

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

      // Volumes, areas, and lengths
      volume_j_(j) = cos_m - cos_p;
      face1_area_j_(j) = volume_j_(j);
      face2_area_j_(j) = sin_m;
      face3_area_j_(j) = volume_j_(j);
      edge1_length_j_(j) = face2_area_j_(j);
      edge2_length_j_(j) = volume_j_(j);
      edge3_length_j_(j) = face2_area_j_(j);

      // Source terms
      src_terms_j1_(j) = 1.0/6.0 * (4.0 - std::cos(2.0*theta_m)
          - 2.0 * cos_m * cos_p - std::cos(2.0*theta_p));
      src_terms_j2_(j) = 1.0/3.0 * (sin_m*sin_m*sin_m - sin_p*sin_p*sin_p)
          / (cos_m - cos_p);
      src_terms_j3_(j) = (sin_p - sin_m) / (cos_m - cos_p);  // cot((theta_p+theta_m)/2)

      // Cell-centered metric
      metric_cell_j1_(j) = sin_c*sin_c;
      metric_cell_j2_(j) = 1.0 / metric_cell_j1_(j);

      // Face-centered metric
      metric_face1_j1_(j) = metric_cell_j1_(j);
      metric_face1_j2_(j) = 1.0 / metric_face1_j1_(j);
      metric_face2_j1_(j) = sin_m*sin_m;
      metric_face2_j2_(j) = 1.0 / metric_face2_j1_(j);
      metric_face3_j1_(j) = metric_cell_j1_(j);
      metric_face3_j2_(j) = 1.0 / metric_face3_j1_(j);

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
    Real theta_c = pb->x2v(pb->js);
    Real theta_m = pb->x2f(pb->js);
    Real theta_p = pb->x2f(pb->js+1);
    Real sin_c = std::sin(theta_c);
    Real sin_m = std::sin(theta_m);
    Real sin_p = std::sin(theta_p);
    Real cos_m = std::cos(theta_m);
    Real cos_p = std::cos(theta_p);

    // Volumes and areas
    volume_j_(pb->js) = cos_m - cos_p;
    face1_area_j_(pb->js) = volume_j_(pb->js);
    face2_area_j_(pb->js) = sin_m;
    face3_area_j_(pb->js) = volume_j_(pb->js);
    edge1_length_j_(pb->js) = face2_area_j_(pb->js);
    edge2_length_j_(pb->js) = volume_j_(pb->js);
    edge3_length_j_(pb->js) = face2_area_j_(pb->js);

    // Source terms
    src_terms_j1_(pb->js) = 1.0/6.0 * (4.0 - std::cos(2.0*theta_m)
        - 2.0 * cos_m * cos_p - std::cos(2.0*theta_p));
    src_terms_j2_(pb->js) = 1.0/3.0 * (sin_m*sin_m*sin_m - sin_p*sin_p*sin_p)
        / (cos_m - cos_p);
    src_terms_j3_(pb->js) = tan(PI/2.0 - 0.5 * (theta_m + theta_p));

    // Cell-centered metric
    metric_cell_j1_(pb->js) = sin_c*sin_c;
    metric_cell_j2_(pb->js) = 1.0 / metric_cell_j1_(pb->js);

    // Face-centered metric
    metric_face1_j1_(pb->js) = metric_cell_j1_(pb->js);
    metric_face1_j2_(pb->js) = 1.0 / metric_face1_j1_(pb->js);
    metric_face2_j1_(pb->js) = sin_m*sin_m;
    metric_face2_j2_(pb->js) = 1.0 / metric_face2_j1_(pb->js);
    metric_face3_j1_(pb->js) = metric_cell_j1_(pb->js);
    metric_face3_j2_(pb->js) = 1.0 / metric_face3_j1_(pb->js);

    // Coordinate transformations
    trans_face1_j1_(pb->js) = sin_c;
    trans_face1_j2_(pb->js) = 1.0 / sin_c;
    trans_face2_j1_(pb->js) = sin_m;
    trans_face2_j2_(pb->js) = 1.0 / sin_m;
    trans_face3_j1_(pb->js) = sin_m;
    trans_face3_j2_(pb->js) = trans_face2_j2_(pb->js);
  }
}

// Destructor
Coordinates::~Coordinates()
{
  volume_i_.DeleteAthenaArray();
  face1_area_i_.DeleteAthenaArray();
  face2_area_i_.DeleteAthenaArray();
  face3_area_i_.DeleteAthenaArray();
  edge1_length_i_.DeleteAthenaArray();
  edge2_length_i_.DeleteAthenaArray();
  edge3_length_i_.DeleteAthenaArray();
  src_terms_i1_.DeleteAthenaArray();
  src_terms_i2_.DeleteAthenaArray();
  src_terms_i3_.DeleteAthenaArray();
  src_terms_i4_.DeleteAthenaArray();
  volume_j_.DeleteAthenaArray();
  face1_area_j_.DeleteAthenaArray();
  face2_area_j_.DeleteAthenaArray();
  face3_area_j_.DeleteAthenaArray();
  edge1_length_j_.DeleteAthenaArray();
  edge2_length_j_.DeleteAthenaArray();
  edge3_length_j_.DeleteAthenaArray();
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
  metric_face1_i4_.DeleteAthenaArray();
  metric_face1_i5_.DeleteAthenaArray();
  metric_face1_i6_.DeleteAthenaArray();
  metric_face1_j1_.DeleteAthenaArray();
  metric_face1_j2_.DeleteAthenaArray();
  metric_face2_i1_.DeleteAthenaArray();
  metric_face2_i2_.DeleteAthenaArray();
  metric_face2_i3_.DeleteAthenaArray();
  metric_face2_i4_.DeleteAthenaArray();
  metric_face2_i5_.DeleteAthenaArray();
  metric_face2_i6_.DeleteAthenaArray();
  metric_face2_j1_.DeleteAthenaArray();
  metric_face2_j2_.DeleteAthenaArray();
  metric_face3_i1_.DeleteAthenaArray();
  metric_face3_i2_.DeleteAthenaArray();
  metric_face3_i3_.DeleteAthenaArray();
  metric_face3_i4_.DeleteAthenaArray();
  metric_face3_i5_.DeleteAthenaArray();
  metric_face3_i6_.DeleteAthenaArray();
  metric_face3_j1_.DeleteAthenaArray();
  metric_face3_j2_.DeleteAthenaArray();
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

// Function for computing lengths of r-edges
void Coordinates::Edge1Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &lengths)
{
  Real &sin_theta = edge1_length_j_(j);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    Real &third_delta_r_cb = edge1_length_i_(i);
    Real &length = lengths(i);
    length = third_delta_r_cb * sin_theta;
  }
  return;
}

// Function for computing lengths of theta-edges
void Coordinates::Edge2Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &lengths)
{
  Real &neg_delta_cos_theta = edge2_length_j_(j);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    Real &r_sq = edge2_length_i_(i);
    Real &length = lengths(i);
    length = r_sq * neg_delta_cos_theta;
  }
  return;
}

// Function for computing lengths of phi-edges
void Coordinates::Edge3Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &lengths)
{
  Real &sin_theta = edge3_length_j_(j);
  Real &delta_phi = pmy_block->dx3f(k);
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    Real &r_sq = edge3_length_i_(i);
    Real &length = lengths(i);
    length = r_sq * sin_theta * delta_phi;
  }
  return;
}

//--------------------------------------------------------------------------------------
// TODO: write functions
// Cell-center Width functions: returns physical width at cell-center
Real Coordinates::CenterWidth1(const int k, const int j, const int i)
{
  return (pmy_block->dx1f(i));
}
Real Coordinates::CenterWidth2(const int k, const int j, const int i)
{
  return (pmy_block->dx2f(j));
}
Real Coordinates::CenterWidth3(const int k, const int j, const int i)
{
  return (pmy_block->dx3f(k));
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
void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
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
void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
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
void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
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
    Real &volume = volumes(i);
    Real &third_delta_r_cb = volume_i_(i);
    volume = third_delta_r_cb * neg_delta_cos_theta * delta_phi;
  }
  return;
}

// Function for computing source terms
// Inputs:
//   dt: size of timestep
//   prim: full grid of primitive values at beginning of half timestep
//   cons: full grid of conserved variables at end of half timestep
// Outputs:
//   cons: source terms added
void Coordinates::CoordinateSourceTerms(Real dt, const AthenaArray<Real> &prim,
    AthenaArray<Real> &cons)
{
  // Extract ratio of specific heats
  const Real gamma_adi = pmy_block->pfluid->pf_eos->GetGamma();
  const Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Go through cells
  for (int k = pmy_block->ks; k <= pmy_block->ke; k++)
    for (int j = pmy_block->js; j <= pmy_block->je; j++)
    {
      // Extract geometric quantities that do not depend on r
      Real &gamma_233 = src_terms_j2_(j);
      Real &gamma_323 = src_terms_j3_(j);
      Real &gamma_332 = gamma_323;

      // Go through cells in r-direction
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

        // Calculate source terms
        Real s0 = gamma_010 * t10 + gamma_100 * t01;
        Real s1 = gamma_001 * t00 + gamma_111 * t11 + gamma_221 * t22 + gamma_331 * t33;
        Real s2 = gamma_122 * t21 + gamma_212 * t12 + gamma_332 * t33;
        Real s3 = gamma_133 * t31 + gamma_233 * t32 + gamma_313 * t13 + gamma_323 * t23;

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
    }
  return;
}

// Function for computing cell-centered metric coefficients
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
  // Extract geometric quantities that do not depend on r
  Real &sin_sq_theta = metric_face1_j1_(j);
  Real &csc_sq_theta = metric_face1_j2_(j);

  // Go through 1D block of faces
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie+1; i++)
  {
    // Extract remaining geometric quantities
    Real &t_factor = metric_face1_i1_(i);
    Real &r_factor = metric_face1_i2_(i);
    Real &r_sq = metric_face1_i3_(i);
    Real &t_inv_factor = metric_face1_i4_(i);
    Real &r_inv_factor = metric_face1_i5_(i);
    Real &r_inv_sq = metric_face1_i6_(i);

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
  // Extract geometric quantities that do not depend on r
  Real &sin_sq_theta = metric_face2_j1_(j);
  Real &csc_sq_theta = metric_face2_j2_(j);

  // Go through 1D block of faces
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
  {
    // Extract remaining geometric quantities
    Real &t_factor = metric_face2_i1_(i);
    Real &r_factor = metric_face2_i2_(i);
    Real &r_sq = metric_face2_i3_(i);
    Real &t_inv_factor = metric_face2_i4_(i);
    Real &r_inv_factor = metric_face2_i5_(i);
    Real &r_inv_sq = metric_face2_i6_(i);

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
  // Extract geometric quantities that do not depend on r
  Real &sin_sq_theta = metric_face3_j1_(j);
  Real &csc_sq_theta = metric_face3_j2_(j);

  // Go through 1D block of faces
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
  {
    // Extract remaining geometric quantities
    Real &t_factor = metric_face3_i1_(i);
    Real &r_factor = metric_face3_i2_(i);
    Real &r_sq = metric_face3_i3_(i);
    Real &t_inv_factor = metric_face3_i4_(i);
    Real &r_inv_factor = metric_face3_i5_(i);
    Real &r_inv_sq = metric_face3_i6_(i);

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
//   b1_vals: 3D array of normal components B^1 of magnetic field, in global coordinates
//   prim_left: 1D array of left primitives, using global coordinates
//   prim_right: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_left: values overwritten in local coordinates
//   prim_right: values overwritten in local coordinates
//   bx: 1D array of longitudinal magnetic fields, in local coordinates
void Coordinates::PrimToLocal1(const int k, const int j,
    const AthenaArray<Real> &b1_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bx)
{
  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie+1; i++)
  {
    // Extract geometric quantities
    Real &g00 = metric_face1_i1_(i);
    Real &g11 = metric_face1_i2_(i);
    Real &g22 = metric_face1_i3_(i);
    Real g33 = metric_face1_i3_(i) * metric_face1_j1_(j);
    Real &mt0 = trans_face1_i1_(i);
    Real &mx1 = trans_face1_i2_(i);
    Real &my2 = trans_face1_i3_(i);
    Real mz3 = trans_face1_i3_(i) * trans_face1_j1_(j);

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
    Real utl = mt0 * u0l;
    Real uxl = mx1 * u1l;
    Real uyl = my2 * u2l;
    Real uzl = mz3 * u3l;
    Real utr = mt0 * u0r;
    Real uxr = mx1 * u1r;
    Real uyr = my2 * u2r;
    Real uzr = mz3 * u3r;

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

      // Construct global covariant magnetic fields
      Real bcov0l = g11*b1*u1l + g22*b2l*u2l + g33*b3l*u3l;
      Real bcov1l = (b1 + bcov0l * u1l) / u0l;
      Real bcov2l = (b2l + bcov0l * u2l) / u0l;
      Real bcov3l = (b3l + bcov0l * u3l) / u0l;
      Real bcov0r = g11*b1*u1r + g22*b2r*u2r + g33*b3r*u3r;
      Real bcov1r = (b1 + bcov0r * u1r) / u0r;
      Real bcov2r = (b2r + bcov0r * u2r) / u0r;
      Real bcov3r = (b3r + bcov0r * u3r) / u0r;

      // Transform covariant magnetic fields
      Real bcovtl = mt0 * bcov0l;
      Real bcovxl = mx1 * bcov1l;
      Real bcovyl = my2 * bcov2l;
      Real bcovzl = mz3 * bcov3l;
      Real bcovtr = mt0 * bcov0r;
      Real bcovxr = mx1 * bcov1r;
      Real bcovyr = my2 * bcov2r;
      Real bcovzr = mz3 * bcov3r;

      // Set local magnetic fields
      // TODO: deal with possible inconsistent B^x (shouldn't happen here)
      Real bxl = utl * bcovxl - uxl * bcovtl;
      Real bxr = utr * bcovxr - uxr * bcovtr;
      bx(i) = 0.5 * (bxl + bxr);
      b2l = utl * bcovyl - uyl * bcovtl;
      b3l = utl * bcovzl - uzl * bcovtl;
      b2r = utr * bcovyr - uyr * bcovtr;
      b3r = utr * bcovzr - uzr * bcovtr;
    }
  }
  return;
}

// Function for transforming primitives to locally flat frame: theta-interface
// Inputs:
//   k: phi-index
//   j: theta-index
//   b2_vals: 3D array of normal components B^2 of magnetic field, in global coordinates
//   prim_left: 1D array of left primitives, using global coordinates
//   prim_right: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_left: values overwritten in local coordinates
//   prim_right: values overwritten in local coordinates
//   by: 1D array of longitudinal magnetic fields, in local coordinates
void Coordinates::PrimToLocal2(const int k, const int j,
    const AthenaArray<Real> &b2_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &by)
{
  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
  {
    // Extract geometric quantities
    Real &g00 = metric_face2_i1_(i);
    Real &g11 = metric_face2_i2_(i);
    Real &g22 = metric_face2_i3_(i);
    Real g33 = metric_face2_i3_(i) * metric_face2_j1_(j);
    Real &mt0 = trans_face2_i1_(i);
    Real &mx1 = trans_face2_i2_(i);
    Real &my2 = trans_face2_i3_(i);
    Real mz3 = trans_face2_i3_(i) * trans_face2_j1_(j);

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
    Real utl = mt0 * u0l;
    Real uxl = mx1 * u1l;
    Real uyl = my2 * u2l;
    Real uzl = mz3 * u3l;
    Real utr = mt0 * u0r;
    Real uxr = mx1 * u1r;
    Real uyr = my2 * u2r;
    Real uzr = mz3 * u3r;

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
      const Real &b2 = b2_vals(k,j,i);
      Real &b3l = prim_left(IBY,i);
      Real &b1l = prim_left(IBZ,i);
      Real &b3r = prim_right(IBY,i);
      Real &b1r = prim_right(IBZ,i);

      // Construct global covariant magnetic fields
      Real bcov0l = g11*b1l*u1l + g22*b2*u2l + g33*b3l*u3l;
      Real bcov1l = (b1l + bcov0l * u1l) / u0l;
      Real bcov2l = (b2 + bcov0l * u2l) / u0l;
      Real bcov3l = (b3l + bcov0l * u3l) / u0l;
      Real bcov0r = g11*b1r*u1r + g22*b2*u2r + g33*b3r*u3r;
      Real bcov1r = (b1r + bcov0r * u1r) / u0r;
      Real bcov2r = (b2 + bcov0r * u2r) / u0r;
      Real bcov3r = (b3r + bcov0r * u3r) / u0r;

      // Transform covariant magnetic fields
      Real bcovtl = mt0 * bcov0l;
      Real bcovxl = mx1 * bcov1l;
      Real bcovyl = my2 * bcov2l;
      Real bcovzl = mz3 * bcov3l;
      Real bcovtr = mt0 * bcov0r;
      Real bcovxr = mx1 * bcov1r;
      Real bcovyr = my2 * bcov2r;
      Real bcovzr = mz3 * bcov3r;

      // Set local magnetic fields
      // TODO: deal with possible inconsistent B^y (shouldn't happen here)
      Real byl = utl * bcovyl - uyl * bcovtl;
      Real byr = utr * bcovyr - uyr * bcovtr;
      by(i) = 0.5 * (byl + byr);
      b3l = utl * bcovzl - uzl * bcovtl;
      b1l = utl * bcovxl - uxl * bcovtl;
      b3r = utr * bcovzr - uzr * bcovtr;
      b1r = utr * bcovxr - uxr * bcovtr;
    }
  }
  return;
}

// Function for transforming primitives to locally flat frame: phi-interface
// Inputs:
//   k: phi-index
//   j: theta-index
//   b3_vals: 3D array of normal components B^3 of magnetic field, in global coordinates
//   prim_left: 1D array of left primitives, using global coordinates
//   prim_right: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_left: values overwritten in local coordinates
//   prim_right: values overwritten in local coordinates
//   bz: 1D array of longitudinal magnetic fields, in local coordinates
void Coordinates::PrimToLocal3(const int k, const int j,
    const AthenaArray<Real> &b3_vals, AthenaArray<Real> &prim_left,
    AthenaArray<Real> &prim_right, AthenaArray<Real> &bz)
{
  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
  {
    // Extract geometric quantities
    Real &g00 = metric_face3_i1_(i);
    Real &g11 = metric_face3_i2_(i);
    Real &g22 = metric_face3_i3_(i);
    Real g33 = metric_face3_i3_(i) * metric_face3_j1_(j);
    Real &mt0 = trans_face3_i1_(i);
    Real &mx1 = trans_face3_i2_(i);
    Real &my2 = trans_face3_i3_(i);
    Real mz3 = trans_face3_i3_(i) * trans_face3_j1_(j);

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
    Real utl = mt0 * u0l;
    Real uxl = mx1 * u1l;
    Real uyl = my2 * u2l;
    Real uzl = mz3 * u3l;
    Real utr = mt0 * u0r;
    Real uxr = mx1 * u1r;
    Real uyr = my2 * u2r;
    Real uzr = mz3 * u3r;

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
      const Real &b3 = b3_vals(k,j,i);
      Real &b1l = prim_left(IBY,i);
      Real &b2l = prim_left(IBZ,i);
      Real &b1r = prim_right(IBY,i);
      Real &b2r = prim_right(IBZ,i);

      // Construct global covariant magnetic fields
      Real bcov0l = g11*b1l*u1l + g22*b2l*u2l + g33*b3*u3l;
      Real bcov1l = (b1l + bcov0l * u1l) / u0l;
      Real bcov2l = (b2l + bcov0l * u2l) / u0l;
      Real bcov3l = (b3 + bcov0l * u3l) / u0l;
      Real bcov0r = g11*b1r*u1r + g22*b2r*u2r + g33*b3*u3r;
      Real bcov1r = (b1r + bcov0r * u1r) / u0r;
      Real bcov2r = (b2r + bcov0r * u2r) / u0r;
      Real bcov3r = (b3 + bcov0r * u3r) / u0r;

      // Transform covariant magnetic fields
      Real bcovtl = mt0 * bcov0l;
      Real bcovxl = mx1 * bcov1l;
      Real bcovyl = my2 * bcov2l;
      Real bcovzl = mz3 * bcov3l;
      Real bcovtr = mt0 * bcov0r;
      Real bcovxr = mx1 * bcov1r;
      Real bcovyr = my2 * bcov2r;
      Real bcovzr = mz3 * bcov3r;

      // Set local magnetic fields
      // TODO: deal with possible inconsistent B^z (shouldn't happen here)
      Real bzl = utl * bcovzl - uzl * bcovtl;
      Real bzr = utr * bcovzr - uzr * bcovtr;
      bz(i) = 0.5 * (bzl + bzr);
      b1l = utl * bcovxl - uxl * bcovtl;
      b2l = utl * bcovyl - uyl * bcovtl;
      b1r = utr * bcovxr - uxr * bcovtr;
      b2r = utr * bcovyr - uyr * bcovtr;
    }
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

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      Real fby = flux(IBY,i);
      Real fbz = flux(IBZ,i);
      flux(IBY,i) = m1x * m2y * fby;
      flux(IBZ,i) = m1x * m3z * fbz;
    }
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

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      Real fbz = flux(IBY,i);
      Real fbx = flux(IBZ,i);
      flux(IBY,i) = m2y * m3z * fbz;
      flux(IBZ,i) = m2y * m1x * fbx;
    }
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

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      Real fbx = flux(IBY,i);
      Real fby = flux(IBZ,i);
      flux(IBY,i) = m3z * m1x * fbx;
      flux(IBZ,i) = m3z * m2y * fby;
    }
  }
  return;
}

// Function for converting all primitives to conserved variables
// Inputs:
//   prim: 3D array of primitives
//   b: 3D array of cell-centered magnetic fields
// Outputs:
//   cons: 3D array of conserved variables
void Coordinates::PrimToCons(const AthenaArray<Real> &prim, const AthenaArray<Real> &b,
    AthenaArray<Real> &cons)
{
  // Extract ratio of specific heats
  const Real gamma_adi = pmy_block->pfluid->pf_eos->GetGamma();
  const Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

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
#pragma simd
      for (int i = il; i <= iu; i++)
      {
        // Extract geometric quantities
        Real &g00 = metric_cell_i1_(i);
        Real &g11 = metric_cell_i2_(i);
        Real &g22 = metric_cell_i3_(i);
        Real g33 = metric_cell_i3_(i) * metric_cell_j1_(j);

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
        Real u_0 = g00 * u0;
        Real u_1 = g11 * u1;
        Real u_2 = g22 * u2;
        Real u_3 = g33 * u3;

        // Calculate covariant magnetic field
        Real bcov0 = g11*b1*u1 + g22*b2*u2 + g33*b3*u3;
        Real bcov1 = 1.0/u0 * (b1 + bcov0 * u1);
        Real bcov2 = 1.0/u0 * (b2 + bcov0 * u2);
        Real bcov3 = 1.0/u0 * (b3 + bcov0 * u3);
        Real bcov_0 = g00 * bcov0;
        Real bcov_1 = g11 * bcov1;
        Real bcov_2 = g22 * bcov2;
        Real bcov_3 = g33 * bcov3;
        Real bcov_sq = bcov0*bcov_0 + bcov1*bcov_1 + bcov2*bcov_2 + bcov3*bcov_3;

        // Extract conserved quantities
        Real &rho_u0 = cons(IDN,k,j,i);
        Real &t0_0 = cons(IEN,k,j,i);
        Real &t0_1 = cons(IVX,k,j,i);
        Real &t0_2 = cons(IVY,k,j,i);
        Real &t0_3 = cons(IVZ,k,j,i);

        // Set conserved quantities
        rho_u0 = rho * u0;
        Real rho_h = rho + gamma_adi_red * pgas;
        t0_0 = (rho_h + bcov_sq) * u0 * u_0 - bcov0 * bcov_0 + pgas + 0.5*bcov_sq;
        t0_1 = (rho_h + bcov_sq) * u0 * u_1 - bcov0 * bcov_1;
        t0_2 = (rho_h + bcov_sq) * u0 * u_2 - bcov0 * bcov_2;
        t0_3 = (rho_h + bcov_sq) * u0 * u_3 - bcov0 * bcov_3;
      }
  return;
}
