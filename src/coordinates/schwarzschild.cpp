//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file schwarzschild.cpp
//  \brief implements functions for Schwarzschild spacetime and spherical (t,r,theta,phi)
//  coordinates in a derived class of the Coordinates abstract base class.
//  Original implementation by CJ White.
//
// Notes:
//   coordinates: t, r, theta, phi
//   parameters: M (mass)
//   metric:
//     ds^2 = -\alpha^2 dt^2 + 1/\alpha^2 * dr^2 + r^2 (d\theta^2 + \sin^2\theta d\phi^2)
//     where \alpha = \sqrt(1 - 2M/r)

// C++ headers
#include <cmath>  // abs(), acos(), cos(), log(), pow(), sin(), sqrt()

// Athena++ headers
#include "coordinates.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../eos/eos.hpp"
#include "../mesh/mesh.hpp"

//----------------------------------------------------------------------------------------
// Schwarzschild Constructor
// Inputs:
//   pmb: pointer to MeshBlock containing this grid
//   pin: pointer to runtime inputs
//   flag: true if object is for coarse grid only in an AMR calculation

Schwarzschild::Schwarzschild(MeshBlock *pmb, ParameterInput *pin, bool flag)
  : Coordinates(pmb, pin, flag) {
  // Set indices
  pmy_block = pmb;
  coarse_flag = flag;
  int il, iu, jl, ju, kl, ku, ng;
  if (coarse_flag == true) {
    il = pmb->cis;
    iu = pmb->cie;
    jl = pmb->cjs;
    ju = pmb->cje;
    kl = pmb->cks;
    ku = pmb->cke;
    ng = pmb->cnghost;
  } else {
    il = pmb->is;
    iu = pmb->ie;
    jl = pmb->js;
    ju = pmb->je;
    kl = pmb->ks;
    ku = pmb->ke;
    ng = NGHOST;
  }
  Mesh *pm = pmy_block->pmy_mesh;
  RegionSize& mesh_size = pmy_block->pmy_mesh->mesh_size;
  RegionSize& block_size = pmy_block->block_size;

  // Allocate arrays for volume-centered coordinates and positions of cells
  int ncells1 = (iu-il+1) + 2*ng;
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = (ju-jl+1) + 2*ng;
  if (block_size.nx3 > 1) ncells3 = (ku-kl+1) + 2*ng;
  dx1v.NewAthenaArray(ncells1);
  dx2v.NewAthenaArray(ncells2);
  dx3v.NewAthenaArray(ncells3);
  x1v.NewAthenaArray(ncells1);
  x2v.NewAthenaArray(ncells2);
  x3v.NewAthenaArray(ncells3);

  // Allocate arrays for area weighted positions for AMR/SMR MHD
  if (pm->multilevel && MAGNETIC_FIELDS_ENABLED) {
    x1s2.NewAthenaArray(ncells1);
    x1s3.NewAthenaArray(ncells1);
    x2s1.NewAthenaArray(ncells2);
    x2s3.NewAthenaArray(ncells2);
    x3s1.NewAthenaArray(ncells3);
    x3s2.NewAthenaArray(ncells3);
  }

  // Set parameters
  bh_mass_ = pin->GetReal("coord", "m");
  const Real &m = bh_mass_;

  // Initialize volume-averaged coordinates and spacings: r-direction
  for (int i = il-ng; i <= iu+ng; ++i) {
    Real r_m = x1f(i);
    Real r_p = x1f(i+1);
    x1v(i) = std::pow(0.5 * (r_m*r_m*r_m + r_p*r_p*r_p), ONE_3RD);
  }
  for (int i = il-ng; i <= iu+ng-1; ++i) {
    dx1v(i) = x1v(i+1) - x1v(i);
  }

  // Initialize volume-averaged coordinates and spacings: theta-direction
  if (pmb->block_size.nx2 == 1) {
    Real theta_m = x2f(jl);
    Real theta_p = x2f(jl+1);
    x2v(jl) = std::acos(0.5 * (std::cos(theta_m) + std::cos(theta_p)));
    dx2v(jl) = dx2f(jl);
  } else {
    for (int j = jl-ng; j <= ju+ng; ++j) {
      Real theta_m = x2f(j);
      Real theta_p = x2f(j+1);
      x2v(j) = std::acos(0.5 * (std::cos(theta_m) + std::cos(theta_p)));
    }
    for (int j = jl-ng; j <= ju+ng-1; ++j) {
      dx2v(j) = x2v(j+1) - x2v(j);
    }
  }

  // Initialize volume-averaged coordinates and spacings: phi-direction
  if (pmb->block_size.nx3 == 1) {
    Real phi_m = x3f(kl);
    Real phi_p = x3f(kl+1);
    x3v(kl) = 0.5 * (phi_m + phi_p);
    dx3v(kl) = dx3f(kl);
  } else {
    for (int k = kl-ng; k <= ku+ng; ++k) {
      Real phi_m = x3f(k);
      Real phi_p = x3f(k+1);
      x3v(k) = 0.5 * (phi_m + phi_p);
    }
    for (int k = kl-ng; k <= ku+ng-1; ++k) {
      dx3v(k) = x3v(k+1) - x3v(k);
    }
  }

  // Initialize area-averaged coordinates used with MHD AMR
  if (pmb->pmy_mesh->multilevel && MAGNETIC_FIELDS_ENABLED) {
    for (int i = il-ng; i <= iu+ng; ++i) {
      x1s2(i) = x1s3(i) = x1v(i);
    }
    if (pmb->block_size.nx2 == 1) {
      x2s1(jl) = x2s3(jl) = x2v(jl);
    } else {
      for (int j = jl-ng; j <= ju+ng; ++j) {
        x2s1(j) = x2s3(j) = x2v(j);
      }
    }
    if (pmb->block_size.nx3 == 1) {
      x3s1(kl) = x3s2(kl) = x3v(kl);
    } else {
      for (int k = kl-ng; k <= ku+ng; ++k) {
        x3s1(k) = x3s2(k) = x3v(k);
      }
    }
  }

  // Allocate and compute arrays for intermediate geometric quantities always needed
  metric_cell_i1_.NewAthenaArray(ncells1);
  metric_cell_j1_.NewAthenaArray(ncells2);
  for (int i = il-ng; i <= iu+ng; ++i) {
    Real r_c = x1v(i);
    Real alpha_c = std::sqrt(1.0 - 2.0*m/r_c);
    metric_cell_i1_(i) = SQR(alpha_c);
  }

  int jll, juu;
  if (pmb->block_size.nx2 > 1) {
    jll = jl - ng; juu = ju + ng;
  } else {
    jll = jl; juu = ju;
  }
  for (int j = jll; j <= juu; ++j) {
    Real sin_c = std::sin(x2v(j));
    Real sin_c_sq = SQR(sin_c);
    metric_cell_j1_(j) = sin_c_sq;
  }

  // Allocate and compute arrays for intermediate geometric quantities that are only
  // needed if object is NOT a coarse mesh
  if (coarse_flag == false) {

    // Allocate arrays for intermediate geometric quantities: r-direction
    coord_vol_i1_.NewAthenaArray(ncells1);
    coord_area1_i1_.NewAthenaArray(ncells1+1);
    coord_area2_i1_.NewAthenaArray(ncells1);
    coord_area3_i1_.NewAthenaArray(ncells1);
    coord_len1_i1_.NewAthenaArray(ncells1);
    coord_len2_i1_.NewAthenaArray(ncells1+1);
    coord_len3_i1_.NewAthenaArray(ncells1+1);
    coord_width1_i1_.NewAthenaArray(ncells1);
    metric_face1_i1_.NewAthenaArray(ncells1+1);
    metric_face2_i1_.NewAthenaArray(ncells1);
    metric_face3_i1_.NewAthenaArray(ncells1);
    trans_face1_i1_.NewAthenaArray(ncells1+1);
    trans_face2_i1_.NewAthenaArray(ncells1);
    trans_face3_i1_.NewAthenaArray(ncells1);
    g_.NewAthenaArray(NMETRIC, ncells1+1);
    gi_.NewAthenaArray(NMETRIC, ncells1+1);

    // Allocate arrays for intermediate geometric quantities: theta-direction
    coord_vol_j1_.NewAthenaArray(ncells2);
    coord_area1_j1_.NewAthenaArray(ncells2);
    coord_area2_j1_.NewAthenaArray(ncells2+1);
    coord_area3_j1_.NewAthenaArray(ncells2);
    coord_len1_j1_.NewAthenaArray(ncells2+1);
    coord_len2_j1_.NewAthenaArray(ncells2);
    coord_len3_j1_.NewAthenaArray(ncells2+1);
    coord_width3_j1_.NewAthenaArray(ncells2);
    coord_src_j1_.NewAthenaArray(ncells2);
    coord_src_j2_.NewAthenaArray(ncells2);
    metric_face1_j1_.NewAthenaArray(ncells2);
    metric_face2_j1_.NewAthenaArray(ncells2+1);
    metric_face3_j1_.NewAthenaArray(ncells2);
    trans_face1_j1_.NewAthenaArray(ncells2);
    trans_face2_j1_.NewAthenaArray(ncells2+1);
    trans_face3_j1_.NewAthenaArray(ncells2);

    // Calculate intermediate geometric quantities: r-direction
    for (int i = il-ng; i <= iu+ng; ++i) {

      // Useful quantities
      Real r_c = x1v(i);
      Real r_m = x1f(i);
      Real r_p = x1f(i+1);
      Real alpha_c = std::sqrt(1.0 - 2.0*m/r_c);
      Real alpha_m = std::sqrt(1.0 - 2.0*m/r_m);
      Real alpha_p = std::sqrt(1.0 - 2.0*m/r_p);
      Real r_p_cu = r_p*r_p*r_p;
      Real r_m_cu = r_m*r_m*r_m;

      // Volumes, areas, lengths, and widths
      coord_vol_i1_(i) = ONE_3RD * (r_p_cu - r_m_cu);
      coord_area1_i1_(i) = SQR(r_m);
      if (i == (iu+ng)) {
        coord_area1_i1_(i+1) = SQR(r_p);
      }
      coord_area2_i1_(i) = coord_vol_i1_(i);
      coord_area3_i1_(i) = coord_vol_i1_(i);
      coord_len1_i1_(i) = coord_vol_i1_(i);
      coord_len2_i1_(i) = coord_area1_i1_(i);
      coord_len3_i1_(i) = coord_area1_i1_(i);
      if (i == (iu+ng)) {
        coord_len2_i1_(i+1) = coord_area1_i1_(i+1);
        coord_len3_i1_(i+1) = coord_area1_i1_(i+1);
      }
      coord_width1_i1_(i) = r_p*alpha_p - r_m*alpha_m
          + m * std::log((r_p*(1.0+alpha_p)-m) / (r_m*(1.0+alpha_m)-m));

      // Metric coefficients
      metric_face1_i1_(i) = SQR(alpha_m);
      if (i == (iu+ng)) {
        metric_face1_i1_(i+1) = SQR(alpha_p);
      }
      metric_face2_i1_(i) = SQR(alpha_c);
      metric_face3_i1_(i) = SQR(alpha_c);

      // Coordinate transformations
      trans_face1_i1_(i) = alpha_m;
      if (i == (iu+ng)) {
        trans_face1_i1_(i+1) = alpha_p;
      }
      trans_face2_i1_(i) = alpha_c;
      trans_face3_i1_(i) = alpha_c;
    }

    // Calculate intermediate geometric quantities: theta-direction
    for (int j = jll; j <= juu; ++j) {

      // Useful quantities
      Real theta_c = x2v(j);
      Real theta_m = x2f(j);
      Real theta_p = x2f(j+1);
      Real sin_c = std::sin(theta_c);
      Real sin_m = std::sin(theta_m);
      Real sin_p = std::sin(theta_p);
      Real cos_c = std::cos(theta_c);
      Real cos_m = std::cos(theta_m);
      Real cos_p = std::cos(theta_p);
      Real sin_c_sq = SQR(sin_c);
      Real sin_m_sq = SQR(sin_m);
      Real sin_p_sq = SQR(sin_p);

      // Volumes, areas, lengths, and widths
      coord_vol_j1_(j) = std::abs(cos_m - cos_p);
      coord_area1_j1_(j) = coord_vol_j1_(j);
      coord_area2_j1_(j) = std::abs(sin_m);
      if (j == juu) {
        coord_area2_j1_(j+1) = std::abs(sin_p);
      }
      coord_area3_j1_(j) = coord_vol_j1_(j);
      coord_len1_j1_(j) = coord_area2_j1_(j);
      if (j == juu) {
        coord_len1_j1_(j+1) = coord_area2_j1_(j+1);
      }
      coord_len2_j1_(j) = coord_vol_j1_(j);
      coord_len3_j1_(j) = coord_area2_j1_(j);
      if (j == juu) {
        coord_len3_j1_(j+1) = coord_area2_j1_(j+1);
      }
      coord_width3_j1_(j) = std::abs(sin_c);

      // Source terms
      coord_src_j1_(j) = sin_c;
      coord_src_j2_(j) = cos_c;

      // Metric coefficients
      metric_face1_j1_(j) = sin_c_sq;
      metric_face2_j1_(j) = sin_m_sq;
      if (j == juu) {
        metric_face2_j1_(j+1) = sin_p_sq;
      }
      metric_face3_j1_(j) = sin_c_sq;

      // Coordinate transformations
      trans_face1_j1_(j) = std::abs(sin_c);
      trans_face2_j1_(j) = std::abs(sin_m);
      if (j == juu) {
        trans_face2_j1_(j+1) = std::abs(sin_p);
      }
      trans_face3_j1_(j) = std::abs(sin_c);
    }
  }
}

//----------------------------------------------------------------------------------------
// Destructor

Schwarzschild::~Schwarzschild() {
  dx1v.DeleteAthenaArray();
  dx2v.DeleteAthenaArray();
  dx3v.DeleteAthenaArray();
  x1v.DeleteAthenaArray();
  x2v.DeleteAthenaArray();
  x3v.DeleteAthenaArray();
  if (pmy_block->pmy_mesh->multilevel && MAGNETIC_FIELDS_ENABLED) {
    x1s2.DeleteAthenaArray();
    x1s3.DeleteAthenaArray();
    x2s1.DeleteAthenaArray();
    x2s3.DeleteAthenaArray();
    x3s1.DeleteAthenaArray();
    x3s2.DeleteAthenaArray();
  }
  metric_cell_i1_.DeleteAthenaArray();
  metric_cell_j1_.DeleteAthenaArray();
  if (coarse_flag == false) {
    coord_vol_i1_.DeleteAthenaArray();
    coord_area1_i1_.DeleteAthenaArray();
    coord_area2_i1_.DeleteAthenaArray();
    coord_area3_i1_.DeleteAthenaArray();
    coord_len1_i1_.DeleteAthenaArray();
    coord_len2_i1_.DeleteAthenaArray();
    coord_len3_i1_.DeleteAthenaArray();
    coord_width1_i1_.DeleteAthenaArray();
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
    g_.DeleteAthenaArray();
    gi_.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
// EdgeXLength functions: compute physical length at cell edge-X as vector
// Edge1(i,j,k) located at (i,j-1/2,k-1/2), i.e. (x1v(i), x2f(j), x3f(k))
// Edge2(i,j,k) located at (i-1/2,j,k-1/2), i.e. (x1f(i), x2v(j), x3f(k))
// Edge3(i,j,k) located at (i-1/2,j-1/2,k), i.e. (x1f(i), x2f(j), x3v(k))

void Schwarzschild::Edge1Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &lengths) {
  // \Delta L = 1/3 \Delta(r^3) \sin\theta
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {
    lengths(i) = coord_len1_i1_(i) * coord_len1_j1_(j);
  }
  return;
}

void Schwarzschild::Edge2Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &lengths) {
  // \Delta L = r^2 (-\Delta\cos\theta)
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {
    lengths(i) = coord_len2_i1_(i) * coord_len2_j1_(j);
  }
  return;
}

void Schwarzschild::Edge3Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &lengths) {
  // \Delta L = r^2 \sin\theta \Delta\phi
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {
    lengths(i) = coord_len3_i1_(i) * coord_len3_j1_(j) * dx3f(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetEdgeXLength functions: return length of edge-X at (i,j,k)

Real Schwarzschild::GetEdge1Length(const int k, const int j, const int i) {
  // \Delta L = 1/3 \Delta(r^3) \sin\theta
  return coord_len1_i1_(i) * coord_len1_j1_(j);
}

Real Schwarzschild::GetEdge2Length(const int k, const int j, const int i) {
  // \Delta L = r^2 (-\Delta\cos\theta)
  return coord_len2_i1_(i) * coord_len2_j1_(j);
}

Real Schwarzschild::GetEdge3Length(const int k, const int j, const int i) {
  // \Delta L = r^2 \sin\theta \Delta\phi
  return coord_len3_i1_(i) * coord_len3_j1_(j) * dx3f(k);
}

//----------------------------------------------------------------------------------------
// CenterWidthX functions: return physical width in X-dir at (i,j,k) cell-center

void Schwarzschild::CenterWidth1(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &dx1) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // \Delta W = \Delta(r \alpha) + M \Delta\log(r(1+\alpha)-M)
    dx1(i) = coord_width1_i1_(i);
  }
  return;
}

void Schwarzschild::CenterWidth2(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &dx2) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // \Delta W = r \Delta\theta
    dx2(i) = x1v(i) * dx1f(j);
  }
  return;
}

void Schwarzschild::CenterWidth3(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &dx3) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // \Delta W = r \sin\theta \Delta\phi
    dx3(i) = x1v(i) * coord_width3_j1_(j) * dx3f(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// FaceXArea functions: compute area of face with normal in X-dir as vector
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to X-face

void Schwarzschild::Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas) {
  //  \Delta A = r^2 (-\Delta\cos\theta) \Delta\phi
  #pragma omp simd
  for (int i = il; i <= iu; ++i)
    areas(i) = coord_area1_i1_(i) * coord_area1_j1_(j) * dx3f(k);
  return;
}

void Schwarzschild::Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas) {
  // \Delta A = 1/3 \Delta(r^3) \sin\theta \Delta\phi
  #pragma omp simd
  for (int i = il; i <= iu; ++i)
    areas(i) = coord_area2_i1_(i) * coord_area2_j1_(j) * dx3f(k);
  return;
}

void Schwarzschild::Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas) {
  // \Delta A = 1/3 \Delta(r^3) (-\Delta\cos\theta)
  #pragma omp simd
  for (int i = il; i <= iu; ++i)
    areas(i) = coord_area3_i1_(i) * coord_area3_j1_(j);
  return;
}

//----------------------------------------------------------------------------------------
// GetFaceXArea functions: return area of face with normal in X-dir at (i,j,k)
// Inputs:
//   k,j,i: phi- theta- and r-indices
// return:
//   interface area orthogonal to X-face

Real Schwarzschild::GetFace1Area(const int k, const int j, const int i) {
  // \Delta A = r^2 (-\Delta\cos\theta) \Delta\phi
  return coord_area1_i1_(i) * coord_area1_j1_(j) * dx3f(k);
}

Real Schwarzschild::GetFace2Area(const int k, const int j, const int i) {
  // \Delta A = 1/3 \Delta(r^3) \sin\theta \Delta\phi
  return coord_area2_i1_(i) * coord_area2_j1_(j) * dx3f(k);
}

Real Schwarzschild::GetFace3Area(const int k, const int j, const int i) {
  // \Delta A = 1/3 \Delta(r^3) (-\Delta\cos\theta)
  return coord_area3_i1_(i) * coord_area3_j1_(j);
}

//----------------------------------------------------------------------------------------
// Cell Volume function: compute volume of cell as vector
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   volumes: 1D array of cell volumes
// Notes:
//   \Delta V = 1/3 * \Delta(r^3) (-\Delta\cos\theta) \Delta\phi

void Schwarzschild::CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &volumes) {
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {
    volumes(i) = coord_vol_i1_(i) * coord_vol_j1_(j) * dx3f(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetCellVolume: returns cell volume at (i,j,k)
// Inputs:
//   k,j,i: phi-, theta-, and r-indices
// Outputs:
//   returned value: cell volume
// Notes:
//   \Delta V = 1/3 * \Delta(r^3) (-\Delta\cos\theta) \Delta\phi

Real Schwarzschild::GetCellVolume(const int k, const int j, const int i) {
  return coord_vol_i1_(i) * coord_vol_j1_(j) * dx3f(k);
}

//----------------------------------------------------------------------------------------
// Coordinate (geometric) source term function
// Inputs:
//   dt: size of timestep
//   flux: 3D array of fluxes
//   prim: 3D array of primitive values at beginning of half timestep
//   bb_cc: 3D array of cell-centered magnetic fields
// Outputs:
//   cons: source terms added to 3D array of conserved variables

void Schwarzschild::CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bb_cc,
    AthenaArray<Real> &cons) {
  // Extract ratio of specific heats
  const Real gamma_adi = pmy_block->peos->GetGamma();

  // Extract geometric quantities that do not depend on location
  const Real &m = bh_mass_;

  // Go through cells
  for (int k = pmy_block->ks; k <= pmy_block->ke; ++k) {
    for (int j = pmy_block->js; j <= pmy_block->je; ++j) {

      // Extract geometric quantities that do not depend on r
      const Real &sin = coord_src_j1_(j);
      const Real &cos = coord_src_j2_(j);
      Real sin2 = SQR(sin);
      Real sincos = sin * cos;

      // Calculate metric coefficients
      CellMetric(k, j, pmy_block->is, pmy_block->ie, g_, gi_);

      // Go through 1D slice
      #pragma omp simd
      for (int i = pmy_block->is; i <= pmy_block->ie; ++i) {

        // Extract geometric quantities
        const Real &g_00 = g_(I00,i);
        const Real &g_11 = g_(I11,i);
        const Real &g_22 = g_(I22,i);
        const Real &g_33 = g_(I33,i);
        const Real &g00 = gi_(I00,i);
        const Real &g11 = gi_(I11,i);
        const Real &g22 = gi_(I22,i);
        const Real &g33 = gi_(I33,i);
        Real alpha = std::sqrt(-1.0/g00);
        const Real &r = x1v(i);
        Real r2 = SQR(r);
        Real d1_g_00 = -2.0*m / r2;
        Real d1_g_11 = -2.0*m / r2 * SQR(g_11);
        Real d1_g_22 = 2.0 * r;
        Real d1_g_33 = 2.0 * r * sin2;
        Real d2_g_33 = 2.0 * r2 * sincos;

        // Extract primitives
        const Real &rho = prim(IDN,k,j,i);
        const Real &pgas = prim(IEN,k,j,i);
        const Real &uu1 = prim(IVX,k,j,i);
        const Real &uu2 = prim(IVY,k,j,i);
        const Real &uu3 = prim(IVZ,k,j,i);

        // Calculate 4-velocity
        Real uu_sq = g_11*uu1*uu1 + g_22*uu2*uu2 + g_33*uu3*uu3;
        Real gamma = std::sqrt(1.0 + uu_sq);
        Real u0 = gamma / alpha;
        Real u1 = uu1;
        Real u2 = uu2;
        Real u3 = uu3;

        // Extract and calculate magnetic field
        Real b0 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0;
        Real b_sq = 0.0;
        if (MAGNETIC_FIELDS_ENABLED) {
          Real u_1 = g_11*u1;
          Real u_2 = g_22*u2;
          Real u_3 = g_33*u3;
          const Real &bb1 = bb_cc(IB1,k,j,i);
          const Real &bb2 = bb_cc(IB2,k,j,i);
          const Real &bb3 = bb_cc(IB3,k,j,i);
          b0 = u_1*bb1 + u_2*bb2 + u_3*bb3;
          b1 = (bb1 + b0 * u1) / u0;
          b2 = (bb2 + b0 * u2) / u0;
          b3 = (bb3 + b0 * u3) / u0;
          Real b_0 = g_00*b0;
          Real b_1 = g_11*b1;
          Real b_2 = g_22*b2;
          Real b_3 = g_33*b3;
          b_sq = b_0*b0 + b_1*b1 + b_2*b2 + b_3*b3;
        }

        // Calculate stress-energy tensor
        Real wtot = rho + gamma_adi/(gamma_adi-1.0) * pgas + b_sq;
        Real ptot = pgas + 0.5*b_sq;
        Real tt00 = wtot * u0 * u0 + ptot * g00 - b0 * b0;
        Real tt11 = wtot * u1 * u1 + ptot * g11 - b1 * b1;
        Real tt22 = wtot * u2 * u2 + ptot * g22 - b2 * b2;
        Real tt33 = wtot * u3 * u3 + ptot * g33 - b3 * b3;

        // Calculate source terms
        Real s_1 = 0.5 * (d1_g_00*tt00 + d1_g_11*tt11 + d1_g_22*tt22 + d1_g_33*tt33);
        Real s_2 = 0.5 * d2_g_33*tt33;

        // Extract conserved quantities
        Real &m_1 = cons(IM1,k,j,i);
        Real &m_2 = cons(IM2,k,j,i);

        // Add source terms to conserved quantities
        m_1 += dt * s_1;
        m_2 += dt * s_2;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for computing cell-centered metric coefficients
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void Schwarzschild::CellMetric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_cell_j1_(j);

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract remaining geometric quantities
    const Real &alpha_sq = metric_cell_i1_(i);
    const Real &r = x1v(i);
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

//----------------------------------------------------------------------------------------
// Function for computing face-centered metric coefficients: r-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void Schwarzschild::Face1Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face1_j1_(j);

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract remaining geometric quantities
    const Real &alpha_sq = metric_face1_i1_(i);
    const Real &r = x1f(i);
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

//----------------------------------------------------------------------------------------
// Function for computing face-centered metric coefficients: theta-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void Schwarzschild::Face2Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face2_j1_(j);

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract remaining geometric quantities
    const Real &alpha_sq = metric_face2_i1_(i);
    const Real &r = x1v(i);
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

//----------------------------------------------------------------------------------------
// Function for computing face-centered metric coefficients: phi-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void Schwarzschild::Face3Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face3_j1_(j);

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract remaining geometric quantities
    const Real &alpha_sq = metric_face3_i1_(i);
    const Real &r = x1v(i);
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

//----------------------------------------------------------------------------------------
// Function for transforming primitives to locally flat frame: r-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   bb1: 3D array of normal components B^1 of magnetic field, in global coordinates
//   prim_l: 3D array of left primitives, using global coordinates
//   prim_r: 3D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
//   bbx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   expects \tilde{u}^1/\tilde{u}^2/\tilde{u}^3 in IVX/IVY/IVZ slots
//   expects B^1 in bb1
//   expects B^2/B^3 in IBY/IBZ slots
//   puts \tilde{u}^x/\tilde{u}^y/\tilde{u}^z in IVX/IVY/IVZ slots
//   puts B^x in bbx
//   puts B^y/B^z in IBY/IBZ slots
//   u^\hat{i} = M^\hat{i}_j \tilde{u}^j

void Schwarzschild::PrimToLocal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb1, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx) {
  // Calculate metric coefficients
  if (MAGNETIC_FIELDS_ENABLED) {
    Face1Metric(k, j, il, iu, g_, gi_);
  }

  // Extract useful quantities that do not depend on r
  const Real &abs_sin_theta = trans_face1_j1_(j);

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract transformation coefficients
    const Real &r = x1f(i);
    const Real &alpha = trans_face1_i1_(i);
    const Real mt_0 = alpha;
    const Real mx_1 = 1.0/alpha;
    const Real my_2 = r;
    const Real mz_3 = r * abs_sin_theta;

    // Extract global projected 4-velocities
    Real uu1_l = prim_l(IVX,k,j,i);
    Real uu2_l = prim_l(IVY,k,j,i);
    Real uu3_l = prim_l(IVZ,k,j,i);
    Real uu1_r = prim_r(IVX,k,j,i);
    Real uu2_r = prim_r(IVY,k,j,i);
    Real uu3_r = prim_r(IVZ,k,j,i);

    // Transform projected 4-velocities
    Real ux_l = mx_1*uu1_l;
    Real uy_l = my_2*uu2_l;
    Real uz_l = mz_3*uu3_l;
    Real ux_r = mx_1*uu1_r;
    Real uy_r = my_2*uu2_r;
    Real uz_r = mz_3*uu3_r;

    // Set local projected 4-velocities
    prim_l(IVX,k,j,i) = ux_l;
    prim_l(IVY,k,j,i) = uy_l;
    prim_l(IVZ,k,j,i) = uz_l;
    prim_r(IVX,k,j,i) = ux_r;
    prim_r(IVY,k,j,i) = uy_r;
    prim_r(IVZ,k,j,i) = uz_r;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED) {

      // Extract metric coefficients
      const Real &g_00 = g_(I00,i);
      const Real &g_11 = g_(I11,i);
      const Real &g_22 = g_(I22,i);
      const Real &g_33 = g_(I33,i);

      // Calculate global 4-velocities
      Real tmp = g_11*uu1_l*uu1_l + g_22*uu2_l*uu2_l + g_33*uu3_l*uu3_l;
      Real gamma_l = std::sqrt(1.0 + tmp);
      Real u0_l = gamma_l / alpha;
      Real u1_l = uu1_l;
      Real u2_l = uu2_l;
      Real u3_l = uu3_l;
      tmp = g_11*uu1_r*uu1_r + g_22*uu2_r*uu2_r + g_33*uu3_r*uu3_r;
      Real gamma_r = std::sqrt(1.0 + tmp);
      Real u0_r = gamma_r / alpha;
      Real u1_r = uu1_r;
      Real u2_r = uu2_r;
      Real u3_r = uu3_r;

      // Extract global magnetic fields
      const Real &bb1_l = bb1(k,j,i);
      const Real &bb1_r = bb1(k,j,i);
      Real &bb2_l = prim_l(IBY,k,j,i);
      Real &bb3_l = prim_l(IBZ,k,j,i);
      Real &bb2_r = prim_r(IBY,k,j,i);
      Real &bb3_r = prim_r(IBZ,k,j,i);

      // Calculate global 4-magnetic fields
      Real b0_l = g_11*bb1_l*u1_l + g_22*bb2_l*u2_l + g_33*bb3_l*u3_l;
      Real b1_l = (bb1_l + b0_l * u1_l) / u0_l;
      Real b2_l = (bb2_l + b0_l * u2_l) / u0_l;
      Real b3_l = (bb3_l + b0_l * u3_l) / u0_l;
      Real b0_r = g_11*bb1_r*u1_r + g_22*bb2_r*u2_r + g_33*bb3_r*u3_r;
      Real b1_r = (bb1_r + b0_r * u1_r) / u0_r;
      Real b2_r = (bb2_r + b0_r * u2_r) / u0_r;
      Real b3_r = (bb3_r + b0_r * u3_r) / u0_r;

      // Transform 4-velocities
      Real ut_l = gamma_l;
      Real ut_r = gamma_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_1*b1_l;
      Real by_l = my_2*b2_l;
      Real bz_l = mz_3*b3_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_1*b1_r;
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

//----------------------------------------------------------------------------------------
// Function for transforming primitives to locally flat frame: theta-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   bb2: 3D array of normal components B^2 of magnetic field, in global coordinates
//   prim_l: 3D array of left primitives, using global coordinates
//   prim_r: 3D array of right primitives, using global coordinates
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
//   u^\hat{i} = M^\hat{i}_j \tilde{u}^j

void Schwarzschild::PrimToLocal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb2, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx) {
  // Calculate metric coefficients
  if (MAGNETIC_FIELDS_ENABLED) {
    Face2Metric(k, j, il, iu, g_, gi_);
  }

  // Extract useful quantities that do not depend on r
  const Real &abs_sin_theta = trans_face2_j1_(j);

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract transformation coefficients
    const Real &r = x1v(i);
    const Real &alpha = trans_face2_i1_(i);
    const Real mt_0 = alpha;
    const Real mx_2 = 1.0/r;
    const Real my_3 = r * abs_sin_theta;
    const Real &mz_1 = 1.0/alpha;

    // Extract global projected 4-velocities
    Real uu1_l = prim_l(IVX,k,j,i);
    Real uu2_l = prim_l(IVY,k,j,i);
    Real uu3_l = prim_l(IVZ,k,j,i);
    Real uu1_r = prim_r(IVX,k,j,i);
    Real uu2_r = prim_r(IVY,k,j,i);
    Real uu3_r = prim_r(IVZ,k,j,i);

    // Transform projected 4-velocities
    Real ux_l = mx_2*uu2_l;
    Real uy_l = my_3*uu3_l;
    Real uz_l = mz_1*uu1_l;
    Real ux_r = mx_2*uu2_r;
    Real uy_r = my_3*uu3_r;
    Real uz_r = mz_1*uu1_r;

    // Set local projected 4-velocities
    prim_l(IVY,k,j,i) = ux_l;
    prim_l(IVZ,k,j,i) = uy_l;
    prim_l(IVX,k,j,i) = uz_l;
    prim_r(IVY,k,j,i) = ux_r;
    prim_r(IVZ,k,j,i) = uy_r;
    prim_r(IVX,k,j,i) = uz_r;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED) {

      // Extract metric coefficients
      const Real &g_00 = g_(I00,i);
      const Real &g_11 = g_(I11,i);
      const Real &g_22 = g_(I22,i);
      const Real &g_33 = g_(I33,i);

      // Calculate global 4-velocities
      Real tmp = g_11*uu1_l*uu1_l + g_22*uu2_l*uu2_l + g_33*uu3_l*uu3_l;
      Real gamma_l = std::sqrt(1.0 + tmp);
      Real u0_l = gamma_l / alpha;
      Real u1_l = uu1_l;
      Real u2_l = uu2_l;
      Real u3_l = uu3_l;
      tmp = g_11*uu1_r*uu1_r + g_22*uu2_r*uu2_r + g_33*uu3_r*uu3_r;
      Real gamma_r = std::sqrt(1.0 + tmp);
      Real u0_r = gamma_r / alpha;
      Real u1_r = uu1_r;
      Real u2_r = uu2_r;
      Real u3_r = uu3_r;

      // Extract global magnetic fields
      const Real &bb2_l = bb2(k,j,i);
      const Real &bb2_r = bb2(k,j,i);
      Real &bb3_l = prim_l(IBY,k,j,i);
      Real &bb1_l = prim_l(IBZ,k,j,i);
      Real &bb3_r = prim_r(IBY,k,j,i);
      Real &bb1_r = prim_r(IBZ,k,j,i);

      // Calculate global 4-magnetic fields
      Real b0_l = g_11*bb1_l*u1_l + g_22*bb2_l*u2_l + g_33*bb3_l*u3_l;
      Real b1_l = (bb1_l + b0_l * u1_l) / u0_l;
      Real b2_l = (bb2_l + b0_l * u2_l) / u0_l;
      Real b3_l = (bb3_l + b0_l * u3_l) / u0_l;
      Real b0_r = g_11*bb1_r*u1_r + g_22*bb2_r*u2_r + g_33*bb3_r*u3_r;
      Real b1_r = (bb1_r + b0_r * u1_r) / u0_r;
      Real b2_r = (bb2_r + b0_r * u2_r) / u0_r;
      Real b3_r = (bb3_r + b0_r * u3_r) / u0_r;

      // Transform 4-velocities
      Real ut_l = gamma_l;
      Real ut_r = gamma_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_2*b2_l;
      Real by_l = my_3*b3_l;
      Real bz_l = mz_1*b1_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_2*b2_r;
      Real by_r = my_3*b3_r;
      Real bz_r = mz_1*b1_r;

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

//----------------------------------------------------------------------------------------
// Function for transforming primitives to locally flat frame: phi-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   bb3: 3D array of normal components B^3 of magnetic field, in global coordinates
//   prim_l: 3D array of left primitives, using global coordinates
//   prim_r: 3D array of right primitives, using global coordinates
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
//   u^\hat{i} = M^\hat{i}_j \tilde{u}^j

void Schwarzschild::PrimToLocal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb3, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx) {
  // Calculate metric coefficients
  if (MAGNETIC_FIELDS_ENABLED) {
    Face3Metric(k, j, il, iu, g_, gi_);
  }

  // Extract useful quantities that do not depend on r
  const Real &abs_sin_theta = trans_face3_j1_(j);

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract transformation coefficients
    const Real &r = x1v(i);
    const Real &alpha = trans_face3_i1_(i);
    const Real mt_0 = alpha;
    const Real mx_3 = r * abs_sin_theta;
    const Real my_1 = 1.0/alpha;
    const Real mz_2 = r;

    // Extract global projected 4-velocities
    Real uu1_l = prim_l(IVX,k,j,i);
    Real uu2_l = prim_l(IVY,k,j,i);
    Real uu3_l = prim_l(IVZ,k,j,i);
    Real uu1_r = prim_r(IVX,k,j,i);
    Real uu2_r = prim_r(IVY,k,j,i);
    Real uu3_r = prim_r(IVZ,k,j,i);

    // Transform projected 4-velocities
    Real ux_l = mx_3*uu3_l;
    Real uy_l = my_1*uu1_l;
    Real uz_l = mz_2*uu2_l;
    Real ux_r = mx_3*uu3_r;
    Real uy_r = my_1*uu1_r;
    Real uz_r = mz_2*uu2_r;

    // Set local projected 4-velocities
    prim_l(IVZ,k,j,i) = ux_l;
    prim_l(IVX,k,j,i) = uy_l;
    prim_l(IVY,k,j,i) = uz_l;
    prim_r(IVZ,k,j,i) = ux_r;
    prim_r(IVX,k,j,i) = uy_r;
    prim_r(IVY,k,j,i) = uz_r;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED) {

      // Extract metric coefficients
      const Real &g_00 = g_(I00,i);
      const Real &g_11 = g_(I11,i);
      const Real &g_22 = g_(I22,i);
      const Real &g_33 = g_(I33,i);

      // Calculate global 4-velocities
      Real tmp = g_11*uu1_l*uu1_l + g_22*uu2_l*uu2_l + g_33*uu3_l*uu3_l;
      Real gamma_l = std::sqrt(1.0 + tmp);
      Real u0_l = gamma_l / alpha;
      Real u1_l = uu1_l;
      Real u2_l = uu2_l;
      Real u3_l = uu3_l;
      tmp = g_11*uu1_r*uu1_r + g_22*uu2_r*uu2_r + g_33*uu3_r*uu3_r;
      Real gamma_r = std::sqrt(1.0 + tmp);
      Real u0_r = gamma_r / alpha;
      Real u1_r = uu1_r;
      Real u2_r = uu2_r;
      Real u3_r = uu3_r;

      // Extract global magnetic fields
      const Real &bb3_l = bb3(k,j,i);
      const Real &bb3_r = bb3(k,j,i);
      Real &bb1_l = prim_l(IBY,k,j,i);
      Real &bb2_l = prim_l(IBZ,k,j,i);
      Real &bb1_r = prim_r(IBY,k,j,i);
      Real &bb2_r = prim_r(IBZ,k,j,i);

      // Calculate global 4-magnetic fields
      Real b0_l = g_11*bb1_l*u1_l + g_22*bb2_l*u2_l + g_33*bb3_l*u3_l;
      Real b1_l = (bb1_l + b0_l * u1_l) / u0_l;
      Real b2_l = (bb2_l + b0_l * u2_l) / u0_l;
      Real b3_l = (bb3_l + b0_l * u3_l) / u0_l;
      Real b0_r = g_11*bb1_r*u1_r + g_22*bb2_r*u2_r + g_33*bb3_r*u3_r;
      Real b1_r = (bb1_r + b0_r * u1_r) / u0_r;
      Real b2_r = (bb2_r + b0_r * u2_r) / u0_r;
      Real b3_r = (bb3_r + b0_r * u3_r) / u0_r;

      // Transform 4-velocities
      Real ut_l = gamma_l;
      Real ut_r = gamma_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_3*b3_l;
      Real by_l = my_1*b1_l;
      Real bz_l = mz_2*b2_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_3*b3_r;
      Real by_r = my_1*b1_r;
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

//----------------------------------------------------------------------------------------
// Function for transforming fluxes to global frame: r-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   cons: 1D array of conserved quantities, using local coordinates (not used)
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates (not used)
//   flux: 3D array of hydrodynamical fluxes, using local coordinates
//   ey,ez: 3D arrays of magnetic fluxes (electric fields), using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
//   ey,ez: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM1/IM2/IM3 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots and ey/ez
//   puts r-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts r-fluxes of B2/B3 in ey/ez

void Schwarzschild::FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face1_j1_(j);
  const Real &abs_sin_theta = trans_face1_j1_(j);

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract geometric quantities
    const Real &alpha_sq = metric_face1_i1_(i);
    const Real &r = x1f(i);
    const Real r_sq = SQR(r);
    const Real &alpha = trans_face1_i1_(i);
    const Real g00 = -alpha_sq;
    const Real g11 = 1.0/alpha_sq;
    const Real g22 = r_sq;
    const Real g33 = r_sq * sin_sq_theta;
    const Real m0_t = 1.0/alpha;
    const Real m1_x = alpha;
    const Real m2_y = 1.0/r;
    const Real m3_z = 1.0 / (r * abs_sin_theta);

    // Extract local conserved quantities and fluxes
    const Real dx = flux(IDN,k,j,i);
    const Real txt = flux(IEN,k,j,i);
    const Real txx = flux(IM1,k,j,i);
    const Real txy = flux(IM2,k,j,i);
    const Real txz = flux(IM3,k,j,i);

    // Transform stress-energy tensor
    Real t10 = m1_x*m0_t*txt;
    Real t11 = m1_x*m1_x*txx;
    Real t12 = m1_x*m2_y*txy;
    Real t13 = m1_x*m3_z*txz;

    // Extract global fluxes
    Real &d1 = flux(IDN,k,j,i);
    Real &t1_0 = flux(IEN,k,j,i);
    Real &t1_1 = flux(IM1,k,j,i);
    Real &t1_2 = flux(IM2,k,j,i);
    Real &t1_3 = flux(IM3,k,j,i);

    // Set fluxes
    d1 = m1_x*dx;
    t1_0 = g00*t10;
    t1_1 = g11*t11;
    t1_2 = g22*t12;
    t1_3 = g33*t13;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED) {
      Real fyx = -ey(k,j,i);
      Real fzx = ez(k,j,i);
      Real f21 = m2_y*m1_x*fyx;
      Real f31 = m3_z*m1_x*fzx;
      ey(k,j,i) = -f21;
      ez(k,j,i) = f31;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming fluxes to global frame: theta-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   cons: 1D array of conserved quantities, using local coordinates (not used)
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates (not used)
//   flux: 3D array of hydrodynamical fluxes, using local coordinates
//   ey,ez: 3D arrays of magnetic fluxes (electric fields), using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
//   ey,ez: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM2/IM3/IM1 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots and ey/ez
//   puts theta-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts theta-fluxes of B3/B1 in ey/ez

void Schwarzschild::FluxToGlobal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face2_j1_(j);
  const Real &abs_sin_theta = trans_face2_j1_(j);

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract geometric quantities
    const Real &alpha_sq = metric_face2_i1_(i);
    const Real &r = x1v(i);
    const Real r_sq = SQR(r);
    const Real &alpha = trans_face2_i1_(i);
    const Real g00 = -alpha_sq;
    const Real g11 = 1.0/alpha_sq;
    const Real g22 = r_sq;
    const Real g33 = r_sq * sin_sq_theta;
    const Real m0_t = 1.0/alpha;
    const Real m1_z = alpha;
    const Real m2_x = 1.0/r;
    const Real m3_y = 1.0 / (r * abs_sin_theta);

    // Extract local conserved quantities and fluxes
    const Real dx = flux(IDN,k,j,i);
    const Real txt = flux(IEN,k,j,i);
    const Real txx = flux(IM2,k,j,i);
    const Real txy = flux(IM3,k,j,i);
    const Real txz = flux(IM1,k,j,i);

    // Transform stress-energy tensor
    Real t20 = m2_x*m0_t*txt;
    Real t21 = m2_x*m1_z*txz;
    Real t22 = m2_x*m2_x*txx;
    Real t23 = m2_x*m3_y*txy;

    // Extract global fluxes
    Real &d2 = flux(IDN,k,j,i);
    Real &t2_0 = flux(IEN,k,j,i);
    Real &t2_1 = flux(IM1,k,j,i);
    Real &t2_2 = flux(IM2,k,j,i);
    Real &t2_3 = flux(IM3,k,j,i);

    // Set fluxes
    d2 = m2_x*dx;
    t2_0 = g00*t20;
    t2_1 = g11*t21;
    t2_2 = g22*t22;
    t2_3 = g33*t23;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED) {
      Real fyx = -ey(k,j,i);
      Real fzx = ez(k,j,i);
      Real f32 = m3_y*m2_x*fyx;
      Real f12 = m1_z*m2_x*fzx;
      ey(k,j,i) = -f32;
      ez(k,j,i) = f12;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming fluxes to global frame: phi-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   cons: 1D array of conserved quantities, using local coordinates (not used)
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates (not used)
//   flux: 3D array of hydrodynamical fluxes, using local coordinates
//   ey,ez: 3D arrays of magnetic fluxes (electric fields), using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
//   ey,ez: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM3/IM1/IM2 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots and ey/ez
//   puts phi-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts phi-fluxes of B1/B2 in ey/ez

void Schwarzschild::FluxToGlobal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  // Extract geometric quantities that do not depend on r
  const Real &sin_sq_theta = metric_face3_j1_(j);
  const Real &abs_sin_theta = trans_face3_j1_(j);

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract geometric quantities
    const Real &alpha_sq = metric_face3_i1_(i);
    const Real &r = x1v(i);
    const Real r_sq = SQR(r);
    const Real &alpha = trans_face3_i1_(i);
    const Real g00 = -alpha_sq;
    const Real g11 = 1.0/alpha_sq;
    const Real g22 = r_sq;
    const Real g33 = r_sq * sin_sq_theta;
    const Real m0_t = 1.0/alpha;
    const Real m1_y = alpha;
    const Real m2_z = 1.0/r;
    const Real m3_x = 1.0 / (r * abs_sin_theta);

    // Extract local conserved quantities and fluxes
    const Real dx = flux(IDN,k,j,i);
    const Real txt = flux(IEN,k,j,i);
    const Real txx = flux(IM3,k,j,i);
    const Real txy = flux(IM1,k,j,i);
    const Real txz = flux(IM2,k,j,i);

    // Transform stress-energy tensor
    Real t30 = m3_x*m0_t*txt;
    Real t31 = m3_x*m1_y*txy;
    Real t32 = m3_x*m2_z*txz;
    Real t33 = m3_x*m3_x*txx;

    // Extract global fluxes
    Real &d3 = flux(IDN,k,j,i);
    Real &t3_0 = flux(IEN,k,j,i);
    Real &t3_1 = flux(IM1,k,j,i);
    Real &t3_2 = flux(IM2,k,j,i);
    Real &t3_3 = flux(IM3,k,j,i);

    // Set fluxes
    d3 = m3_x*dx;
    t3_0 = g00*t30;
    t3_1 = g11*t31;
    t3_2 = g22*t32;
    t3_3 = g33*t33;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED) {
      Real fyx = -ey(k,j,i);
      Real fzx = ez(k,j,i);
      Real f13 = m1_y*m3_x*fyx;
      Real f23 = m2_z*m3_x*fzx;
      ey(k,j,i) = -f13;
      ez(k,j,i) = f23;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for raising covariant components of a vector
// Inputs:
//   a_0,a_1,a_2,a_3: covariant components of vector
//   k,j,i: indices of cell in which transformation is desired
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to contravariant 4-vector components

void Schwarzschild::RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3) {
  // Extract geometric quantities
  const Real &sin_sq_theta = metric_cell_j1_(j);
  const Real &alpha_sq = metric_cell_i1_(i);
  const Real &r = x1v(i);
  Real r_sq = SQR(r);

  // Calculate metric coefficients
  Real g00 = -1.0/alpha_sq;
  Real g11 = alpha_sq;
  Real g22 = 1.0/r_sq;
  Real g33 = 1.0/(r_sq*sin_sq_theta);

  // Set raised components
  *pa0 = g00 * a_0;
  *pa1 = g11 * a_1;
  *pa2 = g22 * a_2;
  *pa3 = g33 * a_3;
  return;
}

//----------------------------------------------------------------------------------------
// Function for lowering contravariant components of a vector
// Inputs:
//   a0,a1,a2,a3: contravariant components of vector
//   k,j,i: indices of cell in which transformation is desired
// Outputs:
//   pa_0,pa_1,pa_2,pa_3: pointers to covariant 4-vector components

void Schwarzschild::LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j,
    int i, Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3) {
  // Extract geometric quantities
  const Real &sin_sq_theta = metric_cell_j1_(j);
  const Real &alpha_sq = metric_cell_i1_(i);
  const Real &r = x1v(i);
  Real r_sq = SQR(r);

  // Calculate metric coefficients
  Real g_00 = -alpha_sq;
  Real g_11 = 1.0/alpha_sq;
  Real g_22 = r_sq;
  Real g_33 = r_sq * sin_sq_theta;

  // Set lowered components
  *pa_0 = g_00 * a0;
  *pa_1 = g_11 * a1;
  *pa_2 = g_22 * a2;
  *pa_3 = g_33 * a3;
  return;
}
