//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file minkowski_cyl.cpp
//  \brief implements functions for Minkowski spacetime and cylindrical (t,R,phi,z)
//  coordinates in a derived class of the Coordinates abstract base class.
//
// Notes:
//   coordinates: t, R, phi, z
//   metric:
//     ds^2 = -dt^2 + dR^2 + R^2 d\phi^2 + dz^2

// C headers

// C++ headers
#include <cmath>      // cos, hypot, sin, sqrt
#include <cstdlib>    // exit
#include <iostream>   // cout
#include <ostream>    // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "coordinates.hpp"

//----------------------------------------------------------------------------------------
// Cylindrical Minkowski Constructor
// Inputs:
//   pmb: pointer to MeshBlock containing this grid
//   pin: pointer to runtime inputs
//   flag: true if object is for coarse grid only in an AMR calculation

MinkowskiCyl::MinkowskiCyl(MeshBlock *pmb, ParameterInput *pin, bool flag)
    : Coordinates(pmb, pin, flag) {

  // Set parameters
  rad_tetrad_ = pin->GetOrAddString("coord", "rad_tetrad", "none");

  // Initialize volume-averaged coordinates and spacings: R-direction
  for (int i = il-ng; i <= iu+ng; ++i) {
    Real rr_m = x1f(i);
    Real rr_p = x1f(i+1);
    x1v(i) = 2.0/3.0 * (rr_p*SQR(rr_p) - rr_m*SQR(rr_m)) / (SQR(rr_p) - SQR(rr_m));
  }
  for (int i = il-ng; i <= iu+ng-1; ++i) {
    dx1v(i) = x1v(i+1) - x1v(i);
  }

  // Initialize volume-averaged coordinates and spacings: phi-direction
  if (pmb->block_size.nx2 == 1) {
    Real phi_m = x2f(jl);
    Real phi_p = x2f(jl+1);
    x2v(jl) = 0.5 * (phi_m + phi_p);
    dx2v(jl) = dx2f(jl);
  } else {
    for (int j = jl-ng; j <= ju+ng; ++j) {
      Real phi_m = x2f(j);
      Real phi_p = x2f(j+1);
      x2v(j) = 0.5 * (phi_m + phi_p);
    }
    for (int j = jl-ng; j <= ju+ng-1; ++j) {
      dx2v(j) = x2v(j+1) - x2v(j);
    }
  }

  // Initialize volume-averaged coordinates and spacings: z-direction
  if (pmb->block_size.nx3 == 1) {
    Real z_m = x3f(kl);
    Real z_p = x3f(kl+1);
    x3v(kl) = 0.5 * (z_m + z_p);
    dx3v(kl) = dx3f(kl);
  } else {
    for (int k = kl-ng; k <= ku+ng; ++k) {
      Real z_m = x3f(k);
      Real z_p = x3f(k+1);
      x3v(k) = 0.5 * (z_m + z_p);
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

  // Allocate and compute arrays for intermediate geometric quantities that are only
  // needed if object is NOT a coarse mesh
  if (!coarse_flag) {

    // Allocate arrays for intermediate geometric quantities: r-direction
    coord_vol_i1_.NewAthenaArray(nc1);
    coord_area1_i1_.NewAthenaArray(nc1+1);
    coord_area2_i1_.NewAthenaArray(nc1);
    coord_area3_i1_.NewAthenaArray(nc1);
    coord_len1_i1_.NewAthenaArray(nc1);
    coord_len2_i1_.NewAthenaArray(nc1+1);

    // Calculate intermediate geometric quantities: r-direction
    for (int i = il-ng; i <= iu+ng; ++i) {

      // Useful quantities
      Real rr_m = x1f(i);
      Real rr_p = x1f(i+1);

      // Volumes, areas, lengths, and widths
      coord_vol_i1_(i) = 0.5 * (SQR(rr_p) - SQR(rr_m));
      coord_area1_i1_(i) = rr_m;
      if (i == iu+ng) {
        coord_area1_i1_(i+1) = rr_p;
      }
      coord_area2_i1_(i) = coord_vol_i1_(i);
      coord_area3_i1_(i) = coord_vol_i1_(i);
      coord_len1_i1_(i) = coord_vol_i1_(i);
      coord_len2_i1_(i) = coord_area1_i1_(i);
      if (i == iu+ng) {
        coord_len2_i1_(i+1) = coord_area1_i1_(i+1);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
// EdgeXLength functions: compute physical length at cell edge-X as vector
// Edge1(i,j,k) located at (i,j-1/2,k-1/2), i.e. (x1v(i), x2f(j), x3f(k))
// Edge2(i,j,k) located at (i-1/2,j,k-1/2), i.e. (x1f(i), x2v(j), x3f(k))
// Edge3(i,j,k) located at (i-1/2,j-1/2,k), i.e. (x1f(i), x2f(j), x3v(k))

// \Delta L = 1/2 \Delta(R^2)
void MinkowskiCyl::Edge1Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &lengths) {
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {
    lengths(i) = coord_len1_i1_(i);
  }
  return;
}

// \Delta L = R \Delta\phi
void MinkowskiCyl::Edge2Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &lengths) {
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {
    lengths(i) = coord_len2_i1_(i) * dx2f(j);
  }
  return;
}

// \Delta L = \Delta z
void MinkowskiCyl::Edge3Length(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &lengths) {
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {
    lengths(i) = dx3f(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetEdgeXLength functions: return length of edge-X at (i,j,k)

// \Delta L = 1/2 \Delta(R^2)
Real MinkowskiCyl::GetEdge1Length(const int k, const int j, const int i) {
  return coord_len1_i1_(i);
}

// \Delta L = R \Delta\phi
Real MinkowskiCyl::GetEdge2Length(const int k, const int j, const int i) {
  return coord_len2_i1_(i) * dx2f(j);
}

// \Delta L = \Delta z
Real MinkowskiCyl::GetEdge3Length(const int k, const int j, const int i) {
  return dx3f(k);
}

//----------------------------------------------------------------------------------------
// CenterWidthX functions: return physical width in X-dir at (i,j,k) cell-center

// \Delta W = R \Delta\phi
void MinkowskiCyl::CenterWidth2(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &dx2) {
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {
    dx2(i) = x1v(i) * dx2f(j);
  }
  return;
}

//----------------------------------------------------------------------------------------
// FaceXArea functions: compute area of face with normal in X-dir as vector
// Inputs:
//   k,j: z- and phi-indices
//   il,iu: R-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to X-face

// \Delta A = R \Delta\phi \Delta z
void MinkowskiCyl::Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas) {
  #pragma omp simd
  for (int i = il; i <= iu; ++i)
    areas(i) = coord_area1_i1_(i) * dx2f(j) * dx3f(k);
  return;
}

// \Delta A = 1/2 \Delta(R^2) \Delta z
void MinkowskiCyl::Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas) {
  #pragma omp simd
  for (int i = il; i <= iu; ++i)
    areas(i) = coord_area2_i1_(i) * dx3f(k);
  return;
}

// \Delta A = 1/2 \Delta(R^2) \Delta\phi
void MinkowskiCyl::Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas) {
  #pragma omp simd
  for (int i = il; i <= iu; ++i)
    areas(i) = coord_area3_i1_(i) * dx2f(j);
  return;
}

//----------------------------------------------------------------------------------------
// GetFaceXArea functions: return area of face with normal in X-dir at (i,j,k)
// Inputs:
//   k,j,i: z-, phi-, and R-indices
// return:
//   interface area orthogonal to X-face

//  \Delta A = R \Delta\phi \Delta z
Real MinkowskiCyl::GetFace1Area(const int k, const int j, const int i) {
  return coord_area1_i1_(i) * dx2f(j) * dx3f(k);
}

// \Delta A = 1/2 \Delta(R^2) \Delta z
Real MinkowskiCyl::GetFace2Area(const int k, const int j, const int i) {
  return coord_area2_i1_(i) * dx3f(k);
}

// \Delta A = 1/2 \Delta(R^2) \Delta\phi
Real MinkowskiCyl::GetFace3Area(const int k, const int j, const int i) {
  return coord_area3_i1_(i) * dx2f(j);
}

//----------------------------------------------------------------------------------------
// Cell Volume function: compute volume of cell as vector
// Inputs:
//   k,j: z- and phi-indices
//   il,iu: R-index bounds
// Outputs:
//   volumes: 1D array of cell volumes
// Notes:
//   \Delta V = 1/2 * \Delta(R^2) \Delta\phi \Delta z

void MinkowskiCyl::CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &volumes) {
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {
    volumes(i) = coord_vol_i1_(i) * dx2f(j) * dx3f(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetCellVolume: returns cell volume at (i,j,k)
// Inputs:
//   k,j,i: z-, phi-, and R-indices
// Outputs:
//   returned value: cell volume
// Notes:
//   \Delta V = 1/2 * \Delta(R^2) \Delta\phi \Delta z

Real MinkowskiCyl::GetCellVolume(const int k, const int j, const int i) {
  return coord_vol_i1_(i) * dx2f(j) * dx3f(k);
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

void MinkowskiCyl::AddCoordTermsDivergence(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bb_cc,
    AthenaArray<Real> &cons) {

  // Extract ratio of specific heats
  const Real gamma_adi = pmy_block->peos->GetGamma();

  // Go through cells
  for (int k = pmy_block->ks; k <= pmy_block->ke; ++k) {
    for (int j = pmy_block->js; j <= pmy_block->je; ++j) {

      // Go through 1D slice
      #pragma omp simd
      for (int i = pmy_block->is; i <= pmy_block->ie; ++i) {

        // Extract geometric quantities
        Real rr = x1v(i);
        Real rr2 = SQR(rr);
        Real g_00 = -1.0;
        Real g_11 = 1.0;
        Real g_22 = rr2;
        Real g_33 = 1.0;
        Real g22 = 1.0 / rr2;
        Real d1_g_22 = 2.0 * rr;

        // Extract primitives
        const Real &rho = prim(IDN,k,j,i);
        const Real &pgas = prim(IEN,k,j,i);
        const Real &uu1 = prim(IVX,k,j,i);
        const Real &uu2 = prim(IVY,k,j,i);
        const Real &uu3 = prim(IVZ,k,j,i);

        // Calculate 4-velocity
        Real uu_sq = g_11*uu1*uu1 + g_22*uu2*uu2 + g_33*uu3*uu3;
        Real gamma = std::sqrt(1.0 + uu_sq);
        Real u0 = gamma;
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
        Real tt22 = wtot * u2 * u2 + ptot * g22 - b2 * b2;

        // Calculate source terms
        Real s_1 = 0.5 * d1_g_22*tt22;

        // Extract conserved quantities
        Real &m_1 = cons(IM1,k,j,i);

        // Add source terms to conserved quantities
        m_1 += dt * s_1;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for computing cell-centered metric coefficients
// Inputs:
//   k,j: z- and phi-indices
//   il,iu: R-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void MinkowskiCyl::CellMetric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract geometric quantities
    const Real &rr = x1v(i);
    Real rr_sq = SQR(rr);

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
    g00 = -1.0;
    g11 = 1.0;
    g22 = rr_sq;
    g33 = 1.0;
    gi00 = -1.0;
    gi11 = 1.0;
    gi22 = 1.0/rr_sq;
    gi33 = 1.0;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for computing face-centered metric coefficients: R-interface
// Inputs:
//   k,j: z- and phi-indices
//   il,iu: R-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void MinkowskiCyl::Face1Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract geometric quantities
    const Real &rr = x1f(i);
    Real rr_sq = SQR(rr);

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
    g00 = -1.0;
    g11 = 1.0;
    g22 = rr_sq;
    g33 = 1.0;
    gi00 = -1.0;
    gi11 = 1.0;
    gi22 = 1.0/rr_sq;
    gi33 = 1.0;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for computing face-centered metric coefficients: phi-interface
// Inputs:
//   k,j: z- and phi-indices
//   il,iu: R-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void MinkowskiCyl::Face2Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract geometric quantities
    const Real &rr = x1v(i);
    Real rr_sq = SQR(rr);

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
    g00 = -1.0;
    g11 = 1.0;
    g22 = rr_sq;
    g33 = 1.0;
    gi00 = -1.0;
    gi11 = 1.0;
    gi22 = 1.0/rr_sq;
    gi33 = 1.0;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for computing face-centered metric coefficients: z-interface
// Inputs:
//   k,j: z- and phi-indices
//   il,iu: R-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void MinkowskiCyl::Face3Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract geometric quantities
    const Real &rr = x1v(i);
    Real rr_sq = SQR(rr);

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
    g00 = -1.0;
    g11 = 1.0;
    g22 = rr_sq;
    g33 = 1.0;
    gi00 = -1.0;
    gi11 = 1.0;
    gi22 = 1.0/rr_sq;
    gi33 = 1.0;
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming primitives to locally flat frame: R-interface
// Inputs:
//   k,j: z- and phi-indices
//   il,iu: R-index bounds
//   bb1: 1D array of normal components B^1 of magnetic field, in global coordinates
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
//   puts \tilde{u}^x/\tilde{u}^y/\tilde{u}^z in IVX/IVY/IVZ slots
//   puts B^x in bbx
//   puts B^y/B^z in IBY/IBZ slots
//   u^\hat{i} = M^\hat{i}_j \tilde{u}^j

void MinkowskiCyl::PrimToLocal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb1, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx) {

  // Calculate metric coefficients
  if (MAGNETIC_FIELDS_ENABLED) {
    Face1Metric(k, j, il, iu, g_, gi_);
  }

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract transformation coefficients
    const Real &rr = x1f(i);
    const Real mt_0 = 1.0;
    const Real mx_1 = 1.0;
    const Real my_2 = rr;
    const Real mz_3 = 1.0;

    // Extract global projected 4-velocities
    Real uu1_l = prim_l(IVX,i);
    Real uu2_l = prim_l(IVY,i);
    Real uu3_l = prim_l(IVZ,i);
    Real uu1_r = prim_r(IVX,i);
    Real uu2_r = prim_r(IVY,i);
    Real uu3_r = prim_r(IVZ,i);

    // Transform projected 4-velocities
    Real ux_l = mx_1*uu1_l;
    Real uy_l = my_2*uu2_l;
    Real uz_l = mz_3*uu3_l;
    Real ux_r = mx_1*uu1_r;
    Real uy_r = my_2*uu2_r;
    Real uz_r = mz_3*uu3_r;

    // Set local projected 4-velocities
    prim_l(IVX,i) = ux_l;
    prim_l(IVY,i) = uy_l;
    prim_l(IVZ,i) = uz_l;
    prim_r(IVX,i) = ux_r;
    prim_r(IVY,i) = uy_r;
    prim_r(IVZ,i) = uz_r;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED) {

      // Extract metric coefficients
      const Real &g_11 = g_(I11,i);
      const Real &g_22 = g_(I22,i);
      const Real &g_33 = g_(I33,i);

      // Calculate global 4-velocities
      Real tmp = g_11*uu1_l*uu1_l + g_22*uu2_l*uu2_l + g_33*uu3_l*uu3_l;
      Real gamma_l = std::sqrt(1.0 + tmp);
      Real u0_l = gamma_l;
      Real u1_l = uu1_l;
      Real u2_l = uu2_l;
      Real u3_l = uu3_l;
      tmp = g_11*uu1_r*uu1_r + g_22*uu2_r*uu2_r + g_33*uu3_r*uu3_r;
      Real gamma_r = std::sqrt(1.0 + tmp);
      Real u0_r = gamma_r;
      Real u1_r = uu1_r;
      Real u2_r = uu2_r;
      Real u3_r = uu3_r;

      // Extract global magnetic fields
      const Real &bb1_l = bb1(i);
      const Real &bb1_r = bb1(i);
      Real &bb2_l = prim_l(IBY,i);
      Real &bb3_l = prim_l(IBZ,i);
      Real &bb2_r = prim_r(IBY,i);
      Real &bb3_r = prim_r(IBZ,i);

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
// Function for transforming primitives to locally flat frame: phi-interface
// Inputs:
//   k,j: z- and phi-indices
//   il,iu: R-index bounds
//   bb2: 1D array of normal components B^2 of magnetic field, in global coordinates
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
//   u^\hat{i} = M^\hat{i}_j \tilde{u}^j

void MinkowskiCyl::PrimToLocal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb2, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx) {

  // Calculate metric coefficients
  if (MAGNETIC_FIELDS_ENABLED) {
    Face2Metric(k, j, il, iu, g_, gi_);
  }

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract transformation coefficients
    const Real &rr = x1v(i);
    const Real mt_0 = 1.0;
    const Real mx_2 = rr;
    const Real my_3 = 1.0;
    const Real &mz_1 = 1.0;

    // Extract global projected 4-velocities
    Real uu1_l = prim_l(IVX,i);
    Real uu2_l = prim_l(IVY,i);
    Real uu3_l = prim_l(IVZ,i);
    Real uu1_r = prim_r(IVX,i);
    Real uu2_r = prim_r(IVY,i);
    Real uu3_r = prim_r(IVZ,i);

    // Transform projected 4-velocities
    Real ux_l = mx_2*uu2_l;
    Real uy_l = my_3*uu3_l;
    Real uz_l = mz_1*uu1_l;
    Real ux_r = mx_2*uu2_r;
    Real uy_r = my_3*uu3_r;
    Real uz_r = mz_1*uu1_r;

    // Set local projected 4-velocities
    prim_l(IVY,i) = ux_l;
    prim_l(IVZ,i) = uy_l;
    prim_l(IVX,i) = uz_l;
    prim_r(IVY,i) = ux_r;
    prim_r(IVZ,i) = uy_r;
    prim_r(IVX,i) = uz_r;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED) {

      // Extract metric coefficients
      const Real &g_11 = g_(I11,i);
      const Real &g_22 = g_(I22,i);
      const Real &g_33 = g_(I33,i);

      // Calculate global 4-velocities
      Real tmp = g_11*uu1_l*uu1_l + g_22*uu2_l*uu2_l + g_33*uu3_l*uu3_l;
      Real gamma_l = std::sqrt(1.0 + tmp);
      Real u0_l = gamma_l;
      Real u1_l = uu1_l;
      Real u2_l = uu2_l;
      Real u3_l = uu3_l;
      tmp = g_11*uu1_r*uu1_r + g_22*uu2_r*uu2_r + g_33*uu3_r*uu3_r;
      Real gamma_r = std::sqrt(1.0 + tmp);
      Real u0_r = gamma_r;
      Real u1_r = uu1_r;
      Real u2_r = uu2_r;
      Real u3_r = uu3_r;

      // Extract global magnetic fields
      const Real &bb2_l = bb2(i);
      const Real &bb2_r = bb2(i);
      Real &bb3_l = prim_l(IBY,i);
      Real &bb1_l = prim_l(IBZ,i);
      Real &bb3_r = prim_r(IBY,i);
      Real &bb1_r = prim_r(IBZ,i);

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
// Function for transforming primitives to locally flat frame: z-interface
// Inputs:
//   k,j: z- and phi-indices
//   il,iu: R-index bounds
//   bb3: 1D array of normal components B^3 of magnetic field, in global coordinates
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
//   u^\hat{i} = M^\hat{i}_j \tilde{u}^j

void MinkowskiCyl::PrimToLocal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb3, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx) {

  // Calculate metric coefficients
  if (MAGNETIC_FIELDS_ENABLED) {
    Face3Metric(k, j, il, iu, g_, gi_);
  }

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract transformation coefficients
    const Real &rr = x1v(i);
    const Real mt_0 = 1.0;
    const Real mx_3 = 1.0;
    const Real my_1 = 1.0;
    const Real mz_2 = rr;

    // Extract global projected 4-velocities
    Real uu1_l = prim_l(IVX,i);
    Real uu2_l = prim_l(IVY,i);
    Real uu3_l = prim_l(IVZ,i);
    Real uu1_r = prim_r(IVX,i);
    Real uu2_r = prim_r(IVY,i);
    Real uu3_r = prim_r(IVZ,i);

    // Transform projected 4-velocities
    Real ux_l = mx_3*uu3_l;
    Real uy_l = my_1*uu1_l;
    Real uz_l = mz_2*uu2_l;
    Real ux_r = mx_3*uu3_r;
    Real uy_r = my_1*uu1_r;
    Real uz_r = mz_2*uu2_r;

    // Set local projected 4-velocities
    prim_l(IVZ,i) = ux_l;
    prim_l(IVX,i) = uy_l;
    prim_l(IVY,i) = uz_l;
    prim_r(IVZ,i) = ux_r;
    prim_r(IVX,i) = uy_r;
    prim_r(IVY,i) = uz_r;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED) {

      // Extract metric coefficients
      const Real &g_11 = g_(I11,i);
      const Real &g_22 = g_(I22,i);
      const Real &g_33 = g_(I33,i);

      // Calculate global 4-velocities
      Real tmp = g_11*uu1_l*uu1_l + g_22*uu2_l*uu2_l + g_33*uu3_l*uu3_l;
      Real gamma_l = std::sqrt(1.0 + tmp);
      Real u0_l = gamma_l;
      Real u1_l = uu1_l;
      Real u2_l = uu2_l;
      Real u3_l = uu3_l;
      tmp = g_11*uu1_r*uu1_r + g_22*uu2_r*uu2_r + g_33*uu3_r*uu3_r;
      Real gamma_r = std::sqrt(1.0 + tmp);
      Real u0_r = gamma_r;
      Real u1_r = uu1_r;
      Real u2_r = uu2_r;
      Real u3_r = uu3_r;

      // Extract global magnetic fields
      const Real &bb3_l = bb3(i);
      const Real &bb3_r = bb3(i);
      Real &bb1_l = prim_l(IBY,i);
      Real &bb2_l = prim_l(IBZ,i);
      Real &bb1_r = prim_r(IBY,i);
      Real &bb2_r = prim_r(IBZ,i);

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
// Function for transforming fluxes to global frame: R-interface
// Inputs:
//   k,j: z- and phi-indices
//   il,iu: R-index bounds
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

void MinkowskiCyl::FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract geometric quantities
    const Real &rr = x1f(i);
    const Real rr_sq = SQR(rr);
    const Real g00 = -1.0;
    const Real g11 = 1.0;
    const Real g22 = rr_sq;
    const Real g33 = 1.0;
    const Real m0_tm = 1.0;
    const Real m1_x = 1.0;
    const Real m2_y = 1.0/rr;
    const Real m3_z = 1.0;

    // Extract local conserved quantities and fluxes
    const Real dx = flux(IDN,k,j,i);
    const Real txt = flux(IEN,k,j,i);
    const Real txx = flux(IM1,k,j,i);
    const Real txy = flux(IM2,k,j,i);
    const Real txz = flux(IM3,k,j,i);

    // Transform stress-energy tensor
    Real t10 = m1_x*m0_tm*txt;
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
// Function for transforming fluxes to global frame: phi-interface
// Inputs:
//   k,j: z- and phi-indices
//   il,iu: R-index bounds
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

void MinkowskiCyl::FluxToGlobal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract geometric quantities
    const Real &rr = x1v(i);
    const Real rr_sq = SQR(rr);
    const Real g00 = -1.0;
    const Real g11 = 1.0;
    const Real g22 = rr_sq;
    const Real g33 = 1.0;
    const Real m0_tm = 1.0;
    const Real m1_z = 1.0;
    const Real m2_x = 1.0/rr;
    const Real m3_y = 1.0;

    // Extract local conserved quantities and fluxes
    const Real dx = flux(IDN,k,j,i);
    const Real txt = flux(IEN,k,j,i);
    const Real txx = flux(IM2,k,j,i);
    const Real txy = flux(IM3,k,j,i);
    const Real txz = flux(IM1,k,j,i);

    // Transform stress-energy tensor
    Real t20 = m2_x*m0_tm*txt;
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
// Function for transforming fluxes to global frame: z-interface
// Inputs:
//   k,j: z- and phi-indices
//   il,iu: R-index bounds
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

void MinkowskiCyl::FluxToGlobal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {

  // Go through 1D block of cells
  #pragma omp simd
  for (int i = il; i <= iu; ++i) {

    // Extract geometric quantities
    const Real &rr = x1v(i);
    const Real rr_sq = SQR(rr);
    const Real g00 = -1.0;
    const Real g11 = 1.0;
    const Real g22 = rr_sq;
    const Real g33 = 1.0;
    const Real m0_tm = 1.0;
    const Real m1_y = 1.0;
    const Real m2_z = 1.0/rr;
    const Real m3_x = 1.0;

    // Extract local conserved quantities and fluxes
    const Real dx = flux(IDN,k,j,i);
    const Real txt = flux(IEN,k,j,i);
    const Real txx = flux(IM3,k,j,i);
    const Real txy = flux(IM1,k,j,i);
    const Real txz = flux(IM2,k,j,i);

    // Transform stress-energy tensor
    Real t30 = m3_x*m0_tm*txt;
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

void MinkowskiCyl::RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3) {

  // Extract geometric quantities
  const Real &rr = x1v(i);
  Real rr_sq = SQR(rr);

  // Calculate metric coefficients
  Real g00 = -1.0;
  Real g11 = 1.0;
  Real g22 = 1.0/rr_sq;
  Real g33 = 1.0;

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

void MinkowskiCyl::LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j,
    int i, Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3) {

  // Extract geometric quantities
  const Real &rr = x1v(i);
  Real rr_sq = SQR(rr);

  // Calculate metric coefficients
  Real g_00 = -1.0;
  Real g_11 = 1.0;
  Real g_22 = rr_sq;
  Real g_33 = 1.0;

  // Set lowered components
  *pa_0 = g_00 * a0;
  *pa_1 = g_11 * a1;
  *pa_2 = g_22 * a2;
  *pa_3 = g_33 * a3;
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating orthonormal tetrad
// Inputs:
//   rr, ph, z: spatial position
// Outputs:
//   e: 2D array for e_{(\hat{\mu})}^\nu:
//     index 0: covariant orthonormal index
//     index 1: contravariant coordinate index
//   e_cov: 2D array for {e_{(\hat{\mu})}}_\nu:
//     index 0: covariant orthonormal index
//     index 1: covariant coordinate index
//   omega: 3D array for \omega^{\hat{\gamma}}_{\hat{\alpha}\hat{\beta}}:
//     index 0: upper index
//     index 1: first lower index
//     index 2: second lower index
// Notes:
//   tetrad options:
//     "cartesian" (Gram-Schmidt on t, z, x, y)
//     "cylindrical" (Gram-Schmidt on t, z, R, phi)
//     "spherical" (Gram-Schmidt on t, theta, phi, r)

void MinkowskiCyl::Tetrad(Real rr, Real ph, Real z, AthenaArray<Real> &e,
    AthenaArray<Real> &e_cov, AthenaArray<Real> &omega) {

  // Check tetrad
  if (rad_tetrad_ != "cartesian" and rad_tetrad_ != "cylindrical"
      and rad_tetrad_ != "spherical") {
    std::stringstream msg;
    msg << "### FATAL ERROR invalid tetrad choice" << std::endl;
    ATHENA_ERROR(msg);
  }

  // Calculate useful quantities
  Real z2 = SQR(z);
  Real rr2 = SQR(rr);
  Real r = std::hypot(rr, z);
  Real r2 = SQR(r);
  Real r3 = r * r2;
  Real sph = std::sin(ph);
  Real cph = std::cos(ph);

  // Allocate intermediate arrays
  Real eta[4][4] = {};
  Real g[4][4] = {};
  Real gi[4][4] = {};
  Real dg[4][4][4] = {};
  Real ei[4][4] = {};
  Real de[4][4][4] = {};
  Real gamma[4][4][4] = {};

  // Set Minkowski metric
  eta[0][0] = -1.0;
  eta[1][1] = 1.0;
  eta[2][2] = 1.0;
  eta[3][3] = 1.0;

  // Set covariant metric
  g[0][0] = -1.0;
  g[1][1] = 1.0;
  g[2][2] = rr2;
  g[3][3] = 1.0;

  // Set contravariant metric
  gi[0][0] = -1.0;
  gi[1][1] = 1.0;
  gi[2][2] = 1.0 / rr2;
  gi[3][3] = 1.0;

  // Set derivatives of covariant metric
  dg[1][2][2] = 2.0 * r;

  // Set tetrad
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      e(i,j) = 0.0;
    }
  }
  if (rad_tetrad_ == "cartesian") {
    e(0,0) = 1.0;
    e(1,1) = cph;
    e(1,2) = -sph / rr;
    e(2,1) = sph;
    e(2,2) = cph / rr;
    e(3,3) = 1.0;
  } else if (rad_tetrad_ == "cylindrical") {
    e(0,0) = 1.0;
    e(1,1) = 1.0;
    e(2,2) = 1.0 / rr;
    e(3,3) = 1.0;
  } else if (rad_tetrad_ == "spherical") {
    e(0,0) = 1.0;
    e(1,2) = 1.0 / rr;
    e(2,1) = rr / r;
    e(2,3) = z / r;
    e(3,1) = z / r;
    e(3,3) = -rr / r;
  }

  // Calculate covariant tetrad
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      e_cov(i,j) = 0.0;
      for (int k = 0; k < 4; ++k) {
        e_cov(i,j) += g[j][k] * e(i,k)
      }
    }
  }

  // Calculate inverse of tetrad
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        for (int l = 0; l < 4; ++l) {
          ei[i][j] += eta[i][k] * g[j][l] * e(k,l);
        }
      }
    }
  }

  // Set derivatives of tetrad
  if (rad_tetrad_ == "cartesian") {
    de[1][1][2] = sph / rr2;
    de[1][2][2] = -cph / rr2;
    de[2][1][1] = -sph;
    de[2][1][2] = -cph / rr;
    de[2][2][1] = cph;
    de[2][2][2] = -sph / rr;
  } else if (rad_tetrad_ == "cylindrical") {
    de[1][2][2] = -1.0 / rr2;
  } else if (rad_tetrad_ == "spherical") {
    de[1][1][2] = -1.0 / rr2;
    de[1][2][1] = z2 / r3;
    de[1][2][3] = -rr * z / r3;
    de[1][3][1] = -rr * z / r3;
    de[1][3][3] = -z2 / r3;
    de[3][2][1] = -rr * z / r3;
    de[3][2][3] = rr2 / r3;
    de[3][3][1] = rr2 / r3;
    de[3][3][3] = rr * z / r3;
  }

  // Calculate Christoffel connection coefficients
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        for (int l = 0; l < 4; ++l) {
          gamma[i][j][k] += 0.5 * gi[i][l] * (dg[j][k][l] + dg[k][j][l] - dg[l][j][k]);
        }
      }
    }
  }

  // Calculate Ricci rotation coefficients
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        omega(i,j,k) = 0.0;
        for (int l = 0; l < 4; ++l) {
          for (int m = 0; m < 4; ++m) {
            omega(i,j,k) += ei[i][l] * e(k,m) * de[m][j][l];
            for (int n = 0; n < 4; ++n) {
              omega(i,j,k) += ei[i][l] * e(k,m) * gamma[l][m][n] * e(j,n);
            }
          }
        }
      }
    }
  }
  return;
}
