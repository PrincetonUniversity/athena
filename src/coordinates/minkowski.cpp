//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file minkowski.cpp
//  \brief implements functions for Minkowski (flat) spacetime and Cartesian (t,x,y,z)
//  coordinates in a derived class of the Coordinates abstract base class.
//
// Notes:
//   coordinates: t, x, y, z
//   metric: ds^2 = -dt^2 + dx^2 + dy^2 + dz^2

// C headers

// C++ headers
#include <cmath>      // hypot, sqrt
#include <cstdlib>    // exit
#include <iostream>   // cout
#include <ostream>    // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "coordinates.hpp"

//----------------------------------------------------------------------------------------
// Minkowski constructor
// Inputs:
//   pmb: pointer to block containing this grid
//   pin: pointer to runtime inputs (not used)
//   flag: true if object is for coarse grid only in an AMR calculation

Minkowski::Minkowski(MeshBlock *pmb, ParameterInput *pin, bool flag)
    : Coordinates(pmb, pin, flag) {

  // Set parameters
  rad_tetrad_ = pin->GetOrAddString("coord", "rad_tetrad", "none");

  // Initialize volume-averaged coordinates and spacings: x-direction
  for (int i=il-ng; i<=iu+ng; ++i) {
    x1v(i) = 0.5 * (x1f(i) + x1f(i+1));
  }
  for (int i=il-ng; i<=iu+ng-1; ++i) {
    dx1v(i) = x1v(i+1) - x1v(i);
  }

  // Initialize volume-averaged coordinates and spacings: y-direction
  if (pmb->block_size.nx2 == 1) {
    x2v(jl) = 0.5 * (x2f(jl) + x2f(jl+1));
    dx2v(jl) = dx2f(jl);
  } else {
    for (int j=jl-ng; j<=ju+ng; ++j) {
      x2v(j) = 0.5 * (x2f(j) + x2f(j+1));
    }
    for (int j=jl-ng; j<=ju+ng-1; ++j) {
      dx2v(j) = x2v(j+1) - x2v(j);
    }
  }

  // Initialize volume-averaged coordinates and spacings: z-direction
  if (pmb->block_size.nx3 == 1) {
    x3v(kl) = 0.5 * (x3f(kl) + x3f(kl+1));
    dx3v(kl) = dx3f(kl);
  } else {
    for (int k=kl-ng; k<=ku+ng; ++k) {
      x3v(k) = 0.5 * (x3f(k) + x3f(k+1));
    }
    for (int k=kl-ng; k<=ku+ng-1; ++k) {
      dx3v(k) = x3v(k+1) - x3v(k);
    }
  }

  // Initialize area-averaged coordinates used with MHD AMR
  if (pmb->pmy_mesh->multilevel && MAGNETIC_FIELDS_ENABLED) {
    for (int i=il-ng; i<=iu+ng; ++i) {
      x1s2(i) = x1s3(i) = x1v(i);
    }
    if (pmb->block_size.nx2 == 1) {
      x2s1(jl) = x2s3(jl) = x2v(jl);
    } else {
      for (int j=jl-ng; j<=ju+ng; ++j) {
        x2s1(j) = x2s3(j) = x2v(j);
      }
    }
    if (pmb->block_size.nx3 == 1) {
      x3s1(kl) = x3s2(kl) = x3v(kl);
    } else {
      for (int k=kl-ng; k<=ku+ng; ++k) {
        x3s1(k) = x3s2(k) = x3v(k);
      }
    }
  }
}


//----------------------------------------------------------------------------------------
// Function for computing cell-centered metric coefficients
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void Minkowski::CellMetric(const int k, const int j, const int il, const int iu,
                           AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
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

//----------------------------------------------------------------------------------------
// Functions for computing face-centered metric coefficients
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void Minkowski::Face1Metric(const int k, const int j, const int il, const int iu,
                            AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
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

void Minkowski::Face2Metric(const int k, const int j, const int il, const int iu,
                            AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
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

void Minkowski::Face3Metric(const int k, const int j, const int il, const int iu,
                            AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
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

//----------------------------------------------------------------------------------------
// Functions for transforming face-centered primitives to locally flat frame
// Inputs
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   bb1: 1D array of normal components B^1 of magnetic field, in global coordinates
//   prim_l: 1D array of left primitives, using global coordinates
//   prim_r: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
//   bbx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   transformation is trivial

void Minkowski::PrimToLocal1(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb1, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx) {
  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i = il; i <= iu; ++i) {
      bbx(i) = bb1(i);
    }
  }
  return;
}

void Minkowski::PrimToLocal2(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb2, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx) {
  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i = il; i <= iu; ++i) {
      bbx(i) = bb2(i);
    }
  }
  return;
}

void Minkowski::PrimToLocal3(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb3, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx) {
  if (MAGNETIC_FIELDS_ENABLED) {
#pragma omp simd
    for (int i = il; i <= iu; ++i) {
      bbx(i) = bb3(i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming fluxes to global frame: X-interface
// Inputs:
//   k,j: z- and y-indices
//   il,iu: x-index bounds
//   cons: 1D array of conserved quantities, using local coordinates (not used)
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates (not used)
//   flux: 3D array of hydrodynamical fluxes, using local coordinates
//   ey,ez: 3D arrays of magnetic fluxes (electric fields) (not used)
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   transformation is trivial except for sign change from lowering time index

void Minkowski::FluxToGlobal1(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    const Real &txt = flux(IEN,k,j,i);
    Real &t10 = flux(IEN,k,j,i);
    t10 = -txt;
  }
  return;
}

void Minkowski::FluxToGlobal2(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    const Real &tyt = flux(IEN,k,j,i);
    Real &t20 = flux(IEN,k,j,i);
    t20 = -tyt;
  }
  return;
}

void Minkowski::FluxToGlobal3(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    const Real &tzt = flux(IEN,k,j,i);
    Real &t30 = flux(IEN,k,j,i);
    t30 = -tzt;
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

void Minkowski::RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3, int k, int j,
                                int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3) {
  *pa0 = -a_0;
  *pa1 = a_1;
  *pa2 = a_2;
  *pa3 = a_3;
  return;
}

//----------------------------------------------------------------------------------------
// Function for lowering contravariant components of a vector
// Inputs:
//   a0,a1,a2,a3: contravariant components of vector
//   k,j,i: indices of cell in which transformation is desired
// Outputs:
//   pa_0,pa_1,pa_2,pa_3: pointers to covariant 4-vector components

void Minkowski::LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j,
                                int i, Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3) {
  *pa_0 = -a0;
  *pa_1 = a1;
  *pa_2 = a2;
  *pa_3 = a3;
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating orthonormal tetrad
// Inputs:
//   x, y, z: spatial position
// Outputs:
//   e: 2D array for e_{(\hat{\mu})}^\nu:
//     index 0: covariant orthonormal index
//     index 1: contravariant coordinate index
//   e_0: 1D array for {e_{(\hat{\mu})}}_0:
//     index 0: covariant orthonormal index
//   omega: 3D array for \omega^{\hat{\gamma}}_{\hat{\alpha}\hat{\beta}}:
//     index 0: upper index
//     index 1: first lower index
//     index 2: second lower index
// Notes:
//   tetrad options:
//     "cartesian" (Gram-Schmidt on t, x, y, z; trivial)
//     "cylindrical" (Gram-Schmidt on t, R, phi, z)
//     "spherical" (Gram-Schmidt on t, r, theta, phi)

void Minkowski::Tetrad(Real x, Real y, Real z, AthenaArray<Real> &e,
    AthenaArray<Real> &e_0, AthenaArray<Real> &omega) {

  // Check tetrad
  if (rad_tetrad_ != "cartesian" and rad_tetrad_ != "cylindrical"
      and rad_tetrad_ != "spherical") {
    std::stringstream msg;
    msg << "### FATAL ERROR invalid tetrad choice" << std::endl;
    ATHENA_ERROR(msg);
  }

  // Calculate useful quantities
  Real x2 = SQR(x);
  Real y2 = SQR(y);
  Real z2 = SQR(z);
  Real rr = std::hypot(x, y);
  Real rr2 = SQR(rr);
  Real rr3 = rr * rr2;
  Real r = std::hypot(rr, z);
  Real r2 = SQR(r);
  Real r3 = r * r2;

  // Allocate intermediate arrays
  Real eta[4][4] = {};
  Real g[4][4] = {};
  Real ei[4][4] = {};
  Real de[4][4][4] = {};

  // Set Minkowski metric
  eta[0][0] = -1.0;
  eta[1][1] = 1.0;
  eta[2][2] = 1.0;
  eta[3][3] = 1.0;

  // Set covariant metric
  g[0][0] = -1.0;
  g[1][1] = 1.0;
  g[2][2] = 1.0;
  g[3][3] = 1.0;

  // Set tetrad
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      e(i,j) = 0.0;
    }
  }
  if (rad_tetrad_ == "cartesian") {
    for (int i = 0; i < 4; ++i) {
      e(i,i) = 1.0;
    }
  } else if (rad_tetrad_ == "cylindrical") {
    e(0,0) = 1.0;
    e(1,1) = x / rr;
    e(1,2) = y / rr;
    e(2,1) = -y / rr;
    e(2,2) = x / rr;
    e(3,3) = 1.0;
  } else if (rad_tetrad_ == "spherical") {
    e(0,0) = 1.0;
    e(1,1) = x / r;
    e(1,2) = y / r;
    e(1,3) = z / r;
    e(2,1) = x * z / (rr * r);
    e(2,2) = y * z / (rr * r);
    e(2,3) = -rr / r;
    e(3,1) = -y / rr;
    e(3,2) = x / rr;
  }

  // Calculate covariant tetrad
  for (int i = 0; i < 4; ++i) {
    e_0(i) = 0.0;
    for (int j = 0; j < 4; ++j) {
      e_0(i) += g[0][j] * e(i,j);
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
  if (rad_tetrad_ == "cylindrical") {
    de[1][1][1] = y2 / rr3;
    de[1][1][2] = -x * y / rr3;
    de[1][2][1] = x * y / rr3;
    de[1][2][2] = y2 / rr3;
    de[2][1][1] = -x * y / rr3;
    de[2][1][2] = x2 / rr3;
    de[2][2][1] = -x2 / rr3;
    de[2][2][2] = -x * y / rr3;
  } else if (rad_tetrad_ == "spherical") {
    de[1][1][1] = (x2 + z*z) / r3;
    de[1][1][2] = -x * y / r3;
    de[1][1][3] = -x * z / r3;
    de[1][2][1] = (r2 * y2 - rr2 * x2) * z / (rr3 * r3);
    de[1][2][2] = -(rr2 * r2) * x * y * z / (rr3 * r3);
    de[1][2][3] = -x * z*z / (rr * r3);
    de[1][3][1] = x * y / rr3;
    de[1][3][2] = y2 / rr3;
    de[2][1][1] = -x * y / r3;
    de[2][1][2] = (x2 * z*z) / r3;
    de[2][1][3] = -y * z / r3;
    de[2][2][1] = -(rr2 + r2) * x * y * z / (rr3 * r3);
    de[2][2][2] = (r2 * x2 - rr2 * y2) * z / (rr3 * r3);
    de[2][2][3] = -y * z*z / (rr3 * r3);
    de[2][3][1] = -x2 / rr3;
    de[2][3][2] = -x * y / rr3;
  }

  // Calculate Ricci rotation coefficients
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = 0; k < 4; ++k) {
        omega(i,j,k) = 0.0;
        for (int l = 0; l < 4; ++l) {
          for (int m = 0; m < 4; ++m) {
            omega(i,j,k) += ei[i][l] * e(k,m) * de[m][j][l];
          }
        }
      }
    }
  }
  return;
}
