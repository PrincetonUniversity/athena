//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gr_user.cpp
//!  \brief Used for arbitrary stationary coordinates in general relativity, with all
//!  functions evaluated numerically based on the metric.
//!
//!  Original implementation by CJ White.

// C headers

// C++ headers
#include <cmath>  // sqrt()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "coordinates.hpp"

namespace {
// Declarations
Real Determinant(const AthenaArray<Real> &g);
Real Determinant(Real a11, Real a12, Real a13, Real a21, Real a22, Real a23,
                 Real a31, Real a32, Real a33);
Real Determinant(Real a11, Real a12, Real a21, Real a22);
void CalculateTransformation(
    const AthenaArray<Real> &g,
    const AthenaArray<Real> &g_inv, int face, AthenaArray<Real> &transformation);
} // namespace

//----------------------------------------------------------------------------------------
// GRUser Constructor
// Inputs:
//   pmb: pointer to MeshBlock containing this grid
//   pin: pointer to runtime inputs
//   flag: true if object is for coarse grid only in an AMR calculation

GRUser::GRUser(MeshBlock *pmb, ParameterInput *pin, bool flag)
    : Coordinates(pmb, pin, flag) {
  // Set object names
  RegionSize& block_size = pmy_block->block_size;

  // set more indices
  int ill = il - ng;
  int iuu = iu + ng;
  int jll, juu;
  if (block_size.nx2 > 1) {
    jll = jl - ng;
    juu = ju + ng;
  } else {
    jll = jl;
    juu = ju;
  }
  int kll, kuu;
  if (block_size.nx3 > 1) {
    kll = kl - ng;
    kuu = ku + ng;
  } else {
    kll = kl;
    kuu = ku;
  }

  // Set parameters
  bh_mass_ = pin->GetReal("coord", "m");
  bh_spin_ = pin->GetReal("coord", "a");
  // unused:
  // const Real &m = bh_mass_;
  // const Real &a = bh_spin_;

  // Initialize volume-averaged coordinates and spacings: r-direction
  for (int i=ill; i<=iuu; ++i) {
    Real r_m = x1f(i);
    Real r_p = x1f(i+1);
    x1v(i) = 0.5 * (r_m + r_p);  // at least 2nd-order accurate
  }
  for (int i=ill; i<=iuu-1; ++i) {
    dx1v(i) = x1v(i+1) - x1v(i);
  }

  // Initialize volume-averaged coordinates and spacings: theta-direction
  if (block_size.nx2 == 1) {
    Real theta_m = x2f(jl);
    Real theta_p = x2f(jl+1);
    x2v(jl) = 0.5 * (theta_m + theta_p);  // at least 2nd-order accurate
    dx2v(jl) = dx2f(jl);
  } else {
    for (int j=jll; j<=juu; ++j) {
      Real theta_m = x2f(j);
      Real theta_p = x2f(j+1);
      x2v(j) = 0.5 * (theta_m + theta_p);  // at least 2nd-order accurate
    }
    for (int j=jll; j<=juu-1; ++j) {
      dx2v(j) = x2v(j+1) - x2v(j);
    }
  }

  // Initialize volume-averaged coordinates and spacings: phi-direction
  if (block_size.nx3 == 1) {
    Real phi_m = x3f(kl);
    Real phi_p = x3f(kl+1);
    x3v(kl) = 0.5 * (phi_m + phi_p);  // at least 2nd-order accurate
    dx3v(kl) = dx3f(kl);
  } else {
    for (int k=kll; k<=kuu; ++k) {
      Real phi_m = x3f(k);
      Real phi_p = x3f(k+1);
      x3v(k) = 0.5 * (phi_m + phi_p);  // at least 2nd-order accurate
    }
    for (int k=kll; k<=kuu-1; ++k) {
      dx3v(k) = x3v(k+1) - x3v(k);
    }
  }

  // Initialize area-averaged coordinates used with MHD AMR
  if (pm->multilevel && MAGNETIC_FIELDS_ENABLED) {
    for (int i=ill; i<=iuu; ++i) {
      x1s2(i) = x1s3(i) = x1v(i);
    }
    if (block_size.nx2 == 1) {
      x2s1(jl) = x2s3(jl) = x2v(jl);
    } else {
      for (int j=jll; j<=juu; ++j) {
        x2s1(j) = x2s3(j) = x2v(j);
      }
    }
    if (block_size.nx3 == 1) {
      x3s1(kl) = x3s2(kl) = x3v(kl);
    } else {
      for (int k=kll; k<=kuu; ++k) {
        x3s1(k) = x3s2(k) = x3v(k);
      }
    }
  }

  // Allocate arrays for geometric quantities
  metric_cell_kji_.NewAthenaArray(2, NMETRIC, nc3, nc2, nc1);
  if (!coarse_flag) {
    coord_vol_kji_.NewAthenaArray(nc3, nc2, nc1);
    coord_area1_kji_.NewAthenaArray(nc3, nc2, nc1+1);
    coord_area2_kji_.NewAthenaArray(nc3, nc2+1, nc1);
    coord_area3_kji_.NewAthenaArray(nc3+1, nc2, nc1);
    coord_len1_kji_.NewAthenaArray(nc3+1, nc2+1, nc1);
    coord_len2_kji_.NewAthenaArray(nc3+1, nc2, nc1+1);
    coord_len3_kji_.NewAthenaArray(nc3, nc2+1, nc1+1);
    coord_width1_kji_.NewAthenaArray(nc3, nc2, nc1);
    coord_width2_kji_.NewAthenaArray(nc3, nc2, nc1);
    coord_width3_kji_.NewAthenaArray(nc3, nc2, nc1);
    coord_src_kji_.NewAthenaArray(3, NMETRIC, nc3, nc2, nc1);
    metric_face1_kji_.NewAthenaArray(2, NMETRIC, nc3, nc2, nc1+1);
    metric_face2_kji_.NewAthenaArray(2, NMETRIC, nc3, nc2+1, nc1);
    metric_face3_kji_.NewAthenaArray(2, NMETRIC, nc3+1, nc2, nc1);
    trans_face1_kji_.NewAthenaArray(2, NMETRIC, nc3, nc2, nc1+1);
    trans_face2_kji_.NewAthenaArray(2, NMETRIC, nc3, nc2+1, nc1);
    trans_face3_kji_.NewAthenaArray(2, NMETRIC, nc3+1, nc2, nc1);
    g_.NewAthenaArray(NMETRIC, nc1+1);
    gi_.NewAthenaArray(NMETRIC, nc1+1);
  }

  // Allocate scratch arrays
  AthenaArray<Real> g, g_inv, dg_dx1, dg_dx2, dg_dx3, transformation;
  g.NewAthenaArray(NMETRIC);
  g_inv.NewAthenaArray(NMETRIC);
  dg_dx1.NewAthenaArray(NMETRIC);
  dg_dx2.NewAthenaArray(NMETRIC);
  dg_dx3.NewAthenaArray(NMETRIC);
  if (!coarse_flag) {
    transformation.NewAthenaArray(2, NTRIANGULAR);
  }

  // Calculate cell-centered geometric quantities
  for (int k=kll; k<=kuu; ++k) {
    for (int j=jll; j<=juu; ++j) {
      for (int i=ill; i<=iuu; ++i) {
        // Get position and separations
        Real x1 = x1v(i);
        Real x2 = x2v(j);
        Real x3 = x3v(k);
        Real dx1 = dx1f(i);
        Real dx2 = dx2f(j);
        Real dx3 = dx3f(k);

        // Calculate metric coefficients
        Metric(x1, x2, x3, pin, g, g_inv, dg_dx1, dg_dx2, dg_dx3);

        // Calculate volumes
        if (!coarse_flag) {
          Real det = Determinant(g);
          coord_vol_kji_(k,j,i) = std::sqrt(-det) * dx1 * dx2 * dx3;
        }

        // Calculate widths
        if (!coarse_flag) {
          coord_width1_kji_(k,j,i) = std::sqrt(g(I11)) * dx1;
          coord_width2_kji_(k,j,i) = std::sqrt(g(I22)) * dx2;
          coord_width3_kji_(k,j,i) = std::sqrt(g(I33)) * dx3;
        }

        // Store metric derivatives
        if (!coarse_flag) {
          for (int m = 0; m < NMETRIC; ++m) {
            coord_src_kji_(0,m,k,j,i) = dg_dx1(m);
            coord_src_kji_(1,m,k,j,i) = dg_dx2(m);
            coord_src_kji_(2,m,k,j,i) = dg_dx3(m);
          }
        }

        // Set metric coefficients
        for (int n = 0; n < NMETRIC; ++n) {
          metric_cell_kji_(0,n,k,j,i) = g(n);
          metric_cell_kji_(1,n,k,j,i) = g_inv(n);
        }
      }
    }
  }

  // Calculate x1-face-centered geometric quantities
  if (!coarse_flag) {
    for (int k=kll; k<=kuu; ++k) {
      for (int j=jll; j<=juu; ++j) {
        for (int i=ill; i<=iuu+1; ++i) {
          // Get position and separations
          Real x1 = x1f(i);
          Real x2 = x2v(j);
          Real x3 = x3v(k);
          Real dx2 = dx2f(j);
          Real dx3 = dx3f(k);

          // Calculate metric coefficients
          Metric(x1, x2, x3, pin, g, g_inv, dg_dx1, dg_dx2, dg_dx3);

          // Calculate areas
          Real det = Determinant(g);
          coord_area1_kji_(k,j,i) = std::sqrt(-det) * dx2 * dx3;

          // Set metric coefficients
          for (int n = 0; n < NMETRIC; ++n) {
            metric_face1_kji_(0,n,k,j,i) = g(n);
            metric_face1_kji_(1,n,k,j,i) = g_inv(n);
          }

          // Calculate frame transformation
          CalculateTransformation(g, g_inv, 1, transformation);
          for (int n = 0; n < 2; ++n) {
            for (int m = 0; m < NTRIANGULAR; ++m) {
              trans_face1_kji_(n,m,k,j,i) = transformation(n,m);
            }
          }
        }
      }
    }
  }

  // Calculate x2-face-centered geometric quantities
  if (!coarse_flag) {
    for (int k=kll; k<=kuu; ++k) {
      for (int j=jll; j<=juu+1; ++j) {
        for (int i=ill; i<=iuu; ++i) {
          // Get position and separations
          Real x1 = x1v(i);
          Real x2 = x2f(j);
          Real x3 = x3v(k);
          Real dx1 = dx1f(i);
          Real dx3 = dx3f(k);

          // Calculate metric coefficients
          Metric(x1, x2, x3, pin, g, g_inv, dg_dx1, dg_dx2, dg_dx3);

          // Calculate areas
          Real det = Determinant(g);
          coord_area2_kji_(k,j,i) = std::sqrt(-det) * dx1 * dx3;

          // Set metric coefficients
          for (int n = 0; n < NMETRIC; ++n) {
            metric_face2_kji_(0,n,k,j,i) = g(n);
            metric_face2_kji_(1,n,k,j,i) = g_inv(n);
          }

          // Calculate frame transformation
          CalculateTransformation(g, g_inv, 2, transformation);
          for (int n = 0; n < 2; ++n) {
            for (int m = 0; m < NTRIANGULAR; ++m) {
              trans_face2_kji_(n,m,k,j,i) = transformation(n,m);
            }
          }
        }
      }
    }
  }

  // Calculate x3-face-centered geometric quantities
  if (!coarse_flag) {
    for (int k=kll; k<=kuu+1; ++k) {
      for (int j=jll; j<=juu; ++j) {
        for (int i=ill; i<=iuu; ++i) {
          // Get position and separations
          Real x1 = x1v(i);
          Real x2 = x2v(j);
          Real x3 = x3f(k);
          Real dx1 = dx1f(i);
          Real dx2 = dx2f(j);

          // Calculate metric coefficients
          Metric(x1, x2, x3, pin, g, g_inv, dg_dx1, dg_dx2, dg_dx3);

          // Calculate areas
          Real det = Determinant(g);
          coord_area3_kji_(k,j,i) = std::sqrt(-det) * dx1 * dx2;

          // Set metric coefficients
          for (int n = 0; n < NMETRIC; ++n) {
            metric_face3_kji_(0,n,k,j,i) = g(n);
            metric_face3_kji_(1,n,k,j,i) = g_inv(n);
          }

          // Calculate frame transformation
          CalculateTransformation(g, g_inv, 3, transformation);
          for (int n = 0; n < 2; ++n) {
            for (int m = 0; m < NTRIANGULAR; ++m) {
              trans_face3_kji_(n,m,k,j,i) = transformation(n,m);
            }
          }
        }
      }
    }
  }

  // Calculate x1-edge-centered geometric quantities
  if (!coarse_flag) {
    for (int k=kll; k<=kuu+1; ++k) {
      for (int j=jll; j<=juu+1; ++j) {
        for (int i=ill; i<=iuu; ++i) {
          // Get position and separation
          Real x1 = x1v(i);
          Real x2 = x2f(j);
          Real x3 = x3f(k);
          Real dx1 = dx1f(i);

          // Calculate metric coefficients
          Metric(x1, x2, x3, pin, g, g_inv, dg_dx1, dg_dx2, dg_dx3);

          // Calculate lengths
          Real det = Determinant(g);
          coord_len1_kji_(k,j,i) = std::sqrt(-det) * dx1;
        }
      }
    }
  }

  // Calculate x2-edge-centered geometric quantities
  if (!coarse_flag) {
    for (int k=kll; k<=kuu+1; ++k) {
      for (int j=jll; j<=juu; ++j) {
        for (int i=ill; i<=iuu+1; ++i) {
          // Get position and separation
          Real x1 = x1f(i);
          Real x2 = x2v(j);
          Real x3 = x3f(k);
          Real dx2 = dx2f(j);

          // Calculate metric coefficients
          Metric(x1, x2, x3, pin, g, g_inv, dg_dx1, dg_dx2, dg_dx3);

          // Calculate lengths
          Real det = Determinant(g);
          coord_len2_kji_(k,j,i) = std::sqrt(-det) * dx2;
        }
      }
    }
  }

  // Calculate x3-edge-centered geometric quantities
  if (!coarse_flag) {
    for (int k=kll; k<=kuu; ++k) {
      for (int j=jll; j<=juu+1; ++j) {
        for (int i=ill; i<=iuu+1; ++i) {
          // Get position and separation
          Real x1 = x1f(i);
          Real x2 = x2f(j);
          Real x3 = x3v(k);
          Real dx3 = dx3f(k);

          // Calculate metric coefficients
          Metric(x1, x2, x3, pin, g, g_inv, dg_dx1, dg_dx2, dg_dx3);

          // Calculate lengths
          Real det = Determinant(g);
          coord_len3_kji_(k,j,i) = std::sqrt(-det) * dx3;
        }
      }
    }
  }
}


//----------------------------------------------------------------------------------------
// EdgeXLength functions: compute physical length at cell edge-X as vector
// Edge1(i,j,k) located at (i,j-1/2,k-1/2), i.e. (x1v(i), x2f(j), x3f(k))
// Edge2(i,j,k) located at (i-1/2,j,k-1/2), i.e. (x1f(i), x2v(j), x3f(k))
// Edge3(i,j,k) located at (i-1/2,j-1/2,k), i.e. (x1f(i), x2f(j), x3v(k))

void GRUser::Edge1Length(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &lengths) {
  // \Delta L \approx \sqrt{-g} \Delta x^1
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    lengths(i) = coord_len1_kji_(k,j,i);
  }
  return;
}

void GRUser::Edge2Length(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &lengths) {
  // \Delta L \approx \sqrt{-g} \Delta x^2
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    lengths(i) = coord_len2_kji_(k,j,i);
  }
  return;
}

void GRUser::Edge3Length(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &lengths) {
  // \Delta L \approx \sqrt{-g} \Delta x^3
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    lengths(i) = coord_len3_kji_(k,j,i);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetEdgeXLength functions: return length of edge-X at (i,j,k)

Real GRUser::GetEdge1Length(const int k, const int j, const int i) {
  // \Delta L \approx \sqrt{-g} \Delta x^1
  return coord_len1_kji_(k,j,i);
}

Real GRUser::GetEdge2Length(const int k, const int j, const int i) {
  // \Delta L \approx \sqrt{-g} \Delta x^2
  return coord_len2_kji_(k,j,i);
}

Real GRUser::GetEdge3Length(const int k, const int j, const int i) {
  // \Delta L \approx \sqrt{-g} \Delta x^3
  return coord_len3_kji_(k,j,i);
}

//----------------------------------------------------------------------------------------
// CenterWidthX functions: return physical width in X-dir at (i,j,k) cell-center

void GRUser::CenterWidth1(const int k, const int j, const int il, const int iu,
                          AthenaArray<Real> &dx1) {
  // \Delta W \approx \sqrt{g_{11}} \Delta x^1
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    dx1(i) = coord_width1_kji_(k,j,i);
  }
  return;
}

void GRUser::CenterWidth2(const int k, const int j, const int il, const int iu,
                          AthenaArray<Real> &dx2) {
  // \Delta W \approx \sqrt{g_{22}} \Delta x^2
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    dx2(i) = coord_width2_kji_(k,j,i);
  }
  return;
}

void GRUser::CenterWidth3(const int k, const int j, const int il, const int iu,
                          AthenaArray<Real> &dx3) {
  // \Delta W \approx \sqrt{g_{33}} \Delta x^3
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    dx3(i) = coord_width3_kji_(k,j,i);
  }
  return;
}

//----------------------------------------------------------------------------------------
// FaceXArea functions: compute area of face with normal in X-dir as vector
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: x1-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to X-face

void GRUser::Face1Area(const int k, const int j, const int il, const int iu,
                       AthenaArray<Real> &areas) {
  // \Delta A \approx \sqrt{-g} \Delta x^2 \Delta x^3
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    areas(i) = coord_area1_kji_(k,j,i);
  }
  return;
}

void GRUser::Face2Area(const int k, const int j, const int il, const int iu,
                       AthenaArray<Real> &areas) {
  // \Delta A \approx \sqrt{-g} \Delta x^1 \Delta x^3
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    areas(i) = coord_area2_kji_(k,j,i);
  }
  return;
}

void GRUser::Face3Area(const int k, const int j, const int il, const int iu,
                       AthenaArray<Real> &areas) {
  // \Delta A \approx \sqrt{-g} \Delta x^1 \Delta x^2
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    areas(i) = coord_area3_kji_(k,j,i);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetFaceXArea functions: return area of face with normal in X-dir at (i,j,k)
// Inputs:
//   k,j,i: x3-, x2-, and x1-indices
// return:
//   interface area orthogonal to X-face

Real GRUser::GetFace1Area(const int k, const int j, const int i) {
  // \Delta A \approx \sqrt{-g} \Delta x^2 \Delta x^3
  return coord_area1_kji_(k,j,i);
}

Real GRUser::GetFace2Area(const int k, const int j, const int i) {
  // \Delta A \approx \sqrt{-g} \Delta x^1 \Delta x^3
  return coord_area2_kji_(k,j,i);
}

Real GRUser::GetFace3Area(const int k, const int j, const int i) {
  // \Delta A \approx \sqrt{-g} \Delta x^1 \Delta x^2
  return coord_area3_kji_(k,j,i);
}

//----------------------------------------------------------------------------------------
// Cell Volume function: compute volume of cell as vector
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: x1-index bounds
// Outputs:
//   volumes: 1D array of cell volumes

void GRUser::CellVolume(const int k, const int j, const int il, const int iu,
                        AthenaArray<Real> &volumes) {
  // \Delta V \approx \sqrt{-g} \Delta x^1 \Delta x^2 \Delta x^3
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    volumes(i) = coord_vol_kji_(k,j,i);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetCellVolume: returns cell volume at (i,j,k)
// Inputs:
//   k,j,i: phi-, theta-, and r-indices
// Outputs:
//   returned value: cell volume

Real GRUser::GetCellVolume(const int k, const int j, const int i) {
  // \Delta V \approx \sqrt{-g} \Delta x^1 \Delta x^2 \Delta x^3
  return coord_vol_kji_(k,j,i);
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

void GRUser::AddCoordTermsDivergence(const Real dt, const AthenaArray<Real> *flux,
                           const AthenaArray<Real> &prim, const AthenaArray<Real> &bb_cc,
                           AthenaArray<Real> &cons) {
  // Extract indices
  int is = pmy_block->is;
  int ie = pmy_block->ie;
  int js = pmy_block->js;
  int je = pmy_block->je;
  int ks = pmy_block->ks;
  int ke = pmy_block->ke;

  // Extract ratio of specific heats
  Real gamma_adi = pmy_block->peos->GetGamma();

  // Go through cells
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        // Extract metric coefficients
        const Real &g_00 = metric_cell_kji_(0,I00,k,j,i);
        const Real &g_01 = metric_cell_kji_(0,I01,k,j,i);
        const Real &g_02 = metric_cell_kji_(0,I02,k,j,i);
        const Real &g_03 = metric_cell_kji_(0,I03,k,j,i);
        const Real &g_11 = metric_cell_kji_(0,I11,k,j,i);
        const Real &g_12 = metric_cell_kji_(0,I12,k,j,i);
        const Real &g_13 = metric_cell_kji_(0,I13,k,j,i);
        const Real &g_22 = metric_cell_kji_(0,I22,k,j,i);
        const Real &g_23 = metric_cell_kji_(0,I23,k,j,i);
        const Real &g_33 = metric_cell_kji_(0,I33,k,j,i);
        const Real &g00 = metric_cell_kji_(1,I00,k,j,i);
        const Real &g01 = metric_cell_kji_(1,I01,k,j,i);
        const Real &g02 = metric_cell_kji_(1,I02,k,j,i);
        const Real &g03 = metric_cell_kji_(1,I03,k,j,i);
        const Real &g11 = metric_cell_kji_(1,I11,k,j,i);
        const Real &g12 = metric_cell_kji_(1,I12,k,j,i);
        const Real &g13 = metric_cell_kji_(1,I13,k,j,i);
        const Real &g22 = metric_cell_kji_(1,I22,k,j,i);
        const Real &g23 = metric_cell_kji_(1,I23,k,j,i);
        const Real &g33 = metric_cell_kji_(1,I33,k,j,i);
        Real alpha = std::sqrt(-1.0/g00);

        // Extract primitives
        const Real &rho = prim(IDN,k,j,i);
        const Real &pgas = prim(IEN,k,j,i);
        const Real &uu1 = prim(IVX,k,j,i);
        const Real &uu2 = prim(IVY,k,j,i);
        const Real &uu3 = prim(IVZ,k,j,i);

        // Calculate 4-velocity
        Real uu_sq = g_11*uu1*uu1 + 2.0*g_12*uu1*uu2 + 2.0*g_13*uu1*uu3
                     + g_22*uu2*uu2 + 2.0*g_23*uu2*uu3
                     + g_33*uu3*uu3;
        Real gamma = std::sqrt(1.0 + uu_sq);
        Real u0 = gamma / alpha;
        Real u1 = uu1 - alpha * gamma * g01;
        Real u2 = uu2 - alpha * gamma * g02;
        Real u3 = uu3 - alpha * gamma * g03;

        // Extract and calculate magnetic field
        Real b0 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0;
        Real b_sq = 0.0;
        if (MAGNETIC_FIELDS_ENABLED) {
          Real u_1 = g_01*u0 + g_11*u1 + g_12*u2 + g_13*u3;
          Real u_2 = g_02*u0 + g_12*u1 + g_22*u2 + g_23*u3;
          Real u_3 = g_03*u0 + g_13*u1 + g_23*u2 + g_33*u3;
          const Real &bb1 = bb_cc(IB1,k,j,i);
          const Real &bb2 = bb_cc(IB2,k,j,i);
          const Real &bb3 = bb_cc(IB3,k,j,i);
          b0 = u_1*bb1 + u_2*bb2 + u_3*bb3;
          b1 = (bb1 + b0 * u1) / u0;
          b2 = (bb2 + b0 * u2) / u0;
          b3 = (bb3 + b0 * u3) / u0;
          Real b_0 = g_00*b0 + g_01*b1 + g_02*b2 + g_03*b3;
          Real b_1 = g_01*b0 + g_11*b1 + g_12*b2 + g_13*b3;
          Real b_2 = g_02*b0 + g_12*b1 + g_22*b2 + g_23*b3;
          Real b_3 = g_03*b0 + g_13*b1 + g_23*b2 + g_33*b3;
          b_sq = b_0*b0 + b_1*b1 + b_2*b2 + b_3*b3;
        }

        // Calculate stress-energy tensor
        Real wtot = rho + gamma_adi/(gamma_adi-1.0) * pgas + b_sq;
        Real ptot = pgas + 0.5*b_sq;
        Real tt[NMETRIC];
        tt[I00] = wtot * u0 * u0 + ptot * g00 - b0 * b0;
        tt[I01] = wtot * u0 * u1 + ptot * g01 - b0 * b1;
        tt[I02] = wtot * u0 * u2 + ptot * g02 - b0 * b2;
        tt[I03] = wtot * u0 * u3 + ptot * g03 - b0 * b3;
        tt[I11] = wtot * u1 * u1 + ptot * g11 - b1 * b1;
        tt[I12] = wtot * u1 * u2 + ptot * g12 - b1 * b2;
        tt[I13] = wtot * u1 * u3 + ptot * g13 - b1 * b3;
        tt[I22] = wtot * u2 * u2 + ptot * g22 - b2 * b2;
        tt[I23] = wtot * u2 * u3 + ptot * g23 - b2 * b3;
        tt[I33] = wtot * u3 * u3 + ptot * g33 - b3 * b3;

        // Calculate source terms
        Real s_1 = 0.0, s_2 = 0.0, s_3 = 0.0;
        for (int n = 0; n < NMETRIC; ++n) {
          s_1 += coord_src_kji_(0,n,k,j,i) * tt[n];
          s_2 += coord_src_kji_(1,n,k,j,i) * tt[n];
          s_3 += coord_src_kji_(2,n,k,j,i) * tt[n];
        }
        s_1 -= 0.5 * (coord_src_kji_(0,I00,k,j,i) * tt[I00]
                      + coord_src_kji_(0,I11,k,j,i) * tt[I11]
                      + coord_src_kji_(0,I22,k,j,i) * tt[I22]
                      + coord_src_kji_(0,I33,k,j,i) * tt[I33]);
        s_2 -= 0.5 * (coord_src_kji_(1,I00,k,j,i) * tt[I00]
                      + coord_src_kji_(1,I11,k,j,i) * tt[I11]
                      + coord_src_kji_(1,I22,k,j,i) * tt[I22]
                      + coord_src_kji_(1,I33,k,j,i) * tt[I33]);
        s_3 -= 0.5 * (coord_src_kji_(2,I00,k,j,i) * tt[I00]
                      + coord_src_kji_(2,I11,k,j,i) * tt[I11]
                      + coord_src_kji_(2,I22,k,j,i) * tt[I22]
                      + coord_src_kji_(2,I33,k,j,i) * tt[I33]);

        // Extract conserved quantities
        Real &m_1 = cons(IM1,k,j,i);
        Real &m_2 = cons(IM2,k,j,i);
        Real &m_3 = cons(IM3,k,j,i);

        // Add source terms to conserved quantities
        m_1 += dt * s_1;
        m_2 += dt * s_2;
        m_3 += dt * s_3;
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for computing cell-centered and face-centered metric coefficients
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: x1-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D

void GRUser::CellMetric(const int k, const int j, const int il, const int iu,
                        AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  for (int n = 0; n < NMETRIC; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      g(n,i) = metric_cell_kji_(0,n,k,j,i);
      g_inv(n,i) = metric_cell_kji_(1,n,k,j,i);
    }
  }
  return;
}

void GRUser::Face1Metric(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  for (int n = 0; n < NMETRIC; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      g(n,i) = metric_face1_kji_(0,n,k,j,i);
      g_inv(n,i) = metric_face1_kji_(1,n,k,j,i);
    }
  }
  return;
}

void GRUser::Face2Metric(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  for (int n = 0; n < NMETRIC; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      g(n,i) = metric_face2_kji_(0,n,k,j,i);
      g_inv(n,i) = metric_face2_kji_(1,n,k,j,i);
    }
  }
  return;
}

void GRUser::Face3Metric(const int k, const int j, const int il, const int iu,
                         AthenaArray<Real> &g, AthenaArray<Real> &g_inv) {
  for (int n = 0; n < NMETRIC; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      g(n,i) = metric_face3_kji_(0,n,k,j,i);
      g_inv(n,i) = metric_face3_kji_(1,n,k,j,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming primitives to locally flat frame: x1-interface
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: x1-index bounds
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

void GRUser::PrimToLocal1(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb1, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx) {
  // Go through 1D block of cells
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // Extract transformation coefficients
    const Real &mt_0 = trans_face1_kji_(1,T00,k,j,i);
    const Real &mx_0 = trans_face1_kji_(1,T10,k,j,i);
    const Real &mx_1 = trans_face1_kji_(1,T11,k,j,i);
    Real mx_2 = 0.0;
    Real mx_3 = 0.0;
    const Real &my_0 = trans_face1_kji_(1,T20,k,j,i);
    const Real &my_1 = trans_face1_kji_(1,T21,k,j,i);
    const Real &my_2 = trans_face1_kji_(1,T22,k,j,i);
    Real my_3 = 0.0;
    const Real &mz_0 = trans_face1_kji_(1,T30,k,j,i);
    const Real &mz_1 = trans_face1_kji_(1,T31,k,j,i);
    const Real &mz_2 = trans_face1_kji_(1,T32,k,j,i);
    const Real &mz_3 = trans_face1_kji_(1,T33,k,j,i);

    // Extract global projected 4-velocities
    Real uu1_l = prim_l(IVX,i);
    Real uu2_l = prim_l(IVY,i);
    Real uu3_l = prim_l(IVZ,i);
    Real uu1_r = prim_r(IVX,i);
    Real uu2_r = prim_r(IVY,i);
    Real uu3_r = prim_r(IVZ,i);

    // Transform projected 4-velocities
    Real ux_l = mx_1*uu1_l + mx_2*uu2_l + mx_3*uu3_l;
    Real uy_l = my_1*uu1_l + my_2*uu2_l + my_3*uu3_l;
    Real uz_l = mz_1*uu1_l + mz_2*uu2_l + mz_3*uu3_l;
    Real ux_r = mx_1*uu1_r + mx_2*uu2_r + mx_3*uu3_r;
    Real uy_r = my_1*uu1_r + my_2*uu2_r + my_3*uu3_r;
    Real uz_r = mz_1*uu1_r + mz_2*uu2_r + mz_3*uu3_r;

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
      const Real &g_00 = metric_face1_kji_(0,I00,k,j,i);
      const Real &g_01 = metric_face1_kji_(0,I01,k,j,i);
      const Real &g_02 = metric_face1_kji_(0,I02,k,j,i);
      const Real &g_03 = metric_face1_kji_(0,I03,k,j,i);
      const Real &g_10 = metric_face1_kji_(0,I01,k,j,i);
      const Real &g_11 = metric_face1_kji_(0,I11,k,j,i);
      const Real &g_12 = metric_face1_kji_(0,I12,k,j,i);
      const Real &g_13 = metric_face1_kji_(0,I13,k,j,i);
      const Real &g_20 = metric_face1_kji_(0,I02,k,j,i);
      const Real &g_21 = metric_face1_kji_(0,I12,k,j,i);
      const Real &g_22 = metric_face1_kji_(0,I22,k,j,i);
      const Real &g_23 = metric_face1_kji_(0,I23,k,j,i);
      const Real &g_30 = metric_face1_kji_(0,I03,k,j,i);
      const Real &g_31 = metric_face1_kji_(0,I13,k,j,i);
      const Real &g_32 = metric_face1_kji_(0,I23,k,j,i);
      const Real &g_33 = metric_face1_kji_(0,I33,k,j,i);
      const Real &g00 = metric_face1_kji_(1,I00,k,j,i);
      const Real &g01 = metric_face1_kji_(1,I01,k,j,i);
      const Real &g02 = metric_face1_kji_(1,I02,k,j,i);
      const Real &g03 = metric_face1_kji_(1,I03,k,j,i);
      const Real &g10 = metric_face1_kji_(1,I01,k,j,i);
      const Real &g11 = metric_face1_kji_(1,I11,k,j,i);
      const Real &g12 = metric_face1_kji_(1,I12,k,j,i);
      const Real &g13 = metric_face1_kji_(1,I13,k,j,i);
      const Real &g20 = metric_face1_kji_(1,I02,k,j,i);
      const Real &g21 = metric_face1_kji_(1,I12,k,j,i);
      const Real &g22 = metric_face1_kji_(1,I22,k,j,i);
      const Real &g23 = metric_face1_kji_(1,I23,k,j,i);
      const Real &g30 = metric_face1_kji_(1,I03,k,j,i);
      const Real &g31 = metric_face1_kji_(1,I13,k,j,i);
      const Real &g32 = metric_face1_kji_(1,I23,k,j,i);
      const Real &g33 = metric_face1_kji_(1,I33,k,j,i);
      Real alpha = std::sqrt(-1.0/g00);

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
      const Real &bb1_l = bb1(i);
      const Real &bb1_r = bb1(i);
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
      Real ut_l = gamma_l;
      Real ut_r = gamma_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_0*b0_l + mx_1*b1_l + mx_2*b2_l + mx_3*b3_l;
      Real by_l = my_0*b0_l + my_1*b1_l + my_2*b2_l + my_3*b3_l;
      Real bz_l = mz_0*b0_l + mz_1*b1_l + mz_2*b2_l + mz_3*b3_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_0*b0_r + mx_1*b1_r + mx_2*b2_r + mx_3*b3_r;
      Real by_r = my_0*b0_r + my_1*b1_r + my_2*b2_r + my_3*b3_r;
      Real bz_r = mz_0*b0_r + mz_1*b1_r + mz_2*b2_r + mz_3*b3_r;

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
// Function for transforming primitives to locally flat frame: x2-interface
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: x1-index bounds
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

void GRUser::PrimToLocal2(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb2, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx) {
  // Go through 1D block of cells
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // Extract transformation coefficients
    const Real &mt_0 = trans_face2_kji_(1,T00,k,j,i);
    const Real &mx_0 = trans_face2_kji_(1,T10,k,j,i);
    Real mx_1 = 0.0;
    const Real &mx_2 = trans_face2_kji_(1,T11,k,j,i);
    Real mx_3 = 0.0;
    const Real &my_0 = trans_face2_kji_(1,T20,k,j,i);
    Real my_1 = 0.0;
    const Real &my_2 = trans_face2_kji_(1,T21,k,j,i);
    const Real &my_3 = trans_face2_kji_(1,T22,k,j,i);
    const Real &mz_0 = trans_face2_kji_(1,T30,k,j,i);
    const Real &mz_1 = trans_face2_kji_(1,T33,k,j,i);
    const Real &mz_2 = trans_face2_kji_(1,T31,k,j,i);
    const Real &mz_3 = trans_face2_kji_(1,T32,k,j,i);

    // Extract global projected 4-velocities
    Real uu1_l = prim_l(IVX,i);
    Real uu2_l = prim_l(IVY,i);
    Real uu3_l = prim_l(IVZ,i);
    Real uu1_r = prim_r(IVX,i);
    Real uu2_r = prim_r(IVY,i);
    Real uu3_r = prim_r(IVZ,i);

    // Transform projected 4-velocities
    Real ux_l = mx_1*uu1_l + mx_2*uu2_l + mx_3*uu3_l;
    Real uy_l = my_1*uu1_l + my_2*uu2_l + my_3*uu3_l;
    Real uz_l = mz_1*uu1_l + mz_2*uu2_l + mz_3*uu3_l;
    Real ux_r = mx_1*uu1_r + mx_2*uu2_r + mx_3*uu3_r;
    Real uy_r = my_1*uu1_r + my_2*uu2_r + my_3*uu3_r;
    Real uz_r = mz_1*uu1_r + mz_2*uu2_r + mz_3*uu3_r;

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
      const Real &g_00 = metric_face2_kji_(0,I00,k,j,i);
      const Real &g_01 = metric_face2_kji_(0,I01,k,j,i);
      const Real &g_02 = metric_face2_kji_(0,I02,k,j,i);
      const Real &g_03 = metric_face2_kji_(0,I03,k,j,i);
      const Real &g_10 = metric_face2_kji_(0,I01,k,j,i);
      const Real &g_11 = metric_face2_kji_(0,I11,k,j,i);
      const Real &g_12 = metric_face2_kji_(0,I12,k,j,i);
      const Real &g_13 = metric_face2_kji_(0,I13,k,j,i);
      const Real &g_20 = metric_face2_kji_(0,I02,k,j,i);
      const Real &g_21 = metric_face2_kji_(0,I12,k,j,i);
      const Real &g_22 = metric_face2_kji_(0,I22,k,j,i);
      const Real &g_23 = metric_face2_kji_(0,I23,k,j,i);
      const Real &g_30 = metric_face2_kji_(0,I03,k,j,i);
      const Real &g_31 = metric_face2_kji_(0,I13,k,j,i);
      const Real &g_32 = metric_face2_kji_(0,I23,k,j,i);
      const Real &g_33 = metric_face2_kji_(0,I33,k,j,i);
      const Real &g00 = metric_face2_kji_(1,I00,k,j,i);
      const Real &g01 = metric_face2_kji_(1,I01,k,j,i);
      const Real &g02 = metric_face2_kji_(1,I02,k,j,i);
      const Real &g03 = metric_face2_kji_(1,I03,k,j,i);
      const Real &g10 = metric_face2_kji_(1,I01,k,j,i);
      const Real &g11 = metric_face2_kji_(1,I11,k,j,i);
      const Real &g12 = metric_face2_kji_(1,I12,k,j,i);
      const Real &g13 = metric_face2_kji_(1,I13,k,j,i);
      const Real &g20 = metric_face2_kji_(1,I02,k,j,i);
      const Real &g21 = metric_face2_kji_(1,I12,k,j,i);
      const Real &g22 = metric_face2_kji_(1,I22,k,j,i);
      const Real &g23 = metric_face2_kji_(1,I23,k,j,i);
      const Real &g30 = metric_face2_kji_(1,I03,k,j,i);
      const Real &g31 = metric_face2_kji_(1,I13,k,j,i);
      const Real &g32 = metric_face2_kji_(1,I23,k,j,i);
      const Real &g33 = metric_face2_kji_(1,I33,k,j,i);
      Real alpha = std::sqrt(-1.0/g00);

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
      const Real &bb2_l = bb2(i);
      const Real &bb2_r = bb2(i);
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
      Real ut_l = gamma_l;
      Real ut_r = gamma_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_0*b0_l + mx_1*b1_l + mx_2*b2_l + mx_3*b3_l;
      Real by_l = my_0*b0_l + my_1*b1_l + my_2*b2_l + my_3*b3_l;
      Real bz_l = mz_0*b0_l + mz_1*b1_l + mz_2*b2_l + mz_3*b3_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_0*b0_r + mx_1*b1_r + mx_2*b2_r + mx_3*b3_r;
      Real by_r = my_0*b0_r + my_1*b1_r + my_2*b2_r + my_3*b3_r;
      Real bz_r = mz_0*b0_r + mz_1*b1_r + mz_2*b2_r + mz_3*b3_r;

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
// Function for transforming primitives to locally flat frame: x3-interface
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: x1-index bounds
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

void GRUser::PrimToLocal3(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb3, AthenaArray<Real> &prim_l, AthenaArray<Real> &prim_r,
    AthenaArray<Real> &bbx) {
  // Go through 1D block of cells
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // Extract transformation coefficients
    const Real &mt_0 = trans_face3_kji_(1,T00,k,j,i);
    const Real &mx_0 = trans_face3_kji_(1,T10,k,j,i);
    Real mx_1 = 0.0;
    Real mx_2 = 0.0;
    const Real &mx_3 = trans_face3_kji_(1,T11,k,j,i);
    const Real &my_0 = trans_face3_kji_(1,T20,k,j,i);
    const Real &my_1 = trans_face3_kji_(1,T22,k,j,i);
    Real my_2 = 0.0;
    const Real &my_3 = trans_face3_kji_(1,T21,k,j,i);
    const Real &mz_0 = trans_face3_kji_(1,T30,k,j,i);
    const Real &mz_1 = trans_face3_kji_(1,T32,k,j,i);
    const Real &mz_2 = trans_face3_kji_(1,T33,k,j,i);
    const Real &mz_3 = trans_face3_kji_(1,T31,k,j,i);

    // Extract global projected 4-velocities
    Real uu1_l = prim_l(IVX,i);
    Real uu2_l = prim_l(IVY,i);
    Real uu3_l = prim_l(IVZ,i);
    Real uu1_r = prim_r(IVX,i);
    Real uu2_r = prim_r(IVY,i);
    Real uu3_r = prim_r(IVZ,i);

    // Transform projected 4-velocities
    Real ux_l = mx_1*uu1_l + mx_2*uu2_l + mx_3*uu3_l;
    Real uy_l = my_1*uu1_l + my_2*uu2_l + my_3*uu3_l;
    Real uz_l = mz_1*uu1_l + mz_2*uu2_l + mz_3*uu3_l;
    Real ux_r = mx_1*uu1_r + mx_2*uu2_r + mx_3*uu3_r;
    Real uy_r = my_1*uu1_r + my_2*uu2_r + my_3*uu3_r;
    Real uz_r = mz_1*uu1_r + mz_2*uu2_r + mz_3*uu3_r;

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
      const Real &g_00 = metric_face3_kji_(0,I00,k,j,i);
      const Real &g_01 = metric_face3_kji_(0,I01,k,j,i);
      const Real &g_02 = metric_face3_kji_(0,I02,k,j,i);
      const Real &g_03 = metric_face3_kji_(0,I03,k,j,i);
      const Real &g_10 = metric_face3_kji_(0,I01,k,j,i);
      const Real &g_11 = metric_face3_kji_(0,I11,k,j,i);
      const Real &g_12 = metric_face3_kji_(0,I12,k,j,i);
      const Real &g_13 = metric_face3_kji_(0,I13,k,j,i);
      const Real &g_20 = metric_face3_kji_(0,I02,k,j,i);
      const Real &g_21 = metric_face3_kji_(0,I12,k,j,i);
      const Real &g_22 = metric_face3_kji_(0,I22,k,j,i);
      const Real &g_23 = metric_face3_kji_(0,I23,k,j,i);
      const Real &g_30 = metric_face3_kji_(0,I03,k,j,i);
      const Real &g_31 = metric_face3_kji_(0,I13,k,j,i);
      const Real &g_32 = metric_face3_kji_(0,I23,k,j,i);
      const Real &g_33 = metric_face3_kji_(0,I33,k,j,i);
      const Real &g00 = metric_face3_kji_(1,I00,k,j,i);
      const Real &g01 = metric_face3_kji_(1,I01,k,j,i);
      const Real &g02 = metric_face3_kji_(1,I02,k,j,i);
      const Real &g03 = metric_face3_kji_(1,I03,k,j,i);
      const Real &g10 = metric_face3_kji_(1,I01,k,j,i);
      const Real &g11 = metric_face3_kji_(1,I11,k,j,i);
      const Real &g12 = metric_face3_kji_(1,I12,k,j,i);
      const Real &g13 = metric_face3_kji_(1,I13,k,j,i);
      const Real &g20 = metric_face3_kji_(1,I02,k,j,i);
      const Real &g21 = metric_face3_kji_(1,I12,k,j,i);
      const Real &g22 = metric_face3_kji_(1,I22,k,j,i);
      const Real &g23 = metric_face3_kji_(1,I23,k,j,i);
      const Real &g30 = metric_face3_kji_(1,I03,k,j,i);
      const Real &g31 = metric_face3_kji_(1,I13,k,j,i);
      const Real &g32 = metric_face3_kji_(1,I23,k,j,i);
      const Real &g33 = metric_face3_kji_(1,I33,k,j,i);
      Real alpha = std::sqrt(-1.0/g00);

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
      const Real &bb3_l = bb3(i);
      const Real &bb3_r = bb3(i);
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
      Real ut_l = gamma_l;
      Real ut_r = gamma_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_0*b0_l + mx_1*b1_l + mx_2*b2_l + mx_3*b3_l;
      Real by_l = my_0*b0_l + my_1*b1_l + my_2*b2_l + my_3*b3_l;
      Real bz_l = mz_0*b0_l + mz_1*b1_l + mz_2*b2_l + mz_3*b3_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_0*b0_r + mx_1*b1_r + mx_2*b2_r + mx_3*b3_r;
      Real by_r = my_0*b0_r + my_1*b1_r + my_2*b2_r + my_3*b3_r;
      Real bz_r = mz_0*b0_r + mz_1*b1_r + mz_2*b2_r + mz_3*b3_r;

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
// Function for transforming fluxes to global frame: x1-interface
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: x1-index bounds
//   cons: 1D array of conserved quantities, using local coordinates
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates
//   flux: 3D array of hydrodynamical fluxes, using local coordinates
//   ey,ez: 3D arrays of magnetic fluxes (electric fields), using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
//   ey,ez: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM1/IM2/IM3 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots and ey/ez
//   puts x1-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts x1-fluxes of B2/B3 in ey/ez

void GRUser::FluxToGlobal1(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  // Go through 1D block of cells
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // Extract transformation coefficients
    const Real &m0_tm = trans_face1_kji_(0,T00,k,j,i);
    const Real &m1_tm = trans_face1_kji_(0,T10,k,j,i);
    const Real &m1_x = trans_face1_kji_(0,T11,k,j,i);
    const Real &m2_tm = trans_face1_kji_(0,T20,k,j,i);
    const Real &m2_x = trans_face1_kji_(0,T21,k,j,i);
    const Real &m2_y = trans_face1_kji_(0,T22,k,j,i);
    const Real &m3_tm = trans_face1_kji_(0,T30,k,j,i);
    const Real &m3_x = trans_face1_kji_(0,T31,k,j,i);
    const Real &m3_y = trans_face1_kji_(0,T32,k,j,i);
    const Real &m3_z = trans_face1_kji_(0,T33,k,j,i);

    // Extract local conserved quantities and fluxes
    Real jt = cons(IDN,i);
    Real ttt = cons(IEN,i);
    Real ttx = cons(IM1,i);
    Real tty = cons(IM2,i);
    Real ttz = cons(IM3,i);
    Real jx = flux(IDN,k,j,i);
    Real txt = flux(IEN,k,j,i);
    Real txx = flux(IM1,k,j,i);
    Real txy = flux(IM2,k,j,i);
    Real txz = flux(IM3,k,j,i);

    // Transform stress-energy tensor
    Real t10 = m1_tm * (m0_tm*ttt)
               + m1_x * (m0_tm*txt);
    Real t11 = m1_tm * (m1_tm*ttt + m1_x*ttx)
               + m1_x * (m1_tm*txt + m1_x*txx);
    Real t12 = m1_tm * (m2_tm*ttt + m2_x*ttx + m2_y*tty)
               + m1_x * (m2_tm*txt + m2_x*txx + m2_y*txy);
    Real t13 = m1_tm * (m3_tm*ttt + m3_x*ttx + m3_y*tty + m3_z*ttz)
               + m1_x * (m3_tm*txt + m3_x*txx + m3_y*txy + m3_z*txz);

    // Extract metric coefficients
    const Real &g_00 = metric_face1_kji_(0,I00,k,j,i);
    const Real &g_01 = metric_face1_kji_(0,I01,k,j,i);
    const Real &g_02 = metric_face1_kji_(0,I02,k,j,i);
    const Real &g_03 = metric_face1_kji_(0,I03,k,j,i);
    const Real &g_10 = metric_face1_kji_(0,I01,k,j,i);
    const Real &g_11 = metric_face1_kji_(0,I11,k,j,i);
    const Real &g_12 = metric_face1_kji_(0,I12,k,j,i);
    const Real &g_13 = metric_face1_kji_(0,I13,k,j,i);
    const Real &g_20 = metric_face1_kji_(0,I02,k,j,i);
    const Real &g_21 = metric_face1_kji_(0,I12,k,j,i);
    const Real &g_22 = metric_face1_kji_(0,I22,k,j,i);
    const Real &g_23 = metric_face1_kji_(0,I23,k,j,i);
    const Real &g_30 = metric_face1_kji_(0,I03,k,j,i);
    const Real &g_31 = metric_face1_kji_(0,I13,k,j,i);
    const Real &g_32 = metric_face1_kji_(0,I23,k,j,i);
    const Real &g_33 = metric_face1_kji_(0,I33,k,j,i);

    // Extract global fluxes
    Real &j1 = flux(IDN,k,j,i);
    Real &t1_0 = flux(IEN,k,j,i);
    Real &t1_1 = flux(IM1,k,j,i);
    Real &t1_2 = flux(IM2,k,j,i);
    Real &t1_3 = flux(IM3,k,j,i);

    // Set fluxes
    j1 = m1_tm*jt + m1_x*jx;
    t1_0 = g_00*t10 + g_01*t11 + g_02*t12 + g_03*t13;
    t1_1 = g_10*t10 + g_11*t11 + g_12*t12 + g_13*t13;
    t1_2 = g_20*t10 + g_21*t11 + g_22*t12 + g_23*t13;
    t1_3 = g_30*t10 + g_31*t11 + g_32*t12 + g_33*t13;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED) {
      Real ftt = 0.0;
      Real fxt = bbx(i);
      Real fyt = cons(IBY,i);
      Real fzt = cons(IBZ,i);
      Real ftx = -bbx(i);
      Real fxx = 0.0;
      Real fyx = -ey(k,j,i);
      Real fzx = ez(k,j,i);
      Real f21 = m1_tm * (m2_tm*ftt + m2_x*fxt + m2_y*fyt)
                 + m1_x * (m2_tm*ftx + m2_x*fxx + m2_y*fyx);
      Real f31 = m1_tm * (m3_tm*ftt + m3_x*fxt + m3_y*fyt + m3_z*fzt)
                 + m1_x * (m3_tm*ftx + m3_x*fxx + m3_y*fyx + m3_z*fzx);
      ey(k,j,i) = -f21;
      ez(k,j,i) = f31;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming fluxes to global frame: x2-interface
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: x1-index bounds
//   cons: 1D array of conserved quantities, using local coordinates
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates
//   flux: 3D array of hydrodynamical fluxes, using local coordinates
//   ey,ez: 3D arrays of magnetic fluxes (electric fields), using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
//   ey,ez: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM2/IM3/IM1 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots and ey/ez
//   puts x2-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts x2-fluxes of B3/B1 in ey/ez

void GRUser::FluxToGlobal2(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  // Go through 1D block of cells
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // Extract transformation coefficients
    const Real &m0_tm = trans_face2_kji_(0,T00,k,j,i);
    const Real &m1_tm = trans_face2_kji_(0,T30,k,j,i);
    const Real &m1_x = trans_face2_kji_(0,T31,k,j,i);
    const Real &m1_y = trans_face2_kji_(0,T32,k,j,i);
    const Real &m1_z = trans_face2_kji_(0,T33,k,j,i);
    const Real &m2_tm = trans_face2_kji_(0,T10,k,j,i);
    const Real &m2_x = trans_face2_kji_(0,T11,k,j,i);
    const Real &m3_tm = trans_face2_kji_(0,T20,k,j,i);
    const Real &m3_x = trans_face2_kji_(0,T21,k,j,i);
    const Real &m3_y = trans_face2_kji_(0,T22,k,j,i);

    // Extract local conserved quantities and fluxes
    Real jt = cons(IDN,i);
    Real ttt = cons(IEN,i);
    Real ttx = cons(IM1,i);
    Real tty = cons(IM2,i);
    Real ttz = cons(IM3,i);
    Real jx = flux(IDN,k,j,i);
    Real txt = flux(IEN,k,j,i);
    Real txx = flux(IM2,k,j,i);
    Real txy = flux(IM3,k,j,i);
    Real txz = flux(IM1,k,j,i);

    // Transform stress-energy tensor
    Real t20 = m2_tm * (m0_tm*ttt)
               + m2_x * (m0_tm*txt);
    Real t21 = m2_tm * (m1_tm*ttt + m1_x*ttx + m1_y*tty + m1_z*ttz)
               + m2_x * (m1_tm*txt + m1_x*txx + m1_y*txy + m1_z*txz);
    Real t22 = m2_tm * (m2_tm*ttt + m2_x*ttx)
               + m2_x * (m2_tm*txt + m2_x*txx);
    Real t23 = m2_tm * (m3_tm*ttt + m3_x*ttx + m3_y*tty)
               + m2_x * (m3_tm*txt + m3_x*txx + m3_y*txy);

    // Extract metric coefficients
    const Real &g_00 = metric_face2_kji_(0,I00,k,j,i);
    const Real &g_01 = metric_face2_kji_(0,I01,k,j,i);
    const Real &g_02 = metric_face2_kji_(0,I02,k,j,i);
    const Real &g_03 = metric_face2_kji_(0,I03,k,j,i);
    const Real &g_10 = metric_face2_kji_(0,I01,k,j,i);
    const Real &g_11 = metric_face2_kji_(0,I11,k,j,i);
    const Real &g_12 = metric_face2_kji_(0,I12,k,j,i);
    const Real &g_13 = metric_face2_kji_(0,I13,k,j,i);
    const Real &g_20 = metric_face2_kji_(0,I02,k,j,i);
    const Real &g_21 = metric_face2_kji_(0,I12,k,j,i);
    const Real &g_22 = metric_face2_kji_(0,I22,k,j,i);
    const Real &g_23 = metric_face2_kji_(0,I23,k,j,i);
    const Real &g_30 = metric_face2_kji_(0,I03,k,j,i);
    const Real &g_31 = metric_face2_kji_(0,I13,k,j,i);
    const Real &g_32 = metric_face2_kji_(0,I23,k,j,i);
    const Real &g_33 = metric_face2_kji_(0,I33,k,j,i);

    // Extract global fluxes
    Real &j2 = flux(IDN,k,j,i);
    Real &t2_0 = flux(IEN,k,j,i);
    Real &t2_1 = flux(IM1,k,j,i);
    Real &t2_2 = flux(IM2,k,j,i);
    Real &t2_3 = flux(IM3,k,j,i);

    // Set fluxes
    j2 = m2_tm*jt + m2_x*jx;
    t2_0 = g_00*t20 + g_01*t21 + g_02*t22 + g_03*t23;
    t2_1 = g_10*t20 + g_11*t21 + g_12*t22 + g_13*t23;
    t2_2 = g_20*t20 + g_21*t21 + g_22*t22 + g_23*t23;
    t2_3 = g_30*t20 + g_31*t21 + g_32*t22 + g_33*t23;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED) {
      Real ftt = 0.0;
      Real fxt = bbx(i);
      Real fyt = cons(IBY,i);
      Real fzt = cons(IBZ,i);
      Real ftx = -bbx(i);
      Real fxx = 0.0;
      Real fyx = -ey(k,j,i);
      Real fzx = ez(k,j,i);
      Real f32 = m2_tm * (m3_tm*ftt + m3_x*fxt + m3_y*fyt)
                 + m2_x * (m3_tm*ftx + m3_x*fxx + m3_y*fyx);
      Real f12 = m2_tm * (m1_tm*ftt + m1_x*fxt + m1_y*fyt + m1_z*fzt)
                 + m2_x * (m1_tm*ftx + m1_x*fxx + m1_y*fyx + m1_z*fzx);
      ey(k,j,i) = -f32;
      ez(k,j,i) = f12;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for transforming fluxes to global frame: x3-interface
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: x1-index bounds
//   cons: 1D array of conserved quantities, using local coordinates
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates
//   flux: 3D array of hydrodynamical fluxes, using local coordinates
//   ey,ez: 3D arrays of magnetic fluxes (electric fields), using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
//   ey,ez: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM3/IM1/IM2 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots and ey/ez
//   puts x3-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts x3-fluxes of B1/B2 in ey/ez

void GRUser::FluxToGlobal3(
    const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx, AthenaArray<Real> &flux,
    AthenaArray<Real> &ey, AthenaArray<Real> &ez) {
  // Go through 1D block of cells
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // Extract transformation coefficients
    const Real &m0_tm = trans_face3_kji_(0,T00,k,j,i);
    const Real &m1_tm = trans_face3_kji_(0,T20,k,j,i);
    const Real &m1_x = trans_face3_kji_(0,T21,k,j,i);
    const Real &m1_y = trans_face3_kji_(0,T22,k,j,i);
    const Real &m2_tm = trans_face3_kji_(0,T30,k,j,i);
    const Real &m2_x = trans_face3_kji_(0,T31,k,j,i);
    const Real &m2_y = trans_face3_kji_(0,T32,k,j,i);
    const Real &m2_z = trans_face3_kji_(0,T33,k,j,i);
    const Real &m3_tm = trans_face3_kji_(0,T10,k,j,i);
    const Real &m3_x = trans_face3_kji_(0,T11,k,j,i);

    // Extract local conserved quantities and fluxes
    Real jt = cons(IDN,i);
    Real ttt = cons(IEN,i);
    Real ttx = cons(IM1,i);
    Real tty = cons(IM2,i);
    Real ttz = cons(IM3,i);
    Real jx = flux(IDN,k,j,i);
    Real txt = flux(IEN,k,j,i);
    Real txx = flux(IM3,k,j,i);
    Real txy = flux(IM1,k,j,i);
    Real txz = flux(IM2,k,j,i);

    // Transform stress-energy tensor
    Real t30 = m3_tm * (m0_tm*ttt)
               + m3_x * (m0_tm*txt);
    Real t31 = m3_tm * (m1_tm*ttt + m1_x*ttx + m1_y*tty)
               + m3_x * (m1_tm*txt + m1_x*txx + m1_y*txy);
    Real t32 = m3_tm * (m2_tm*ttt + m2_x*ttx + m2_y*tty + m2_z*ttz)
               + m3_x * (m2_tm*txt + m2_x*txx + m2_y*txy + m2_z*txz);
    Real t33 = m3_tm * (m3_tm*ttt + m3_x*ttx)
               + m3_x * (m3_tm*txt + m3_x*txx);

    // Extract metric coefficients
    const Real &g_00 = metric_face3_kji_(0,I00,k,j,i);
    const Real &g_01 = metric_face3_kji_(0,I01,k,j,i);
    const Real &g_02 = metric_face3_kji_(0,I02,k,j,i);
    const Real &g_03 = metric_face3_kji_(0,I03,k,j,i);
    const Real &g_10 = metric_face3_kji_(0,I01,k,j,i);
    const Real &g_11 = metric_face3_kji_(0,I11,k,j,i);
    const Real &g_12 = metric_face3_kji_(0,I12,k,j,i);
    const Real &g_13 = metric_face3_kji_(0,I13,k,j,i);
    const Real &g_20 = metric_face3_kji_(0,I02,k,j,i);
    const Real &g_21 = metric_face3_kji_(0,I12,k,j,i);
    const Real &g_22 = metric_face3_kji_(0,I22,k,j,i);
    const Real &g_23 = metric_face3_kji_(0,I23,k,j,i);
    const Real &g_30 = metric_face3_kji_(0,I03,k,j,i);
    const Real &g_31 = metric_face3_kji_(0,I13,k,j,i);
    const Real &g_32 = metric_face3_kji_(0,I23,k,j,i);
    const Real &g_33 = metric_face3_kji_(0,I33,k,j,i);

    // Extract global fluxes
    Real &j3 = flux(IDN,k,j,i);
    Real &t3_0 = flux(IEN,k,j,i);
    Real &t3_1 = flux(IM1,k,j,i);
    Real &t3_2 = flux(IM2,k,j,i);
    Real &t3_3 = flux(IM3,k,j,i);

    // Set fluxes
    j3 = m3_tm*jt + m3_x*jx;
    t3_0 = g_00*t30 + g_01*t31 + g_02*t32 + g_03*t33;
    t3_1 = g_10*t30 + g_11*t31 + g_12*t32 + g_13*t33;
    t3_2 = g_20*t30 + g_21*t31 + g_22*t32 + g_23*t33;
    t3_3 = g_30*t30 + g_31*t31 + g_32*t32 + g_33*t33;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED) {
      Real ftt = 0.0;
      Real fxt = bbx(i);
      Real fyt = cons(IBY,i);
      Real fzt = cons(IBZ,i);
      Real ftx = -bbx(i);
      Real fxx = 0.0;
      Real fyx = -ey(k,j,i);
      Real fzx = ez(k,j,i);
      Real f13 = m3_tm * (m1_tm*ftt + m1_x*fxt + m1_y*fyt)
                 + m3_x * (m1_tm*ftx + m1_x*fxx + m1_y*fyx);
      Real f23 = m3_tm * (m2_tm*ftt + m2_x*fxt + m2_y*fyt + m2_z*fzt)
                 + m3_x * (m2_tm*ftx + m2_x*fxx + m2_y*fyx + m2_z*fzx);
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

void GRUser::RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3, int k, int j, int i,
                             Real *pa0, Real *pa1, Real *pa2, Real *pa3) {
  // Extract metric coefficients
  const Real &g00 = metric_cell_kji_(1,I00,k,j,i);
  const Real &g01 = metric_cell_kji_(1,I01,k,j,i);
  const Real &g02 = metric_cell_kji_(1,I02,k,j,i);
  const Real &g03 = metric_cell_kji_(1,I03,k,j,i);
  const Real &g10 = metric_cell_kji_(1,I01,k,j,i);
  const Real &g11 = metric_cell_kji_(1,I11,k,j,i);
  const Real &g12 = metric_cell_kji_(1,I12,k,j,i);
  const Real &g13 = metric_cell_kji_(1,I13,k,j,i);
  const Real &g20 = metric_cell_kji_(1,I02,k,j,i);
  const Real &g21 = metric_cell_kji_(1,I12,k,j,i);
  const Real &g22 = metric_cell_kji_(1,I22,k,j,i);
  const Real &g23 = metric_cell_kji_(1,I23,k,j,i);
  const Real &g30 = metric_cell_kji_(1,I03,k,j,i);
  const Real &g31 = metric_cell_kji_(1,I13,k,j,i);
  const Real &g32 = metric_cell_kji_(1,I23,k,j,i);
  const Real &g33 = metric_cell_kji_(1,I33,k,j,i);

  // Set raised components
  *pa0 = g00*a_0 + g01*a_1 + g02*a_2 + g03*a_3;
  *pa1 = g10*a_0 + g11*a_1 + g12*a_2 + g13*a_3;
  *pa2 = g20*a_0 + g21*a_1 + g22*a_2 + g23*a_3;
  *pa3 = g30*a_0 + g31*a_1 + g32*a_2 + g33*a_3;
  return;
}

//----------------------------------------------------------------------------------------
// Function for lowering contravariant components of a vector
// Inputs:
//   a0,a1,a2,a3: contravariant components of vector
//   k,j,i: indices of cell in which transformation is desired
// Outputs:
//   pa_0,pa_1,pa_2,pa_3: pointers to covariant 4-vector components

void GRUser::LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j, int i,
                             Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3) {
  // Extract metric coefficients
  const Real &g_00 = metric_cell_kji_(0,I00,k,j,i);
  const Real &g_01 = metric_cell_kji_(0,I01,k,j,i);
  const Real &g_02 = metric_cell_kji_(0,I02,k,j,i);
  const Real &g_03 = metric_cell_kji_(0,I03,k,j,i);
  const Real &g_10 = metric_cell_kji_(0,I01,k,j,i);
  const Real &g_11 = metric_cell_kji_(0,I11,k,j,i);
  const Real &g_12 = metric_cell_kji_(0,I12,k,j,i);
  const Real &g_13 = metric_cell_kji_(0,I13,k,j,i);
  const Real &g_20 = metric_cell_kji_(0,I02,k,j,i);
  const Real &g_21 = metric_cell_kji_(0,I12,k,j,i);
  const Real &g_22 = metric_cell_kji_(0,I22,k,j,i);
  const Real &g_23 = metric_cell_kji_(0,I23,k,j,i);
  const Real &g_30 = metric_cell_kji_(0,I03,k,j,i);
  const Real &g_31 = metric_cell_kji_(0,I13,k,j,i);
  const Real &g_32 = metric_cell_kji_(0,I23,k,j,i);
  const Real &g_33 = metric_cell_kji_(0,I33,k,j,i);

  // Set lowered components
  *pa_0 = g_00*a0 + g_01*a1 + g_02*a2 + g_03*a3;
  *pa_1 = g_10*a0 + g_11*a1 + g_12*a2 + g_13*a3;
  *pa_2 = g_20*a0 + g_21*a1 + g_22*a2 + g_23*a3;
  *pa_3 = g_30*a0 + g_31*a1 + g_32*a2 + g_33*a3;
  return;
}

namespace {
//----------------------------------------------------------------------------------------
// Functions for calculating determinant
// Inputs:
//   g: array of covariant metric coefficients
//   a11,a12,a13,a21,a22,a23,a31,a32,a33: elements of matrix
//   a11,a12,a21,a22: elements of matrix
// Outputs:
//   returned value: determinant

Real Determinant(const AthenaArray<Real> &g) {
  const Real &a11 = g(I00);
  const Real &a12 = g(I01);
  const Real &a13 = g(I02);
  const Real &a14 = g(I03);
  const Real &a21 = g(I01);
  const Real &a22 = g(I11);
  const Real &a23 = g(I12);
  const Real &a24 = g(I13);
  const Real &a31 = g(I02);
  const Real &a32 = g(I12);
  const Real &a33 = g(I22);
  const Real &a34 = g(I23);
  const Real &a41 = g(I03);
  const Real &a42 = g(I13);
  const Real &a43 = g(I23);
  const Real &a44 = g(I33);
  Real det = a11 * Determinant(a22, a23, a24, a32, a33, a34, a42, a43, a44)
             - a12 * Determinant(a21, a23, a24, a31, a33, a34, a41, a43, a44)
             + a13 * Determinant(a21, a22, a24, a31, a32, a34, a41, a42, a44)
             - a14 * Determinant(a21, a22, a23, a31, a32, a33, a41, a42, a43);
  return det;
}

Real Determinant(Real a11, Real a12, Real a13, Real a21, Real a22, Real a23,
                 Real a31, Real a32, Real a33) {
  Real det = a11 * Determinant(a22, a23, a32, a33)
             - a12 * Determinant(a21, a23, a31, a33)
             + a13 * Determinant(a21, a22, a31, a32);
  return det;
}

Real Determinant(Real a11, Real a12, Real a21, Real a22) {
  return a11 * a22 - a12 * a21;
}

//----------------------------------------------------------------------------------------
// Function for calculating frame transformation coefficients
// Inputs:
//   g,g_inv: arrays of covariant and contravariant metric coefficients
//   face: 1, 2, or 3 depending on which face is considered
// Outputs:
//   transformation: array of transformation coefficients

void CalculateTransformation(
    const AthenaArray<Real> &g,
    const AthenaArray<Real> &g_inv, int face, AthenaArray<Real> &transformation) {
  // Prepare indices
  int index[4][4];
  index[0][0] = I00;
  index[0][1] = I01; index[0][2] = I02; index[0][3] = I03;
  index[1][1] = I11; index[1][2] = I12; index[1][3] = I13;
  index[2][1] = I12; index[2][2] = I22; index[2][3] = I23;
  index[3][1] = I13; index[3][2] = I23; index[3][3] = I33;

  // Shift indices according to face
  int i0 = 0;
  int i1 = face;
  int i2 = 1 + face%3;
  int i3 = 1 + (face+1)%3;

  // Extract metric coefficients
  const Real &g_22 = g(index[i2][i2]);
  const Real &g_23 = g(index[i2][i3]);
  const Real &g_33 = g(index[i3][i3]);
  const Real &g00 = g_inv(index[i0][i0]);
  const Real &g01 = g_inv(index[i0][i1]);
  const Real &g02 = g_inv(index[i0][i2]);
  const Real &g03 = g_inv(index[i0][i3]);
  const Real &g11 = g_inv(index[i1][i1]);
  const Real &g12 = g_inv(index[i1][i2]);
  const Real &g13 = g_inv(index[i1][i3]);

  // Calculate intermediate quantities
  Real aa = -1.0 / std::sqrt(-g00);
  Real bb = 1.0 / std::sqrt(g00 * (g00 * g11 - SQR(g01)));
  Real cc = 1.0 / std::sqrt(g_33);
  Real dd = 1.0 / std::sqrt(g_33 * (g_22 * g_33 - SQR(g_23)));
  Real ee = g01 * g12 - g11 * g02;
  Real ff = g01 * g02 - g00 * g12;
  Real gg = g01 * g13 - g11 * g03;
  Real hh = g01 * g03 - g00 * g13;
  Real ii = SQR(bb)/cc * g00 * (gg + ee * g_23/g_33);
  Real jj = SQR(bb)/cc * g00 * (hh + ff * g_23/g_33);

  // Set local-to-global transformation coefficients
  transformation(0,T00) = aa * g00;
  transformation(0,T10) = aa * g01;
  transformation(0,T20) = aa * g02;
  transformation(0,T30) = aa * g03;
  transformation(0,T11) = bb * (g01 * g01 - g00 * g11);
  transformation(0,T21) = bb * (g01 * g02 - g00 * g12);
  transformation(0,T31) = bb * (g01 * g03 - g00 * g13);
  transformation(0,T22) = dd * g_33;
  transformation(0,T32) = -dd * g_23;
  transformation(0,T33) = cc;

  // Set global-to-local transformation coefficients
  transformation(1,T00) = -aa;
  transformation(1,T10) = bb * g01;
  transformation(1,T11) = -bb * g00;
  transformation(1,T20) = SQR(bb)*ee/dd * g00/g_33;
  transformation(1,T21) = SQR(bb)*ff/dd * g00/g_33;
  transformation(1,T22) = 1.0 / (dd * g_33);
  transformation(1,T30) = ii;
  transformation(1,T31) = jj;
  transformation(1,T32) = 1.0/cc * g_23/g_33;
  transformation(1,T33) = 1.0/cc;
  return;
}
} // namespace
