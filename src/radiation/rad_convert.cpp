//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of conversion/inversion-related functions in class Radiation

// C++ headers
#include <cmath>  // sqrt

// Athena++ headers
#include "radiation.hpp"
#include "../athena.hpp"                   // Real, indices
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../mesh/mesh.hpp"                // MeshBlock

//----------------------------------------------------------------------------------------
// Radiation conversion from primitive to conserved variables
// Inputs:
//   prim_in: primitives
//   pcoord: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   cons_out: conserved quantities

void Radiation::PrimitiveToConserved(const AthenaArray<Real> &prim_in,
    AthenaArray<Real> &cons_out, Coordinates *pcoord, int il, int iu, int jl, int ju,
    int kl, int ku) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            cons_out(lm,k,j,i) = n0_n_mu_(0,l,m,k,j,i) * prim_in(lm,k,j,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Radiation inversion from conserved to primitive variables
// Inputs:
//   cons_in: conserved quantities
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   prim_out: primitives
// Notes:
//   Primitives are floored at 0. Conserved quantities are adjusted to match.
//   This should be the only place where angular ghost zones need to be set.

void Radiation::ConservedToPrimitive(AthenaArray<Real> &cons_in,
    AthenaArray<Real> &prim_out, int il, int iu, int jl, int ju, int kl, int ku) {

  // Calculate primitive intensities
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim_out(lm,k,j,i) = cons_in(lm,k,j,i) / n0_n_mu_(0,l,m,k,j,i);
            if (prim_out(lm,k,j,i) < 0.0) {
              prim_out(lm,k,j,i) = 0.0;
              cons_in(lm,k,j,i) = 0.0;
            }
          }
        }
      }
    }
  }

  // Populate angular ghost zones in azimuthal angle
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps-NGHOST; m <= ps-1; ++m) {
      int m_src = pe - ps + 1 + m;
      int lm = AngleInd(l, m);
      int lm_src = AngleInd(l, m_src);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim_out(lm,k,j,i) = prim_out(lm_src,k,j,i);
          }
        }
      }
    }
    for (int m = pe+1; m <= pe+NGHOST; ++m) {
      int m_src = ps - pe - 1 + m;
      int lm = AngleInd(l, m);
      int lm_src = AngleInd(l, m_src);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim_out(lm,k,j,i) = prim_out(lm_src,k,j,i);
          }
        }
      }
    }
  }

  // Populate angular ghost zones in polar angle
  for (int l = zs-NGHOST; l <= zs-1; ++l) {
    for (int m = ps-NGHOST; m <= pe+NGHOST; ++m) {
      int l_src = 2*zs - 1 - l;
      int m_src = (m + npsi/2) % (npsi + 2*NGHOST);
      int lm = AngleInd(l, m);
      int lm_src = AngleInd(l_src, m_src);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim_out(lm,k,j,i) = prim_out(lm_src,k,j,i);
          }
        }
      }
    }
  }
  for (int l = ze+1; l <= ze+NGHOST; ++l) {
    for (int m = ps-NGHOST; m <= pe+NGHOST; ++m) {
      int l_src = 2*ze + 1 - l;
      int m_src = (m + npsi/2) % (npsi + 2*NGHOST);
      int lm = AngleInd(l, m);
      int lm_src = AngleInd(l_src, m_src);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim_out(lm,k,j,i) = prim_out(lm_src,k,j,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Radiation inversion from conserved to primitive variables, including setting moments
// Inputs:
//   cons_in: conserved quantities
//   prim_hydro: up-to-date primitive hydro quantities
//   pcoord: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   prim_out: primitives
// Notes:
//   Defers to ConservedToPrimitive() for actual inversion.
//   Updates all relevant moments to match updated radiation variables.

void Radiation::ConservedToPrimitiveWithMoments(AthenaArray<Real> &cons_in,
    AthenaArray<Real> &prim_out, const AthenaArray<Real> &prim_hydro, Coordinates *pcoord,
    int il, int iu, int jl, int ju, int kl, int ku) {
  ConservedToPrimitive(cons_in, prim_out, il, iu, jl, ju, kl, ku);
  SetMoments(prim_hydro, pcoord, il, iu, jl, ju, kl, ku);
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating moments of radiation field
// Inputs:
//   prim_hydro: up-to-date primitive hydro quantities
//   pcoord: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region for which moments should be calculated
// Outputs: (none)
// Notes:
//   Populates moments_coord array with 10 components of radiation stress tensor in
//       coordinate frame.
//   Populates moments_tetrad array with 10 components of radiation stress tensor in
//       orthonormal tetrad frame.
//   Populates moments_fluid array with 10 components of radiation stress tensor in
//       orthonormal fluid frame if radiation is coupled to matter.

void Radiation::SetMoments(const AthenaArray<Real> &prim_hydro, Coordinates *pcoord,
    int il, int iu, int jl, int ju, int kl, int ku) {

  // Zero moments
  for (int n = 0; n < 10; ++n) {
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = il; i <= iu; ++i) {
          moments_coord(n,k,j,i) = 0.0;
          moments_tetrad(n,k,j,i) = 0.0;
          if (coupled_to_matter) {
            moments_fluid(n,k,j,i) = 0.0;
          }
        }
      }
    }
  }

  // Set coordinate-frame components
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int n1 = 0, n12 = 0; n1 < 4; ++n1) {
        for (int n2 = n1; n2 < 4; ++n2, ++n12) {
          for (int k = kl; k <= ku; ++k) {
            for (int j = jl; j <= ju; ++j) {
              for (int i = il; i <= iu; ++i) {
                moments_coord(n12,k,j,i) += nmu_(n1,l,m,k,j,i) * nmu_(n2,l,m,k,j,i)
                    * prim(lm,k,j,i) * solid_angle(l,m);
              }
            }
          }
        }
      }
    }
  }

  // Set tetrad-frame components
  AthenaArray<Real> e, e_cov, omega;
  e.NewAthenaArray(4, 4);
  e_cov.NewAthenaArray(4, 4);
  omega.NewAthenaArray(4, 4, 4);
  Real moments_coord_full[4][4];
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);
        moments_coord_full[0][0] = moments_coord(0,k,j,i);
        moments_coord_full[0][1] = moments_coord_full[1][0] = moments_coord(1,k,j,i);
        moments_coord_full[0][2] = moments_coord_full[2][0] = moments_coord(2,k,j,i);
        moments_coord_full[0][3] = moments_coord_full[3][0] = moments_coord(3,k,j,i);
        moments_coord_full[1][1] = moments_coord(4,k,j,i);
        moments_coord_full[1][2] = moments_coord_full[2][1] = moments_coord(5,k,j,i);
        moments_coord_full[1][3] = moments_coord_full[3][1] = moments_coord(6,k,j,i);
        moments_coord_full[2][2] = moments_coord(7,k,j,i);
        moments_coord_full[2][3] = moments_coord_full[3][2] = moments_coord(8,k,j,i);
        moments_coord_full[3][3] = moments_coord(9,k,j,i);
        for (int n1 = 0, n12 = 0; n1 < 4; ++n1) {
          for (int n2 = n1; n2 < 4; ++n2, ++n12) {
            for (int m1 = 0; m1 < 4; ++m1) {
              for (int m2 = 0; m2 < 4; ++m2) {
                moments_tetrad(n12,k,j,i) +=
                    e_cov(n1,m1) * e_cov(n2,m2) * moments_coord_full[m1][m2];
              }
            }
          }
        }
        moments_tetrad(1,k,j,i) *= -1.0;
        moments_tetrad(2,k,j,i) *= -1.0;
        moments_tetrad(3,k,j,i) *= -1.0;
      }
    }
  }

  // Set fluid-frame components
  if (coupled_to_matter) {
    Real tet_to_fluid[4][4];
    Real moments_tetrad_full[4][4];
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        pmy_block->pcoord->CellMetric(k, j, il, iu, g_, gi_);
        for (int i = il; i <= iu; ++i) {

          // Calculate fluid velocity in tetrad frame
          Real uu1 = prim_hydro(IVX,k,j,i);
          Real uu2 = prim_hydro(IVY,k,j,i);
          Real uu3 = prim_hydro(IVZ,k,j,i);
          Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
              + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22) * SQR(uu2)
              + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
          Real uu0 = std::sqrt(1.0 + temp_var);
          Real utet0 = norm_to_tet_(0,0,k,j,i) * uu0 + norm_to_tet_(0,1,k,j,i) * uu1
              + norm_to_tet_(0,2,k,j,i) * uu2 + norm_to_tet_(0,3,k,j,i) * uu3;
          Real utet1 = norm_to_tet_(1,0,k,j,i) * uu0 + norm_to_tet_(1,1,k,j,i) * uu1
              + norm_to_tet_(1,2,k,j,i) * uu2 + norm_to_tet_(1,3,k,j,i) * uu3;
          Real utet2 = norm_to_tet_(2,0,k,j,i) * uu0 + norm_to_tet_(2,1,k,j,i) * uu1
              + norm_to_tet_(2,2,k,j,i) * uu2 + norm_to_tet_(2,3,k,j,i) * uu3;
          Real utet3 = norm_to_tet_(3,0,k,j,i) * uu0 + norm_to_tet_(3,1,k,j,i) * uu1
              + norm_to_tet_(3,2,k,j,i) * uu2 + norm_to_tet_(3,3,k,j,i) * uu3;

          // Construct Lorentz boost from tetrad frame to orthonormal fluid frame
          tet_to_fluid[0][0] = utet0;
          tet_to_fluid[0][1] = tet_to_fluid[1][0] = -utet1;
          tet_to_fluid[0][2] = tet_to_fluid[2][0] = -utet2;
          tet_to_fluid[0][3] = tet_to_fluid[3][0] = -utet3;
          tet_to_fluid[1][1] = SQR(utet1) / (1.0 + utet0) + 1.0;
          tet_to_fluid[1][2] = tet_to_fluid[2][1] = utet1 * utet2 / (1.0 + utet0);
          tet_to_fluid[1][3] = tet_to_fluid[3][1] = utet1 * utet3 / (1.0 + utet0);
          tet_to_fluid[2][2] = SQR(utet2) / (1.0 + utet0) + 1.0;
          tet_to_fluid[2][3] = tet_to_fluid[3][2] = utet2 * utet3 / (1.0 + utet0);
          tet_to_fluid[3][3] = SQR(utet3) / (1.0 + utet0) + 1.0;

          // Transform moments
          moments_tetrad_full[0][0] = moments_tetrad(0,k,j,i);
          moments_tetrad_full[0][1] = moments_tetrad_full[1][0] = moments_tetrad(1,k,j,i);
          moments_tetrad_full[0][2] = moments_tetrad_full[2][0] = moments_tetrad(2,k,j,i);
          moments_tetrad_full[0][3] = moments_tetrad_full[3][0] = moments_tetrad(3,k,j,i);
          moments_tetrad_full[1][1] = moments_tetrad(4,k,j,i);
          moments_tetrad_full[1][2] = moments_tetrad_full[2][1] = moments_tetrad(5,k,j,i);
          moments_tetrad_full[1][3] = moments_tetrad_full[3][1] = moments_tetrad(6,k,j,i);
          moments_tetrad_full[2][2] = moments_tetrad(7,k,j,i);
          moments_tetrad_full[2][3] = moments_tetrad_full[3][2] = moments_tetrad(8,k,j,i);
          moments_tetrad_full[3][3] = moments_tetrad(9,k,j,i);
          for (int n1 = 0, n12 = 0; n1 < 4; ++n1) {
            for (int n2 = n1; n2 < 4; ++n2, ++n12) {
              for (int m1 = 0; m1 < 4; ++m1) {
                for (int m2 = 0; m2 < 4; ++m2) {
                  moments_fluid(n12,k,j,i) += tet_to_fluid[n1][m1] * tet_to_fluid[n2][m2]
                      * moments_tetrad_full[m1][m2];
                }
              }
            }
          }
        }
      }
    }
  }
  return;
}
