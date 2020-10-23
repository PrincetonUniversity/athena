//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of setup-related functions in class Radiation

// C++ headers
#include <algorithm>  // max
#include <cmath>      // cos, hypot, sin, sqrt

// Athena++ headers
#include "radiation.hpp"
#include "../athena.hpp"                   // Real, indices
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../mesh/mesh.hpp"                // MeshBlock

//----------------------------------------------------------------------------------------
// Function for calculating conserved intensity corresponding to beam source
// Inputs:
//   pos_1, pos_2, pos_3: coordinates of beam origin
//   width: full proper diameter of beam
//   dir_1, dir_2, dir_3: relative directions of beam center
//   spread: full spread of beam in direction, in degrees
//   dii_dt: injected I per unit time
//   cylindrical: flag indicating coordinates are cylindrical
//   spherical: flag indicating coordinates are spherical
// Outputs:
//   dcons_dt: conserved values (n^0 n_0 I) per unit time set
// Notes:
//   Arrays should be 4D, with first index holding both zeta and psi.
//   Cylindrical coordinates:
//     phi (x2) will be mapped to [-pi, pi]:
//       Beams near phi = 0 should work fine.
//       Beams near phi = pi will likely not be initialized correctly.
//   Spherical coordinates:
//     theta (x2) will be mapped to [0, pi], adjusting phi (x3) as necessary.
//     phi (x3) will be mapped to [-pi, pi]:
//       Beams near phi = 0 should work fine.
//       Beams near phi = pi will likely not be initialized correctly.

void Radiation::CalculateBeamSource(Real pos_1, Real pos_2, Real pos_3, Real width,
    Real dir_1, Real dir_2, Real dir_3, Real spread, Real dii_dt,
    AthenaArray<Real> &dcons_dt, bool cylindrical, bool spherical) {

  // Account for cylindrical/spherical coordinates in beam origin
  if (cylindrical) {
    if (pos_2 > PI) {
      pos_2 -= 2.0*PI;
    }
  }
  if (spherical) {
    if (pos_2 < 0.0) {
      pos_2 = -pos_2;
      pos_3 -= PI;
    }
    if (pos_2 > PI) {
      pos_2 = 2.0*PI - pos_2;
      pos_3 -= PI;
    }
    if (pos_3 > PI) {
      pos_3 -= 2.0*PI;
    }
  }

  // Allocate scratch arrays
  AthenaArray<Real> g, gi, e, e_cov, omega, nh;
  g.NewAthenaArray(NMETRIC, ie + 1);
  gi.NewAthenaArray(NMETRIC, ie + 1);
  e.NewAthenaArray(4, 4);
  e_cov.NewAthenaArray(4, 4);
  omega.NewAthenaArray(4, 4, 4);
  nh.NewAthenaArray(ze + 1, pe + 1, 4);

  // Calculate unit normal components in orthonormal frame
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      nh(l,m,0) = 1.0;
      nh(l,m,1) = std::sin(zetav(l)) * std::cos(psiv(m));
      nh(l,m,2) = std::sin(zetav(l)) * std::sin(psiv(m));
      nh(l,m,3) = std::cos(zetav(l));
    }
  }

  // Calculate minimum angle between directions
  Real mu_min = std::cos(spread/2.0 * PI/180.0);

  // Go through cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pmy_block->pcoord->CellMetric(k, j, is, ie, g, gi);
      for (int i = is; i <= ie; ++i) {

        // Extract position, accounting for cylindrical/spherical coordinates
        Real x1 = pmy_block->pcoord->x1v(i);
        Real x2 = pmy_block->pcoord->x2v(j);
        Real x3 = pmy_block->pcoord->x3v(k);
        if (cylindrical) {
          if (x2 > PI) {
            x2 -= 2.0*PI;
          }
        }
        if (spherical) {
          if (x2 < 0.0) {
            x2 = -x2;
            x3 -= PI;
          }
          if (x2 > PI) {
            x2 = 2.0*PI - x2;
            x3 -= PI;
          }
          if (x3 > PI) {
            x3 -= 2.0*PI;
          }
        }

        // Calculate proper distance to beam origin
        Real dx1 = x1 - pos_1;
        Real dx2 = x2 - pos_2;
        Real dx3 = x3 - pos_3;
        Real dx_sq = g(I11,i) * SQR(dx1) + 2.0 * g(I12,i) * dx1 * dx2
            + 2.0 * g(I13,i) * dx1 * dx3 + g(I22,i) * SQR(dx2)
            + 2.0 * g(I23,i) * dx2 * dx3 + g(I33,i) * SQR(dx3);

        // Set to 0 if too far from beam in space
        if (dx_sq >= SQR(width/2.0)) {
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              dcons_dt(lm,k,j,i) = 0.0;
            }
          }
          continue;
        }

        // Calculate tetrad
        pmy_block->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);

        // Calculate contravariant time component of direction
        Real temp_a = g(I00,i);
        Real temp_b = 2.0 * (g(I01,i) * dir_1 + g(I02,i) * dir_2 + g(I03,i) * dir_3);
        Real temp_c = g(I11,i) * SQR(dir_1) + 2.0 * g(I12,i) * dir_1 * dir_2
            + 2.0 * g(I13,i) * dir_1 * dir_3 + g(I22,i) * SQR(dir_2)
            + 2.0 * g(I23,i) * dir_2 * dir_3 + g(I33,i) * SQR(dir_3);
        Real dir_0 =
            (-temp_b - std::sqrt(SQR(temp_b) - 4.0 * temp_a * temp_c)) / (2.0 * temp_a);

        // Calculate covariant direction
        Real dir_cov_0, dir_cov_1, dir_cov_2, dir_cov_3;
        pmy_block->pcoord->LowerVectorCell(dir_0, dir_1, dir_2, dir_3, k, j, i,
            &dir_cov_0, &dir_cov_1, &dir_cov_2, &dir_cov_3);

        // Calculate covariant direction in tetrad frame
        Real dir_tet_cov_0 = e(0,0) * dir_cov_0 + e(0,1) * dir_cov_1 + e(0,2) * dir_cov_2
            + e(0,3) * dir_cov_3;
        Real dir_tet_cov_1 = e(1,0) * dir_cov_0 + e(1,1) * dir_cov_1 + e(1,2) * dir_cov_2
            + e(1,3) * dir_cov_3;
        Real dir_tet_cov_2 = e(2,0) * dir_cov_0 + e(2,1) * dir_cov_1 + e(2,2) * dir_cov_2
            + e(2,3) * dir_cov_3;
        Real dir_tet_cov_3 = e(3,0) * dir_cov_0 + e(3,1) * dir_cov_1 + e(3,2) * dir_cov_2
            + e(3,3) * dir_cov_3;

        // Normalize covariant spatial direction in tetrad frame
        dir_tet_cov_1 /= -dir_tet_cov_0;
        dir_tet_cov_2 /= -dir_tet_cov_0;
        dir_tet_cov_3 /= -dir_tet_cov_0;

        // Go through angles
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {

            // Calculate combined angle index
            int lm = AngleInd(l, m);

            // Calculate angle to beam direction
            Real mu = nh(l,m,1) * dir_tet_cov_1 + nh(l,m,2) * dir_tet_cov_2
                + nh(l,m,3) * dir_tet_cov_3;

            // Set to 0 if too far from beam in angle
            if (mu <= mu_min) {
              dcons_dt(lm,k,j,i) = 0.0;
              continue;
            }

            // Set to nonzero value
            dcons_dt(lm,k,j,i) = dii_dt * n0_n_mu_(0,l,m,k,j,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating conserved intensity corresponding to spatially constant source
// Inputs:
//   energy: coordinate-frame energy density
//   u1, u2, u3: contravariant 4-velocity components of isotropic radiation frame
// Outputs:
//   cons_out: conserved values (n^0 n_0 I) set

void Radiation::CalculateConstantRadiation(Real energy, Real u1, Real u2, Real u3,
    AthenaArray<Real> &cons_out) {

  // Allocate scratch arrays
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, ie + 1);
  gi.NewAthenaArray(NMETRIC, ie + 1);

  // Go through cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pmy_block->pcoord->CellMetric(k, j, is, ie, g, gi);
      for (int i = is; i <= ie; ++i) {
        CalculateRadiationInCellM1(energy, u1, u2, u3, k, j, i, g, cons_out);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating conserved intensity from energy density and velocity
// Inputs:
//   energy: coordinate-frame energy density
//   u1, u2, u3: contravariant 4-velocity components of isotropic radiation frame
//   k, j, i: indices for cell to set
//   g: covariant metric
// Outputs:
//   cons_out: conserved values (n^0 n_0 I) set in given cell
// Notes:
//   Assumes intensity field associated with M1 closure.

void Radiation::CalculateRadiationInCellM1(Real energy, Real u1, Real u2, Real u3, int k,
    int j, int i, const AthenaArray<Real> &g, AthenaArray<Real> &cons_out) {

  // Calculate contravariant time component of isotropic radiation frame velocity
  Real temp_a = g(I00,i);
  Real temp_b = 2.0 * (g(I01,i) * u1 + g(I02,i) * u2 + g(I03,i) * u3);
  Real temp_c = g(I11,i) * SQR(u1) + 2.0 * g(I12,i) * u1 * u2
    + 2.0 * g(I13,i) * u1 * u3 + g(I22,i) * SQR(u2) + 2.0 * g(I23,i) * u2 * u3
    + g(I33,i) * SQR(u3) + 1.0;
  Real temp_d = std::max(SQR(temp_b) - 4.0 * temp_a * temp_c, 0.0);
  Real u0 =
      (-temp_b - std::sqrt(temp_d)) / (2.0 * temp_a);

  // Calculate covariant radiation velocity
  Real u_0, u_1, u_2, u_3;
  pmy_block->pcoord->LowerVectorCell(u0, u1, u2, u3, k, j, i, &u_0, &u_1, &u_2,
      &u_3);

  // Set conserved quantity at each angle, tracking energy density
  Real energy_sum = 0.0;
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      Real u_n = u_0 * nmu_(0,l,m,k,j,i) + u_1 * nmu_(1,l,m,k,j,i)
          + u_2 * nmu_(2,l,m,k,j,i) + u_3 * nmu_(3,l,m,k,j,i);
      Real ii = 1.0 / SQR(SQR(-u_n));
      cons_out(lm,k,j,i) = n0_n_mu_(0,l,m,k,j,i) * ii;
      energy_sum += SQR(nmu_(0,l,m,k,j,i)) * ii * solid_angle(l,m);
    }
  }

  // Normalize conserved quantities
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      cons_out(lm,k,j,i) *= energy / energy_sum;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating conserved intensity from energy density and velocity
// Inputs:
//   ee_f: fluid-frame energy density E_f
//   ff1_f, ff2_f, ff3_f: fluid-frame flux F_f^i
//   uu1, uu2, uu3: normal-frame fluid 4-velocity u^{i'}
//   k, j, i: indices for cell to set
//   g: covariant metric
// Outputs:
//   cons_out: conserved values (n^0 n_0 I) set in given cell
// Notes:
//   Assumes intensity field from Minerbo that is linear in cosine of angle with respect
//       to preferred direction.
//   Satisfies Eddington approximation in fluid frame so long as magnitude of F_f is less
//       than E_f/3.

void Radiation::CalculateRadiationInCellLinear(Real ee_f, Real ff1_f, Real ff2_f,
    Real ff3_f, Real uu1, Real uu2, Real uu3, int k, int j, int i,
    const AthenaArray<Real> &g, AthenaArray<Real> &cons_out) {

  // Calculate normalized flux in fluid frame
  Real ff_f = std::sqrt(SQR(ff1_f) + SQR(ff2_f) + SQR(ff3_f));
  Real f_f = 0.0;
  Real f1_f = 0.0;
  Real f2_f = 0.0;
  Real f3_f = 0.0;
  if (ee_f > 0.0 and ff_f > 0.0) {
    f_f = ff_f / ee_f;
    f1_f = ff1_f / ff_f;
    f2_f = ff2_f / ff_f;
    f3_f = ff3_f / ff_f;
  }

  // Calculate fluid velocity in tetrad frame
  Real temp_var = g(I11,i) * SQR(uu1) + 2.0 * g(I12,i) * uu1 * uu2
      + 2.0 * g(I13,i) * uu1 * uu3 + g(I22) * SQR(uu2)
      + 2.0 * g(I23,i) * uu2 * uu3 + g(I33,i) * SQR(uu3);
  Real uu0 = std::sqrt(1.0 + temp_var);
  Real u0_t = norm_to_tet_(0,0,k,j,i) * uu0 + norm_to_tet_(0,1,k,j,i) * uu1
      + norm_to_tet_(0,2,k,j,i) * uu2 + norm_to_tet_(0,3,k,j,i) * uu3;
  Real u1_t = norm_to_tet_(1,0,k,j,i) * uu0 + norm_to_tet_(1,1,k,j,i) * uu1
      + norm_to_tet_(1,2,k,j,i) * uu2 + norm_to_tet_(1,3,k,j,i) * uu3;
  Real u2_t = norm_to_tet_(2,0,k,j,i) * uu0 + norm_to_tet_(2,1,k,j,i) * uu1
      + norm_to_tet_(2,2,k,j,i) * uu2 + norm_to_tet_(2,3,k,j,i) * uu3;
  Real u3_t = norm_to_tet_(3,0,k,j,i) * uu0 + norm_to_tet_(3,1,k,j,i) * uu1
      + norm_to_tet_(3,2,k,j,i) * uu2 + norm_to_tet_(3,3,k,j,i) * uu3;

  // Go through each angle
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);

      // Calculate direction in fluid frame
      Real un_t = u1_t * nh_cc_(1,l,m) + u2_t * nh_cc_(2,l,m) + u3_t * nh_cc_(3,l,m);
      Real n0_f = u0_t * nh_cc_(0,l,m) - un_t;
      Real n1_f = -u1_t * nh_cc_(0,l,m) + u1_t / (u0_t + 1.0) * un_t + nh_cc_(1,l,m);
      Real n2_f = -u2_t * nh_cc_(0,l,m) + u2_t / (u0_t + 1.0) * un_t + nh_cc_(2,l,m);
      Real n3_f = -u3_t * nh_cc_(0,l,m) + u3_t / (u0_t + 1.0) * un_t + nh_cc_(3,l,m);

      // Calculate intensity in fluid frame
      Real fn_f = f1_f * n1_f + f2_f * n2_f + f3_f * n3_f;
      Real ii_f = 0.0;
      if (f_f <= 1.0/3.0) {
        ii_f = ee_f / (4.0*PI) * (1.0 + 3.0 * f_f * fn_f);
      } else {
        ii_f = ee_f / (9.0*PI) * (fn_f - 3.0 * f_f + 2.0) / SQR(1.0 - f_f);
      }

      // Calculate intensity in tetrad frame
      Real ii = ii_f / SQR(SQR(n0_f));

      // Store conserved quantity in output
      cons_out(lm,k,j,i) = n0_n_mu_(0,l,m,k,j,i) * ii;
    }
  }
  return;
}
