//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file emission.cpp
//  \brief coupling of radiation to matter via optically thin emission

// C++ headers
#include <cmath>  // isnan, pow, sqrt

// Athena++ headers
#include "../radiation.hpp"
#include "../../athena.hpp"         // Real
#include "../../athena_arrays.hpp"  // AthenaArray
#include "../../defs.hpp"           // macros
#include "../../eos/eos.hpp"        // EquationofState
#include "../../mesh/mesh.hpp"      // MeshBlock

// Declarations
bool FourthPolyRoot(const Real coef4, const Real tconst, Real &root);

//----------------------------------------------------------------------------------------
// Function for calculating effect on radiation from optically thin emission
// Inputs:
//   prim_hydro: primitive hydro variables:
//     index 0: density (IDN) or pressure (IPR)
//     indices 1-3: k, j, i
//   normal: unit normal null vector n:
//     index 0: fluid-frame component (0 through 3)
//     index 1: angle (0 through nzeta * npsi - 1)
//     index 2: i
//   n0: unit normal null vector n, time component in coordinate frame:
//     index 1: angle (0 through nzeta * npsi - 1)
//     indices 1-3: k, j, i
//   omega: fractional solid angle (normalized to 1) in fluid frame:
//     index 0: angle (0 through nzeta * npsi - 1)
//     index 1: i
//   dt: coordinate time interval over which coupling should be applied:
//     index 0: i
//   dtau: fluid proper time interval over which coupling should be applied:
//     index 0: i
//   k, j: x3- and x2-indices
//   intensity: fluid-frame I:
//     index 0: angle (0 through nzeta * npsi - 1)
//     index 1: i
// Outputs:
//   intensity: fluid-frame I updated
// Notes:
//   Implements emission with no absorption or scattering.
//   n^0 is -u_\mu n^\mu, which is equal to gamma * (1 - v \dot n) in SR.
//   n^1/n^0 is the cosine of the angle the direction makes with the fluid-frame
//       x-direction, and similarly for n^2/n^0 and n^3/n^0.
//   dt and dtau are related according to dt = u^0 * dtau, where u^0 is the time component
//       of the fluid 4-velocity in the coordinate frame.

void Radiation::Coupling(const AthenaArray<Real> &prim_hydro,
    const AthenaArray<Real> &normal, const AthenaArray<Real> &n0,
    const AthenaArray<Real> &omega, const AthenaArray<Real> &dt,
    const AthenaArray<Real> &dtau, int k, int j, AthenaArray<Real> &intensity) {

  // Prepare to go through cells
  Real gamma = pmy_block->peos->GetGamma();
  int nang_local = nzeta * npsi;
  Real coef[2] = {};

  // Go through cells
  for (int i = is; i <= ie; ++i) {

    // Extract temperature
    Real rho = prim_hydro(IDN,k,j,i);
    Real pgas = prim_hydro(IPR,k,j,i);
    Real tgas = pgas / rho;
    Real tgasnew = tgas;

    // Extract angle-dependent terms
    for (int n = 0; n < nang_local; ++n) {
      intensity_scr_(n) = intensity(n,i);
      tran_coef_(n) = normal(0,n,i);
      weight_(n) = omega(n,i);
    }

    // Extract opacity
    Real sigma_s = opacity(OPAS,k,j,i);
    Real sigma_a = opacity(OPAA,k,j,i);
    Real sigma_p = 0.0;
    if (using_planck_mean) {
      sigma_p = opacity(OPAP,k,j,i);
    }
    Real dtaucsigmaa = dtau(i) * sigma_a;
    Real dtaucsigmas = dtau(i) * sigma_s;
    Real dtaucsigmap = dtau(i) * sigma_p;

    Real dtcsigmaa = dt(i) * sigma_a;
    Real dtcsigmas = dt(i) * sigma_s;
    Real dtcsigmap = dt(i) * sigma_p;

    // Calculate polynomial coefficients
    Real jr_cm = 0.0;
    Real suma1 = 0.0;
    Real suma2 = 0.0;
    for (int n = 0; n < nang_local; n++) {
       Real n0_local = n0(n,k,j,i);
       Real vncsigma = 1.0 / (n0_local + (dtcsigmaa + dtcsigmas) * tran_coef_(n));
       vncsigma2_(n) = tran_coef_(n) * vncsigma;
       Real ir_weight = intensity_scr_(n) * weight_(n);
       jr_cm += ir_weight;
       suma1 += weight_(n) * vncsigma2_(n);
       suma2 += ir_weight * n0_local * vncsigma;
    }
    Real suma3 = suma1 * (dtcsigmas - dtcsigmap);
    suma1 *= (dtcsigmaa + dtcsigmap);
    coef[1] = (dtaucsigmaa + dtaucsigmap - (dtaucsigmaa + dtaucsigmap) * suma1
        / (1.0 - suma3)) * arad * (gamma - 1.0) / rho;
    coef[0] = -tgas - (dtaucsigmaa + dtaucsigmap) * suma2 * (gamma - 1.0)
        / (rho * (1.0 - suma3));

    // Calculate new gas temperature
    bool badcell = false;
    if (fabs(coef[1]) > TINY_NUMBER) {
      bool flag = FourthPolyRoot(coef[1], coef[0], tgasnew);
      if (not flag or std::isnan(tgasnew)) {
        badcell = true;
        tgasnew = tgas;
      }
    } else {
      tgasnew = -coef[0];
    }

    // Apply changes to intensity
    if (not badcell) {

      // Calculate emission coefficient
      Real emission = arad * SQR(SQR(tgasnew));

      // Calculate updated jr_cm
      jr_cm = (suma1 * emission + suma2) / (1.0 - suma3);

      // Update the comoving frame specific intensity
      // Note: even if tr == told, intensity can change
      for (int n = 0; n < nang_local; n++) {
        intensity_scr_(n) += ((dtcsigmas - dtcsigmap) * jr_cm + (dtcsigmaa + dtcsigmap)
            * emission - (dtcsigmas + dtcsigmaa) * intensity_scr_(n)) * vncsigma2_(n);
      }

      // Copy intensity back 
      for (int n = 0; n < nang_local; ++n) {
        intensity(n,i) = intensity_scr_(n);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Exact solution for fourth order polynomial
// Inputs:
//   coef4: quartic coefficient
//   tconst: constant coefficient
// Outputs:
//   root: solution to equation
//   returned value: flag indicating success
// Notes:
//   Polynomial has the form coef4 * x^4 + x + tconst = 0.

bool FourthPolyRoot(const Real coef4, const Real tconst, Real &root) {

  // Calculate real root of z^3 - 4*tconst/coef4 * z - 1/coef4^2 = 0
  Real asquar = coef4 * coef4;
  Real acubic = coef4 * asquar;
  Real ccubic = tconst * tconst * tconst;
  Real delta1 = 0.25 - 64.0 * ccubic * coef4 / 27.0;
  if (delta1 < 0.0) {
    return false;
  }
  delta1 = std::sqrt(delta1);
  if (delta1 < 0.5) {
    return false;
  }
  Real zroot;
  if (delta1 > 1.0e11) {  // to avoid small number cancellation
    zroot = std::pow(delta1, -2.0/3.0) / 3.0;
  } else {
    zroot = std::pow(0.5 + delta1, 1.0/3.0) - std::pow(-0.5 + delta1, 1.0/3.0);
  }
  if (zroot < 0.0) {
    return false;
  }
  zroot *= std::pow(coef4, -2.0/3.0);

  // Calculate quartic root using cubic root
  Real rcoef = std::sqrt(zroot);
  Real delta2 = -zroot + 2.0 / (coef4 * rcoef);
  if (delta2 < 0.0) {
    return false;
  }
  delta2 = std::sqrt(delta2);
  root = 0.5 * (delta2 - rcoef);
  if (root < 0.0) {
    return false;
  }
  return true;
}
