//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file emission.cpp
//  \brief coupling of radiation to matter via optically thin emission

// Athena++ headers
#include "../radiation.hpp"
#include "../../athena_arrays.hpp"  // AthenaArray

//----------------------------------------------------------------------------------------
// Function for calculating effect on radiation from optically thin emission
// Inputs:
//   prim_hydro: primitive hydro variables:
//     index 0: density (IDN) or pressure (IPR)
//     indices 1-3: k, j, i
//   n: unit normal null vector:
//     index 0: fluid-frame component (0 through 3)
//     index 1: angle (0 through nzeta * npsi - 1)
//     index 2: i
//   omega: fractional solid angle (normalized to 1) in fluid frame:
//     index 0: angle (0 through nzeta * npsi - 1)
//     index 1: i
//   dt: time interval over which coupling should be applied
//   k, j: x3- and x2-indices
//   intensity: fluid-frame I:
//     index 0: angle (0 through nzeta * npsi - 1)
//     index 1: i
// Outputs:
//   intensity: fluid-frame I updated
// Notes:
//   Implements emission with no absorption or scattering.
void Radiation::Coupling(const AthenaArray<Real> &prim_hydro, const AthenaArray<Real> &n,
    const AthenaArray<Real> &omega, Real dt, int k, int j, AthenaArray<Real> &intensity) {
  Real factor = 1.0/PI * sigma_sb * kappa * dt;
  for (int n = 0; n < nzeta * npsi; ++n) {
    for (int i = is; i <= ie; ++i) {
      Real rho = prim_hydro(IDN,k,j,i);
      Real pgas = prim_hydro(IPR,k,j,i);
      Real temperature = pgas / rho;
      intensity(n,i) += factor * rho * SQR(SQR(temperature));
    }
  }
  return;
}
