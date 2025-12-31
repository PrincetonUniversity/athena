//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file orbital_remapping.cpp
//! \brief functions for remapping in the orbital direction

// C/C++ headers
#include <algorithm>  // min()
#include <cmath>      // abs
#include <iostream>   // cout, endl
#include <sstream>    //
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

// this class header
#include "./orbital_advection.hpp"

//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::RemapFluxPlm(AthenaArray<Real> &pflux_,
//!                                         const AthenaArray<Real> &pbuf_,
//!                                         const Real eps_, const int osgn_,
//!                                         const int k, const int j, const int il,
//!                                         const int iu, const int shift_)
//! \brief Remap & Calculate flux with plm

void OrbitalAdvection::RemapFluxPlm(AthenaArray<Real> &pflux_,
                                    const AthenaArray<Real> &pbuf_,
                                    const Real eps_, const int osgn_, const int k,
                                    const int j, const int il, const int iu,
                                    const int shift_) {
  const Real odi     = osgn_-0.5;
  const Real coeff   = odi*(1.0-2.0*odi*eps_);
#pragma omp simd simdlen(SIMD_WIDTH)
  for(int i = il; i <= iu; i++) {
    Real dul  = pbuf_(k,j,i+shift_)-pbuf_(k,j,i+shift_-1);
    Real dur  = pbuf_(k,j,i+shift_+1)-pbuf_(k,j,i+shift_);
    Real du2  = dul*dur;
    Real dum  = 2.0*du2/(dul+dur);
    if (du2 <= 0.0) dum = 0.0;
    pflux_(i) = eps_*(pbuf_(k,j,i+shift_)+coeff*dum);
  }
  return;
}

void OrbitalAdvection::RemapFluxPlm(AthenaArray<Real> &pflux_,
                                    const AthenaArray<Real> &pbuf_,
                                    const Real eps_, const int osgn_, const int k,
                                    const int j, const int il, const int iu,
                                    const int nl, const int nu, const int shift_) {
  const Real odi     = osgn_-0.5;
  const Real coeff   = odi*(1.0-2.0*odi*eps_);

  for (int i = il; i <= iu; i++) {
    for (int n=nl; n<=nu; ++n) {
      Real dul  = pbuf_(k,j,i+shift_,n)-pbuf_(k,j,i+shift_-1,n);
      Real dur  = pbuf_(k,j,i+shift_+1,n)-pbuf_(k,j,i+shift_,n);
      Real du2  = dul*dur;
      Real dum  = 2.0*du2/(dul+dur);
      if (du2 <= 0.0) dum = 0.0;
      pflux_(i,n) = eps_*(pbuf_(k,j,i+shift_,n)+coeff*dum);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::RemapFluxPpm(AthenaArray<Real> &pflux_,
//!                                         AthenaArray<Real> &pbuf_,
//!                                         const Real eps_, const int osgn_,
//!                                         const int k, const int j, const int il,
//!                                         const int iu, const int shift_)
//! \brief Remap & Calculate flux with ppm

void OrbitalAdvection::RemapFluxPpm(AthenaArray<Real> &pflux_,
                                    AthenaArray<Real> &pbuf_,
                                    const Real eps_, const int osgn_,
                                    const int k, const int j, const int il,
                                    const int iu, const int shift_) {
  const Real odi = 2.0*(osgn_-0.5);
  const Real qx = odi*TWO_3RD*eps_;

//--- Step 1. -------------------------------------------------------------------------
#pragma omp simd
  for (int i = il; i <= iu; i++) {
    Real q_m2_i = pbuf_(k, j, i+shift_-2);
    Real q_m1_i = pbuf_(k, j, i+shift_-1);
    Real q_i    = pbuf_(k, j, i+shift_);
    Real q_p1_i = pbuf_(k, j, i+shift_+1);
    Real q_p2_i = pbuf_(k, j, i+shift_+2);

    Real qa = q_i - q_m1_i;
    Real qb = q_p1_i - q_i;
    Real dd_m1_i = 0.5*(qa + q_m1_i - q_m2_i);
    Real dd_i    = 0.5*(qb + qa);
    Real dd_p1_i = 0.5*(q_p2_i - q_p1_i + qb);

    Real dph_i    = 0.5*(q_m1_i + q_i) + (dd_m1_i - dd_i)/6.0;
    Real dph_p1_i = 0.5*(q_i + q_p1_i) + (dd_i - dd_p1_i)/6.0;
    
    // limiting - strict monotonicity
    dph_i    = std::min(dph_i,    std::max(q_i, q_m1_i));
    dph_p1_i = std::min(dph_p1_i, std::max(q_i, q_p1_i));
    dph_i    = std::max(dph_i,    std::min(q_i, q_m1_i));
    dph_p1_i = std::max(dph_p1_i, std::min(q_i, q_p1_i));

    Real qminus_i = dph_i;
    Real qplus_i  = dph_p1_i;

    Real dqf_minus_i = q_i - qminus_i;
    Real dqf_plus_i  = qplus_i - q_i;

    if (dqf_minus_i*dqf_plus_i <= 0.0) { // extrema
      qminus_i = qplus_i = q_i;
    } else { // no extrema
      qminus_i = q_i - 2.0*dqf_plus_i;
      qplus_i  = q_i + 2.0*dqf_minus_i;
    }

    Real dU = qplus_i-qminus_i;
    Real U6 = 6.0*(q_i - 0.5*(qplus_i + qminus_i));
    pflux_(i) = eps_*(osgn_*qplus_i + (1 - osgn_)*qminus_i
               -0.5*eps_*(dU - odi*(1.0 - qx)*U6));
  }
  return;
}

void OrbitalAdvection::RemapFluxPpm(AthenaArray<Real> &pflux_,
                                    AthenaArray<Real> &pbuf_,
                                    const Real eps_, const int osgn_,
                                    const int k, const int j, const int il,
                                    const int iu, const int nl, const int nu,
                                    const int shift_) {
  // will implement later
  return;
}
