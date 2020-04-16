//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file puncture_z4c.cpp
//  \brief implementation of functions in the Z4c class for initializing puntures evolution

// C++ standard headers
#include <cmath> // pow

// Athena++ headers
#include "z4c.hpp"
#include "z4c_macro.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"

//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMOnePuncture(AthenaArray<Real> & u)
// \brief Initialize ADM vars to single puncture (no spin)

void Z4c::ADMOnePuncture(AthenaArray<Real> & u_adm) {
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  // Flat spacetime
  ADMMinkowski(u_adm);

  GLOOP2(k,j) {
    // Isotropic radius
    GLOOP1(i) {
      r(i) = std::sqrt(SQR(mbi.x1(i)) + SQR(mbi.x2(j)) + SQR(mbi.x3(k)));
    }
    // psi4
    GLOOP1(i) {
      adm.psi4(k,j,i) = std::pow(1.0+0.5*opt.punc_ADM_mass/r(i),4);
    }
    // g_ab
    for(int a = 0; a < NDIM; ++a)
      for(int b = a; b < NDIM; ++b) {
        GLOOP1(i) {
          adm.g_dd(a,b,k,j,i) *= adm.psi4(k,j,i);
        }
      }
  }

}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugePreCollapsedLapse(AthenaArray<Real> & u)
// \brief Initialize precollapsed lapse and zero shift for single punture evolution

void Z4c::GaugePreCollapsedLapse(AthenaArray<Real> & u) {
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  GLOOP2(k,j) {
    // Isotropic radius
    GLOOP1(i) {
      r(i) = std::sqrt(SQR(mbi.x1(i)) + SQR(mbi.x2(j)) + SQR(mbi.x3(k)));
    }
    // lapse
    GLOOP1(i) {
      z4c.alpha(k,j,i) = 1.0/std::pow(1.0+0.5*opt.punc_ADM_mass/r(i),2);
    }
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMOnePunctureSpin(AthenaArray<Real> & u)
// \brief Initialize ADM vars to single puncture with spin

void Z4c::ADMOnePunctureSpin(AthenaArray<Real> & u_adm) {
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  // TODO: here we want to readin from input file and interp.

}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMTwoPunctures(AthenaArray<Real> & u)
// \brief Initialize ADM vars to two punctures

void Z4c::ADMTwoPunctures(AthenaArray<Real> & u_adm) {
  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  // TODO: here we want to readin from input file and interp.

}
