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
// \!fn void Z4c::GaugePreCollapsedLapse(AthenaArray<Real> & u)
// \brief Initialize precollapsed lapse and zero shift for single punture evolution

void Z4c::GaugePreCollapsedLapse(AthenaArray<Real> & u_adm, AthenaArray<Real> & u)
{
  ADM_vars adm;
  SetADMAliases(u_adm, adm);
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.Fill(0.);

  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  GLOOP2(k,j) {
    GLOOP1(i) {
      z4c.alpha(k,j,i) = std::pow(adm.psi4(k,j,i),-0.5);
    }
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugeGeodesic(AthenaArray<Real> & u)
// \brief Initialize lapse to 1 and shift to 0

void Z4c::GaugeGeodesic(AthenaArray<Real> & u)
{
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.ZeroClear();
}

