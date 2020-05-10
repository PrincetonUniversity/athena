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

void Z4c::ADMOnePuncture(ParameterInput *pin, AthenaArray<Real> & u_adm)
{
  ADM_vars adm;
  SetADMAliases(u_adm, adm);
  Real ADM_mass = pin->GetOrAddReal("problem", "punc_ADM_mass", 1.);
   
  MeshBlock * pmb = pmy_block;
  Coordinates * pco = pmb->pcoord;

  // Flat spacetime
  ADMMinkowski(u_adm);

  GLOOP2(k,j) {
    // Isotropic radius
    GLOOP1(i) {
      r(i) = std::sqrt(SQR(pco->x1v(i)) + SQR(pco->x2v(j)) + SQR(pco->x3v(k)));
    }
    // psi4
    GLOOP1(i) {
      adm.psi4(k,j,i) = std::pow(1.0+0.5*ADM_mass/r(i),4);
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
