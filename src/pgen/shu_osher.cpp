//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file shu_osher.cpp
//! \brief Problem generator for Shu-Osher shocktube test, involving
//!  interaction of a Mach 3 shock with a sine wave density distribution.
//!
//! REFERENCE: C.W. Shu & S. Osher, "Efficient implementation of essentially
//!   non-oscillatory shock-capturing schemes, II", JCP, 83, 32 (1998)
//========================================================================================

// C headers

// C++ headers
#include <cmath>  // sin()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Shu-Osher test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // setup dependent variables
  Real dl = 3.857143;
  Real pl = 10.33333;
  Real ul = 2.629369;
  Real vl = 0.0;
  Real wl = 0.0;

  Real gm1 = (peos->GetGamma()) - 1.0;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        if (pcoord->x1v(i) < -0.8) {
          phydro->u(IDN,k,j,i) = dl;
          phydro->u(IM1,k,j,i) = ul*dl;
          phydro->u(IM2,k,j,i) = vl*dl;
          phydro->u(IM3,k,j,i) = wl*dl;
          phydro->u(IEN,k,j,i) = pl/gm1 + 0.5*dl*(ul*ul + vl*vl + wl*wl);
        } else {
          phydro->u(IDN,k,j,i) = 1.0 + 0.2*std::sin(5.0*PI*(pcoord->x1v(i)));
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          phydro->u(IEN,k,j,i) = 1.0/gm1;
        }
      }
    }
  }

  return;
}
