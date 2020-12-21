//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file twoibw.cpp
//! \brief Problem generator for two interacting blast waves test.
//!
//! REFERENCE: P. Woodward & P. Colella, "The numerical simulation of two-dimensional
//!   fluid flow with strong shocks", JCP, 54, 115, sect. IVa
//========================================================================================

// C headers

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the two interacting blast waves
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  std::stringstream msg;

  // parse shock direction: {1,2,3} -> {x1,x2,x3}
  int shk_dir = pin->GetInteger("problem","shock_dir");
  if (shk_dir < 1 || shk_dir > 3) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "shk_dir = " << shk_dir << " must be either 1,2, or 3" << std::endl;
    ATHENA_ERROR(msg);
  }

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if ((shk_dir==1 && pcoord->x1v(i) < 0.1) ||
            (shk_dir==2 && pcoord->x2v(j) < 0.1) ||
            (shk_dir==3 && pcoord->x3v(k) < 0.1)) {
          phydro->u(IEN,k,j,i)= 1.0e3/(peos->GetGamma() - 1.0);
        } else if ((shk_dir==1 && pcoord->x1v(i) > 0.9) ||
                   (shk_dir==2 && pcoord->x2v(j) > 0.9) ||
                   (shk_dir==3 && pcoord->x3v(k) > 0.9)) {
          phydro->u(IEN,k,j,i)= 1.0e2/(peos->GetGamma() - 1.0);
        } else {
          phydro->u(IEN,k,j,i)= 0.01/(peos->GetGamma() - 1.0);
        }
      }
    }
  }

  return;
}
