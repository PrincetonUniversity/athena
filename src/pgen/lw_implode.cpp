//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file lw_implode.cpp
//! \brief Problem generator for square implosion problem
//!
//! REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)
//========================================================================================

// C headers

// C++ headers

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
//  \brief Liska & Wendroff implosion test problem generator
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real d_in = pin->GetReal("problem","d_in");
  Real p_in = pin->GetReal("problem","p_in");

  Real d_out = pin->GetReal("problem","d_out");
  Real p_out = pin->GetReal("problem","p_out");

  Real gm1 = peos->GetGamma() - 1.0;

  // to make sure the ICs are symmetric, set y0 to be in between cell centers
  Real y0 = 0.5*(pmy_mesh->mesh_size.x2max + pmy_mesh->mesh_size.x2min);
  for (int j=js; j<=je; j++) {
    if (pcoord->x2v(j) > y0) {
      // TODO(felker): check this condition for multi-meshblock setups
      // further adjust y0 to be between cell center and lower x2 face
      y0 = pcoord->x2f(j) + 0.5*pcoord->dx2f(j);
      break;
    }
  }

  // Set initial conditions
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (pcoord->x2v(j) > (y0 - pcoord->x1v(i))) {
          phydro->u(IDN,k,j,i) = d_out;
          phydro->u(IEN,k,j,i) = p_out/gm1;
        } else {
          phydro->u(IDN,k,j,i) = d_in;
          phydro->u(IEN,k,j,i) = p_in/gm1;
        }
      }
    }
  }

  return;
}
