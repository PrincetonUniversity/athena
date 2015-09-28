//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file two_ibw.c
//  \brief Problem generator for two interacting blast waves test.
//
// REFERENCE: P. Woodward & P. Colella, "The numerical simulation of two-dimensional
//   fluid flow with strong shocks", JCP, 54, 115, sect. IVa
//======================================================================================

// C/C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

void Mesh::ProblemGenerator(Hydro *pfl, Field *pfd, ParameterInput *pin)
{
  MeshBlock *pmb = pfl->pmy_block;
  Coordinates *pco = pmb->pcoord;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  std::stringstream msg;

// parse shock direction: {1,2,3} -> {x1,x2,x3}

  int shk_dir = pin->GetInteger("problem","shock_dir");
  if (shk_dir < 1 || shk_dir > 3) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "shk_dir = " << shk_dir << " must be either 1,2, or 3" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
#pragma simd
    for (int i=is; i<=ie; ++i) {
      pfl->u(IDN,k,j,i) = 1.0;
      pfl->u(IM1,k,j,i) = 0.0;
      pfl->u(IM2,k,j,i) = 0.0;
      pfl->u(IM3,k,j,i) = 0.0;
      if ((shk_dir==1 && pco->x1v(i) < 0.1) ||
          (shk_dir==2 && pco->x2v(j) < 0.1) ||
          (shk_dir==3 && pco->x3v(k) < 0.1)) {
        pfl->u(IEN,k,j,i)= 1.0e3/(pfl->pf_eos->GetGamma() - 1.0);
      }
      else if ((shk_dir==1 && pco->x1v(i) > 0.9) ||
               (shk_dir==2 && pco->x2v(j) > 0.9) ||
               (shk_dir==3 && pco->x3v(k) > 0.9)) {
        pfl->u(IEN,k,j,i)= 1.0e2/(pfl->pf_eos->GetGamma() - 1.0);
      }
      else {
        pfl->u(IEN,k,j,i)= 0.01/(pfl->pf_eos->GetGamma() - 1.0);
      }
    }
  }}

  return;
}
