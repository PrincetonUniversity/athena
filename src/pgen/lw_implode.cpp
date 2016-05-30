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
//! \file lw_implode.cpp
//  \brief Problem generator for square implosion problem
//
// REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)    */
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/eos/eos.hpp"
#include "../coordinates/coordinates.hpp"


#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Liska & Wendroff implosion test problem generator
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real d_in = pin->GetReal("problem","d_in");
  Real p_in = pin->GetReal("problem","p_in");

  Real d_out = pin->GetReal("problem","d_out");
  Real p_out = pin->GetReal("problem","p_out");

  Real gm1 = (phydro->peos->GetGamma() - 1.0);
  Real y0 = 0.5*(pmy_mesh->mesh_size.x2max + pmy_mesh->mesh_size.x2min);

// Set initial conditions
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if(pcoord->x2v(j) > (y0 - pcoord->x1v(i))) {
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
