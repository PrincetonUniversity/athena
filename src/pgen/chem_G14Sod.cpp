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
//! \file chem_Z16Sod.cpp
//  \brief problem generator, uniform mesh with chemistry
//  \ Only work with single core cpu
//  \ Ref: https://www.aanda.org/articles/aa/pdf/2016/02/aa27262-15.pdf
//======================================================================================

// C++ headers
#include <string>     // c_str()
#include <iostream>   // endl
#include <vector>     // vector container
#include <sstream>    // stringstream
#include <stdio.h>    // c style file
#include <string.h>   // strcmp()
#include <algorithm>  // std::find()
#include <stdexcept>  // std::runtime_error()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../scalars/scalars.hpp"
#include "../chemistry/utils/thermo.hpp"
#include "../radiation/radiation.hpp"
#include "../radiation/integrators/rad_integrators.hpp"
#include "../field/field.hpp"
#include "../eos/eos.hpp"
#include "../coordinates/coordinates.hpp"
//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief initialize problem
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Define the constant
  const Real mu = pin->GetReal("hydro", "mu");
  const Real k_b = 1.380658e-16;
  const Real m_h = 1.6733e-24;
  const Real unit_E_cgs = punit->code_energydensity_cgs;
  const Real unit_d_cgs = punit->code_density_cgs;

  //read density and radiation field strength
  //const Real nH = pin->GetReal("problem", "nH");
  const Real dl = pin->GetReal("problem", "LHS_rho");
  const Real dr = pin->GetReal("problem", "RHS_rho");
  const Real vxl = pin->GetOrAddReal("problem","vxl", 0.);
  const Real vyl = pin->GetOrAddReal("problem","vyl", 0.);
  const Real vzl = pin->GetOrAddReal("problem","vzl", 0.);
  const Real vxr = pin->GetOrAddReal("problem","vxr", 0.);
  const Real vyr = pin->GetOrAddReal("problem","vyr", 0.);
  const Real vzr = pin->GetOrAddReal("problem","vzr", 0.);

  const Real TL   = pin->GetReal("problem", "LHS_T");  // unit in K
  const Real TR   = pin->GetReal("problem", "RHS_T");
  const Real s_init = pin->GetOrAddReal("problem", "s_init", 0.);
  const Real gm1  = peos->GetGamma() - 1.0;; 
  const Real Lpres = (dl*unit_d_cgs)*TL*k_b/mu/m_h;
  const Real Rpres = (dr*unit_d_cgs)*TR*k_b/mu/m_h;
  // parse shock location (must be inside grid)
  Real xshock = pin->GetReal("problem","xshock");

  // Initialize the discontinuity   ---------------------------------
  for (int k=ks; k<=ke; ++k) {
  	for (int j=js; j<=je; ++j) {
	  for (int i=is; i<=ie; ++i) {
		if (pcoord->x1v(i) < xshock){
			phydro->u(IDN,k,j,i) = dl;     
			phydro->u(IM1,k,j,i) = dl*vxl;
      phydro->u(IM2,k,j,i) = dl*vyl;
      phydro->u(IM3,k,j,i) = dl*vzl;

			if (NON_BAROTROPIC_EOS) {
        		phydro->u(IEN,k,j,i) = Lpres/gm1/unit_E_cgs;
          }
		}else{
			phydro->u(IDN, k, j, i) = dr; 
			phydro->u(IM1,k,j,i) = dr*vxr;
      phydro->u(IM2,k,j,i) = dr*vyr;
      phydro->u(IM3,k,j,i) = dr*vzr;
			if (NON_BAROTROPIC_EOS) {
	    		phydro->u(IEN,k,j,i) = Rpres/gm1/unit_E_cgs;
      }
		}
	  }
    }
  }
  //intialize chemical species
  if (NSCALARS > 0) {
    for (int k=ks; k<=ke; ++k) {
		  for (int j=js; j<=je; ++j) {
			  for (int i=is; i<=ie; ++i) {
		      for (int ispec=0; ispec < NSCALARS; ++ispec) {
		        pscalars->s(ispec, k, j, i) = s_init*phydro->u(IDN, k, j, i);
#ifdef INCLUDE_CHEMISTRY
		        Real s_ispec = pin->GetOrAddReal("problem",
		          "s_init_"+pscalars->chemnet.species_names[ispec], -1);
		        if (s_ispec >= 0.) {
		          pscalars->s(ispec, k, j, i) = s_ispec*phydro->u(IDN, k, j, i);
            }
#endif
        }
      }
  	  }
    }
  }

  return;
}

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  return;
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  return;
}