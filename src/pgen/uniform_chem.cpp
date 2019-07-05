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
//! \file uniform_chem.cpp
//  \brief problem generator, uniform mesh with chemistry
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
#include "../chemistry/thermo.hpp"
#include "../radiation/radiation.hpp"
#include "../field/field.hpp"
#include "../eos/eos.hpp"
#include "../coordinates/coordinates.hpp"


//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================
void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  FILE *pf = fopen("chem_network.dat", "w");
  pblock->pscalars->pchemnet->OutputProperties(pf);
  fclose(pf);
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief initialize problem by reading in vtk file.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  //dimensions of meshblock
  const int Nx = ie - is + 1;
  const int Ny = je - js + 1;
  const int Nz = ke - ks + 1;
	//read density and radiation field strength
	const Real nH = pin->GetReal("problem", "nH");
	const Real G0 = pin->GetReal("problem", "G0");
	const Real s_init = pin->GetReal("problem", "s_init");

	//set density
	for (int k=ks; k<=ke; ++k) {
		for (int j=js; j<=je; ++j) {
			for (int i=is; i<=ie; ++i) {
				phydro->u(IDN, k, j, i) = nH;
			}
		}
	}

	//intialize radiation field
	if (RADIATION_ENABLED) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
						for (int ifreq=0; ifreq < prad->nfreq; ++ifreq) {
							for (int iang=0; iang < prad->nang; ++iang) {
								prad->ir(k, j, i, ifreq * prad->nang + iang) = G0;
							}
						}
          }
        }
      }
	}

	//intialize chemical species
#ifdef INCLUDE_CHEMISTRY
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        for (int ispec=0; ispec < NSCALARS; ++ispec) {
          pscalars->s(ispec, k, j, i) = s_init;
        }
        //temperature
        pscalars->s(pscalars->pchemnet->iE_, k, j, i) = 20. *
          Thermo::CvCold(0., 0.1, 0.);
      }
    }
  }
#endif

  return;
}
