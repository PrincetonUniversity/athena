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
//! \file G14Sod.cpp
//  \brief problem generator, uniform mesh with chemistry
//  \ Only work with single core cpu
//  \ Ref: https://www.aanda.org/articles/aa/pdf/2016/02/aa27262-15.pdf
//======================================================================================

// C headers

// C++ headers
#include <algorithm>  // std::find()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // std::runtime_error()
#include <string>     // c_str()
#include <vector>     // vector container

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../chem_rad/chem_rad.hpp"
#include "../chem_rad/integrators/rad_integrators.hpp"
#include "../chemistry/utils/thermo.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief initialize problem
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Define the constants
  const Real mu = pin->GetReal("problem", "mu");
  const Real k_b = Constants::k_boltzmann_cgs;
  const Real m_h = Constants::hydrogen_mass_cgs;

  // read density and temperature
  const Real dl_cgs = pin->GetReal("problem", "LHS_rho_cgs");
  const Real dr_cgs = pin->GetReal("problem", "RHS_rho_cgs");
  const Real TL   = pin->GetReal("problem", "LHS_T_cgs");  // unit in K
  const Real TR   = pin->GetReal("problem", "RHS_T_cgs");

  const Real r_init = pin->GetOrAddReal("problem", "r_init", 0.);
  const Real gm1  = peos->GetGamma() - 1.0;;

  // calculate density and pressure in code units
  const Real dl = dl_cgs / pmy_mesh->punit->code_density_cgs;
  const Real dr = dr_cgs / pmy_mesh->punit->code_density_cgs;
  const Real Lpres = dl_cgs*TL*k_b/mu/m_h / pmy_mesh->punit->code_pressure_cgs;
  const Real Rpres = dr_cgs*TR*k_b/mu/m_h / pmy_mesh->punit->code_pressure_cgs;
  // parse shock location (must be inside grid)
  Real xshock = pin->GetReal("problem","xshock");

  // Initialize the discontinuity   ---------------------------------
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        if (pcoord->x1v(i) < xshock) {
          phydro->u(IDN,k,j,i) = dl;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;

          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = Lpres/gm1;
          }
        } else {
          phydro->u(IDN, k, j, i) = dr;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = Rpres/gm1;
          }
        }
      }
    }
  }
  // intialize chemical species
  if (NSCALARS > 0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ispec=0; ispec < NSCALARS; ++ispec) {
            pscalars->s(ispec, k, j, i) = r_init*phydro->u(IDN, k, j, i);
            if (CHEMISTRY_ENABLED) {
              Real s_ispec = pin->GetOrAddReal("problem",
                  "r_init_"+pscalars->chemnet.species_names[ispec], -1);
              if (s_ispec >= 0.) {
                pscalars->s(ispec, k, j, i) = s_ispec*phydro->u(IDN, k, j, i);
              }
            }
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
