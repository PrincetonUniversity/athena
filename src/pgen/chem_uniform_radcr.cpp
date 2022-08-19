//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file chem_uniform_radcr.cpp
//! \brief problem generator, uniform chemistry and radiation
//======================================================================================

// c headers
#include <stdio.h>    // c style file
#include <string.h>   // strcmp()

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
#include "../chemistry/utils/thermo.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../radiation/integrators/rad_integrators.hpp"
#include "../radiation/radiation.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/units.hpp"


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief initialize problem of uniform chemistry and radiation
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  //dimensions of meshblock
  const int Nx = ie - is + 1;
  const int Ny = je - js + 1;
  const int Nz = ke - ks + 1;
  //read density and radiation field strength
  const Real nH = pin->GetReal("problem", "nH");
  const Real vx = pin->GetOrAddReal("problem", "vx_kms", 0);
  const Real s_init = pin->GetReal("problem", "s_init");
  const Real iso_cs = pin->GetReal("hydro", "iso_sound_speed");
  const Real pres = nH*SQR(iso_cs);
  const Real gm1  = peos->GetGamma() - 1.0;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        //density
        phydro->u(IDN, k, j, i) = nH;
        //velocity, x direction
        phydro->u(IM1, k, j, i) = nH*vx;
        //energy
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN, k, j, i) = pres/gm1;
        }
      }
    }
  }

  //intialize radiation field
  if (RADIATION_ENABLED) {
    const Real G0 = pin->GetReal("radiation", "G0");
    const Real cr_rate = pin->GetOrAddReal("radiation", "CR", 2e-16);
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ifreq=0; ifreq < prad->nfreq; ++ifreq) {
            for (int iang=0; iang < prad->nang; ++iang) {
              prad->ir(k, j, i, ifreq * prad->nang + iang) = G0;
            }
          }
#ifdef INCLUDE_CHEMISTRY
          for (int iang=0; iang < prad->nang; ++iang) {
            //cr rate
            prad->ir(k, j, i,
                pscalars->chemnet.index_cr_ * prad->nang + iang) = cr_rate;
          }
#endif
        }
      }
    }
    //calculate the average radiation field for output of the initial condition
    prad->pradintegrator->CopyToOutput();
  }

  //intialize chemical species
  if (NSCALARS > 0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ispec=0; ispec < NSCALARS; ++ispec) {
            pscalars->s(ispec, k, j, i) = s_init*nH;
#ifdef INCLUDE_CHEMISTRY
            Real s_ispec = pin->GetOrAddReal("problem",
                "s_init_"+pscalars->chemnet.species_names[ispec], -1);
            if (s_ispec >= 0.) {
              pscalars->s(ispec, k, j, i) = s_ispec*nH;
            }
#endif
          }
        }
      }
    }
  }

  return;
}
