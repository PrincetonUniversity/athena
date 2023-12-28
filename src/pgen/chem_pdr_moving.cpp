//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file chem_uniform_sixray.cpp
//! \brief problem generator, uniform chemistry and radiation with six-ray
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
#include "../bvals/bvals.hpp"
#include "../chem_rad/chem_rad.hpp"
#include "../chem_rad/integrators/rad_integrators.hpp"
#include "../chemistry/utils/thermo.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "../units/units.hpp"

// Radiation boundary
namespace {
  AthenaArray<Real> G0_iang;
  Real G0, cr_rate;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  G0 = pin->GetOrAddReal("chem_radiation", "G0", 0.);
  G0_iang.NewAthenaArray(6);
  G0_iang(BoundaryFace::inner_x1) = pin->GetOrAddReal("chem_radiation","G0_inner_x1",G0);
  G0_iang(BoundaryFace::inner_x2) = pin->GetOrAddReal("chem_radiation","G0_inner_x2",G0);
  G0_iang(BoundaryFace::inner_x3) = pin->GetOrAddReal("chem_radiation","G0_inner_x3",G0);
  G0_iang(BoundaryFace::outer_x1) = pin->GetOrAddReal("chem_radiation","G0_outer_x1",G0);
  G0_iang(BoundaryFace::outer_x2) = pin->GetOrAddReal("chem_radiation","G0_outer_x2",G0);
  G0_iang(BoundaryFace::outer_x3) = pin->GetOrAddReal("chem_radiation","G0_outer_x3",G0);
  cr_rate = pin->GetOrAddReal("chem_radiation", "CR", 2e-16);
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief initialize problem of uniform chemistry and radiation
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // read density and radiation field strength
  const Real vx = pin->GetReal("problem", "vx_kms");
  const Real r_init = pin->GetOrAddReal("problem", "r_init", 0.);
  // 2 phase initial condition
  const Real nc = 100.;
  const Real Tc = 40.;
  const Real nw = pin->GetOrAddReal("problem", "nw", 1e-1);
  //const Real Tw = pin->GetOrAddReal("problem", "Tw", 4e4);
  const Real cv = Thermo::CvCold(0.5, 0.1, 0.);
  const Real Eunit = pmy_mesh->punit->code_energydensity_cgs;
  const Real Eth = nc * Tc * cv / Eunit;
  const Real xc_start = 5.;
  const Real xc_end = 45.;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x1 = pcoord->x1v(i);
        if (x1 >= xc_start && x1 <= xc_end) {
          // density, cold
          phydro->u(IDN, k, j, i) = nc;
          // velocity, x direction, cold
          phydro->u(IM1, k, j, i) = nc*vx;
          // energy, cold
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN, k, j, i) = Eth + 0.5*nc*SQR(vx);
          }
        } else {
          // density, warm
          phydro->u(IDN, k, j, i) = nw;
          // velocity, x direction, warm
          phydro->u(IM1, k, j, i) = nw*vx;
          // energy, warm
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN, k, j, i) = Eth + 0.5*nw*SQR(vx);
          }
        }
      }
    }
  }

  // intialize radiation field
  if (CHEMRADIATION_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ifreq=0; ifreq < pchemrad->nfreq; ++ifreq) {
            for (int iang=0; iang < pchemrad->nang; ++iang) {
              pchemrad->ir(k, j, i, ifreq * pchemrad->nang + iang) = G0_iang(iang);
            }
          }
          if (CHEMISTRY_ENABLED) {
            for (int iang=0; iang < pchemrad->nang; ++iang) {
              // cr rate
              pchemrad->ir(k, j, i,
                  pscalars->chemnet.index_cr_ * pchemrad->nang + iang) = cr_rate;
            }
          }
        }
      }
    }
    // calculate the average radiation field for output of the initial condition
    pchemrad->pchemradintegrator->CopyToOutput();
  }

  // intialize chemical species
  if (NSPECIES > 0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ispec=0; ispec < NSPECIES; ++ispec) {
            Real x1 = pcoord->x1v(i);
            Real nH;
            if (x1 >= xc_start && x1 <= xc_end) {
              nH = nc;
            } else {
              nH = nw;
            }
            pscalars->s(ispec, k, j, i) = r_init*nH;
            if (CHEMISTRY_ENABLED) {
              Real s_ispec = pin->GetOrAddReal("problem",
                  "r_init_"+pscalars->chemnet.species_names[ispec], -1);
              if (s_ispec >= 0.) {
                pscalars->s(ispec, k, j, i) = s_ispec*nH;
              }
            }
          }
        }
      }
    }
  }
  return;
}
