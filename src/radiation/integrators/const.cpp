//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file const.cpp
//! \brief implementation of radiation integrators: constant radiation

// this class header
#include "rad_integrators.hpp"

// Athena++ headers
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../radiation.hpp"

//----------------------------------------------------------------------------------------
//! constructor, for constant radiation integrator
RadIntegrator::RadIntegrator(Radiation *prad, ParameterInput *pin)
#ifdef INCLUDE_CHEMISTRY
    : col(0, 0, 0, 0),
    col_bvar(prad->pmy_block, &col)
#endif //INCLUDE_CHEMISTRY
{
  pmy_mb = prad->pmy_block;
  pmy_rad = prad;
}

//----------------------------------------------------------------------------------------
//! destructor
RadIntegrator::~RadIntegrator() {}

//----------------------------------------------------------------------------------------
//! \fn void RadIntegrator::CopyToOutput()
//! \brief average radiation field over all angles and copy values to output
void RadIntegrator::CopyToOutput() {
  int is = pmy_mb->is;
  int js = pmy_mb->js;
  int ks = pmy_mb->ks;
  int ie = pmy_mb->ie;
  int je = pmy_mb->je;
  int ke = pmy_mb->ke;
  for (int k=ks-NGHOST; k<=ke+NGHOST; ++k) {
    for (int j=js-NGHOST; j<=je+NGHOST; ++j) {
      for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
        for (int ifreq=0; ifreq < pmy_rad->nfreq; ++ifreq) {
          pmy_rad->ir_avg(ifreq, k, j, i) =
            pmy_rad->ir(k, j, i, ifreq * pmy_rad->nang);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadIntegrator::UpdateRadiation()
//! \brief update radiation field
void RadIntegrator::UpdateRadiation() {}

#ifdef INCLUDE_CHEMISTRY
void RadIntegrator::GetColMB(BoundaryFace direction) {}
void RadIntegrator::UpdateCol(BoundaryFace direction) {}
#endif //INCLUDE_CHEMISTRY
