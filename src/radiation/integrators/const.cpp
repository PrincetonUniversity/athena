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
//! \file const.cpp
//  \brief implementation of radiation integrators: constant radiation
//======================================================================================


// Athena++ headers
#include "../radiation.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"

// Class header
#include "rad_integrators.hpp"

RadIntegrator::RadIntegrator(Radiation *prad, ParameterInput *pin)
{
  pmy_mb = prad->pmy_block;
  pmy_rad = prad;
#ifdef INCLUDE_CHEMISTRY
  pmy_chemnet = &(pmy_mb->pscalars->chemnet);
  ncol = pmy_chemnet->n_cols_;
  col.NewAthenaArray(pmy_rad->nang,
                     pmy_mb->ncells3, pmy_mb->ncells2, pmy_mb->ncells1, ncol);
  col_avg.NewAthenaArray(ncol, pmy_mb->ncells3, pmy_mb->ncells2, pmy_mb->ncells1);
#endif
}

RadIntegrator::~RadIntegrator() {}

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

void RadIntegrator::UpdateRadiation(int direction){}

