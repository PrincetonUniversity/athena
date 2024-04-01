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
//! \file cr_integrators.cpp
//  \brief implementation of cosmic ray integrators
//======================================================================================

// C headers

// C++ headers
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../cr.hpp"
#include "cr_integrators.hpp"


CRIntegrator::CRIntegrator(CosmicRay *pcr, ParameterInput *pin) {
  pmy_cr = pcr;
  MeshBlock *pmb=pcr->pmy_block;

  cr_xorder = pin->GetOrAddInteger("time","cr_xorder",2);
  if (cr_xorder == 3) {
    if (NGHOST < 3) {
      std::stringstream msg;
      msg << "### FATAL ERROR in cosmic ray reconstruction constructor" << std::endl
          << "cr_xorder=" << cr_xorder <<
          " (PPM) reconstruction selected, but nghost=" << NGHOST << std::endl
          << "Reconfigure with --nghost=3  " <<std::endl;
      ATHENA_ERROR(msg);
    }
  }

  int ncells1 = pmb->ncells1, ncells2 = pmb->ncells2,
  ncells3 = pmb->ncells3;

  x1face_area_.NewAthenaArray(ncells1+1);
  if(ncells2 > 1) {
    x2face_area_.NewAthenaArray(ncells1);
    x2face_area_p1_.NewAthenaArray(ncells1);
  }
  if(ncells3 > 1) {
    x3face_area_.NewAthenaArray(ncells1);
    x3face_area_p1_.NewAthenaArray(ncells1);
  }
  cell_volume_.NewAthenaArray(ncells1);


  cwidth2_.NewAthenaArray(ncells1);
  cwidth3_.NewAthenaArray(ncells1);

  dflx_.NewAthenaArray(NCR,ncells1);

  // arrays for spatial recontruction
  //NCR+1 represent Ecr, Fcr1, Fcr2, Fcr3 and vel
  ucr_l_.NewAthenaArray(NCR+1,ncells1);
  ucr_lb_.NewAthenaArray(NCR+1,ncells1);

  ucr_r_.NewAthenaArray(NCR+1,ncells1);
  ucr_vel_.NewAthenaArray(NCR+1,ncells3,ncells2,ncells1);

  vdiff_l_.NewAthenaArray(ncells1);
  vdiff_r_.NewAthenaArray(ncells1);


  taufact_ = pin->GetOrAddReal("cr","taucell",1.0);
  vel_flx_flag_ = pin->GetOrAddInteger("cr","vflx",0);

  new_sol_.NewAthenaArray(NCR,ncells1);

  grad_pc_.NewAthenaArray(3,ncells3,ncells2,ncells1);
  ec_source_.NewAthenaArray(ncells3,ncells2,ncells1);
  coord_source_.NewAthenaArray(NCR,ncells3,ncells2,ncells1);
}
