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
//! \file rad_integrators.cpp
//  \brief implementation of radiation integrators
//======================================================================================

#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../tc.hpp"
#include "tc_integrators.hpp"
#include "../../coordinates/coordinates.hpp"



TCIntegrator::TCIntegrator(ThermalConduction *ptc, ParameterInput *pin)
{

  pmy_tc = ptc;

  MeshBlock *pmb=ptc->pmy_block;  

  tc_xorder = pin->GetOrAddInteger("time","tc_xorder",2);
  if (tc_xorder == 3) {
    if (NGHOST < 3){ 
      std::stringstream msg;
      msg << "### FATAL ERROR in thermal conduction reconstruction constructor" 
          << std::endl
          << "tc_xorder=" << tc_xorder <<
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

  dflx_.NewAthenaArray(NTC,ncells1);

  // arrays for spatial recontruction 
  //NCR+1 represent Ecr, Fcr1, Fcr2, Fcr3 and vel
  utc_l_.NewAthenaArray(NTC+2,ncells1);
  utc_lb_.NewAthenaArray(NTC+2,ncells1);

  utc_r_.NewAthenaArray(NTC+2,ncells1);
  utc_rho_t_.NewAthenaArray(NTC+2,ncells3,ncells2,ncells1);

  vdiff_.NewAthenaArray(3,ncells3,ncells2,ncells1);
  vdiff_l_.NewAthenaArray(ncells1);
  vdiff_r_.NewAthenaArray(ncells1);
    

  taufact_ = pin->GetOrAddReal("tc","taucell",5.0);

  tc_esource_.NewAthenaArray(ncells3,ncells2,ncells1);
  coord_source_.NewAthenaArray(NCR,ncells3,ncells2,ncells1);

}
// destructor

TCIntegrator::~TCIntegrator()
{

  x1face_area_.DeleteAthenaArray();
  if(pmy_tc->pmy_block->ncells2 > 1) {
    x2face_area_.DeleteAthenaArray();
    x2face_area_p1_.DeleteAthenaArray();
  }
  if(pmy_tc->pmy_block->ncells3 > 1) {
    x3face_area_.DeleteAthenaArray();
    x3face_area_p1_.DeleteAthenaArray();
  }
  cell_volume_.DeleteAthenaArray();


  cwidth2_.DeleteAthenaArray();
  cwidth3_.DeleteAthenaArray();
  dflx_.DeleteAthenaArray();
  utc_l_.DeleteAthenaArray();
  utc_r_.DeleteAthenaArray();
  utc_lb_.DeleteAthenaArray();
  vdiff_l_.DeleteAthenaArray();
  vdiff_r_.DeleteAthenaArray();
  vdiff_.DeleteAthenaArray();

  tc_esource_.DeleteAthenaArray();
  coord_source_.DeleteAthenaArray();
  utc_rho_t_.DeleteAthenaArray();

}




