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
//! \file radiation.cpp
//  \brief implementation of functions in class Radiation
//======================================================================================

//C++ headers
#include <string>
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "radiation.hpp"
#include "../mesh/mesh.hpp"
#include "integrators/rad_integrators.hpp"

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin)
{
  // read in the parameters
  integrator = RADIATION_INTEGRATOR;
	nfreq = pin->GetOrAddInteger("radiation","n_frequency",1);
  
  pmy_block = pmb;
  
	if (integrator == "six_ray") {
		nang = 6;
	} else if (integrator == "loc_jeans" or integrator == "const") {
		nang = 1;
	} else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Radiation constructor" << std::endl
        << "integrator=" << integrator << " not valid radiation integrator, " << std::endl
        << "choose from {jeans, six_ray, const}" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  
  n_fre_ang = nang * nfreq;
  
  
  // allocate arrays
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);
  // store frequency and angles as [nfre][ang]
  ir.NewAthenaArray(ncells3, ncells2, ncells1, n_fre_ang);
  ir_avg.NewAthenaArray(nfreq, ncells3, ncells2, ncells1);

  //radiation integrator
  pradintegrator = new RadIntegrator(this, pin);
}


