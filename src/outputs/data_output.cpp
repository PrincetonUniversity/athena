//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in
 * the code distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "data_output.hpp"

//======================================================================================
/*! \file data_output.cpp
 *  \brief implements functions for fluid data outputs
 *====================================================================================*/

// constructor

DataOutput::DataOutput(ParameterInput *pin)
{
  pfirst_block_ = NULL;
  InputBlock *pin_block = pin->pfirst_block_;
  OutputBlock *pob = NULL;
  
// loop over input block names.  Find those that start with "output", and create new
// node in OutputBlock list.  Read input paramters in block into node.

  while (pin_block != NULL) {
    OutputBlock* plast = pob;
    if (pin_block->block_name.compare(0,6,"output") == 0) { 
      pob = new OutputBlock;
// extract integer number of output block.  Save name and number 
      std::string outn = pin_block->block_name.substr(6); // 6 because starts at 0!
      pob->block_number = atoi(outn.c_str());
      pob->block_name.assign(pin_block->block_name);
      pob->pnext = NULL;
// set time of last output, time between outputs
      pob->last_time = pin->GetOrAddReal(pob->block_name,"last_time",0.0);
      pob->dt = pin->GetReal(pob->block_name,"dt");
// Create filename

      std::cout << pob->block_number << std::endl;
      std::cout << std::setprecision(6) << pob->last_time << std::endl;
      std::cout << std::setprecision(6) << pob->dt << std::endl;

// if this is the first output block in list, save pointer to it in class

     if (pfirst_block_ == NULL) {
       pfirst_block_ = pob;
     } else {
       plast->pnext = pob;  // add to end of list
     }
    }
    pin_block = pin_block->pnext;  // move to next input block name
  }
  
}

// destructor

DataOutput::~DataOutput()
{
}

//--------------------------------------------------------------------------------------
//! \fn void
//  \brief

void DataOutput::CheckForOutputs(Mesh *pm)
{
  OutputBlock* pob = pfirst_block_;

  while (pob != NULL) {
    if ( (pm->time == pm->start_time) || (pm->time >= (pob->last_time + pob->dt)) ||
         (pm->time == pm->tlim) ) {
//      pob->OutputData.ComputeData();
//      pob->OutputData.WriteFile();
//      IncremntFilename(pob->filename);
      pob->last_time += pob->dt;
      std::cout << "Making output " << pob->block_number << std::endl;
    }
    pob = pob->pnext;

  }
}
