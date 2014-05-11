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
#include "../parameter_input.hpp"
#include "data_output.hpp"

//======================================================================================
/*! \file data_output.cpp
 *  \brief implements functions for fluid data outputs
 *====================================================================================*/

// constructor

DataOutput::DataOutput(ParameterInput *pin)
{
  InputBlock *pin_block = pin->pfirst_block_;
  OutputBlock *pb = pfirst_block_;
  
// loop over input block names.  Find those that start with "output", and create new
// node in OutputBlock list.  Read input paramters in block into node.

  while (pin_block != NULL) {
    if (pin_block->block_name.compare(0,6,"output") == 0) { 
      pb = new OutputBlock;
      pb->pnext = NULL;
// extract integer number of output block.  Save name and number 
      std::string outn = pin_block->block_name.substr(6); // 6 because starts at 0!
      pb->block_number = atoi(outn.c_str());
      pb->block_name   = pin_block->block_name;
// extract time of last output, time between outputs
      pb->last_time = pin->GetOrAddReal(pb->block_name,"last_time",0.0);
      pb->dt = pin->GetReal(pb->block_name,"dt");

      std::cout << pb->block_number << std::endl;
      std::cout << std::setprecision(6) << pb->last_time << std::endl;
      std::cout << std::setprecision(6) << pb->dt << std::endl;

      pb = pb->pnext;
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

void DataOutput::CheckForOutputs()
{
}
