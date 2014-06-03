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
#include <stdlib.h>
#include <stdio.h>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../fluid.hpp"
#include "outputs.hpp"

//======================================================================================
/*! \file history.cpp
 *  \brief writes history output data.  History data are volume-averaged quantities that
 *  are output frequently in time to trace their history.
 *====================================================================================*/

//--------------------------------------------------------------------------------------
// HistoryOutput constructor

HistoryOutput::HistoryOutput(OutputBlock out_blk, Block *pb)
  : OutputType(out_blk,pb)
{
}

//--------------------------------------------------------------------------------------
/*! \fn void HistoryOutput::ComputeDataList()
 *  \brief
 */

void HistoryOutput::ComputeOutputData(OutputData *pod)
{
}

//--------------------------------------------------------------------------------------
/*! \fn void HistoryOutput:::WriteOutputData()
 *  \brief writes DataBlock to file in history format in c style printf for the moment
 */

void HistoryOutput::WriteOutputData()
{
  std::stringstream msg;

// create filename
  std::string fname;
  fname.assign(output_block.file_basename);
  fname.append(".hst");

// open file for output
  FILE *pfile;
  if((pfile = fopen(fname.c_str(),"a")) == NULL){
    msg << "### FATAL ERROR in function [HistoryOutput::WriteOutputData]" << std::endl
        << "Output file '" << fname << "' could not be opened";
    throw std::runtime_error(msg.str().c_str());
  }
 
  fclose(pfile);
  return;
}
