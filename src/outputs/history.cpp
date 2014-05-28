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
#include "../datablock.hpp"
#include "output.hpp"

//======================================================================================
/*! \file history.cpp
 *  \brief writes history output data.  History data are volume-averaged quantities that
 *  are output frequently in time to trace their history.
 *====================================================================================*/

//--------------------------------------------------------------------------------------
// HistoryOutput constructor

HistoryOutput::HistoryOutput(InputBlock *pin_block, Mesh *pm)
  : Output(pin_block,pm)
{
  output_block.out = "cons";
}

//--------------------------------------------------------------------------------------
/*! \fn void HistoryOutput::ComputeFromMesh(Mesh *pm)
 *  \brief Performs any operations on Mesh and stores results in the output's DataBlock
 */

void HistoryOutput::ComputeFromMesh(Mesh *pm)
{

// Read data from Mesh
  pdata = LoadFluidFromMesh(pm);

// Loop over transformation list
  DataBlockTransform* pdbt = ptrans;
  while (pdbt != NULL) {
    pdbt->Transform(pdata);
    pdbt = pdbt->GetNext();
  }

  SumOverAll sum_trans;
  sum_trans.Transform(pdata);

}

//--------------------------------------------------------------------------------------
/*! \fn void HistoryOutput:::Write()
 *  \brief writes DataBlock to file in history format
 *         uses c style printf for the moment
 */

void HistoryOutput::Write()
{
  std::stringstream msg;

// create filename
  std::stringstream ftmp;
  ftmp << "Sod.hst";//hardwire for now
  std::string filename=ftmp.str();

// open file for output
  FILE *pfile;
  if((pfile = fopen(filename.c_str(),"a")) == NULL){
    msg << "### FATAL ERROR in function [HistoryOutput::Write]" << std::endl
        << "Output file '" << filename << "' could not be opened";
    throw std::runtime_error(msg.str().c_str());
  }
 
  Real dVol = 1.0; // Temporary!!!
  DataNode* pdn;
  if (Number() == 0) {
    fprintf(pfile,"# Athena history dump for volume=%e\n",dVol);
// write header description
    if(pdata->GetDescriptor().compare("") != 0)
      fprintf(pfile,"# %s\n",pdata->GetDescriptor().c_str());

// Loop over headers and write row headings
    int col_cnt = 0;
    fprintf(pfile,"#");
    fprintf(pfile," [%d]=%s",col_cnt++,"time");
    fprintf(pfile," [%d]=%s",col_cnt++,"dt");
    pdn = pdata->GetNode(2);
    while (pdn->GetNext() != NULL){
      pdn = pdn->GetNext();
      fprintf(pfile," [%d]=%s",col_cnt,pdn->GetVariableTypes().c_str());
      col_cnt++;
    }
    fprintf(pfile,"\n");
  }
// Write time and timestep
  fprintf( pfile, output_block.dat_fmt.c_str(), pdata->GetTime() );
  fprintf( pfile, output_block.dat_fmt.c_str(), pdata->GetTimeStep() );
// Loop over DataNodes and write sums
  pdn = pdata->GetNode(2);
  while (pdn->GetNext() != NULL){
    pdn = pdn->GetNext();
    fprintf( pfile, output_block.dat_fmt.c_str(), (*pdn->GetData())(0) );
  }
  fprintf(pfile,"\n");
  fclose(pfile);
  delete pdata;  // Release any memory allocated this output cycle
}
