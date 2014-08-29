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
#include "../coordinates/coordinates.hpp"
#include "../fluid.hpp"
#include "../mesh.hpp"
#include "outputs.hpp"

//======================================================================================
/*! \file history.cpp
 *  \brief writes history output data.  History data are volume-averaged quantities that
 *  are output frequently in time to trace their history.
 *====================================================================================*/

//--------------------------------------------------------------------------------------
// HistoryOutput constructor

HistoryOutput::HistoryOutput(OutputParameters oparams, MeshBlock *pb)
  : OutputType(oparams,pb)
{
}

// destructor - not required for this derived class

//--------------------------------------------------------------------------------------
/*! \fn void HistoryOutput::LoadOutputData(OutputData *pod)
 *  \brief computes data to be included in output data container (OutputData).  This
 *  version over-rides the default defined in the base class.
 */

void HistoryOutput::LoadOutputData(OutputData *pod)
{
  Fluid *pf = pmy_block->pfluid;;
  AthenaArray<Real> vol = pmy_block->pcoord->cell_volume.ShallowCopy();

// add header

  std::stringstream str;
  str << "# Athena++ history data" << std::endl;
  pod->data_header.descriptor.append(str.str());
  pod->data_header.il = pmy_block->is; pod->data_header.iu = pmy_block->ie;
  pod->data_header.jl = pmy_block->js; pod->data_header.ju = pmy_block->je;
  pod->data_header.kl = pmy_block->ks; pod->data_header.ku = pmy_block->ke;

// Add text for column headers to var_header

  OutputVariableHeader var_header;
  AthenaArray<Real> *phistory_datum;
  phistory_datum = new AthenaArray<Real>;
  phistory_datum->NewAthenaArray(10,1,1,1);
  var_header.type = "SCALARS";
  var_header.name.assign("[1]=time     ");
  var_header.name.append("[2]=dt       ");
  var_header.name.append("[3]=mass     ");
  var_header.name.append("[4]=1-mom    ");
  var_header.name.append("[5]=2-mom    ");
  var_header.name.append("[6]=3-mom    ");
  var_header.name.append("[7]=1-KE     ");
  var_header.name.append("[8]=2-KE     ");
  var_header.name.append("[9]=3-KE     ");
  var_header.name.append("[10]=tot-E   ");

// Add time, time step

  (*phistory_datum)(0,0,0,0) = pf->pmy_block->pmy_domain->pmy_mesh->time;
  (*phistory_datum)(1,0,0,0) = pf->pmy_block->pmy_domain->pmy_mesh->dt;

// Sum over cells, add mass, mom, KE, and total-E

  for (int n=2; n<10; ++n) (*phistory_datum)(n,0,0,0) = 0.0;

  for (int k=(pod->data_header.kl); k<=(pod->data_header.ku); ++k) {
  for (int j=(pod->data_header.jl); j<=(pod->data_header.ju); ++j) {
    Real partial_sum[8];
    for (int i=0; i<8; ++i) partial_sum[i] = 0.0;
    pmy_block->pcoord->CellVolume(k,j,(pod->data_header.il),(pod->data_header.iu),vol);

    for (int i=(pod->data_header.il); i<=(pod->data_header.iu); ++i) {
      partial_sum[0] += vol(i)*pf->u(IDN,k,j,i);
      partial_sum[1] += vol(i)*pf->u(IM1,k,j,i);
      partial_sum[2] += vol(i)*pf->u(IM2,k,j,i);
      partial_sum[3] += vol(i)*pf->u(IM3,k,j,i);
      partial_sum[4] += vol(i)*pf->w(IDN,k,j,i)*pf->w(IM1,k,j,i)*pf->w(IM1,k,j,i);
      partial_sum[5] += vol(i)*pf->w(IDN,k,j,i)*pf->w(IM2,k,j,i)*pf->w(IM2,k,j,i);
      partial_sum[6] += vol(i)*pf->w(IDN,k,j,i)*pf->w(IM3,k,j,i)*pf->w(IM3,k,j,i);
      partial_sum[7] += vol(i)*pf->u(IEN,k,j,i);
    }
    for (int n=0; n<8; ++n) (*phistory_datum)(n+2,0,0,0) += partial_sum[n];

  }}

// Append node to linked list in OutputData

  pod->AppendNode(phistory_datum,var_header);

  return;
}

//--------------------------------------------------------------------------------------
/*! \fn void HistoryOutput:::WriteOutputFile()
 *  \brief writes OutputData to file in history format using C style fprintf
 */

void HistoryOutput::WriteOutputFile(OutputData *pod)
{
  std::stringstream msg;

// create filename: "file_basename" + ".hst".  There is no file number.

  std::string fname;
  fname.assign(output_params.file_basename);
  fname.append(".hst");

// open file for output

  FILE *pfile;
  if((pfile = fopen(fname.c_str(),"a")) == NULL){
    msg << "### FATAL ERROR in function [HistoryOutput::WriteOutputFile]" << std::endl
        << "Output file '" << fname << "' could not be opened";
    throw std::runtime_error(msg.str().c_str());
  }

// If this is the first output, write header

  if (output_params.file_number == 0) {
    fprintf(pfile,"%s",pod->data_header.descriptor.c_str()); // descriptor is first line
    fprintf(pfile,"# ");                              // start second line with hash
    OutputVariable *pvar = pod->pfirst_var;
    while (pvar != NULL) {
      fprintf(pfile,"%s",pvar->var_header.name.c_str()); // add column headers
      pvar = pvar->pnext;
    }
    fprintf(pfile,"\n");                              // terminate line
  }

// step through linked-list of data nodes and write data on same line

  OutputVariable *pvar = pod->pfirst_var;
  while (pvar != NULL) {
    for (int n=0; n<(pvar->pdata->GetDim4()); ++n) {
      fprintf( pfile, output_params.data_format.c_str(), (*pvar->pdata)(n,0,0,0) );
    }
    pvar = pvar->pnext;
  }
  fprintf(pfile,"\n"); // terminate line

// close output file, increment output counter, update time of last output, clean up

  fclose(pfile);
  output_params.file_number++;
  output_params.next_time += output_params.dt;

  return;
}
