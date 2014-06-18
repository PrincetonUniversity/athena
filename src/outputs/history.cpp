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

HistoryOutput::HistoryOutput(OutputBlock out_blk, Block *pb)
  : OutputType(out_blk,pb)
{
}

//--------------------------------------------------------------------------------------
/*! \fn OutputData* OutputType::LoadOutputData()
 *  \brief computes data to be included in output data container (OuputData).  This
 *  version over-rides the default defined in the base class.
 */

OutputData* HistoryOutput::LoadOutputData()
{
  Fluid *pf = pparent_block->pfluid;;
  AthenaArray<Real> vol = pparent_block->pcoord->cell_volume.ShallowCopy();

// Allocate OutputData, add header

  OutputData *pod = new OutputData;
  std::stringstream str;
  str << "# Athena++ history data" << std::endl;
  pod->header.descriptor.append(str.str());
  pod->header.il = pparent_block->is; pod->header.iu = pparent_block->ie;
  pod->header.jl = pparent_block->js; pod->header.ju = pparent_block->je;
  pod->header.kl = pparent_block->ks; pod->header.ku = pparent_block->ke;

// Add text for column headers to node_header

  OutputDataNodeHeader node_header;
  AthenaArray<Real> *phistory_datum;
  phistory_datum = new AthenaArray<Real>;
  phistory_datum->NewAthenaArray(10,1,1,1);
  node_header.type = "SCALARS";
  node_header.name.assign("[1]=time     ");
  node_header.name.append("[2]=dt       ");
  node_header.name.append("[3]=mass     ");
  node_header.name.append("[4]=1-mom    ");
  node_header.name.append("[5]=2-mom    ");
  node_header.name.append("[6]=3-mom    ");
  node_header.name.append("[7]=1-KE     ");
  node_header.name.append("[8]=2-KE     ");
  node_header.name.append("[9]=3-KE     ");
  node_header.name.append("[10]=tot-E   ");

// Add time, time step

  (*phistory_datum)(0,0,0,0) = pf->pparent_block->pparent_domain->pparent_mesh->time;
  (*phistory_datum)(1,0,0,0) = pf->pparent_block->pparent_domain->pparent_mesh->dt;

// Sum over cells, add mass, mom, KE, and total-E

  for (int n=2; n<10; ++n) (*phistory_datum)(n,0,0,0) = 0.0;

  for (int k=(pod->header.kl); k<=(pod->header.ku); ++k) {
  for (int j=(pod->header.jl); j<=(pod->header.ju); ++j) {
    Real partial_sum[8];
    for (int i=0; i<8; ++i) partial_sum[i] = 0.0;
    pparent_block->pcoord->CellVolume(k,j,(pod->header.il),(pod->header.iu),vol);

    for (int i=(pod->header.il); i<=(pod->header.iu); ++i) {
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

  pod->AppendNode(phistory_datum,node_header);

  return pod;
}

//--------------------------------------------------------------------------------------
/*! \fn void HistoryOutput:::WriteOutputData()
 *  \brief writes DataBlock to file in history format using C style fprintf
 */

void HistoryOutput::WriteOutputData()
{
  std::stringstream msg;
  OutputData *pod;

// create OutputData

  pod = LoadOutputData();

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

// If this is the first output, write header

  if (output_block.file_number == 0) {
    fprintf(pfile,"%s",pod->header.descriptor.c_str()); // descriptor is first line
    fprintf(pfile,"# ");                              // start second line with hash
    OutputDataNode *pnode = pod->pfirst_node;
    while (pnode != NULL) {
      fprintf(pfile,"%s",pnode->header.name.c_str()); // add column headers
      pnode = pnode->pnext;
    }
    fprintf(pfile,"\n");                              // terminate line
  }

// step through linked-list of data nodes and write data on same line

  OutputDataNode *pnode = pod->pfirst_node;
  while (pnode != NULL) {
    for (int n=0; n<(pnode->pdata->GetDim4()); ++n) {
      fprintf( pfile, output_block.data_format.c_str(), (*pnode->pdata)(n,0,0,0) );
    }
    pnode = pnode->pnext;
  }
  fprintf(pfile,"\n"); // terminate line

// close output file, increment output counter, update time of last output, clean up

  fclose(pfile);
  output_block.file_number++;
  output_block.next_time += output_block.dt;
  delete pod; // delete OutputData object created in LoadOutputData
 
  fclose(pfile);

  return;
}
