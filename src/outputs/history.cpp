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
//! \file history.cpp
//  \brief writes history output data.  History data are volume-averaged quantities that
//  are output frequently in time to trace their history.
//======================================================================================

// C/C++ headers
#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../fluid/fluid.hpp"
#include "../field/field.hpp"
#include "../mesh.hpp"

//this class header
#include "outputs.hpp"

//--------------------------------------------------------------------------------------
// HistoryOutput constructor

HistoryOutput::HistoryOutput(OutputParameters oparams)
  : OutputType(oparams)
{
}

// destructor - not required for this derived class

//--------------------------------------------------------------------------------------
//! \fn void HistoryOutput::LoadOutputData(OutputData *pod, MeshBlock *pmd)
//  \brief computes data to be included in output data container (OutputData).  This
//  version over-rides the default defined in the base class.

void HistoryOutput::LoadOutputData(OutputData *pod, MeshBlock *pmb)
{
  Hydro *pfl = pmb->pfluid;;
  Field *pfd = pmb->pfield;;
  Mesh *pmm = pmb->pmy_mesh;
  int tid=0;

  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  AthenaArray<Real> cell_volume;
  cell_volume.NewAthenaArray((pmm->GetNumMeshThreads()),ncells1);

  AthenaArray<Real> vol;
  vol.InitWithShallowSlice(cell_volume,2,tid,1);

// add OutputData header

  std::stringstream str;
  str << "# Athena++ history data" << std::endl;
  pod->data_header.descriptor.append(str.str());
  pod->data_header.il = pmb->is; pod->data_header.iu = pmb->ie;
  pod->data_header.jl = pmb->js; pod->data_header.ju = pmb->je;
  pod->data_header.kl = pmb->ks; pod->data_header.ku = pmb->ke;
  pod->data_header.ndata = (pmb->ie - pmb->is + 1)*(pmb->je - pmb->js + 1)
                          *(pmb->ke - pmb->ks + 1);

// Add text for column headers to var_header

  int nvars = 12;
  if (NON_BAROTROPIC_EOS) nvars++;

  OutputVariable *pvar = new OutputVariable;
  pvar->data.NewAthenaArray(nvars,1,1,1);

  pvar->type = "SCALARS";
  pvar->name.assign("[1]=time     ");
  pvar->name.append("[2]=dt       ");
  pvar->name.append("[3]=mass     ");
  pvar->name.append("[4]=1-mom    ");
  pvar->name.append("[5]=2-mom    ");
  pvar->name.append("[6]=3-mom    ");
  pvar->name.append("[7]=1-KE     ");
  pvar->name.append("[8]=2-KE     ");
  pvar->name.append("[9]=3-KE     ");
  if (NON_BAROTROPIC_EOS) pvar->name.append("[10]=tot-E   ");
  if (MAGNETIC_FIELDS_ENABLED) {
    pvar->name.append("[11]=1-ME    ");
    pvar->name.append("[12]=2-ME    ");
    pvar->name.append("[13]=3-ME    ");
  }

// Add time, time step

  pvar->data(0,0,0,0) = pmb->pmy_mesh->time;
  pvar->data(1,0,0,0) = pmb->pmy_mesh->dt;

// Sum over cells, add mass, mom, KE, and total-E

  for (int n=2; n<nvars; ++n) pvar->data(n,0,0,0) = 0.0;

  for (int k=(pod->data_header.kl); k<=(pod->data_header.ku); ++k) {
  for (int j=(pod->data_header.jl); j<=(pod->data_header.ju); ++j) {
    Real partial_sum[(nvars-2)];
    for (int i=0; i<(nvars-2); ++i) partial_sum[i] = 0.0;
    pmb->pcoord->CellVolume(k,j,(pod->data_header.il),(pod->data_header.iu),vol);

#pragma simd
    for (int i=(pod->data_header.il); i<=(pod->data_header.iu); ++i) {
      Real& u_d  = pfl->u(IDN,k,j,i);
      Real& u_mx = pfl->u(IM1,k,j,i);
      Real& u_my = pfl->u(IM2,k,j,i);
      Real& u_mz = pfl->u(IM3,k,j,i);

      partial_sum[0] += vol(i)*u_d;
      partial_sum[1] += vol(i)*u_mx;
      partial_sum[2] += vol(i)*u_my;
      partial_sum[3] += vol(i)*u_mz;
      partial_sum[4] += vol(i)*0.5*SQR(u_mx)/u_d;
      partial_sum[5] += vol(i)*0.5*SQR(u_my)/u_d;
      partial_sum[6] += vol(i)*0.5*SQR(u_mz)/u_d;

      if (NON_BAROTROPIC_EOS) {
        Real& u_e = pfl->u(IEN,k,j,i);;
        partial_sum[7] += vol(i)*u_e;
      }
      if (MAGNETIC_FIELDS_ENABLED) {
        Real& bcc1 = pfd->bcc(IB1,k,j,i);
        Real& bcc2 = pfd->bcc(IB2,k,j,i);
        Real& bcc3 = pfd->bcc(IB3,k,j,i);
        partial_sum[8]  += vol(i)*0.5*bcc1*bcc1;
        partial_sum[9]  += vol(i)*0.5*bcc2*bcc2;
        partial_sum[10] += vol(i)*0.5*bcc3*bcc3;
      }
    }
    for (int n=0; n<(nvars-2); ++n) pvar->data(n+2,0,0,0) += partial_sum[n];

  }}

// Append node to linked list in OutputData

  pod->AppendNode(pvar);

// Modify OutputData header

  pod->data_header.il = 1; pod->data_header.iu = 1;
  pod->data_header.jl = 1; pod->data_header.ju = 1;
  pod->data_header.kl = 1; pod->data_header.ku = 1;
  pod->data_header.ndata = 1;

  cell_volume.DeleteAthenaArray();

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void HistoryOutput:::WriteOutputFile()
//  \brief writes OutputData to file in history format using C style fprintf

void HistoryOutput::WriteOutputFile(OutputData *pod, MeshBlock *pmb)
{
  std::stringstream msg;
  if (pod->data_header.ndata == 0) return;  // slice out of range, etc.

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
      fprintf(pfile,"%s",pvar->name.c_str()); // add column headers
      pvar = pvar->pnext;
    }
    fprintf(pfile,"\n");                              // terminate line
  }

// step through linked-list of data nodes and write data on same line

  OutputVariable *pvar = pod->pfirst_var;
  while (pvar != NULL) {
    for (int n=0; n<(pvar->data.GetDim4()); ++n) {
      fprintf( pfile, output_params.data_format.c_str(), pvar->data(n,0,0,0) );
    }
    pvar = pvar->pnext;
  }
  fprintf(pfile,"\n"); // terminate line

// close output file, increment output counter, update time of last output, clean up

  fclose(pfile);

  return;
}
