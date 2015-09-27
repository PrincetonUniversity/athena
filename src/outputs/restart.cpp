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
//! \file restart.cpp
//  \brief writes restart dump files
//======================================================================================

// C/C++ headers
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <stdio.h>
#include <fstream>

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh.hpp"
#include "../parameter_input.hpp"
#include "../fluid/fluid.hpp"
#include "../field/field.hpp"

// This class header
#include "outputs.hpp"

RestartOutput::RestartOutput(OutputParameters oparams)
  : OutputType(oparams)
{
}

//--------------------------------------------------------------------------------------
//! \fn void RestartOutput::Initialize(Mesh *pM, ParameterInput *pin)
//  \brief open the restarting file, output the parameter and header blocks

void RestartOutput::Initialize(Mesh *pM, ParameterInput *pin)
{
  std::stringstream msg;
  std::string fname;
  std::stringstream ost;
  FILE *fp;
  MeshBlock *pmb;
  int *nblocks, *displ;
  IOWrapperSize_t *myblocksize;
  int i, level;
  IOWrapperSize_t listsize;

  // create single output, filename:"file_basename"+"."+"file_id"+"."+XXXXX+".rst",
  // where XXXXX = 5-digit file_number
  char number[6]; // array to store 4-digit number and end-of-string char
  sprintf(number,"%05d",output_params.file_number);
  fname.assign(output_params.file_basename);
  fname.append(".");
  fname.append(output_params.file_id);
  fname.append(".");
  fname.append(number);
  fname.append(".rst");

  // count up here for the restarting file.
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);
  pin->ParameterDump(ost);

  resfile.Open(fname.c_str(),WRAPPER_WRITE_MODE);

  if(Globals::my_rank==0) {
    // output the input parameters; this part is serial
    std::string sbuf=ost.str();
    resfile.Write(sbuf.c_str(),sizeof(char),sbuf.size());

    // output Mesh information; this part is serial
    resfile.Write(&(pM->nbtotal), sizeof(int), 1);
    resfile.Write(&(pM->root_level), sizeof(int), 1);
    resfile.Write(&(pM->mesh_size), sizeof(RegionSize), 1);
    resfile.Write(pM->mesh_bcs, sizeof(int), 6);
    resfile.Write(&(pM->time), sizeof(Real), 1);
    resfile.Write(&(pM->dt), sizeof(Real), 1);
    resfile.Write(&(pM->ncycle), sizeof(int), 1);
  }

  // the size of an element of the ID list
  listsize=sizeof(int)+sizeof(LogicalLocation)+sizeof(Real)+sizeof(IOWrapperSize_t);

  int mynb=pM->nbend-pM->nbstart+1;
  blocksize=new IOWrapperSize_t[pM->nbtotal+1];
  offset=new IOWrapperSize_t[mynb];

#ifdef MPI_PARALLEL
  displ = new int[Globals::nranks];
  if(Globals::my_rank==0) mynb++; // the first process includes the information block
  myblocksize=new IOWrapperSize_t[mynb];

  displ[0]=0;
  for(i=1;i<Globals::nranks;i++)
    displ[i]=pM->nslist[i]+1;

  i=0;
  if(Globals::my_rank==0) { // the information block
    myblocksize[0]=resfile.GetPosition();
    i=1;
  }
  pmb=pM->pblock;
  while(pmb!=NULL) {
    myblocksize[i]=pmb->GetBlockSizeInBytes();
    i++;
    pmb=pmb->next;
  }

  // distribute the size of each block + header size
  pM->nblist[0]++; // include the information block
  MPI_Allgatherv(myblocksize,mynb,MPI_LONG,
     blocksize,pM->nblist,displ,MPI_LONG,MPI_COMM_WORLD);
  pM->nblist[0]--; // recover the original list

  // clean up
  delete [] myblocksize;
  delete [] displ;

#else // serial
  pmb=pM->pblock;
  blocksize[0]=resfile.GetPosition();
  i=1;
  while(pmb!=NULL) {
    blocksize[i]=pmb->GetBlockSizeInBytes();
    i++;
    pmb=pmb->next;
  }
#endif

  // blocksize[0] = offset to the end of the header
  IOWrapperSize_t firstoffset=blocksize[0]+listsize*pM->nbtotal; // end of the id list
  for(i=0;i<pM->nbstart;i++)
    firstoffset+=blocksize[i+1];
  offset[0]=firstoffset;
//  offset[0]=blocksize[0]+listsize*pM->nbtotal;
  for(i=pM->nbstart;i<pM->nbend;i++) // calculate the offsets
    offset[i-pM->nbstart+1]=offset[i-pM->nbstart]+blocksize[i+1]; // blocksize[i]=the size of i-1th block

  // seek to the head of the ID list of this process
  resfile.Seek(blocksize[0]+pM->nbstart*listsize);

  pmb=pM->pblock;
  i=0;
  while(pmb!=NULL) {
    // perhaps these data should be packed and dumped at once
    resfile.Write(&(pmb->gid),sizeof(int),1);
    resfile.Write(&(pmb->loc),sizeof(LogicalLocation),1);
    resfile.Write(&(pmb->cost),sizeof(Real),1);
    resfile.Write(&(offset[i]),sizeof(IOWrapperSize_t),1);
    i++;
    pmb=pmb->next;
  }
  // If any additional physics is implemented, modify here accordingly.
  // Especially, the offset from the top of the file must be recalculated.
  // For MPI-IO, it is important to know the absolute location in advance.

  // leave the file open; it will be closed in Finalize()
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RestartOutput::Finalize(ParameterInput *pin)
//  \brief close the file

void RestartOutput::Finalize(ParameterInput *pin)
{
  resfile.Close();
  delete [] blocksize;
  delete [] offset;
}

//--------------------------------------------------------------------------------------
//! \fn void RestartOutput:::WriteOutputFile(OutputData *pod, MeshBlock *pmb)
//  \brief writes OutputData to file in Restart format

void RestartOutput::WriteOutputFile(OutputData *pod, MeshBlock *pmb)
{
  resfile.Seek(offset[pmb->gid - pmb->pmy_mesh->nbstart]);
  resfile.Write(&(pmb->block_size), sizeof(RegionSize), 1);
  resfile.Write(pmb->block_bcs, sizeof(int), 6);
  resfile.Write(pmb->pfluid->u.GetArrayPointer(),sizeof(Real),
                       pmb->pfluid->u.GetSize());
  if (GENERAL_RELATIVITY) {
    resfile.Write(pmb->pfluid->w.GetArrayPointer(),sizeof(Real),
                         pmb->pfluid->w.GetSize());
    resfile.Write(pmb->pfluid->w1.GetArrayPointer(),sizeof(Real),
                         pmb->pfluid->w1.GetSize());
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    resfile.Write(pmb->pfield->b.x1f.GetArrayPointer(),sizeof(Real),
                         pmb->pfield->b.x1f.GetSize());
    resfile.Write(pmb->pfield->b.x2f.GetArrayPointer(),sizeof(Real),
                         pmb->pfield->b.x2f.GetSize());
    resfile.Write(pmb->pfield->b.x3f.GetArrayPointer(),sizeof(Real),
                         pmb->pfield->b.x3f.GetSize());
  }
  return;
}
