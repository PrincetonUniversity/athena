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

#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh.hpp"
#include "../parameter_input.hpp"
#include "../fluid/fluid.hpp"
#include "../field/field.hpp"
#include "outputs.hpp"
#include "../blockuid.hpp"

//======================================================================================
//! \file restart.cpp
//  \brief writes restart dump files
//======================================================================================


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
  WrapIOSize_t *myblocksize;
  int i, level;
  int idl=IDLENGTH;
  WrapIOSize_t listsize, idlistoffset;
  ID_t rawid[IDLENGTH];

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

  resfile.Open(fname.c_str(),writemode);

  if(myrank==0) {
    // output the input parameters; this part is serial
    std::string sbuf=ost.str();
    resfile.Write(sbuf.c_str(),sizeof(char),sbuf.size());

    // output Mesh information; this part is serial
    resfile.Write(&(pM->nbtotal), sizeof(int), 1);
    resfile.Write(&idl, sizeof(int), 1); // for extensibility
    resfile.Write(&(pM->root_level), sizeof(int), 1);
    resfile.Write(&(pM->mesh_size), sizeof(RegionSize), 1);
    resfile.Write(pM->mesh_bcs, sizeof(int), 6);
    resfile.Write(&(pM->time), sizeof(Real), 1);
    resfile.Write(&(pM->dt), sizeof(Real), 1);
    resfile.Write(&(pM->ncycle), sizeof(int), 1);
  }

  // the size of an element of the ID list
  listsize=sizeof(int)*2+sizeof(ID_t)*IDLENGTH+sizeof(Real)+sizeof(WrapIOSize_t);

  int mynb=pM->nbend-pM->nbstart+1;
  blocksize=new WrapIOSize_t[pM->nbtotal+1];
  offset=new WrapIOSize_t[mynb];

#ifdef MPI_PARALLEL
  nblocks = new int[nproc];
  displ = new int[nproc];
  if(myrank==0) myblocksize=new WrapIOSize_t[mynb+1];
  else myblocksize=new WrapIOSize_t[mynb];

  // distribute the numbers of the blocks
  MPI_Allgather(&mynb,1,MPI_INTEGER,nblocks,1,MPI_INTEGER,MPI_COMM_WORLD);
  nblocks[0]+=1; // for the first block
  displ[0]=0;
  for(i=1;i<nproc;i++)
    displ[i]=displ[i-1]+nblocks[i-1];

  i=0;
  if(myrank==0) {
    myblocksize[0]=resfile.Tell();
    i=1;
  }
  pmb=pM->pblock;
  while(pmb!=NULL) {
    myblocksize[i]=pmb->GetBlockSizeInBytes();
    i++;
    pmb=pmb->next;
  }

  // distribute the size of each block + header size
  MPI_Allgatherv(myblocksize,nblocks[myrank],MPI_LONG,
     blocksize,nblocks,displ,MPI_LONG,MPI_COMM_WORLD);

  // clean up
  delete [] myblocksize;
  delete [] nblocks;
  delete [] displ;

#else // serial
  pmb=pM->pblock;
  blocksize[0]=resfile.Tell();
  i=1;
  while(pmb!=NULL) {
    blocksize[i]=pmb->GetBlockSizeInBytes();
    i++;
    pmb=pmb->next;
  }
#endif

  // blocksize[0] = offset to the end of the header
  WrapIOSize_t firstoffset=blocksize[0]+listsize*pM->nbtotal; // end of the id list
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
    level=pmb->uid.GetLevel();
    pmb->uid.GetRawUID(rawid);

    // perhaps these data should be packed and dumped at once
    resfile.Write(&(pmb->gid),sizeof(int),1);
    resfile.Write(&level,sizeof(int),1);
    resfile.Write(rawid,sizeof(ID_t),IDLENGTH);
    resfile.Write(&(pmb->cost),sizeof(Real),1);
    resfile.Write(&(offset[i]),sizeof(WrapIOSize_t),1);
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


