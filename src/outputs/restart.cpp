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
#include "../blockuid/blockuid.hpp"

//======================================================================================
//! \file restart.cpp
//  \brief writes restart dump files
//======================================================================================


RestartOutput::RestartOutput(OutputParameters oparams)
  : OutputType(oparams)
{
}


//--------------------------------------------------------------------------------------
//! \fn void RestartOutput::Initialize(Mesh *pM, ParameterInput *pi)
//  \brief open the restarting file, output the parameter and header blocks
void RestartOutput::Initialize(Mesh *pM, ParameterInput *pi)
{
  std::stringstream msg;
  std::string fname;
  std::ofstream ost;
  FILE *fp;
  MeshBlock *pmb;
  int i, level;
  int idl=IDLENGTH;
  ResSize_t listsize, idlistoffset;
  ID_t rawid[IDLENGTH];

  // create single output, filename:"file_basename"+"."+"file_id"+"."+XXXX+".rst",
  // where XXXX = 4-digit file_number
  char number[8]; // array to store 4-digit number and end-of-string char
  sprintf(number,"%04d",output_params.file_number);
  fname.assign(output_params.file_basename);
  fname.append(".");
  fname.append(output_params.file_id);
  fname.append(".");
  fname.append(number);
  fname.append(".rst");
  // open the stream and output the athinput block
  ost.open(fname.c_str());
  if (!ost){
    msg << "### FATAL ERROR in function [RestartOutput::WriteOutputFile]"
        << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  pi->ParameterDump(ost);
  ost.close();
  ost.clear();

  // reopen it for MPI compatibility
  resfile.ResFileOpen(fname.c_str());
  resfile.ResFileSeekToEnd();

  // output Mesh information; this part is serial
  resfile.ResFileWrite(&(pM->nbtotal), sizeof(int), 1);
  resfile.ResFileWrite(&idl, sizeof(int), 1); // for extensibility
  resfile.ResFileWrite(&(pM->root_level), sizeof(int), 1);
  resfile.ResFileWrite(&(pM->max_level), sizeof(int), 1);
  resfile.ResFileWrite(&(pM->mesh_size), sizeof(RegionSize), 1);
  resfile.ResFileWrite(&(pM->mesh_bcs), sizeof(RegionBCs), 1);
  resfile.ResFileWrite(&(pM->time), sizeof(Real), 1);
  resfile.ResFileWrite(&(pM->dt), sizeof(Real), 1);
  resfile.ResFileWrite(&(pM->ncycle), sizeof(int), 1);

  // seek the file pointer to the end of the file
  resfile.ResFileSeekToEnd();

  // output ID list: gid,uid(level+id),offset
  // This part must be here in the "common" area
  // in order to allow direct access through MPI.
  // This part must be parallelized.
  listsize=sizeof(int)*2+sizeof(ID_t)*IDLENGTH+sizeof(Real)+sizeof(ResSize_t);

  offset=new size_t[pM->nbtotal+1];
  blocksize=new size_t[pM->nbtotal];
  i=0;
  pmb=pM->pblock;
  while(pmb!=NULL) // must be parallelized for MPI
  {
    blocksize[i]=pmb->GetBlockSize();
    i++;
    pmb=pmb->next;
  }

  idlistoffset=resfile.ResFileTell();
  offset[0]=idlistoffset+listsize*pM->nbtotal;
  for(i=1;i<pM->nbtotal;i++) // serial
    offset[i]=offset[i-1]+blocksize[i-1];

  // pack the offset to the head of the ID list in the last element of the offset array
  offset[pM->nbtotal]=idlistoffset;
  // the offset values must be distributed for MPI

  // seek to the head of the ID list of this process
  resfile.ResFileSeek(idlistoffset+pM->nbstart*listsize);

  pmb=pM->pblock;
  i=0;
  while(pmb!=NULL)
  {
    level=pmb->uid.GetLevel();
    pmb->uid.GetRawUID(rawid);

    resfile.ResFileWrite(&(pmb->gid),sizeof(int),1);
    resfile.ResFileWrite(&level,sizeof(int),1);
    resfile.ResFileWrite(rawid,sizeof(ID_t),IDLENGTH);
    resfile.ResFileWrite(&(pmb->cost),sizeof(Real),1);
    resfile.ResFileWrite(&(offset[i]),sizeof(ResSize_t),1);
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
//! \fn void RestartOutput::Finalize(void)
//  \brief close the file
void RestartOutput::Finalize(void)
{
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  resfile.ResFileClose();
  delete [] blocksize;
  delete [] offset;
}
//--------------------------------------------------------------------------------------
//! \fn void RestartOutput:::WriteOutputFile(OutputData *pod, MeshBlock *pmb)
//  \brief writes OutputData to file in Restart format

void RestartOutput::WriteOutputFile(OutputData *pod, MeshBlock *pmb)
{
  resfile.ResFileSeek(offset[pmb->gid]);
  resfile.ResFileWrite(&(pmb->block_size), sizeof(RegionSize), 1);
  resfile.ResFileWrite(&(pmb->block_bcs), sizeof(RegionBCs), 1);
  resfile.ResFileWrite(&(pmb->neighbor), sizeof(NeighborBlock), 6*2*2);
  resfile.ResFileWrite(pmb->x1f.GetArrayPointer(),sizeof(Real),pmb->x1f.GetSize());
  resfile.ResFileWrite(pmb->x2f.GetArrayPointer(),sizeof(Real),pmb->x2f.GetSize());
  resfile.ResFileWrite(pmb->x3f.GetArrayPointer(),sizeof(Real),pmb->x3f.GetSize());
  resfile.ResFileWrite(pmb->pfluid->u.GetArrayPointer(),sizeof(Real),
                       pmb->pfluid->u.GetSize());
  if (MAGNETIC_FIELDS_ENABLED) {
    resfile.ResFileWrite(pmb->pfield->b.x1f.GetArrayPointer(),sizeof(Real),
                         pmb->pfield->b.x1f.GetSize());
    resfile.ResFileWrite(pmb->pfield->b.x2f.GetArrayPointer(),sizeof(Real),
                         pmb->pfield->b.x2f.GetSize());
    resfile.ResFileWrite(pmb->pfield->b.x3f.GetArrayPointer(),sizeof(Real),
                         pmb->pfield->b.x3f.GetSize());
  }
  return;
}


