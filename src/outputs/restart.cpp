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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <fstream>
#include <string.h>

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"

// This class header
#include "outputs.hpp"

RestartOutput::RestartOutput(OutputParameters oparams)
  : OutputType(oparams)
{
}

//--------------------------------------------------------------------------------------
//! \fn void RestartOutput::Initialize(Mesh *pM, ParameterInput *pin, bool wtflag)
//  \brief open the restarting file, output the parameter and header blocks

void RestartOutput::Initialize(Mesh *pM, ParameterInput *pin, bool wtflag)
{
/*
  std::string fname;
  std::stringstream ost;
  MeshBlock *pmb;

  // create single output, filename:"file_basename"+"."+"file_id"+"."+XXXXX+".rst",
  // where XXXXX = 5-digit file_number
  char number[6]; // array to store 4-digit number and end-of-string char
  sprintf(number,"%05d",output_params.file_number);
  fname.assign(output_params.file_basename);
  fname.append(".");
  fname.append(output_params.file_id);
  fname.append(".");
  if(wtflag==false) 
    fname.append(number);
  else 
    fname.append("last");
  fname.append(".rst");

  // count up here for the restarting file.
  if(wtflag==false) {
    output_params.file_number++;
    output_params.next_time += output_params.dt;
    pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
    pin->SetReal(output_params.block_name, "next_time", output_params.next_time);
  }
  resfile.Open(fname.c_str(),IO_WRAPPER_WRITE_MODE);

  // prepare the input parameters
  pin->ParameterDump(ost);
  std::string sbuf=ost.str();

  // calculate the header size
  IOWrapperSize_t udsize=0;
  for(int n=0; n<pM->nint_user_mesh_data_; n++)
    udsize+=pM->iusermeshdata[n].GetSizeInBytes();
  for(int n=0; n<pM->nreal_user_mesh_data_; n++)
    udsize+=pM->rusermeshdata[n].GetSizeInBytes();

  headeroffset=sbuf.size()*sizeof(char)+3*sizeof(int)+sizeof(RegionSize)
              +2*sizeof(Real)+sizeof(IOWrapperSize_t)+udsize;
  // the size of an element of the ID list
  listsize=sizeof(LogicalLocation)+sizeof(Real);
  // the size of each MeshBlock
  datasize = pM->pblock->GetBlockSizeInBytes();
  nbtotal=pM->nbtotal;
  myns=pM->nslist[Globals::my_rank];
  mynb=pM->nblist[Globals::my_rank];

  // write the header
  if(Globals::my_rank==0) {
    // output the input parameters; this part is serial
    resfile.Write(sbuf.c_str(),sizeof(char),sbuf.size());

    // output Mesh information; this part is serial
    resfile.Write(&(pM->nbtotal), sizeof(int), 1);
    resfile.Write(&(pM->root_level), sizeof(int), 1);
    resfile.Write(&(pM->mesh_size), sizeof(RegionSize), 1);
    resfile.Write(&(pM->time), sizeof(Real), 1);
    resfile.Write(&(pM->dt), sizeof(Real), 1);
    resfile.Write(&(pM->ncycle), sizeof(int), 1);
    resfile.Write(&(datasize), sizeof(IOWrapperSize_t), 1);

    // collect and write user Mesh data
    if(udsize!=0) {
      char *ud = new char[udsize];
      IOWrapperSize_t udoffset = 0;
      for(int n=0; n<pM->nint_user_mesh_data_; n++) {
        memcpy(&(ud[udoffset]), pM->iusermeshdata[n].data(),
               pM->iusermeshdata[n].GetSizeInBytes());
        udoffset+=pM->iusermeshdata[n].GetSizeInBytes();
      }
      for(int n=0; n<pM->nreal_user_mesh_data_; n++) {
        memcpy(&(ud[udoffset]), pM->rusermeshdata[n].data(),
               pM->rusermeshdata[n].GetSizeInBytes());
        udoffset+=pM->rusermeshdata[n].GetSizeInBytes();
      }
      resfile.Write(ud, 1, udsize);
      delete [] ud;
    }
  }

  // allocate memory for the ID list and the data
  char *idlist=new char [listsize*mynb];
  data = new char [mynb*datasize];
  pmb=pM->pblock;
  int os=0;
  while(pmb!=NULL) {
    // pack the meta data
    memcpy(&(idlist[os]), &(pmb->loc), sizeof(LogicalLocation));
    os+=sizeof(LogicalLocation);
    memcpy(&(idlist[os]), &(pmb->cost), sizeof(Real));
    os+=sizeof(Real);
    pmb=pmb->next;
  }
  // write the ID list collectively
  IOWrapperSize_t myoffset=headeroffset+listsize*myns;
  resfile.Write_at_all(idlist,listsize,mynb,myoffset);

  // deallocate the idlist array
  delete [] idlist;

  // leave the file open; it will be closed in Finalize()
*/
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RestartOutput::LoadOutputData(OutputData *pout_data, MeshBlock *pblock)
//  \brief Load the data array from all the MeshBlocks
void RestartOutput::LoadOutputData(OutputData *pout_data, MeshBlock *pblock)
{
/*
  // pack the data
  char *pdata=&(data[pblock->lid*datasize]);
  memcpy(pdata,pblock->phydro->u.data(), pblock->phydro->u.GetSizeInBytes());
  pdata+=pblock->phydro->u.GetSizeInBytes();
  if (GENERAL_RELATIVITY) {
    memcpy(pdata,pblock->phydro->w.data(), pblock->phydro->w.GetSizeInBytes());
    pdata+=pblock->phydro->w.GetSizeInBytes();
    memcpy(pdata,pblock->phydro->w1.data(), pblock->phydro->w1.GetSizeInBytes());
    pdata+=pblock->phydro->w1.GetSizeInBytes();
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    memcpy(pdata,pblock->pfield->b.x1f.data(), pblock->pfield->b.x1f.GetSizeInBytes());
    pdata+=pblock->pfield->b.x1f.GetSizeInBytes();
    memcpy(pdata,pblock->pfield->b.x2f.data(), pblock->pfield->b.x2f.GetSizeInBytes());
    pdata+=pblock->pfield->b.x2f.GetSizeInBytes();
    memcpy(pdata,pblock->pfield->b.x3f.data(), pblock->pfield->b.x3f.GetSizeInBytes());
    pdata+=pblock->pfield->b.x3f.GetSizeInBytes();
  }

  // add additional physics here
  // also update MeshBlock::GetBlockSizeInBytes accordingly


  // pack the user MeshBlock data
  for(int n=0; n<pblock->nint_user_meshblock_data_; n++) {
    memcpy(pdata, pblock->iusermeshblockdata[n].data(),
           pblock->iusermeshblockdata[n].GetSizeInBytes());
    pdata+=pblock->iusermeshblockdata[n].GetSizeInBytes();
  }
  for(int n=0; n<pblock->nreal_user_meshblock_data_; n++) {
    memcpy(pdata, pblock->rusermeshblockdata[n].data(),
           pblock->rusermeshblockdata[n].GetSizeInBytes());
    pdata+=pblock->rusermeshblockdata[n].GetSizeInBytes();
  }
*/
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RestartOutput::Finalize(ParameterInput *pin)
//  \brief perform collective data output and clean up
/*
void RestartOutput::Finalize(ParameterInput *pin)
{
  // call actual write here
  IOWrapperSize_t myoffset=headeroffset+listsize*nbtotal+datasize*myns;
  resfile.Write_at_all(data,datasize,mynb,myoffset);
  resfile.Close();
  delete [] data;
}
*/

