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
/*! \file vtk.cpp
 *  \brief writes output data in (legacy) vtk format.
 *  Data is written in RECTILINEAR_GRID geometry, in BINARY format, and in FLOAT type
 *====================================================================================*/

//--------------------------------------------------------------------------------------
// VTKOutput constructor

VTKOutput::VTKOutput(OutputBlock out_blk, Block *pb)
  : OutputType(out_blk,pb)
{
}

//--------------------------------------------------------------------------------------
// Functions to detect big endian machine, and to byte-swap 32-bit words.

int IsBigEndian(void)
{
  short int n = 1;
  char *ep = (char *)&n;
  return (*ep == 0); // Returns 1 on a big endian machine
}

static inline void Swap4Bytes(void *vdat) {
  char tmp, *dat = (char *) vdat;
  tmp = dat[0];  dat[0] = dat[3];  dat[3] = tmp;
  tmp = dat[1];  dat[1] = dat[2];  dat[2] = tmp;
}

//--------------------------------------------------------------------------------------
/*! \fn void VTKOutput:::WriteOutputData()
 *  \brief writes DataBlock to file in (legacy) vtk format  */

void VTKOutput::WriteOutputData()
{
  std::stringstream msg;
  OutputData *pod;
  int big_end = IsBigEndian(); // =1 on big endian machine

// create OutputData, apply transforms (slices, sums, etc)

  pod = LoadOutputData();
  ComputeOutputData(pod);

// create filename

  std::string fname;
  fname.assign(output_block.file_basename);
  fname.append(".");
  char number[5];
  sprintf(number,"%04d",output_block.file_number);
  fname.append(number);
  fname.append(".vtk");

// open file for output

  FILE *pfile;
  if ((pfile = fopen(fname.c_str(),"w")) == NULL){
    msg << "### FATAL ERROR in function [VTKOutput::WriteOutputData]"
        << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// There are five basic parts to the VTK "legacy" file format.
//  1. Write file version and identifier

  fprintf(pfile,"# vtk DataFile Version 2.0\n");

//  2. Header

  fprintf(pfile,"%s",pod->header.descriptor.c_str());

//  3. File format

  fprintf(pfile,"BINARY\n");

//  4. Dataset structure


  int ncells1 = pod->header.iu - pod->header.il + 1;
  int ncells2 = pod->header.ju - pod->header.jl + 1;
  int ncells3 = pod->header.ku - pod->header.kl + 1;
  int ncoord1 = ncells1 + 1;
  int ncoord2 = ncells2; if (ncells2 > 1) ncoord2++;
  int ncoord3 = ncells3; if (ncells3 > 1) ncoord3++;

  float *data;
  int ndata = std::max(ncoord1,ncoord2);
  ndata = std::max(ndata,ncoord3);
  data = new float[3*ndata];

// Specify the type of data, dimensions, and coordinates.  If N>1, then write N+1 cell
// faces as binary floats.  If N=1, then write 1 cell center position.

  fprintf(pfile,"DATASET RECTILINEAR_GRID\n");
  fprintf(pfile,"DIMENSIONS %d %d %d\n",ncoord1,ncoord2,ncoord3);

// write x1-coordinates as binary float in big endian order

  fprintf(pfile,"X_COORDINATES %d float\n",ncoord1);
  for (int i=(pod->header.il); i<=(pod->header.iu)+1; ++i) {
    data[i-(pod->header.il)] = (float)pparent_block->x1f(i);
  }
  if (!big_end) {for (int i=0; i<ncoord1; ++i) Swap4Bytes(&data[i]);}
  fwrite(data,sizeof(float),(size_t)ncoord1,pfile);

// write x2-coordinates as binary float in big endian order

  fprintf(pfile,"\nY_COORDINATES %d float\n",ncoord2);
  if (ncells2 == 1) {
      data[0] = (float)pparent_block->x2v(pod->header.jl);
  } else {
    for (int j=(pod->header.jl); j<=(pod->header.ju)+1; ++j) {
      data[j-(pod->header.jl)] = (float)pparent_block->x2f(j);
    }
  }
  if (!big_end) {for (int i=0; i<ncoord2; ++i) Swap4Bytes(&data[i]);}
  fwrite(data,sizeof(float),(size_t)ncoord2,pfile);

// write x3-coordinates as binary float in big endian order

  fprintf(pfile,"\nZ_COORDINATES %d float\n",ncoord3);
  if (ncells3 == 1) {
      data[0] = (float)pparent_block->x3v(pod->header.kl);
  } else {
    for (int k=(pod->header.kl); k<=(pod->header.ku)+1; ++k) {
      data[k-(pod->header.kl)] = (float)pparent_block->x3f(k);
    }
  }
  if (!big_end) {for (int i=0; i<ncoord3; ++i) Swap4Bytes(&data[i]);}
  fwrite(data,sizeof(float),(size_t)ncoord3,pfile);

//  5. Data.  An arbitrary number of scalars and vectors can be written (every node in
//  in the OutputData linked lists), all in binary floats format

  fprintf(pfile,"\nCELL_DATA %d ", (ncells1)*(ncells2)*(ncells3));

// step through linked-list of data nodes and write out each data array

  OutputDataNode *pnode = pod->pfirst_node;
  while (pnode != NULL) {

// write data type (SCALARS or VECTORS) and name

    fprintf(pfile,"\n%s %s float\n",pnode->header.type.c_str(),
      pnode->header.name.c_str());

    int nvar = pnode->pdata->GetDim4();
    if (nvar == 1) fprintf(pfile,"LOOKUP_TABLE default\n");
    for (int k=(pod->header.kl); k<=(pod->header.ku); ++k) {
    for (int j=(pod->header.jl); j<=(pod->header.ju); ++j) {

      for (int i=(pod->header.il); i<=(pod->header.iu); ++i) {
      for (int n=0; n<nvar; ++n) {
        data[nvar*(i-(pod->header.il))+n] = (float)(*pnode->pdata)(n,k,j,i);
      }}

// write data in big endian order

      if (!big_end) {for (int i=0; i<(nvar*ncells1); ++i) Swap4Bytes(&data[i]);}
      fwrite(data,sizeof(float),(size_t)(nvar*ncells1),pfile);
     
    }}

    pnode = pnode->pnext;
  }

// close output file, increment file number, update time of last output, clean up

  fclose(pfile);
  output_block.file_number++;
  output_block.next_time += output_block.dt;
  delete pod; // delete OutputData object created in LoadOutputData
  delete data;

  return;
}
