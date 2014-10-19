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

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh.hpp"
#include "../fluid/fluid.hpp"
#include "outputs.hpp"

//--------------------------------------------------------------------------------------
// Functions to detect big endian machine, and to byte-swap 32-bit words.  The vtk
// legacy format requires data to be stored as big-endian.

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

//======================================================================================
//! \file vtk.cpp
//  \brief writes output data in (legacy) vtk format.
//  Data is written in RECTILINEAR_GRID geometry, in BINARY format, and in FLOAT type
//======================================================================================

//--------------------------------------------------------------------------------------
// VTKOutput constructor

VTKOutput::VTKOutput(OutputParameters oparams)
  : OutputType(oparams)
{
}

// destructor - not needed for this derived class

//--------------------------------------------------------------------------------------
//! \fn void VTKOutput:::WriteOutputFile(OutputData *pod, MeshBlock *pmb)
//  \brief writes OutputData to file in (legacy) vtk format

void VTKOutput::WriteOutputFile(OutputData *pod, MeshBlock *pmb)
{
  std::stringstream msg;
  int big_end = IsBigEndian(); // =1 on big endian machine
  if (pod->data_header.ndata == 0) return;  // slice out of range, etc.

// create filename: "file_basename" + "." + "file_id" + "." + XXXX + ".vtk",
// where XXXX = 4-digit file_number

  std::string fname;
  char number[5]; // array to store 4-digit number and end-of-string char
  sprintf(number,"%04d",output_params.file_number);

  fname.assign(output_params.file_basename);
  fname.append(".");
  fname.append(output_params.file_id);
  fname.append(".");
  fname.append(number);
  fname.append(".vtk");

// open file for output

  FILE *pfile;
  if ((pfile = fopen(fname.c_str(),"w")) == NULL){
    msg << "### FATAL ERROR in function [VTKOutput::WriteOutputFile]"
        << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// There are five basic parts to the VTK "legacy" file format.
//  1. Write file version and identifier

  fprintf(pfile,"# vtk DataFile Version 2.0\n");

//  2. Header

  fprintf(pfile,"%s",pod->data_header.descriptor.c_str());

//  3. File format

  fprintf(pfile,"BINARY\n");

//  4. Dataset structure


  int ncells1 = pod->data_header.iu - pod->data_header.il + 1;
  int ncells2 = pod->data_header.ju - pod->data_header.jl + 1;
  int ncells3 = pod->data_header.ku - pod->data_header.kl + 1;
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
  for (int i=(pod->data_header.il); i<=(pod->data_header.iu)+1; ++i) {
    data[i-(pod->data_header.il)] = (float)pmb->x1f(i);
  }
  if (!big_end) {for (int i=0; i<ncoord1; ++i) Swap4Bytes(&data[i]);}
  fwrite(data,sizeof(float),(size_t)ncoord1,pfile);

// write x2-coordinates as binary float in big endian order

  fprintf(pfile,"\nY_COORDINATES %d float\n",ncoord2);
  if (ncells2 == 1) {
      data[0] = (float)pmb->x2v(pod->data_header.jl);
  } else {
    for (int j=(pod->data_header.jl); j<=(pod->data_header.ju)+1; ++j) {
      data[j-(pod->data_header.jl)] = (float)pmb->x2f(j);
    }
  }
  if (!big_end) {for (int i=0; i<ncoord2; ++i) Swap4Bytes(&data[i]);}
  fwrite(data,sizeof(float),(size_t)ncoord2,pfile);

// write x3-coordinates as binary float in big endian order

  fprintf(pfile,"\nZ_COORDINATES %d float\n",ncoord3);
  if (ncells3 == 1) {
      data[0] = (float)pmb->x3v(pod->data_header.kl);
  } else {
    for (int k=(pod->data_header.kl); k<=(pod->data_header.ku)+1; ++k) {
      data[k-(pod->data_header.kl)] = (float)pmb->x3f(k);
    }
  }
  if (!big_end) {for (int i=0; i<ncoord3; ++i) Swap4Bytes(&data[i]);}
  fwrite(data,sizeof(float),(size_t)ncoord3,pfile);

//  5. Data.  An arbitrary number of scalars and vectors can be written (every node in
//  in the OutputData linked lists), all in binary floats format

  fprintf(pfile,"\nCELL_DATA %d ", (ncells1)*(ncells2)*(ncells3));

// step through linked-list of data nodes and write out each data array

  OutputVariable *pvar = pod->pfirst_var;
  while (pvar != NULL) {

// write data type (SCALARS or VECTORS) and name

    fprintf(pfile,"\n%s %s float\n",pvar->var_header.type.c_str(),
      pvar->var_header.name.c_str());

    int nvar = pvar->pdata->GetDim4();
    if (nvar == 1) fprintf(pfile,"LOOKUP_TABLE default\n");
    for (int k=(pod->data_header.kl); k<=(pod->data_header.ku); ++k) {
    for (int j=(pod->data_header.jl); j<=(pod->data_header.ju); ++j) {

      for (int i=(pod->data_header.il); i<=(pod->data_header.iu); ++i) {
      for (int n=0; n<nvar; ++n) {
        data[nvar*(i-(pod->data_header.il))+n] = (float)(*pvar->pdata)(n,k,j,i);
      }}

// write data in big endian order

      if (!big_end) {for (int i=0; i<(nvar*ncells1); ++i) Swap4Bytes(&data[i]);}
      fwrite(data,sizeof(float),(size_t)(nvar*ncells1),pfile);
     
    }}

    pvar = pvar->pnext;
  }

// close output file, increment file number, update time of last output, clean up

  fclose(pfile);
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  delete data;

  return;
}
