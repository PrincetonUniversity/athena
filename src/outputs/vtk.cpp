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
//! \file vtk.cpp
//  \brief writes output data in (legacy) vtk format.
//  Data is written in RECTILINEAR_GRID geometry, in BINARY format, and in FLOAT type
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
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "outputs.hpp"
#include "../coordinates/coordinates.hpp"

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

//--------------------------------------------------------------------------------------
// VTKOutput constructor
// destructor - not needed for this derived class

VTKOutput::VTKOutput(OutputParameters oparams)
  : OutputType(oparams)
{
}

//--------------------------------------------------------------------------------------
//! \fn void VTKOutput:::WriteOutputFile(Mesh *pm)
//  \brief writes OutputData to file in (legacy) vtk format

void VTKOutput::WriteOutputFile(Mesh *pm)
{
  MeshBlock *pmb=pm->pblock;
  int big_end = IsBigEndian(); // =1 on big endian machine

  // Loop over MeshBlocks
  while (pmb != NULL) {
    oil=pmb->is; oiu=pmb->ie;
    ojl=pmb->js; oju=pmb->je;
    okl=pmb->ks; oku=pmb->ke;
    if (output_params.include_ghost_zones) {
      oil -= NGHOST; oiu += NGHOST;
      if (ojl != oju) {ojl -= NGHOST; oju += NGHOST;}
      if (okl != oku) {okl -= NGHOST; oku += NGHOST;}
    }

    // set ptrs to data in OutputData linked list
    LoadOutputData(pmb);
    if (TransformOutputData(pmb) == false) {continue;} // skip if slice out of range

    // create filename: "file_basename"+ "."+"bloclid"+"."+"file_id"+"."+XXXXX+".vtk",
    // where XXXXX = 5-digit file_number
    std::string fname;
    char number[6];
    sprintf(number,"%05d",output_params.file_number);
    char blockid[12];
    sprintf(blockid,"block%d",pmb->gid);

    fname.assign(output_params.file_basename);
    fname.append(".");
    fname.append(blockid);
    fname.append(".");
    fname.append(output_params.file_id);
    fname.append(".");
    fname.append(number);
    fname.append(".vtk");

    // open file for output
    FILE *pfile;
    std::stringstream msg;
    if ((pfile = fopen(fname.c_str(),"w")) == NULL){
      msg << "### FATAL ERROR in function [VTKOutput::WriteOutputFile]"
          <<std::endl<< "Output file '" <<fname<< "' could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    // There are five basic parts to the VTK "legacy" file format.
    //  1. Write file version and identifier

    fprintf(pfile,"# vtk DataFile Version 2.0\n");

    //  2. Header

    // print file header
    fprintf(pfile,"# Athena++ data at time=%e",pm->time);
    fprintf(pfile,"  cycle=%d",pmb->pmy_mesh->ncycle);
    fprintf(pfile,"  variables=%s \n",output_params.variable.c_str());

    //  3. File format

    fprintf(pfile,"BINARY\n");

    //  4. Dataset structure

    int ncells1 = oiu - oil + 1;
    int ncells2 = oju - ojl + 1;
    int ncells3 = oku - okl + 1;
    int ncoord1 = ncells1 + 1;
    int ncoord2 = ncells2; if (ncells2 > 1) ncoord2++;
    int ncoord3 = ncells3; if (ncells3 > 1) ncoord3++;

    float *data;
    int ndata = std::max(ncoord1,ncoord2);
    ndata = std::max(ndata,ncoord3);
    data = new float[3*ndata];

    // Specify the type of data, dimensions, and coordinates.  If N>1, then write N+1
    // cell faces as binary floats.  If N=1, then write 1 cell center position.

    fprintf(pfile,"DATASET RECTILINEAR_GRID\n");
    fprintf(pfile,"DIMENSIONS %d %d %d\n",ncoord1,ncoord2,ncoord3);

    // write x1-coordinates as binary float in big endian order

    fprintf(pfile,"X_COORDINATES %d float\n",ncoord1);
    for (int i=oil; i<=oiu+1; ++i) {
      data[i-oil] = (float)(pmb->pcoord->x1f(i));
    }
    if (!big_end) {for (int i=0; i<ncoord1; ++i) Swap4Bytes(&data[i]);}
    fwrite(data,sizeof(float),(size_t)ncoord1,pfile);

    // write x2-coordinates as binary float in big endian order

    fprintf(pfile,"\nY_COORDINATES %d float\n",ncoord2);
    if (ncells2 == 1) {
        data[0] = (float)(pmb->pcoord->x2v(ojl));
    } else {
      for (int j=ojl; j<=oju+1; ++j) {
        data[j-ojl] = (float)(pmb->pcoord->x2f(j));
      }
    }
    if (!big_end) {for (int i=0; i<ncoord2; ++i) Swap4Bytes(&data[i]);}
    fwrite(data,sizeof(float),(size_t)ncoord2,pfile);

    // write x3-coordinates as binary float in big endian order

    fprintf(pfile,"\nZ_COORDINATES %d float\n",ncoord3);
    if (ncells3 == 1) {
        data[0] = (float)(pmb->pcoord->x3v(okl));
    } else {
      for (int k=okl; k<=oku+1; ++k) {
        data[k-okl] = (float)(pmb->pcoord->x3f(k));
      }
    }
    if (!big_end) {for (int i=0; i<ncoord3; ++i) Swap4Bytes(&data[i]);}
    fwrite(data,sizeof(float),(size_t)ncoord3,pfile);

    //  5. Data.  An arbitrary number of scalars and vectors can be written (every node
    //  in the OutputData linked lists), all in binary floats format

    fprintf(pfile,"\nCELL_DATA %d", (ncells1)*(ncells2)*(ncells3));

    OutputData *pdata = pfirst_data_;
    while (pdata != NULL) {

      // write data type (SCALARS or VECTORS) and name
      fprintf(pfile,"\n%s %s float\n",pdata->type.c_str(), pdata->name.c_str());

      int nvar = pdata->data.GetDim4();
      if (nvar == 1) fprintf(pfile,"LOOKUP_TABLE default\n");
      for (int k=okl; k<=oku; ++k) {
      for (int j=ojl; j<=oju; ++j) {

        for (int i=oil; i<=oiu; ++i) {
        for (int n=0; n<nvar; ++n) {
          data[nvar*(i-oil)+n] = (float)pdata->data(n,k,j,i);
        }}

        // write data in big endian order
        if (!big_end) {for (int i=0; i<(nvar*ncells1); ++i) Swap4Bytes(&data[i]);}
        fwrite(data,sizeof(float),(size_t)(nvar*ncells1),pfile);
     
      }}

      pdata = pdata->pnext;
    }

    // don't forget to close the output file and clean up ptrs to data in OutputData
    fclose(pfile);
    ClearOutputData();
    delete data;

    pmb=pmb->next;

  }  // end loop over MeshBlocks

  return;
}
