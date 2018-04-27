//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file vtk.cpp
//  \brief writes output data in (legacy) vtk format.
//  Data is written in RECTILINEAR_GRID geometry, in BINARY format, and in FLOAT type
//  Writes one file per MeshBlock.

// C headers
#include <stdio.h>
#include <stdlib.h>

// C++ headers
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "outputs.hpp"

//----------------------------------------------------------------------------------------
// Functions to detect big endian machine, and to byte-swap 32-bit words.  The vtk
// legacy format requires data to be stored as big-endian.

int IsBigEndian(void) {
  int32_t n = 1;
  // careful! although int -> char * -> int round-trip conversion is safe,
  // an arbitrary char* may not be converted to int*
  char *ep = reinterpret_cast<char *>(&n);
  return (*ep == 0); // Returns 1 (true) on a big endian machine
}

static inline void Swap4Bytes(void *vdat) {
  char tmp, *dat = static_cast<char *>(vdat);
  tmp = dat[0];  dat[0] = dat[3];  dat[3] = tmp;
  tmp = dat[1];  dat[1] = dat[2];  dat[2] = tmp;
}

//----------------------------------------------------------------------------------------
// VTKOutput constructor
// destructor - not needed for this derived class

VTKOutput::VTKOutput(OutputParameters oparams)
  : OutputType(oparams) {
}

//----------------------------------------------------------------------------------------
//! \fn void VTKOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes OutputData in (legacy) vtk format, one
//         MeshBlock per file

void VTKOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin, bool flag) {
  MeshBlock *pmb=pm->pblock;
  int big_end = IsBigEndian(); // =1 on big endian machine

  // Loop over MeshBlocks
  while (pmb != NULL) {
    // set start/end array indices depending on whether ghost zones are included
    out_is=pmb->is; out_ie=pmb->ie;
    out_js=pmb->js; out_je=pmb->je;
    out_ks=pmb->ks; out_ke=pmb->ke;
    if (output_params.include_ghost_zones) {
      out_is -= NGHOST; out_ie += NGHOST;
      if (out_js != out_je) {out_js -= NGHOST; out_je += NGHOST;}
      if (out_ks != out_ke) {out_ks -= NGHOST; out_ke += NGHOST;}
    }

    // set ptrs to data in OutputData linked list, then slice/sum as needed
    LoadOutputData(pmb);
    if (TransformOutputData(pmb) == false) {
      ClearOutputData();  // required when LoadOutputData() is used.
      pmb=pmb->next;
      continue;
    } // skip if slice was out of range

    // create filename: "file_basename"+ "."+"blockid"+"."+"file_id"+"."+XXXXX+".vtk",
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
    if ((pfile = fopen(fname.c_str(),"w")) == NULL) {
      msg << "### FATAL ERROR in function [VTKOutput::WriteOutputFile]"
          <<std::endl<< "Output file '" <<fname<< "' could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    // There are five basic parts to the VTK "legacy" file format.
    //  1. Write file version and identifier
    fprintf(pfile,"# vtk DataFile Version 2.0\n");

    //  2. Header
    fprintf(pfile,"# Athena++ data at time=%e",pm->time);
    fprintf(pfile,"  cycle=%d",pmb->pmy_mesh->ncycle);
    fprintf(pfile,"  variables=%s \n",output_params.variable.c_str());

    //  3. File format
    fprintf(pfile,"BINARY\n");

    //  4. Dataset structure
    int ncells1 = out_ie - out_is + 1;
    int ncells2 = out_je - out_js + 1;
    int ncells3 = out_ke - out_ks + 1;

    int ncoord1 = ncells1;
    if (ncells1 > 1) ncoord1++;
    int ncoord2 = ncells2;
    if (ncells2 > 1) ncoord2++;
    int ncoord3 = ncells3;
    if (ncells3 > 1) ncoord3++;

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
    if (ncells1 == 1) {
        data[0] = static_cast<float>(pmb->pcoord->x1v(out_is));
    } else {
      for (int i=out_is; i<=out_ie+1; ++i) {
        data[i-out_is] = static_cast<float>(pmb->pcoord->x1f(i));
      }
    }
    if (!big_end) {for (int i=0; i<ncoord1; ++i) Swap4Bytes(&data[i]);}
    fwrite(data,sizeof(float),static_cast<size_t>(ncoord1),pfile);

    // write x2-coordinates as binary float in big endian order
    fprintf(pfile,"\nY_COORDINATES %d float\n",ncoord2);
    if (ncells2 == 1) {
        data[0] = static_cast<float>(pmb->pcoord->x2v(out_js));
    } else {
      for (int j=out_js; j<=out_je+1; ++j) {
        data[j-out_js] = static_cast<float>(pmb->pcoord->x2f(j));
      }
    }
    if (!big_end) {for (int i=0; i<ncoord2; ++i) Swap4Bytes(&data[i]);}
    fwrite(data,sizeof(float),static_cast<size_t>(ncoord2),pfile);

    // write x3-coordinates as binary float in big endian order
    fprintf(pfile,"\nZ_COORDINATES %d float\n",ncoord3);
    if (ncells3 == 1) {
        data[0] = static_cast<float>(pmb->pcoord->x3v(out_ks));
    } else {
      for (int k=out_ks; k<=out_ke+1; ++k) {
        data[k-out_ks] = static_cast<float>(pmb->pcoord->x3f(k));
      }
    }
    if (!big_end) {for (int i=0; i<ncoord3; ++i) Swap4Bytes(&data[i]);}
    fwrite(data,sizeof(float),static_cast<size_t>(ncoord3),pfile);

    //  5. Data.  An arbitrary number of scalars and vectors can be written (every node
    //  in the OutputData linked lists), all in binary floats format
    fprintf(pfile,"\nCELL_DATA %d", (ncells1)*(ncells2)*(ncells3));

    OutputData *pdata = pfirst_data_;
    while (pdata != NULL) {
      // write data type (SCALARS or VECTORS) and name
      fprintf(pfile,"\n%s %s float\n",pdata->type.c_str(), pdata->name.c_str());

      int nvar = pdata->data.GetDim4();
      if (nvar == 1) fprintf(pfile,"LOOKUP_TABLE default\n");
      for (int k=out_ks; k<=out_ke; ++k) {
      for (int j=out_js; j<=out_je; ++j) {

        for (int i=out_is; i<=out_ie; ++i) {
        for (int n=0; n<nvar; ++n) {
          data[nvar*(i-out_is)+n] = static_cast<float>(pdata->data(n,k,j,i));
        }}

        // write data in big endian order
        if (!big_end) {for (int i=0; i<(nvar*ncells1); ++i) Swap4Bytes(&data[i]);}
        fwrite(data,sizeof(float),static_cast<size_t>(nvar*ncells1),pfile);

      }}

      pdata = pdata->pnext;
    }

    // don't forget to close the output file and clean up ptrs to data in OutputData
    fclose(pfile);
    ClearOutputData();  // required when LoadOutputData() is used.
    delete [] data;

    pmb=pmb->next;
  }  // end loop over MeshBlocks

  // increment counters
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);

  return;
}
