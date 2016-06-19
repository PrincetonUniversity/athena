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
#include "../globals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"

//this class header
#include "outputs.hpp"

#define NHISTORY_VARS ((NHYDRO)+(NFIELD)+3)

//--------------------------------------------------------------------------------------
// HistoryOutput constructor

HistoryOutput::HistoryOutput(OutputParameters oparams)
  : OutputType(oparams)
{
}

// destructor - not needed for this derived class

//--------------------------------------------------------------------------------------
//! \fn void OutputType::HistoryFile()
//  \brief Writes a history file

void HistoryOutput::WriteOutputFile(Mesh *pm)
{
  MeshBlock *pmb=pm->pblock;
  AthenaArray<Real> vol;

  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  vol.NewAthenaArray(ncells1);

  Real data_sum[(NHISTORY_VARS)];
  for (int n=0; n<NHISTORY_VARS; ++n) data_sum[n]=0.0;

  // Loop over MeshBlocks
  while (pmb != NULL) {
    Hydro *phyd = pmb->phydro;;
    Field *pfld = pmb->pfield;;

    // Sum history variables over cells
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,vol);

#pragma simd
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real& u_d  = phyd->u(IDN,k,j,i);
        Real& u_mx = phyd->u(IM1,k,j,i);
        Real& u_my = phyd->u(IM2,k,j,i);
        Real& u_mz = phyd->u(IM3,k,j,i);

        data_sum[0] += vol(i)*u_d;
        data_sum[1] += vol(i)*u_mx;
        data_sum[2] += vol(i)*u_my;
        data_sum[3] += vol(i)*u_mz;
        data_sum[4] += vol(i)*0.5*SQR(u_mx)/u_d;
        data_sum[5] += vol(i)*0.5*SQR(u_my)/u_d;
        data_sum[6] += vol(i)*0.5*SQR(u_mz)/u_d;

        if (NON_BAROTROPIC_EOS) {
          Real& u_e = phyd->u(IEN,k,j,i);;
          data_sum[7] += vol(i)*u_e;
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          Real& bcc1 = pfld->bcc(IB1,k,j,i);
          Real& bcc2 = pfld->bcc(IB2,k,j,i);
          Real& bcc3 = pfld->bcc(IB3,k,j,i);
          data_sum[NHYDRO + 3] += vol(i)*0.5*bcc1*bcc1;
          data_sum[NHYDRO + 4] += vol(i)*0.5*bcc2*bcc2;
          data_sum[NHYDRO + 5] += vol(i)*0.5*bcc3*bcc3;
        }
      }
    }}
    pmb=pmb->next;

  }  // end loop over MeshBlocks

#ifdef MPI_PARALLEL
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &data_sum, (NHISTORY_VARS), MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&data_sum, &data_sum, (NHISTORY_VARS), MPI_ATHENA_REAL, MPI_SUM, 0,
               MPI_COMM_WORLD);
  }
#endif

  // only the master rank writes the file
  // create filename: "file_basename" + ".hst".  There is no file number.
  if (Globals::my_rank == 0) {
    std::string fname;
    fname.assign(output_params.file_basename);
    fname.append(".hst");

  // open file for output
    FILE *pfile;
    std::stringstream msg;
    if((pfile = fopen(fname.c_str(),"a")) == NULL){
      msg << "### FATAL ERROR in function [OutputType::HistoryFile]" << std::endl
          << "Output file '" << fname << "' could not be opened";
      throw std::runtime_error(msg.str().c_str());
    }

  // If this is the first output, write header
    if (output_params.file_number == 0) {
      fprintf(pfile,"# Athena++ history data\n"); // descriptor is first line
      fprintf(pfile,"# [1]=time     ");
      fprintf(pfile,"[2]=dt       ");
      fprintf(pfile,"[3]=mass     ");
      fprintf(pfile,"[4]=1-mom    ");
      fprintf(pfile,"[5]=2-mom    ");
      fprintf(pfile,"[6]=3-mom    ");
      fprintf(pfile,"[7]=1-KE     ");
      fprintf(pfile,"[8]=2-KE     ");
      fprintf(pfile,"[9]=3-KE     ");
      if (NON_BAROTROPIC_EOS) fprintf(pfile,"[10]=tot-E   ");
      if (MAGNETIC_FIELDS_ENABLED) {
        fprintf(pfile,"[11]=1-ME    ");
        fprintf(pfile,"[12]=2-ME    ");
        fprintf(pfile,"[13]=3-ME    ");
      }
      fprintf(pfile,"\n");                              // terminate line
    }

  // write history variables
    fprintf(pfile, output_params.data_format.c_str(), pm->time);
    fprintf(pfile, output_params.data_format.c_str(), pm->dt);
    for (int n=0; n<NHISTORY_VARS; ++n) {
      fprintf(pfile, output_params.data_format.c_str(), data_sum[n]);
    }
    fprintf(pfile,"\n"); // terminate line
    fclose(pfile);
  }

  vol.DeleteAthenaArray();
  return;
}
