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
//! \file formatted_table.cpp
//  \brief writes output data as a formatted table.  Should not be used to output large
//  3D data sets as this format is very slow and memory intensive.  Most useful for 1D
//  slices and/or sums.
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
#include "../mesh/mesh.hpp"
#include "outputs.hpp"
#include "../coordinates/coordinates.hpp"

//--------------------------------------------------------------------------------------
// FormattedTableOutput constructor
// destructor not required for this derived class

FormattedTableOutput::FormattedTableOutput(OutputParameters oparams)
  : OutputType(oparams) 
{
}

//--------------------------------------------------------------------------------------
//! \fn void FormattedTableOutput:::WriteOutputFile(Mesh *pm)
//  \brief writes OutputData to file in tabular format using C style fprintf

void FormattedTableOutput::WriteOutputFile(Mesh *pm)
{
  MeshBlock *pmb=pm->pblock;

  // Loop over MeshBlocks
  while (pmb != NULL) {
    int il=pmb->is; int iu=pmb->ie;
    int jl=pmb->js; int ju=pmb->je;
    int kl=pmb->ks; int ku=pmb->ke;

    // set ptrs to data in OutputData linked list
    LoadOutputData(pmb);

//  if (pod->data_header.ndata == 0) return;  // slice out of range, etc.

    // create filename: "file_basename"+ "."+"bloclid"+"."+"file_id"+"."+XXXXX+".tab",
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
    fname.append(".tab");

    // open file for output
    FILE *pfile;
    std::stringstream msg;
    if ((pfile = fopen(fname.c_str(),"w")) == NULL){
      msg << "### FATAL ERROR in function [FormattedTableOutput::WriteOutputFile]"
          <<std::endl<< "Output file '" <<fname<< "' could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    // print header
    fprintf(pfile,"# Athena++ data at time=%e",pm->time);
    fprintf(pfile,"  cycle=%d",pmb->pmy_mesh->ncycle);
    fprintf(pfile,"  variables=%s \n",output_params.variable.c_str());

//  fprintf(pfile,"%s",pod->data_header.transforms.c_str());

    // loop over all cells in data arrays
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
    for (int i=il; i<=iu; ++i) {

      // write x1, x2, x3 indices and coordinates on start of new line
      if (il != iu) {
        fprintf(pfile,"%04d",i);
        fprintf(pfile,output_params.data_format.c_str(),pmb->pcoord->x1v(i));
      }
      if (jl != ju) {
        fprintf(pfile," %04d",j);  // note extra space for formatting
        fprintf(pfile,output_params.data_format.c_str(),pmb->pcoord->x2v(j));
      }
      if (kl != ku) {
        fprintf(pfile," %04d",k);  // note extra space for formatting
        fprintf(pfile,output_params.data_format.c_str(),pmb->pcoord->x3v(k));
      }

      // step through linked-list of OutputData's and write each on same line
      OutputData *pdata = pfirst_data_;
      while (pdata != NULL) {
        for (int n=0; n<(pdata->data.GetDim4()); ++n) {
          fprintf(pfile, output_params.data_format.c_str(), pdata->data(n,k,j,i));
        }
        pdata = pdata->pnext;
      }

      fprintf(pfile,"\n"); // terminate line
    }}}

    // don't forget to close the output file and clean up ptrs to data in OutputData
    fclose(pfile);
    ClearOutputData();

    pmb=pmb->next;

  }  // end loop over MeshBlocks

  return;
}
