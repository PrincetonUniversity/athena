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
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

#include "../athena.hpp"
#include "../mesh.hpp"
#include "outputs.hpp"

//======================================================================================
/*! \file formatted_table.cpp
 *  \brief writes output data as a formatted table.  Should not be used to output large
 *  3D data sets as this format is very memory intensive.  Most useful for 1D slices.
 *====================================================================================*/

//--------------------------------------------------------------------------------------
// FormattedTableOutput constructor

FormattedTableOutput::FormattedTableOutput(OutputParameters oparams, MeshBlock *pb)
  : OutputType(oparams,pb)
{
}

// destructor - not required for this derived class

//--------------------------------------------------------------------------------------
/*! \fn void FormattedTableOutput:::WriteOutputFile(OutputData *pod)
 *  \brief writes OutputData to file in tabular format using C style fprintf
 */

void FormattedTableOutput::WriteOutputFile(OutputData *pod)
{
  std::stringstream msg;

// create filename: "file_basename" + XXXX + ".tab", where XXXX = 4-digit file_number

  std::string fname;
  fname.assign(output_params.file_basename);
  fname.append(".");
  char number[5];
  sprintf(number,"%04d",output_params.file_number);
  fname.append(number);
  fname.append(".tab");

// open file for output
  FILE *pfile;
  if ((pfile = fopen(fname.c_str(),"w")) == NULL){
    msg << "### FATAL ERROR in function [FormattedTableOutput::WriteOutputFile]"
        << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// print header (data descriptor and transforms)

  fprintf(pfile,"%s",pod->data_header.descriptor.c_str());
  fprintf(pfile,"%s",pod->data_header.transforms.c_str());

// loop over all cells in data arrays

  for (int k=(pod->data_header.kl); k<=(pod->data_header.ku); ++k) {
  for (int j=(pod->data_header.jl); j<=(pod->data_header.ju); ++j) {
  for (int i=(pod->data_header.il); i<=(pod->data_header.iu); ++i) {

// write x1, x2, x3 indices and coordinates on start of new line

    if (pod->data_header.il != pod->data_header.iu) {
      fprintf(pfile,"%04d",i);
      fprintf(pfile,output_params.data_format.c_str(),pmy_block->x1v(i));
    }

    if (pod->data_header.jl != pod->data_header.ju) {
      fprintf(pfile,"%04d",j);
      fprintf(pfile,output_params.data_format.c_str(),pmy_block->x2v(j));
    }

    if (pod->data_header.kl != pod->data_header.ku) {
      fprintf(pfile,"%04d",k);
      fprintf(pfile,output_params.data_format.c_str(),pmy_block->x3v(k));
    }

// step through linked-list of variables and write data on same line

    OutputVariable *pvar = pod->pfirst_var;
    while (pvar != NULL) {
      for (int n=0; n<(pvar->pdata->GetDim4()); ++n) {
        fprintf( pfile, output_params.data_format.c_str(), (*pvar->pdata)(n,k,j,i) );
      }
      pvar = pvar->pnext;
    }

    fprintf(pfile,"\n"); // terminate line

  }}}

// close output file, increment file number, update time of last output, clean up

  fclose(pfile);
  output_params.file_number++;
  output_params.next_time += output_params.dt;

  return;
}
