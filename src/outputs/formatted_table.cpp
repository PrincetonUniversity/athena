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
#include "../mesh.hpp"
#include "outputs.hpp"

//======================================================================================
/*! \file formatted_table.cpp
 *  \brief writes output data as a formatted table.  Should not be used to output large
 *  3D data sets as this format is very memory intensive.  Most useful for 1D slices.
 *====================================================================================*/

//--------------------------------------------------------------------------------------
// FormattedTableOutput constructor

FormattedTableOutput::FormattedTableOutput(OutputBlock out_blk, Block *pb)
  : OutputType(out_blk,pb)
{
}

//--------------------------------------------------------------------------------------
/*! \fn void FormattedTableOutput:::WriteOutputData()
 *  \brief writes DataBlock to file in tabular format using C style fprintf
 */

void FormattedTableOutput::WriteOutputData()
{
  std::stringstream msg;
  OutputData *pod;

// create OutputData, apply transforms (slices, sums, etc)

  pod = LoadOutputData();
  TransformOutputData(pod);

// create filename
  std::string fname;
  fname.assign(output_block.file_basename);
  fname.append(".");
  char number[5];
  sprintf(number,"%04d",output_block.file_number);
  fname.append(number);
  fname.append(".tab");

// open file for output
  FILE *pfile;
  if ((pfile = fopen(fname.c_str(),"w")) == NULL){
    msg << "### FATAL ERROR in function [FormattedTableOutput::WriteOutputData]"
        << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

// print header (data descriptor and transforms)

  fprintf(pfile,"%s",pod->header.descriptor.c_str());
  fprintf(pfile,"%s",pod->header.transforms.c_str());

// loop over all cells in data arrays

  for (int k=(pod->header.kl); k<=(pod->header.ku); ++k) {
  for (int j=(pod->header.jl); j<=(pod->header.ju); ++j) {
  for (int i=(pod->header.il); i<=(pod->header.iu); ++i) {

// write x1, x2, x3 indices and coordinates on start of new line

    if (pod->header.il != pod->header.iu) {
      fprintf(pfile,"%04d",i);
      fprintf(pfile,output_block.data_format.c_str(),pparent_block->x1v(i));
    }

    if (pod->header.jl != pod->header.ju) {
      fprintf(pfile,"%04d",j);
      fprintf(pfile,output_block.data_format.c_str(),pparent_block->x2v(j));
    }

    if (pod->header.kl != pod->header.ku) {
      fprintf(pfile,"%04d",k);
      fprintf(pfile,output_block.data_format.c_str(),pparent_block->x3v(k));
    }

// step through linked-list of data nodes and write data on same line

    OutputDataNode *pnode = pod->pfirst_node;
    while (pnode != NULL) {
      for (int n=0; n<(pnode->pdata->GetDim4()); ++n) {
        fprintf( pfile, output_block.data_format.c_str(), (*pnode->pdata)(n,k,j,i) );
      }
      pnode = pnode->pnext;
    }

    fprintf(pfile,"\n"); // terminate line

  }}}

// close output file, increment file number, update time of last output, clean up

  fclose(pfile);
  output_block.file_number++;
  output_block.next_time += output_block.dt;
  delete pod;

  return;
}
