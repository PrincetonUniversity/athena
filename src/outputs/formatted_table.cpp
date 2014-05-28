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
#include "../datablock.hpp"
#include "output.hpp"

//======================================================================================
/*! \file formatted_table.cpp
 *  \brief writes output data as a formatted table.  Should not be used to output large
 *  3D data sets as this format is very memory intensive.  Most useful for 1D slices.
 *====================================================================================*/

//--------------------------------------------------------------------------------------
// FormattedTableOutput constructor

FormattedTableOutput::FormattedTableOutput(InputBlock *pin_block, Mesh *pm)
  : Output(pin_block,pm)
{
}

//--------------------------------------------------------------------------------------
/*! \fn void TabularOutput::ComputeFromMesh(Mesh *pm)
 *  \brief Performs any operations on Mesh and stores results in the output's DataBlock
 */

void FormattedTableOutput::ComputeFromMesh(Mesh *pm)
{

// Read data from Mesh
  pdata = LoadFluidFromMesh(pm);

// Determine initial dimensions of output
  pdata->GetNode(3)->GetDimensions(nx3,nx2,nx1);

// Loop over transformation list
  DataBlockTransform* pdbt = ptrans;
  while (pdbt != NULL) {
    pdbt->Transform(pdata);
    pdbt = pdbt->GetNext();
  }

//Determine final dimensions of output
  pdata->GetNode(3)->GetDimensions(nx3,nx2,nx1);

}

//--------------------------------------------------------------------------------------
/*! \fn
 *  \brief
 */

void FormattedTableOutput::WriteNodeData(DataNode *pdn, FILE *pfile, std::string fmt)
{
  int is,ie,js,je,ks,ke;
  pdn->GetRanges(is,ie,js,je,ks,ke);

  for (int n=0; n<pdn->GetData()->GetDim4(); ++n){
    for (int k=ks; k<=ke; ++k){
      for (int j=js; j<=je; ++j){
        for (int i=is; i<=ie; ++i){
          // Not certain this wil autovectorize well ???
          fprintf( pfile, fmt.c_str(), (*pdn->GetData())(n,k,j,i) );
        }
        fprintf(pfile,"\n");
      }}}

}

//--------------------------------------------------------------------------------------
/*! \fn void TabularOutput:::Write()
 *  \brief writes DataBlock to file in tabular format
 *         uses c style printf for the moment
 */

void FormattedTableOutput::Write()
{
  std::stringstream msg;

// create filename
  std::stringstream ftmp;
  ftmp << "Sod." << Number() << ".tab";//hardwire for now
  std::string filename=ftmp.str();

// open file for output
  FILE *pfile;
  if((pfile = fopen(filename.c_str(),"w")) == NULL){
    msg << "### FATAL ERROR in function [TabularOutput::Write]" << std::endl
        << "Output file '" << filename << "' could not be opened";
    throw std::runtime_error(msg.str().c_str());
  }

// write header description
  if(pdata->GetDescriptor().compare("") != 0)
    fprintf(pfile,"# %s\n",pdata->GetDescriptor().c_str());
// write dimensions
  if (nx1 > 1) fprintf(pfile,"# Nx1 = %d\n",nx1);
  if (nx2 > 1) fprintf(pfile,"# Nx2 = %d\n",nx2);
  if (nx3 > 1) fprintf(pfile,"# Nx3 = %d\n",nx3);

// Loop over headers and write row headings
  int col_cnt = 0;
  fprintf(pfile,"#");
  if (nx1 > 1) {
    fprintf(pfile," [%d]=%s",col_cnt,pdata->GetNode(0)->GetVariableTypes().c_str());
    col_cnt++;
  }
  if (nx2 > 1) {
    fprintf(pfile," [%d]=%s",col_cnt,pdata->GetNode(1)->GetVariableTypes().c_str());
    col_cnt++;
  }
  if (nx3 > 1) {
    fprintf(pfile," [%d]=%s",col_cnt,pdata->GetNode(2)->GetVariableTypes().c_str());
    col_cnt++;
  }
  DataNode* pdn = pdata->GetNode(2);
  while (pdn->GetNext() != NULL){
    pdn = pdn->GetNext();
    fprintf(pfile," [%d]=%s",col_cnt,pdn->GetVariableTypes().c_str());
    col_cnt++;
  }
  fprintf(pfile,"\n");

// First three nodes always contain mesh data: write if more than one element
  if (nx1 > 1) WriteNodeData(pdata->GetNode(0),pfile,output_block.dat_fmt);
  if (nx2 > 1) WriteNodeData(pdata->GetNode(1),pfile,output_block.dat_fmt);
  if (nx3 > 1) WriteNodeData(pdata->GetNode(2),pfile,output_block.dat_fmt);

// Loop over remaining data nodes and write output
  pdn = pdata->GetNode(2);
  while (pdn->GetNext() != NULL){
    pdn = pdn->GetNext();
    WriteNodeData(pdn,pfile,output_block.dat_fmt);
  }
  fclose(pfile);
  delete pdata;  // Release any memory allocated this output cycle
}
