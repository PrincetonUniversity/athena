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
/*! \file outputs.cpp
 *  \brief implements functions for fluid data outputs
 *====================================================================================*/

//--------------------------------------------------------------------------------------
// DataNode constructor

DataNode::DataNode(AthenaArray<Real> *pinit_data, DataNodeHeader init_head)
{
  header = init_head;
  pdata = pinit_data;
  pnext = NULL;
}

//--------------------------------------------------------------------------------------
// DataNode destructor

DataNode::~DataNode()
{
  pdata->DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
// DataList constructor

DataList::DataList()
{
  pfirst_node = NULL;
  plast_node_ = NULL;
}

//--------------------------------------------------------------------------------------
// DataList destructor

DataList::~DataList()
{
}

//--------------------------------------------------------------------------------------
// OutputType constructor

OutputType::OutputType(OutputBlock out_blk, Block *pb)
{
  output_block = out_blk;
  pparent_block = pb;
  pnext = NULL; // Terminate linked list with NULL ptr
}

//--------------------------------------------------------------------------------------
// OutputType destructor

OutputType::~OutputType()
{
}

//--------------------------------------------------------------------------------------
// OutputList constructor

OutputList::OutputList(Block *pb)
{
  pparent_block = pb;
  pfirst_out_ = NULL;
}

// OutputList destructor - iterates through linked list of OutputTypes and deletes nodes

OutputList::~OutputList()
{
  OutputType *pout=pfirst_out_;
  while(pout != NULL) {
    OutputType *pold_out = pout;
    pout = pout->pnext;
    delete pold_out;
  }
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputList::InitOutputs()
 *  \brief Creates linked list of OutputTypes based on <ouput> blocks in input file.
 */

void OutputList::InitOutputs(ParameterInput *pin)
{
  std::stringstream msg;
  InputBlock *pib = pin->pfirst_block;
  OutputType *pnew_out;
  OutputType *plast = pfirst_out_;

// loop over input block names.  Find those that start with "output", read parameters,
// and construct linked list of OutputTypes.

  while (pib != NULL) {
    if (pib->block_name.compare(0,6,"output") == 0) {
      OutputBlock ob;  // define temporary OutputBlock struct

// extract integer number of output block.  Save name and number 

      std::string outn = pib->block_name.substr(6); // 6 because starts at 0!
      ob.block_number = atoi(outn.c_str());
      ob.block_name.assign(pib->block_name);

// set time of last output, time between outputs

      ob.last_time = pin->GetOrAddReal(ob.block_name,"last_time",0.0);
      ob.dt = pin->GetReal(ob.block_name,"dt");

// set file number, basename, id, and format

      ob.file_number = pin->GetOrAddInteger(ob.block_name,"file_number",0);
      ob.file_basename = pin->GetString("job","problem_id");
      char define_id[10];
      sprintf(define_id,"out%d",ob.block_number);  // default id="outN"
      ob.file_id = pin->GetOrAddString(ob.block_name,"id",define_id);
      ob.file_format = pin->GetString(ob.block_name,"file_format");

// set output variable and optional data format string used in formatted writes

      if (ob.file_format.compare("hst") != 0)
        ob.variable = pin->GetString(ob.block_name,"variable");
      ob.data_format = pin->GetOrAddString(ob.block_name,"data_format","%12e.5");

// Construct new OutputType according to file format

      if (ob.file_format.compare("tab") == 0) {
        pnew_out = new FormattedTableOutput(ob, pparent_block);
      } else if (ob.file_format.compare("hst") == 0) {
        pnew_out = new HistoryOutput(ob, pparent_block);
      } else {
        msg << "### FATAL ERROR in function [OutputList::InitOutputs]"
            << std::endl << "Unrecognized file format = '" << ob.file_format 
            << "' in output block '" << pib->block_name << "'";
        throw std::runtime_error(msg.str().c_str());
      }

// Add new type as node in linked list 

      if (pfirst_out_ == NULL) {
        pfirst_out_ = pnew_out;
      } else {
        plast->pnext = pnew_out;
      }
      plast = pnew_out;
    }
    pib = pib->pnext;  // move to next input block name
  }
}

//--------------------------------------------------------------------------------------
/*! \fn void DataList::InsertNode()
 *  \brief
 */

void DataList::InsertNode(AthenaArray<Real> *pdata, DataNodeHeader init_head)
{
  DataNode *pnew_node = new DataNode(pdata, init_head);

  if (pfirst_node == NULL)
    pfirst_node = plast_node_ = pnew_node;
  else {
    plast_node_->pnext = pnew_node;
    plast_node_ = pnew_node;
  }
}

//--------------------------------------------------------------------------------------
/*! \fn DataList* OutputType::LoadDataList()
 *  \brief selects data to be included in output data container (DataList)
 */

DataList* OutputType::LoadDataList()
{
  DataNodeHeader node_header;
  DataList *pdl = new DataList;
  Fluid *pf = pparent_block->pfluid;;

  pdl->header.il = pparent_block->is;
  pdl->header.iu = pparent_block->ie;
  pdl->header.jl = pparent_block->js;
  pdl->header.ju = pparent_block->je;
  pdl->header.kl = pparent_block->ks;
  pdl->header.ku = pparent_block->ke;

  node_header.type = "SCALARS";
  node_header.type = "dens";
  pdl->InsertNode(pf->u.ShallowCopy(IDN,1),node_header);

  node_header.type = "SCALARS";
  node_header.type = "Etot";
  pdl->InsertNode(pf->u.ShallowCopy(IEN,1),node_header);

  node_header.type = "VECTORS";
  node_header.type = "mom";
  pdl->InsertNode(pf->u.ShallowCopy(IM1,3),node_header);

  return pdl;
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputType::ComputeOutputData()
 *  \brief 
 */

void OutputType::ComputeDataList()
{
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputList::MakeOutputs()
 *  \brief scans through OutputList and makes any outputs needed.
 */

void OutputList::MakeOutputs()
{
  OutputType* pout=pfirst_out_;
  Mesh* pm = pparent_block->pparent_domain->pparent_mesh;

  while (pout != NULL) {
    if ((pm->time == pm->start_time) || 
        (pm->time >= (pout->output_block.last_time + pout->output_block.dt)) ||
        (pm->time >= pm->tlim)) {

//      pout->ComputeFromMesh(pm); // Load data and perform transformations
 //     pout->WriteOutputData();  // Write data to file
  //    pout->Update();  // Update number and time since last output
      // std::cout << "Making output " << pout->BlockNumber() << std::endl;
    
/*
      std::cout << std::setprecision(6) << pout->output_block.last_time << std::endl;
      std::cout << std::setprecision(6) << pout->output_block.dt << std::endl;
      std::cout << pout->output_block.block_number << std::endl;
      std::cout << pout->output_block.file_number << std::endl;
      std::cout << pout->output_block.block_name << std::endl;
      std::cout << pout->output_block.file_basename << std::endl;
      std::cout << pout->output_block.file_id << std::endl;
      std::cout << pout->output_block.variable << std::endl;
      std::cout << pout->output_block.file_format << std::endl;
      std::cout << pout->output_block.data_format << std::endl;
*/

      pout->WriteOutputData();

    }
    pout = pout->pnext;
  }
}
