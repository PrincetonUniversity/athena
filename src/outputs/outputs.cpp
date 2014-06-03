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
// OutputDataNode constructor

OutputDataNode::OutputDataNode(AthenaArray<Real> *parray, OutputDataNodeHeader head)
{
  header = head;
  pdata = parray;
  pnext = NULL;
  pprev = NULL;
}

//--------------------------------------------------------------------------------------
// OutputDataNode destructor

OutputDataNode::~OutputDataNode()
{
  if (!pdata->IsShallowCopy()) pdata->DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
// DataList constructor

OutputData::OutputData()
{
  pfirst_node = NULL;
  plast_node = NULL;
}

//--------------------------------------------------------------------------------------
// DataList destructor

OutputData::~OutputData()
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
  Block *pb = pparent_block;

// loop over input block names.  Find those that start with "output", read parameters,
// and construct linked list of OutputTypes.

  while (pib != NULL) {
    if (pib->block_name.compare(0,6,"output") == 0) {
      int create_output = 1;
      OutputBlock ob;  // define temporary OutputBlock struct

// extract integer number of output block.  Save name and number 

      std::string outn = pib->block_name.substr(6); // 6 because starts at 0!
      ob.block_number = atoi(outn.c_str());
      ob.block_name.assign(pib->block_name);

// set time of last output, time between outputs

      ob.next_time = pin->GetOrAddReal(ob.block_name,"next_time",0.0);
      ob.dt = pin->GetReal(ob.block_name,"dt");

// set file number, basename, id, and format

      ob.file_number = pin->GetOrAddInteger(ob.block_name,"file_number",0);
      ob.file_basename = pin->GetString("job","problem_id");
      char define_id[10];
      sprintf(define_id,"out%d",ob.block_number);  // default id="outN"
      ob.file_id = pin->GetOrAddString(ob.block_name,"id",define_id);
      ob.file_format = pin->GetString(ob.block_name,"file_format");

// read slicing options.  Check that slice is within range

      if (pin->ParameterExists(ob.block_name,"x1_slice")) {
        Real x1 = pin->GetReal(ob.block_name,"x1_slice");
        if (x1 >= pb->block_size.x1min && x1 < pb->block_size.x1max) {
          for (int i=pb->is+1; i<=pb->ie+1; ++i) {
            if (pb->x1f(i) > x1) {
              ob.islice = i-1;
              break;
            }
          }
        } else {
          create_output=0;;
        }
      } else {ob.islice = -999;}

      if (pin->ParameterExists(ob.block_name,"x2_slice")) {
        Real x2 = pin->GetReal(ob.block_name,"x2_slice");
        if (x2 >= pb->block_size.x2min && x2 < pb->block_size.x2max) {
          for (int j=pb->js+1; j<=pb->je+1; ++j) {
            if (pb->x2f(j) > x2) {
              ob.jslice = j-1;
              break;
            }
          }
        } else {
          create_output=0;;
        }
      } else {ob.jslice = -999;}

      if (pin->ParameterExists(ob.block_name,"x3_slice")) {
        Real x3 = pin->GetReal(ob.block_name,"x3_slice");
        if (x3 >= pb->block_size.x3min && x3 < pb->block_size.x3max) {
          for (int k=pb->ks+1; k<=pb->ke+1; ++k) {
            if (pb->x3f(k) > x3) {
              ob.kslice = k-1;
              break;
            }
          }
        } else {
          create_output=0;;
        }
      } else {ob.kslice = -999;}

// read sum options.  Check for conflicts with slicing.

      if (pin->ParameterExists(ob.block_name,"x1_sum")) {
        if (pin->ParameterExists(ob.block_name,"x1_slice")) {
          msg << "### FATAL ERROR in function [OutputList::InitOutputs]"
              << std::endl << "Cannot request both slice and sum along x1-direction"
              << " in output block '" << ob.block_name << "'" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        } else {
          ob.isum = pin->GetInteger(ob.block_name,"x1_sum");;
        }
      } else {ob.isum = 0;}

      if (pin->ParameterExists(ob.block_name,"x2_sum")) {
        if (pin->ParameterExists(ob.block_name,"x2_slice")) {
          msg << "### FATAL ERROR in function [OutputList::InitOutputs]"
              << std::endl << "Cannot request both slice and sum along x2-direction"
              << " in output block '" << ob.block_name << "'" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        } else {
          ob.jsum = pin->GetInteger(ob.block_name,"x2_sum");;
        }
      } else {ob.jsum = 0;}

      if (pin->ParameterExists(ob.block_name,"x3_sum")) {
        if (pin->ParameterExists(ob.block_name,"x3_slice")) {
          msg << "### FATAL ERROR in function [OutputList::InitOutputs]"
              << std::endl << "Cannot request both slice and sum along x3-direction"
              << " in output block '" << ob.block_name << "'" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        } else {
          ob.ksum = pin->GetInteger(ob.block_name,"x3_sum");;
        }
      } else {ob.ksum = 0;}

      if (create_output) {  // skip output if slice not in range
// set output variable and optional data format string used in formatted writes

        if (ob.file_format.compare("hst") != 0) {
          ob.variable = pin->GetString(ob.block_name,"variable");
        }
        ob.data_format = pin->GetOrAddString(ob.block_name,"data_format","%12e.5");

// Construct new OutputType according to file format

        if (ob.file_format.compare("tab") == 0) {
          pnew_out = new FormattedTableOutput(ob, pparent_block);
        } else if (ob.file_format.compare("hst") == 0) {
          pnew_out = new HistoryOutput(ob, pparent_block);
        } else {
          msg << "### FATAL ERROR in function [OutputList::InitOutputs]"
              << std::endl << "Unrecognized file format = '" << ob.file_format 
              << "' in output block '" << ob.block_name << "'" << std::endl;
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
    }
    pib = pib->pnext;  // move to next input block name
  }
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputData::AppendNode()
 *  \brief
 */

void OutputData::AppendNode(AthenaArray<Real> *parray, OutputDataNodeHeader head)
{
  OutputDataNode *pnew_node = new OutputDataNode(parray, head);

  if (pfirst_node == NULL)
    pfirst_node = plast_node = pnew_node;
  else {
    pnew_node->pprev = plast_node;
    plast_node->pnext = pnew_node;
    plast_node = pnew_node;
  }
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputData::ReplaceNode()
 *  \brief
 */

void OutputData::ReplaceNode(OutputDataNode *pold, OutputDataNode *pnew) 
{
  if (pold == pfirst_node) {
    pfirst_node = pnew;
    if (pold->pnext != NULL) {    // there is another node in the list 
      pnew->pnext = pold->pnext;
      pnew->pnext->pprev = pnew;
    } else {                      // there is only one node in the list
      plast_node = pnew;
    }
  } else if (pold == plast_node) {
    plast_node = pnew;
    pnew->pprev = pold->pprev;
    pnew->pprev->pnext = pnew;
  } else {
    pnew->pnext = pold->pnext;
    pnew->pprev = pold->pprev;
    pnew->pprev->pnext = pnew;
    pnew->pnext->pprev = pnew;
  }
  delete pold;
}

//--------------------------------------------------------------------------------------
/*! \fn OutputData* OutputType::LoadOutputData()
 *  \brief selects data to be included in output data container (OuputData)
 */

OutputData* OutputType::LoadOutputData()
{
  OutputDataNodeHeader node_header;
  OutputData *pod = new OutputData;
  Fluid *pf = pparent_block->pfluid;;
  std::stringstream str;

  str << "# Athena++ tabular data at time=" 
      << pparent_block->pparent_domain->pparent_mesh->time << " cycle="
      << pparent_block->pparent_domain->pparent_mesh->ncycle << std::endl;
  pod->header.descriptor.append(str.str());
  pod->header.il = pparent_block->is;
  pod->header.iu = pparent_block->ie;
  pod->header.jl = pparent_block->js;
  pod->header.ju = pparent_block->je;
  pod->header.kl = pparent_block->ks;
  pod->header.ku = pparent_block->ke;

  node_header.type = "SCALARS";
  node_header.type = "dens";
  pod->AppendNode(pf->u.ShallowCopy(IDN,1),node_header);

  node_header.type = "SCALARS";
  node_header.type = "Etot";
  pod->AppendNode(pf->u.ShallowCopy(IEN,1),node_header);

  node_header.type = "VECTORS";
  node_header.type = "mom";
  pod->AppendNode(pf->u.ShallowCopy(IM1,3),node_header);

  return pod;
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputType::ComputeOutputData()
 *  \brief 
 */

void OutputType::ComputeOutputData(OutputData *pod)
{
  if (output_block.kslice != -999) {
    Slice(pod,3);
  }
  if (output_block.jslice != -999) {
    Slice(pod,2);
  }
  if (output_block.islice != -999) {
    Slice(pod,1);
  }
  return;
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputType::Slice(OutputData* pod, int dim)
 *  \brief
 */

void OutputType::Slice(OutputData* pod, int dim)
{
  OutputDataNodeHeader node_header;
  AthenaArray<Real> *pslice;
  std::stringstream str;

// modify OutputData header

  if (dim == 3) {
    str << "# Slice at x3= " << pparent_block->x3v(output_block.kslice) << std::endl;
    pod->header.kl = 0;
    pod->header.ku = 0;
  } else if (dim == 2) {
    str << "# Slice at x2= " << pparent_block->x2v(output_block.jslice) << std::endl;
    pod->header.jl = 0;
    pod->header.ju = 0;
  } else {
    str << "# Slice at x1= " << pparent_block->x1v(output_block.islice) << std::endl;
    pod->header.il = 0;
    pod->header.iu = 0;
  }
  pod->header.descriptor.append(str.str());

// For each node in OutputData linked list, slice arrays containing output data  

  OutputDataNode *pdn;
  pdn = pod->pfirst_node;

  while (pdn != NULL) {
    node_header = pdn->header;
    int nx4 = pdn->pdata->GetDim4();
    int nx3 = pdn->pdata->GetDim3();
    int nx2 = pdn->pdata->GetDim2();
    int nx1 = pdn->pdata->GetDim1();
    pslice = new AthenaArray<Real>;

// Loop over variables and dimensions, extract slice

    if (dim == 3) {
      pslice->NewAthenaArray(nx4,1,nx2,nx1);
      for (int n=0; n<nx4; ++n){
      for (int j=(pod->header.jl); j<=(pod->header.ju); ++j){
        for (int i=(pod->header.il); i<=(pod->header.iu); ++i){
          (*pslice)(n,0,j,i) = (*pdn->pdata)(n,output_block.kslice,j,i);
        }
      }}
    } else if (dim == 2) {
      pslice->NewAthenaArray(nx4,nx3,1,nx1);
      for (int n=0; n<nx4; ++n){
      for (int k=(pod->header.kl); k<=(pod->header.ku); ++k){
        for (int i=(pod->header.il); i<=(pod->header.iu); ++i){
          (*pslice)(n,k,0,i) = (*pdn->pdata)(n,k,output_block.jslice,i);
        }
      }}
    } else {
      pslice->NewAthenaArray(nx4,nx3,nx2,1);
      for (int n=0; n<nx4; ++n){
      for (int k=(pod->header.kl); k<=(pod->header.ku); ++k){
        for (int j=(pod->header.jl); j<=(pod->header.ju); ++j){
          (*pslice)(n,k,j,0) = (*pdn->pdata)(n,k,j,output_block.islice);
        }
      }}
    }

    OutputDataNode *pnew = new OutputDataNode(pslice, node_header);
    pod->ReplaceNode(pdn,pnew);
    pdn = pdn->pnext;
  }
 
  return;
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
        (pm->time >= pout->output_block.next_time) ||
        (pm->time >= pm->tlim)) {

      std::cout << std::setprecision(6) << pout->output_block.next_time << std::endl;
      std::cout << std::setprecision(6) << pout->output_block.dt << std::endl;
      std::cout << pout->output_block.block_number << std::endl;
      std::cout << pout->output_block.file_number << std::endl;
      std::cout << pout->output_block.block_name << std::endl;
      std::cout << pout->output_block.file_basename << std::endl;
      std::cout << pout->output_block.file_id << std::endl;
      std::cout << pout->output_block.variable << std::endl;
      std::cout << pout->output_block.file_format << std::endl;
      std::cout << pout->output_block.data_format << std::endl;
std::cout << "islice=" << pout->output_block.islice
          << "jslice=" << pout->output_block.jslice
          << "kslice=" << pout->output_block.kslice << std::endl;
std::cout << "isum=" << pout->output_block.isum
          << "jsum=" << pout->output_block.jsum
          << "ksum=" << pout->output_block.ksum << std::endl;

      pout->WriteOutputData();

    }
    pout = pout->pnext;
  }
}
