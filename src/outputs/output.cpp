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
/*! \file output.cpp
 *  \brief implements functions for fluid data outputs
 *====================================================================================*/

//--------------------------------------------------------------------------------------
// OutputBlock constructor
// This needs to be made more extensible since new outputs may need new specifiers

OutputBlock::OutputBlock(InputBlock *pin_block)
{

  if (pin_block->block_name.compare(0,6,"output") == 0) {
// extract integer number of output block.  Save name and number 
    std::string outn = pin_block->block_name.substr(6); // 6 because starts at 0!
    block_number = atoi(outn.c_str());
    block_name.assign(pin_block->block_name);
// set time of last output, time between outputs
    last_time = pin_block->GetOrAddRealInThisBlock("last_time",0.0);
    dt = pin_block->GetRealInThisBlock("dt");
// Output cycle # for this block
    number = pin_block->GetOrAddIntInThisBlock("num",0);
// Specifier of variable(s) to be outputted
    out = pin_block->GetOrAddStringInThisBlock("out","prim");
// Read dat_fmt or set defaul, then add white space
    dat_fmt = pin_block->GetOrAddStringInThisBlock("dat_fmt","%12.8e");
// Add white space to dat_fmt
    dat_fmt = " " + dat_fmt;
// Read in optional transform
// Currently allows for single transofrm from preset list -- multiple transforms
// (including user defined) should be added
    trans = pin_block->GetOrAddStringInThisBlock("trans"," ");

  }
}
// OutputBlock destructor
OutputBlock::~OutputBlock()
{

}

//--------------------------------------------------------------------------------------
// Output constructor

Output::Output(InputBlock *pin_block, Mesh *pm)
  : output_block(pin_block)
{

// Initialize Transformation list
// Currently allows for single transofrm from preset list -- multiple transforms
// (including user defined) should be added
  if (output_block.trans.compare("Slice1Dx1") == 0) {
    Real x2 = pin_block->GetOrAddRealInThisBlock("slcx2",0.0);
    Real x3 = pin_block->GetOrAddRealInThisBlock("slcx3",0.0);
    ptrans = new Slice1Dx1(pm,x2,x3);
  } else {
    ptrans = NULL;
  }
  pnext = NULL; // Initialize pnext to NULL
}

//--------------------------------------------------------------------------------------
// OutputList constructor
OutputList::OutputList(ParameterInput *pin, Mesh *pm)
{
  std::stringstream msg;
  InputBlock *pin_block = pin->GetFirstBlock();
  Output *pout, *plast;

  phead = NULL;
// loop over input block names.  Find those that start with "output", and create new
// Output node.
  while (pin_block != NULL) {
    if (pin_block->block_name.compare(0,6,"output") == 0) {
      plast = pout;
      std::string fmt = pin->GetString(pin_block->block_name,"out_fmt");
      if (fmt.compare("hst") == 0) {
        pout = new HistoryOutput(pin_block,pm);
      } else if (fmt.compare("tab") == 0) {
        pout = new FormattedTableOutput(pin_block,pm);
      } else {
        msg << "### FATAL ERROR in function [Output_example::InitializeOutputs]"
            << std::endl << "Unrecognized out_fmt = '" << fmt << "' in output block '"
            << pin_block->block_name << "'";
        throw std::runtime_error(msg.str().c_str());
      }

// if this is the first output in list, save pointer to head node
      if (phead == NULL)
        phead = pout;
      else {
        plast->SetNext(pout); // add to end of list
      }
    }
    pin_block = pin_block->pnext;  // move to next input block name
  }

}

// OutputList destructor
OutputList::~OutputList()
{
  Output *pnext, *pout = phead;

// delete Output nodes, starting with head node
  pnext = pout->GetNext();
  while(pnext != NULL) {
    delete pout;
    pout = pnext;
    pnext = pout->GetNext();
  }
  delete pout;
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputList::CheckForOutputs(Mesh *pm)
 *  \brief Scans through OutputList and checks whether and output is
 *  needed in any Output node.
 */

void OutputList::CheckForOutputs(Mesh *pm)
{
  Output* pout = phead;

  while (pout != NULL) {
    if ( (pm->time == pm->start_time) || (pm->time >= pout->NextTime()) ||
         (pm->time == pm->tlim) ) {

      pout->ComputeFromMesh(pm); // Load data and perform transformations
      pout->Write();  // Write data to file
      pout->Update();  // Update number and time since last output
      // std::cout << "Making output " << pout->BlockNumber() << std::endl;
    }
    pout = pout->GetNext();

  }
}

//--------------------------------------------------------------------------------------
/*! \fn DataBlock* Output::LoadFluidFromMesh(Mesh *pm)
 *  \brief Performs shallow copies from Mesh and stores results in data block
 */

DataBlock* Output::LoadFluidFromMesh(Mesh *pm)
{
  DataBlock *pdb = new DataBlock;
// Set time information
  pdb->SetTime(pm->time);
  pdb->SetTimeStep(pm->dt);
  pdb->SetCycleNumber(pm->ncycle);
// Make shallow copies to the position arrays in Block
  pdb->InsertNode(pm->pdomain->pblock->x1v.ShallowCopy(0));
  pdb->LastNode()->SetVariableTypes("x1");
  pdb->LastNode()->SetRanges(pm->pdomain->pblock->is,pm->pdomain->pblock->ie);
  pdb->InsertNode(pm->pdomain->pblock->x2v.ShallowCopy(0));
  pdb->LastNode()->SetVariableTypes("x2");
  pdb->LastNode()->SetRanges(pm->pdomain->pblock->js,pm->pdomain->pblock->je);
  pdb->InsertNode(pm->pdomain->pblock->x3v.ShallowCopy(0));
  pdb->LastNode()->SetVariableTypes("x3");
  pdb->LastNode()->SetRanges(pm->pdomain->pblock->ks,pm->pdomain->pblock->ke);
// Initialize common header for initializing all DataNodes
  DataNodeHeader hinit;
  hinit.is = pm->pdomain->pblock->is;
  hinit.ie = pm->pdomain->pblock->ie;
  hinit.js = pm->pdomain->pblock->js;
  hinit.je = pm->pdomain->pblock->je;
  hinit.ks = pm->pdomain->pblock->ks;
  hinit.ke = pm->pdomain->pblock->ke;
// Make shallow copies to the arrays of each variable and assign each to a node
  if (output_block.out.compare("prim") == 0) {
    pdb->InsertNode(pm->pdomain->pblock->pfluid->w.ShallowCopy(IDN),hinit);
    pdb->LastNode()->SetVariableTypes("d");
    pdb->InsertNode(pm->pdomain->pblock->pfluid->w.ShallowCopy(IVX),hinit);
    pdb->LastNode()->SetVariableTypes("V1");
    pdb->InsertNode(pm->pdomain->pblock->pfluid->w.ShallowCopy(IVY),hinit);
    pdb->LastNode()->SetVariableTypes("V2");
    pdb->InsertNode(pm->pdomain->pblock->pfluid->w.ShallowCopy(IVZ),hinit);
    pdb->LastNode()->SetVariableTypes("V3");
    pdb->InsertNode(pm->pdomain->pblock->pfluid->w.ShallowCopy(IEN),hinit);
    pdb->LastNode()->SetVariableTypes("P");
  } else if (output_block.out.compare("cons") == 0) {
    pdb->InsertNode(pm->pdomain->pblock->pfluid->u.ShallowCopy(IDN),hinit);
    pdb->LastNode()->SetVariableTypes("d");
    pdb->InsertNode(pm->pdomain->pblock->pfluid->u.ShallowCopy(IM1),hinit);
    pdb->LastNode()->SetVariableTypes("M1");
    pdb->InsertNode(pm->pdomain->pblock->pfluid->u.ShallowCopy(IM2),hinit);
    pdb->LastNode()->SetVariableTypes("M2");
    pdb->InsertNode(pm->pdomain->pblock->pfluid->u.ShallowCopy(IM3),hinit);
    pdb->LastNode()->SetVariableTypes("M3");
    pdb->InsertNode(pm->pdomain->pblock->pfluid->u.ShallowCopy(IEN),hinit);
    pdb->LastNode()->SetVariableTypes("E");
  } else if (output_block.out.compare("d") == 0) {
    pdb->InsertNode(pm->pdomain->pblock->pfluid->w.ShallowCopy(IDN),hinit);
    pdb->LastNode()->SetVariableTypes("d");
  } else if (output_block.out.compare("V1") == 0) {
    pdb->InsertNode(pm->pdomain->pblock->pfluid->w.ShallowCopy(IVX),hinit);
    pdb->LastNode()->SetVariableTypes("V1");
  } else if (output_block.out.compare("V2") == 0) {
    pdb->InsertNode(pm->pdomain->pblock->pfluid->w.ShallowCopy(IVY),hinit);
    pdb->LastNode()->SetVariableTypes("V2");
  } else if (output_block.out.compare("V3") == 0) {
    pdb->InsertNode(pm->pdomain->pblock->pfluid->w.ShallowCopy(IVZ),hinit);
    pdb->LastNode()->SetVariableTypes("V3");
  } else if (output_block.out.compare("P") == 0) {
    pdb->InsertNode(pm->pdomain->pblock->pfluid->w.ShallowCopy(IEN),hinit);
    pdb->LastNode()->SetVariableTypes("P");
  } else if (output_block.out.compare("M1") == 0) {
    pdb->InsertNode(pm->pdomain->pblock->pfluid->u.ShallowCopy(IM1),hinit);
    pdb->LastNode()->SetVariableTypes("M1");
  } else if (output_block.out.compare("M2") == 0) {
    pdb->InsertNode(pm->pdomain->pblock->pfluid->u.ShallowCopy(IM2),hinit);
    pdb->LastNode()->SetVariableTypes("M2");
  } else if (output_block.out.compare("M3") == 0) {
    pdb->InsertNode(pm->pdomain->pblock->pfluid->u.ShallowCopy(IM3),hinit);
    pdb->LastNode()->SetVariableTypes("M3");
  } else if (output_block.out.compare("E") == 0) {
    pdb->InsertNode(pm->pdomain->pblock->pfluid->u.ShallowCopy(IEN),hinit);
    pdb->LastNode()->SetVariableTypes("E");
  } else {
    // Lazy exception handling for now
    std::cout << "Warning: output variable out='" << output_block.out
              << "' in output block <" << output_block.block_name
              << "> is not recognized."  << std::endl;
  }
  return pdb;

}
