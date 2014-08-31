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
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../fluid/fluid.hpp"
#include "outputs.hpp"

//======================================================================================
/*! \file outputs.cpp
 *  \brief implements functions for Athena++ outputs
 *
 * The number and types of outputs are all controlled by the number and values of
 * parameters specified in <outputN> blocks in the input file.  Each output block must
 * be labelled by a unique integer "N".  Following the convention of the parser
 * implemented in the ParameterInput class, a second output block with the same integer
 * "N" of an earlier block will silently overwrite the values read by the first block.
 * The numbering of the output blocks does not need to be consecutive, and blocks may
 * appear in any order in the input file.  Moreover, unlike the C version of Athena, the
 * total number of <outputN> blocks does not need to be specified -- in Athena++ a new
 * output type will be created for each and every <outputN> block in the input file.
 *
 * Required parameters that must be specified in an <outputN> block are:
 *   - variable     = cons,prim,D,d,E,e,m,v
 *   - file_type    = tab,vtk,hst
 *   - dt           = problem time between outputs
 *
 * Optional parameters that may be specified in an <outputN> block are:
 *   - data_format  = format string used in writing data (e.g. %12.5e)
 *   - next_time    = time of next output (useful for restarts)
 *   - id           = any string
 *   - file_number  = any integer with up to 4 digits
 *   - x[123]_slice = specifies data should be a slice at x[123] position
 *   - x[123]_sum   = set to 1 to sum data along specified direction
 *   
 * EXAMPLE of an <outputN> block for a VTK dump:
 *   <output3>
 *   file_type   = tab       # Tabular data dump
 *   variable    = prim      # variables to be output
 *   data_format = %12.5e    # Optional data format string
 *   dt          = 0.01      # time increment between outputs
 *   x2_slice    = 0.0       # slice in x2
 *   x3_slice    = 0.0       # slice in x3
 *
 * Each <outputN> block will result in a new node being created in a linked list of
 * OutputType stored in the Outputs class.  During a simulation, outputs are made when
 * the simulation time satisfies the criteria implemented in the MakeOutputs() function.
 *
 * To implement a new type of output X, write a new derived OutputType class, and
 * construct an object of this class in the Outputs::InitOutputTypes() function.
 *====================================================================================*/

//--------------------------------------------------------------------------------------
// OutputVariable constructor

OutputVariable::OutputVariable(AthenaArray<Real> *parray, OutputVariableHeader vhead)
{
  var_header = vhead;
  pdata = parray;
  pnext = NULL;
  pprev = NULL;
}

// destructor

OutputVariable::~OutputVariable()
{
  if (!pdata->IsShallowCopy()) pdata->DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
// OutputData constructor

OutputData::OutputData()
{
  pfirst_var = NULL;
  plast_var = NULL;
}

// destructor - iterates through linked list of OutputVariables and deletes nodes

OutputData::~OutputData()
{
  OutputVariable *pvar = pfirst_var;
  while(pvar != NULL) {
    OutputVariable *pvar_old = pvar;
    pvar = pvar->pnext;
    delete pvar_old;
  }
}

//--------------------------------------------------------------------------------------
// OutputType constructor

OutputType::OutputType(OutputParameters oparams, MeshBlock *pb)
{
  pmy_block = pb;
  output_params = oparams;
  pnext = NULL; // Terminate linked list with NULL ptr
}

// destructor

OutputType::~OutputType()
{
}

//--------------------------------------------------------------------------------------
// Outputs constructor

Outputs::Outputs(MeshBlock *pb, ParameterInput *pin)
{
  pmy_block = pb;
  pfirst_out_ = NULL;
}

// destructor - iterates through linked list of OutputTypes and deletes nodes

Outputs::~Outputs()
{
  OutputType *pout = pfirst_out_;
  while(pout != NULL) {
    OutputType *pout_old = pout;
    pout = pout->pnext;
    delete pout_old;
  }
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputData::AppendNode()
 *  \brief
 */

void OutputData::AppendNode(AthenaArray<Real> *parray, OutputVariableHeader vhead)
{
  OutputVariable *pnew_var = new OutputVariable(parray, vhead);

  if (pfirst_var == NULL)
    pfirst_var = plast_var = pnew_var;
  else {
    pnew_var->pprev = plast_var;
    plast_var->pnext = pnew_var;
    plast_var = pnew_var;
  }
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputData::ReplaceNode()
 *  \brief
 */

void OutputData::ReplaceNode(OutputVariable *pold, OutputVariable *pnew) 
{
  if (pold == pfirst_var) {
    pfirst_var = pnew;
    if (pold->pnext != NULL) {    // there is another node in the list 
      pnew->pnext = pold->pnext;
      pnew->pnext->pprev = pnew;
    } else {                      // there is only one node in the list
      plast_var = pnew;
    }
  } else if (pold == plast_var) {
    plast_var = pnew;
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
/*! \fn void OutputType::LoadOutputData(OutputData *pod)
 *  \brief initializes output data in OutputData container
 */

void OutputType::LoadOutputData(OutputData *pod)
{
  OutputVariableHeader var_header;
  Fluid *pf = pmy_block->pfluid;;
  std::stringstream str;

// Create OutputData header

  str << "# Athena++ data at time=" << pmy_block->pmy_domain->pmy_mesh->time
      << "  cycle=" << pmy_block->pmy_domain->pmy_mesh->ncycle
      << "  variables=" << output_params.variable << std::endl;
  pod->data_header.descriptor.append(str.str());
  pod->data_header.il = pmy_block->is;
  pod->data_header.iu = pmy_block->ie;
  pod->data_header.jl = pmy_block->js;
  pod->data_header.ju = pmy_block->je;
  pod->data_header.kl = pmy_block->ks;
  pod->data_header.ku = pmy_block->ke;

// Create linked list of OutputVariables containing requested data

  int var_added = 0;
  if (output_params.variable.compare("D") == 0 || 
      output_params.variable.compare("cons") == 0) {
    var_header.type = "SCALARS";
    var_header.name = "dens";
    pod->AppendNode(pf->u.ShallowCopy(IDN,1),var_header); // (lab-frame) density
    var_added = 1;
  }

  if (output_params.variable.compare("d") == 0 || 
      output_params.variable.compare("prim") == 0) {
    var_header.type = "SCALARS";
    var_header.name = "rho";
    pod->AppendNode(pf->w.ShallowCopy(IDN,1),var_header); // (rest-frame) density
    var_added = 1;
  }

  if (output_params.variable.compare("E") == 0 || 
      output_params.variable.compare("cons") == 0) {
    var_header.type = "SCALARS";
    var_header.name = "Etot";
    pod->AppendNode(pf->u.ShallowCopy(IEN,1),var_header); // total energy
    var_added = 1;
  }

  if (output_params.variable.compare("e") == 0 || 
      output_params.variable.compare("prim") == 0) {
    var_header.type = "SCALARS";
    var_header.name = "eint";
    pod->AppendNode(pf->w.ShallowCopy(IEN,1),var_header); // internal energy
    var_added = 1;
  }

  if (output_params.variable.compare("m") == 0 || 
      output_params.variable.compare("cons") == 0) {
    var_header.type = "VECTORS";
    var_header.name = "mom";
    pod->AppendNode(pf->u.ShallowCopy(IM1,3),var_header); // momentum vector
    var_added = 1;
  }

  if (output_params.variable.compare("v") == 0 || 
      output_params.variable.compare("prim") == 0) {
    var_header.type = "VECTORS";
    var_header.name = "vel";
    pod->AppendNode(pf->w.ShallowCopy(IM1,3),var_header); // velocity vector
    var_added = 1;
  }

// throw an error if output variable name not recognized

  if (!var_added) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [OutputType::LoadOutputData]" << std::endl
        << "Output variable '" << output_params.variable << "' not implemented"
        << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return;
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputType::TransformOutputData()
 *  \brief 
 */

void OutputType::TransformOutputData(OutputData *pod)
{
  if (output_params.kslice != -999) {
    Slice(pod,3);
  }
  if (output_params.jslice != -999) {
    Slice(pod,2);
  }
  if (output_params.islice != -999) {
    Slice(pod,1);
  }
  if (output_params.ksum) {
    Sum(pod,3);
  }
  if (output_params.jsum) {
    Sum(pod,2);
  }
  if (output_params.isum) {
    Sum(pod,1);
  }
  return;
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputType::Slice(OutputData* pod, int dim)
 *  \brief
 */

void OutputType::Slice(OutputData* pod, int dim)
{
  OutputVariableHeader var_header;
  AthenaArray<Real> *pslice;
  std::stringstream str;

// For each node in OutputData linked list, slice arrays containing output data  

  OutputVariable *pdn;
  pdn = pod->pfirst_var;

  while (pdn != NULL) {
    var_header = pdn->var_header;
    int nx4 = pdn->pdata->GetDim4();
    int nx3 = pdn->pdata->GetDim3();
    int nx2 = pdn->pdata->GetDim2();
    int nx1 = pdn->pdata->GetDim1();
    pslice = new AthenaArray<Real>;

// Loop over variables and dimensions, extract slice

    if (dim == 3) {
      pslice->NewAthenaArray(nx4,1,nx2,nx1);
      for (int n=0; n<nx4; ++n){
      for (int j=(pod->data_header.jl); j<=(pod->data_header.ju); ++j){
        for (int i=(pod->data_header.il); i<=(pod->data_header.iu); ++i){
          (*pslice)(n,0,j,i) = (*pdn->pdata)(n,output_params.kslice,j,i);
        }
      }}
    } else if (dim == 2) {
      pslice->NewAthenaArray(nx4,nx3,1,nx1);
      for (int n=0; n<nx4; ++n){
      for (int k=(pod->data_header.kl); k<=(pod->data_header.ku); ++k){
        for (int i=(pod->data_header.il); i<=(pod->data_header.iu); ++i){
          (*pslice)(n,k,0,i) = (*pdn->pdata)(n,k,output_params.jslice,i);
        }
      }}
    } else {
      pslice->NewAthenaArray(nx4,nx3,nx2,1);
      for (int n=0; n<nx4; ++n){
      for (int k=(pod->data_header.kl); k<=(pod->data_header.ku); ++k){
        for (int j=(pod->data_header.jl); j<=(pod->data_header.ju); ++j){
          (*pslice)(n,k,j,0) = (*pdn->pdata)(n,k,j,output_params.islice);
        }
      }}
    }

    OutputVariable *pnew = new OutputVariable(pslice, var_header);
    pod->ReplaceNode(pdn,pnew);
    pdn = pdn->pnext;
  }
 
// modify OutputData header

  if (dim == 3) {
    str << "# Slice at x3=" << pmy_block->x3v(output_params.kslice)
        << "  (k-ks)=" << (output_params.kslice - pmy_block->ks) << std::endl;
    pod->data_header.kl = 0;
    pod->data_header.ku = 0;
  } else if (dim == 2) {
    str << "# Slice at x2=" << pmy_block->x2v(output_params.jslice)
        << "  (j-js)=" << (output_params.jslice - pmy_block->js) << std::endl;
    pod->data_header.jl = 0;
    pod->data_header.ju = 0;
  } else {
    str << "# Slice at x1=" << pmy_block->x1v(output_params.islice)
        << "  (i-is)=" << (output_params.islice - pmy_block->is) << std::endl;
    pod->data_header.il = 0;
    pod->data_header.iu = 0;
  }
  pod->data_header.transforms.append(str.str());

  return;
}

//--------------------------------------------------------------------------------------
/*! \fn void OutputType::Sum(OutputData* pod, int dim)
 *  \brief
 */

void OutputType::Sum(OutputData* pod, int dim)
{
  OutputVariableHeader var_header;
  AthenaArray<Real> *psum;
  std::stringstream str;

// For each node in OutputData linked list, sum arrays containing output data  

  OutputVariable *pdn;
  pdn = pod->pfirst_var;

  while (pdn != NULL) {
    var_header = pdn->var_header;
    int nx4 = pdn->pdata->GetDim4();
    int nx3 = pdn->pdata->GetDim3();
    int nx2 = pdn->pdata->GetDim2();
    int nx1 = pdn->pdata->GetDim1();
    psum = new AthenaArray<Real>;

// Loop over variables and dimensions, sum over specified dimension

    if (dim == 3) {
      psum->NewAthenaArray(nx4,1,nx2,nx1);
      for (int n=0; n<nx4; ++n){
      for (int k=(pod->data_header.kl); k<=(pod->data_header.ku); ++k){
      for (int j=(pod->data_header.jl); j<=(pod->data_header.ju); ++j){
        for (int i=(pod->data_header.il); i<=(pod->data_header.iu); ++i){
          (*psum)(n,0,j,i) += (*pdn->pdata)(n,k,j,i);
        }
      }}}
    } else if (dim == 2) {
      psum->NewAthenaArray(nx4,nx3,1,nx1);
      for (int n=0; n<nx4; ++n){
      for (int k=(pod->data_header.kl); k<=(pod->data_header.ku); ++k){
      for (int j=(pod->data_header.jl); j<=(pod->data_header.ju); ++j){
        for (int i=(pod->data_header.il); i<=(pod->data_header.iu); ++i){
          (*psum)(n,k,0,i) += (*pdn->pdata)(n,k,j,i);
        }
      }}}
    } else {
      psum->NewAthenaArray(nx4,nx3,nx2,1);
      for (int n=0; n<nx4; ++n){
      for (int k=(pod->data_header.kl); k<=(pod->data_header.ku); ++k){
      for (int j=(pod->data_header.jl); j<=(pod->data_header.ju); ++j){
        for (int i=(pod->data_header.il); i<=(pod->data_header.iu); ++i){
          (*psum)(n,k,j,0) += (*pdn->pdata)(n,k,j,i);
        }
      }}}
    }

    OutputVariable *pnew = new OutputVariable(psum, var_header);
    pod->ReplaceNode(pdn,pnew);
    pdn = pdn->pnext;
  }
 
// modify OutputData header

  if (dim == 3) {
    str << "# Sum over x3" << std::endl;
    pod->data_header.kl = 0;
    pod->data_header.ku = 0;
  } else if (dim == 2) {
    str << "# Sum over x2" << std::endl;
    pod->data_header.jl = 0;
    pod->data_header.ju = 0;
  } else {
    str << "# Sum over x1" << std::endl;
    pod->data_header.il = 0;
    pod->data_header.iu = 0;
  }
  pod->data_header.transforms.append(str.str());

  return;
}

//--------------------------------------------------------------------------------------
/*! \fn void Outputs::InitOutputTypes()
 *  \brief Creates linked list of OutputTypes based on <ouput> blocks in input file.
 */

void Outputs::InitOutputTypes(ParameterInput *pin)
{
  std::stringstream msg;
  InputBlock *pib = pin->pfirst_block;
  OutputType *pnew_out;
  OutputType *plast = pfirst_out_;
  MeshBlock *pb = pmy_block;

// loop over input block names.  Find those that start with "output", read parameters,
// and construct linked list of OutputTypes.

  while (pib != NULL) {
    if (pib->block_name.compare(0,6,"output") == 0) {
      int create_output = 1;
      OutputParameters op;  // define temporary OutputParameters struct

// extract integer number of output block.  Save name and number 

      std::string outn = pib->block_name.substr(6); // 6 because starts at 0!
      op.block_number = atoi(outn.c_str());
      op.block_name.assign(pib->block_name);

// set time of last output, time between outputs

      op.next_time = pin->GetOrAddReal(op.block_name,"next_time",0.0);
      op.dt = pin->GetReal(op.block_name,"dt");

// set file number, basename, id, and format

      op.file_number = pin->GetOrAddInteger(op.block_name,"file_number",0);
      op.file_basename = pin->GetString("job","problem_id");
      char define_id[10];
      sprintf(define_id,"out%d",op.block_number);  // default id="outN"
      op.file_id = pin->GetOrAddString(op.block_name,"id",define_id);
      op.file_type = pin->GetString(op.block_name,"file_type");

// read slicing options.  Check that slice is within range

      if (pin->DoesParameterExist(op.block_name,"x1_slice")) {
        Real x1 = pin->GetReal(op.block_name,"x1_slice");
        if (x1 >= pb->block_size.x1min && x1 < pb->block_size.x1max) {
          for (int i=pb->is+1; i<=pb->ie+1; ++i) {
            if (pb->x1f(i) > x1) {
              op.islice = i-1;
              break;
            }
          }
        } else {
          create_output=0;;
        }
      } else {op.islice = -999;}

      if (pin->DoesParameterExist(op.block_name,"x2_slice")) {
        Real x2 = pin->GetReal(op.block_name,"x2_slice");
        if (x2 >= pb->block_size.x2min && x2 < pb->block_size.x2max) {
          for (int j=pb->js+1; j<=pb->je+1; ++j) {
            if (pb->x2f(j) > x2) {
              op.jslice = j-1;
              break;
            }
          }
        } else {
          create_output=0;;
        }
      } else {op.jslice = -999;}

      if (pin->DoesParameterExist(op.block_name,"x3_slice")) {
        Real x3 = pin->GetReal(op.block_name,"x3_slice");
        if (x3 >= pb->block_size.x3min && x3 < pb->block_size.x3max) {
          for (int k=pb->ks+1; k<=pb->ke+1; ++k) {
            if (pb->x3f(k) > x3) {
              op.kslice = k-1;
              break;
            }
          }
        } else {
          create_output=0;;
        }
      } else {op.kslice = -999;}

// read sum options.  Check for conflicts with slicing.

      if (pin->DoesParameterExist(op.block_name,"x1_sum")) {
        if (pin->DoesParameterExist(op.block_name,"x1_slice")) {
          msg << "### FATAL ERROR in function [Outputs::InitOutputTypes]"
              << std::endl << "Cannot request both slice and sum along x1-direction"
              << " in output block '" << op.block_name << "'" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        } else {
          op.isum = pin->GetInteger(op.block_name,"x1_sum");;
        }
      } else {op.isum = 0;}

      if (pin->DoesParameterExist(op.block_name,"x2_sum")) {
        if (pin->DoesParameterExist(op.block_name,"x2_slice")) {
          msg << "### FATAL ERROR in function [Outputs::InitOutputTypes]"
              << std::endl << "Cannot request both slice and sum along x2-direction"
              << " in output block '" << op.block_name << "'" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        } else {
          op.jsum = pin->GetInteger(op.block_name,"x2_sum");;
        }
      } else {op.jsum = 0;}

      if (pin->DoesParameterExist(op.block_name,"x3_sum")) {
        if (pin->DoesParameterExist(op.block_name,"x3_slice")) {
          msg << "### FATAL ERROR in function [Outputs::InitOutputTypes]"
              << std::endl << "Cannot request both slice and sum along x3-direction"
              << " in output block '" << op.block_name << "'" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        } else {
          op.ksum = pin->GetInteger(op.block_name,"x3_sum");;
        }
      } else {op.ksum = 0;}

      if (create_output) {  // skip output if slice not in range
// set output variable and optional data format string used in formatted writes

        if (op.file_type.compare("hst") != 0) {
          op.variable = pin->GetString(op.block_name,"variable");
        }
        op.data_format = pin->GetOrAddString(op.block_name,"data_format","%12e.5");
        op.data_format.insert(0," "); // prepend with blank to separate columns

// Construct new OutputType according to file format
// TODO: add any new output types here

        if (op.file_type.compare("tab") == 0) {
          pnew_out = new FormattedTableOutput(op, pmy_block);
        } else if (op.file_type.compare("hst") == 0) {
          pnew_out = new HistoryOutput(op, pmy_block);
        } else if (op.file_type.compare("vtk") == 0) {
          pnew_out = new VTKOutput(op, pmy_block);
        } else {
          msg << "### FATAL ERROR in function [Outputs::InitOutputTypes]"
              << std::endl << "Unrecognized file format = '" << op.file_type 
              << "' in output block '" << op.block_name << "'" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        }

// Add type as node in linked list 

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
/*! \fn void Outputs::MakeOutputs()
 *  \brief scans through linked list of OutputTypes and makes any outputs needed.
 */

void Outputs::MakeOutputs()
{
  OutputType* pout = pfirst_out_;
  Mesh* pm = pmy_block->pmy_domain->pmy_mesh;

  while (pout != NULL) {
    if ((pm->time == pm->start_time) ||
        (pm->time >= pout->output_params.next_time) ||
        (pm->time >= pm->tlim)) {

// Create new OutputData container, load and transform data, then write to file

      OutputData* pod = new OutputData;
      pout->LoadOutputData(pod);
      pout->TransformOutputData(pod);
      pout->WriteOutputFile(pod);
      delete pod;

    }
    pout = pout->pnext; // move to next OutputType in list
  }
}
