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
//! \file outputs.cpp
//  \brief implements functions for Athena++ outputs
//
// The number and types of outputs are all controlled by the number and values of
// parameters specified in <outputN> blocks in the input file.  Each output block must
// be labelled by a unique integer "N".  Following the convention of the parser
// implemented in the ParameterInput class, a second output block with the same integer
// "N" of an earlier block will silently overwrite the values read by the first block.
// The numbering of the output blocks does not need to be consecutive, and blocks may
// appear in any order in the input file.  Moreover, unlike the C version of Athena, the
// total number of <outputN> blocks does not need to be specified -- in Athena++ a new
// output type will be created for each and every <outputN> block in the input file.
//
// Required parameters that must be specified in an <outputN> block are:
//   - variable     = cons,prim,D,d,E,e,m,v
//   - file_type    = rst,tab,vtk,hst
//   - dt           = problem time between outputs
//
// Optional parameters that may be specified in an <outputN> block are:
//   - data_format  = format string used in writing data (e.g. %12.5e)
//   - next_time    = time of next output (useful for restarts)
//   - id           = any string
//   - file_number  = any integer with up to 4 digits
//   - x[123]_slice = specifies data should be a slice at x[123] position
//   - x[123]_sum   = set to "true" to sum data along specified direction
//   
// EXAMPLE of an <outputN> block for a VTK dump:
//   <output3>
//   file_type   = tab       # Tabular data dump
//   variable    = prim      # variables to be output
//   data_format = %12.5e    # Optional data format string
//   dt          = 0.01      # time increment between outputs
//   x2_slice    = 0.0       # slice in x2
//   x3_slice    = 0.0       # slice in x3
//
// Each <outputN> block will result in a new node being created in a linked list of
// OutputType stored in the Outputs class.  During a simulation, outputs are made when
// the simulation time satisfies the criteria implemented in the MakeOutputs() function.
//
// To implement a new output type, write a new derived OutputType class, and construct
// an object of this class in the Outputs constructor at the location indicated by the
// text 'ADD NEW OUTPUT TYPES HERE'.
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
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp" // Coordinates

// this class header
#include "outputs.hpp"

//--------------------------------------------------------------------------------------
// OutputType constructor

OutputType::OutputType(OutputParameters oparams)
{
  output_params = oparams;
  pnext_type = NULL;   // Terminate this node in linked list with NULL ptr

  num_vars_ = 0;
  pfirst_data_ = NULL; // Initialize start of linked list of OutputData's to NULL
  plast_data_ = NULL;  // Initialize end   of linked list of OutputData's to NULL
}

// destructor

OutputType::~OutputType()
{
}

//--------------------------------------------------------------------------------------
// Outputs constructor

Outputs::Outputs(Mesh *pm, ParameterInput *pin)
{
  pfirst_type_ = NULL;
  std::stringstream msg;
  InputBlock *pib = pin->pfirst_block;
  OutputType *pnew_type;
  OutputType *plast = pfirst_type_;

  // loop over input block names.  Find those that start with "output", read parameters,
  // and construct linked list of OutputTypes.
  while (pib != NULL) {
    if (pib->block_name.compare(0,6,"output") == 0) {
      OutputParameters op;  // define temporary OutputParameters struct

      // extract integer number of output block.  Save name and number 
      std::string outn = pib->block_name.substr(6); // 6 because counting starts at 0!
      op.block_number = atoi(outn.c_str());
      op.block_name.assign(pib->block_name);

      // set time of last output, time between outputs
      op.next_time = pin->GetOrAddReal(op.block_name,"next_time",0.0);
      op.dt = pin->GetReal(op.block_name,"dt");

      if (op.dt > 0.0) {  // only add output if dt>0

        // set file number, basename, id, and format
        op.file_number = pin->GetOrAddInteger(op.block_name,"file_number",0);
        op.file_basename = pin->GetString("job","problem_id");
        char define_id[10];
        sprintf(define_id,"out%d",op.block_number);  // default id="outN"
        op.file_id = pin->GetOrAddString(op.block_name,"id",define_id);
        op.file_type = pin->GetString(op.block_name,"file_type");

        // read slicing options.  Check that slice is within mesh
        if (pin->DoesParameterExist(op.block_name,"x1_slice")) {
          Real x1 = pin->GetReal(op.block_name,"x1_slice");
          if (x1 >= pm->mesh_size.x1min && x1 < pm->mesh_size.x1max) {
            op.x1_slice = x1;
            op.output_slicex1 = true;
          } else {
            msg << "### FATAL ERROR in Outputs constructor" << std::endl
                << "Slice at x1=" << x1 << " in output block '" << op.block_name
                << "' is out of range of Mesh" << std::endl;
            throw std::runtime_error(msg.str().c_str());
          }
        }

        if (pin->DoesParameterExist(op.block_name,"x2_slice")) {
          Real x2 = pin->GetReal(op.block_name,"x2_slice");
          if (x2 >= pm->mesh_size.x2min && x2 < pm->mesh_size.x2max) {
            op.x2_slice = x2;
            op.output_slicex2 = true;
          } else {
            msg << "### FATAL ERROR in Outputs constructor" << std::endl
                << "Slice at x2=" << x2 << " in output block '" << op.block_name
                << "' is out of range of Mesh" << std::endl;
            throw std::runtime_error(msg.str().c_str());
          }
        }

        if (pin->DoesParameterExist(op.block_name,"x3_slice")) {
          Real x3 = pin->GetReal(op.block_name,"x3_slice");
          if (x3 >= pm->mesh_size.x3min && x3 < pm->mesh_size.x3max) {
            op.x3_slice = x3;
            op.output_slicex3 = true;
          } else {
            msg << "### FATAL ERROR in Outputs constructor" << std::endl
                << "Slice at x3=" << x3 << " in output block '" << op.block_name
                << "' is out of range of Mesh" << std::endl;
            throw std::runtime_error(msg.str().c_str());
          }
        }

        // read sum options.  Check for conflicts with slicing.
        op.output_sumx1 = pin->GetOrAddBoolean(op.block_name,"x1_sum",false);
        if ((op.output_slicex1) && (op.output_sumx1)) {
          msg << "### FATAL ERROR in Outputs constructor" << std::endl 
              << "Cannot request both slice and sum along x1-direction"
              << " in output block '" << op.block_name << "'" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        }
        op.output_sumx2 = pin->GetOrAddBoolean(op.block_name,"x2_sum",false);
        if ((op.output_slicex2) && (op.output_sumx2)) {
          msg << "### FATAL ERROR in Outputs constructor" << std::endl
              << "Cannot request both slice and sum along x2-direction"
              << " in output block '" << op.block_name << "'" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        }
        op.output_sumx3 = pin->GetOrAddBoolean(op.block_name,"x3_sum",false);
        if ((op.output_slicex3) && (op.output_sumx3)) {
          msg << "### FATAL ERROR in Outputs constructor" << std::endl
              << "Cannot request both slice and sum along x3-direction"
              << " in output block '" << op.block_name << "'" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        }

        // read ghost cell option
        op.include_ghost_zones=pin->GetOrAddBoolean(op.block_name,"ghost_zones",false);

        // set output variable and optional data format string used in formatted writes
        if (op.file_type.compare("hst") != 0 && op.file_type.compare("rst") != 0) {
          op.variable = pin->GetString(op.block_name,"variable");
        }
        op.data_format = pin->GetOrAddString(op.block_name,"data_format","%12.5e");
        op.data_format.insert(0," "); // prepend with blank to separate columns

        // Construct new OutputType according to file format
        // ADD NEW OUTPUT TYPES HERE
        if (op.file_type.compare("hst") == 0) {
          pnew_type = new HistoryOutput(op);
        } else if (op.file_type.compare("tab") == 0) {
          pnew_type = new FormattedTableOutput(op);
        } else if (op.file_type.compare("vtk") == 0) {
          pnew_type = new VTKOutput(op);
        }
        
//        if (op.file_type.compare("rst") == 0) {
//          pnew_type = new RestartOutput(op);
//        }
#ifdef HDF5OUTPUT
//        else if (op.file_type.compare("ath5") == 0 || op.file_type.compare("hdf5") == 0) {
//          pnew_type = new ATHDF5Output(op);
//        }
#endif
        else {
          msg << "### FATAL ERROR in Outputs constructor" << std::endl
              << "Unrecognized file format = '" << op.file_type 
              << "' in output block '" << op.block_name << "'" << std::endl;
          throw std::runtime_error(msg.str().c_str());
        }

        // Add type as node in linked list 
        if (pfirst_type_ == NULL) {
          pfirst_type_ = pnew_type;
        } else {
          plast->pnext_type = pnew_type;
        }
        plast = pnew_type;
      }
    }
    pib = pib->pnext;  // move to next input block name
  }

  // Move restarts to the end of the OutputType list, so file counters for other 
  // output types are up-to-date in restart file
  int pos=0, found=0;
  OutputType *pot=pfirst_type_, *prst;
  while(pot!=NULL) {
    if(pot->output_params.file_type.compare("rst")==0) {
      prst=pot;
      found=1;
      if(pot->pnext_type==NULL) found=2;
      break;
    }
    pos++;
    pot=pot->pnext_type;
  }
  if(found==1) {
    // remove the restarting block
    pot=pfirst_type_;
    if(pos==0) { // first block
      pfirst_type_=pfirst_type_->pnext_type;
    }
    else {
      for(int j=0; j<pos-1; j++) // seek the list
        pot=pot->pnext_type;
      pot->pnext_type=prst->pnext_type; // remove it
    }
    while(pot->pnext_type!=NULL)
      pot=pot->pnext_type; // find the end
    prst->pnext_type=NULL;
    pot->pnext_type=prst;
  }
  // if found==2, do nothing; it's already at the end of the list
}

// destructor - iterates through linked list of OutputTypes and deletes nodes

Outputs::~Outputs()
{
  OutputType *ptype = pfirst_type_;
  while(ptype != NULL) {
    OutputType *ptype_old = ptype;
    ptype = ptype->pnext_type;
    delete ptype_old;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void OutputType::LoadOutputData(MeshBlock *pmb)
//  \brief Create linked list of OutputData's containing requested variables

void OutputType::LoadOutputData(MeshBlock *pmb)
{
  Hydro *phyd = pmb->phydro;
  Field *pfld = pmb->pfield;
  OutputData *pod;
  num_vars_ = 0;

  // (lab-frame) density
  if (output_params.variable.compare("D") == 0 || 
      output_params.variable.compare("cons") == 0) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "dens";
    pod->data.InitWithShallowSlice(phyd->u,4,IDN,1);
    AppendOutputDataNode(pod);
    num_vars_++;
  }

  // (rest-frame) density
  if (output_params.variable.compare("d") == 0 || 
      output_params.variable.compare("prim") == 0) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "rho";
    pod->data.InitWithShallowSlice(phyd->w,4,IDN,1);
    AppendOutputDataNode(pod);
    num_vars_++;
  }

  // total energy
  if (NON_BAROTROPIC_EOS) {
    if (output_params.variable.compare("E") == 0 || 
        output_params.variable.compare("cons") == 0) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "Etot";
      pod->data.InitWithShallowSlice(phyd->u,4,IEN,1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
  }

  // pressure
  if (NON_BAROTROPIC_EOS) {
    if (output_params.variable.compare("p") == 0 || 
        output_params.variable.compare("prim") == 0) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "press";
      pod->data.InitWithShallowSlice(phyd->w,4,IPR,1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
  }

  // momentum vector
  if (output_params.variable.compare("m") == 0 || 
      output_params.variable.compare("cons") == 0) {
    pod = new OutputData;
    pod->type = "VECTORS";
    pod->name = "mom";
    pod->data.InitWithShallowSlice(phyd->u,4,IM1,3);
    AppendOutputDataNode(pod);
    num_vars_+=3;
  }

  // each component of momentum
  if (output_params.variable.compare("m1") == 0) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "mom1";
    pod->data.InitWithShallowSlice(phyd->u,4,IM1,1);
    AppendOutputDataNode(pod);
    num_vars_++;
  }
  if (output_params.variable.compare("m2") == 0) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "mom2";
    pod->data.InitWithShallowSlice(phyd->u,4,IM2,1);
    AppendOutputDataNode(pod);
    num_vars_++;
  }
  if (output_params.variable.compare("m3") == 0) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "mom3";
    pod->data.InitWithShallowSlice(phyd->u,4,IM3,1);
    AppendOutputDataNode(pod);
    num_vars_++;
  }

  // velocity vector
  if (output_params.variable.compare("v") == 0 || 
      output_params.variable.compare("prim") == 0) {
    pod = new OutputData;
    pod->type = "VECTORS";
    pod->name = "vel";
    pod->data.InitWithShallowSlice(phyd->w,4,IVX,3);
    AppendOutputDataNode(pod);
    num_vars_+=3;
  }

  // each component of velocity
  if (output_params.variable.compare("vx") == 0 || 
      output_params.variable.compare("v1") == 0) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "vel1";
    pod->data.InitWithShallowSlice(phyd->w,4,IVX,1);
    AppendOutputDataNode(pod);
    num_vars_++;
  }
  if (output_params.variable.compare("vy") == 0 || 
      output_params.variable.compare("v2") == 0) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "vel2";
    pod->data.InitWithShallowSlice(phyd->w,4,IVY,1);
    AppendOutputDataNode(pod);
    num_vars_++;
  }
  if (output_params.variable.compare("vz") == 0 || 
      output_params.variable.compare("v3") == 0) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "vel3";
    pod->data.InitWithShallowSlice(phyd->w,4,IVZ,1);
    AppendOutputDataNode(pod);
    num_vars_++;
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    // vector of cell-centered magnetic field
    if (output_params.variable.compare("bcc") == 0 || 
        output_params.variable.compare("prim") == 0 ||
        output_params.variable.compare("cons") == 0) {
      pod = new OutputData;
      pod->type = "VECTORS";
      pod->name = "Bcc";
      pod->data.InitWithShallowSlice(pfld->bcc,4,IB1,3);
      AppendOutputDataNode(pod);
      num_vars_+=3;
    }

    // each component of cell-centered magnetic field
    if (output_params.variable.compare("bcc1") == 0) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "Bcc1";
      pod->data.InitWithShallowSlice(pfld->bcc,4,IB1,1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
    if (output_params.variable.compare("bcc2") == 0) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "Bcc2";
      pod->data.InitWithShallowSlice(pfld->bcc,4,IB2,1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
    if (output_params.variable.compare("bcc3") == 0) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "Bcc2";
      pod->data.InitWithShallowSlice(pfld->bcc,4,IB3,1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }

    // each component of face-centered magnetic field
    if (output_params.variable.compare("b1") == 0) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "B1";
      pod->data.InitWithShallowSlice(pfld->b.x1f,4,0,1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
    if (output_params.variable.compare("b2") == 0) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "B2";
      pod->data.InitWithShallowSlice(pfld->b.x2f,4,0,1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
    if (output_params.variable.compare("b3") == 0) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "B3";
      pod->data.InitWithShallowSlice(pfld->b.x3f,4,0,1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
  } // endif (MAGNETIC_FIELDS_ENABLED)

  // throw an error if output variable name not recognized
  if (num_vars_==0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [OutputType::LoadOutputData]" << std::endl
        << "Output variable '" << output_params.variable << "' not implemented"
        << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutputData::AppendOutputDataNode(OutputData *pod)
//  \brief

void OutputType::AppendOutputDataNode(OutputData *pnew_data)
{
  if (pfirst_data_ == NULL) {
    pfirst_data_ = pnew_data;
  } else {
    pnew_data->pprev = plast_data_;
    plast_data_->pnext = pnew_data;
  }
  plast_data_ = pnew_data;
}

//--------------------------------------------------------------------------------------
//! \fn void OutputData::ReplaceOutputDataNode()
//  \brief

void OutputType::ReplaceOutputDataNode(OutputData *pold, OutputData *pnew)
{
  if (pold == pfirst_data_) {
    pfirst_data_ = pnew;
    if (pold->pnext != NULL) {    // there is another node in the list 
      pnew->pnext = pold->pnext;
      pnew->pnext->pprev = pnew;
    } else {                      // there is only one node in the list
      plast_data_ = pnew;
    }
  } else if (pold == plast_data_) {
    plast_data_ = pnew;
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
//! \fn void OutputData::ClearOutputData()
//  \brief

void OutputType::ClearOutputData()
{
  OutputData *pdata = pfirst_data_;
  while (pdata != NULL) {
    OutputData *pdata_old = pdata;
    pdata = pdata->pnext;
    delete pdata_old;
  }
  pfirst_data_ = NULL;
  plast_data_  = NULL;
}

//--------------------------------------------------------------------------------------
//! \fn void Outputs::MakeOutputs(Mesh *pm, ParameterInput *pin, bool wtflag)
//  \brief scans through linked list of OutputTypes and makes any outputs needed.

void Outputs::MakeOutputs(Mesh *pm, ParameterInput *pin, bool wtflag)
{
  OutputType* ptype = pfirst_type_;
  MeshBlock *pmb;

  while (ptype != NULL) {
    if ((pm->time == pm->start_time) ||
        (pm->time >= ptype->output_params.next_time) ||
        (pm->time >= pm->tlim) ||
        (wtflag==true && ptype->output_params.file_type=="rst")) {

      ptype->WriteOutputFile(pm);
      ptype->FinalizeOutput(pin); // increment time, counter, etc.
    }
    ptype = ptype->pnext_type; // move to next OutputType in list
  }

}

//--------------------------------------------------------------------------------------
//! \fn void OutputType::FinalizeOutput(ParameterInput *pin)
//  \brief increment file number, update next output time, store in ParameterInput

void OutputType::FinalizeOutput(ParameterInput *pin)
{
  output_params.file_number++;
  output_params.next_time += output_params.dt;
  pin->SetInteger(output_params.block_name, "file_number", output_params.file_number);
  pin->SetReal(output_params.block_name, "next_time", output_params.next_time);
}

//--------------------------------------------------------------------------------------
//! \fn void OutputType::TransformOutputData(MeshBlock *pmb)
//  \brief Calls sum and slice functions on each direction in turn, in order to allow
//  mulitple operations performed on the same data set

bool OutputType::TransformOutputData(MeshBlock *pmb)
{
  bool flag = true;
  if (output_params.output_slicex3) {
    bool ret = Slice(pmb,3);
    if (ret==false) flag=false;
  }
  if (output_params.output_slicex2) {
    bool ret = Slice(pmb,2);
    if (ret==false) flag=false;
  }
  if (output_params.output_slicex1) {
    bool ret = Slice(pmb,1);
    if (ret==false) flag=false;
  }
  if (output_params.output_sumx1) {
    Sum(pmb,3);
  }
  if (output_params.output_sumx2) {
    Sum(pmb,2);
  }
  if (output_params.output_sumx3) {
    Sum(pmb,1);
  }
  return flag;
}

//--------------------------------------------------------------------------------------
//! \fn bool OutputType::Slice(MeshBlock *pmb, int dim)
//  \brief

bool OutputType::Slice(MeshBlock *pmb, int dim)
{
  int islice, jslice, kslice;

  // Compute i,j,k indices of slice; check if in range of data in this block
  if (dim == 1) {
    if (output_params.x1_slice >= pmb->block_size.x1min && 
        output_params.x1_slice < pmb->block_size.x1max) {
      for (int i=pmb->is+1; i<=pmb->ie+1; ++i) {
        if (pmb->pcoord->x1f(i) > output_params.x1_slice) {
          islice = i-1;
          output_params.islice = islice;
          break;
        }
      }
    } else {
      return false;
    }
  } else if (dim == 2) {
    if (output_params.x2_slice >= pmb->block_size.x2min &&
        output_params.x2_slice < pmb->block_size.x2max) {
      for (int j=pmb->js+1; j<=pmb->je+1; ++j) {
        if (pmb->pcoord->x2f(j) > output_params.x2_slice) {
          jslice = j-1;
          output_params.jslice = jslice;
          break;
        }
      }
    } else {
      return false;
    }
  } else {
    if (output_params.x3_slice >= pmb->block_size.x3min &&
        output_params.x3_slice < pmb->block_size.x3max) {
      for (int k=pmb->ks+1; k<=pmb->ke+1; ++k) {
        if (pmb->pcoord->x3f(k) > output_params.x3_slice) {
          kslice = k-1;
          output_params.kslice = kslice;
          break;
        }
      }
    } else {
      return false;
    }
  }

  // For each node in OutputData linked list, slice arrays containing output data  
  OutputData *pdata,*pnew;
  pdata = pfirst_data_;

  while (pdata != NULL) {
    pnew = new OutputData;
    pnew->type = pdata->type;
    pnew->name = pdata->name;
    int nx4 = pdata->data.GetDim4();
    int nx3 = pdata->data.GetDim3();
    int nx2 = pdata->data.GetDim2();
    int nx1 = pdata->data.GetDim1();

    // Loop over variables and dimensions, extract slice
    if (dim == 3) {
      pnew->data.NewAthenaArray(nx4,1,nx2,nx1);
      for (int n=0; n<nx4; ++n){
      for (int j=ojl; j<=oju; ++j){
        for (int i=oil; i<=oiu; ++i){
          pnew->data(n,0,j,i) = pdata->data(n,kslice,j,i);
        }
      }}
    } else if (dim == 2) {
      pnew->data.NewAthenaArray(nx4,nx3,1,nx1);
      for (int n=0; n<nx4; ++n){
      for (int k=okl; k<=oku; ++k){
        for (int i=oil; i<=oiu; ++i){
          pnew->data(n,k,0,i) = pdata->data(n,k,jslice,i);
        }
      }}
    } else {
      pnew->data.NewAthenaArray(nx4,nx3,nx2,1);
      for (int n=0; n<nx4; ++n){
      for (int k=okl; k<=oku; ++k){
        for (int j=ojl; j<=oju; ++j){
          pnew->data(n,k,j,0) = pdata->data(n,k,j,islice);
        }
      }}
    }

    ReplaceOutputDataNode(pdata,pnew);
    pdata = pdata->pnext;
  }
 
  // modify array indices
  if (dim == 3) {
    okl = 0;
    oku = 0;
  } else if (dim == 2) {
    ojl = 0;
    oju = 0;
  } else {
    oil = 0;
    oiu = 0;
  }

  return true;
}

//--------------------------------------------------------------------------------------
//! \fn void OutputType::Sum(OutputData* pod, int dim)
//  \brief

void OutputType::Sum(MeshBlock* pmb, int dim)
{
  AthenaArray<Real> *psum;
  std::stringstream str;

  // For each node in OutputData linked list, sum arrays containing output data  
  OutputData *pdata,*pnew;
  pdata = pfirst_data_;

  while (pdata != NULL) {
    pnew = new OutputData;
    pnew->type = pdata->type;
    pnew->name = pdata->name;
    int nx4 = pdata->data.GetDim4();
    int nx3 = pdata->data.GetDim3();
    int nx2 = pdata->data.GetDim2();
    int nx1 = pdata->data.GetDim1();
    psum = new AthenaArray<Real>;

    // Loop over variables and dimensions, sum over specified dimension
    if (dim == 3) {
      pnew->data.NewAthenaArray(nx4,1,nx2,nx1);
      for (int n=0; n<nx4; ++n){
      for (int k=okl; k<=oku; ++k){
      for (int j=ojl; j<=oju; ++j){
        for (int i=oil; i<=oiu; ++i){
          pnew->data(n,0,j,i) += pdata->data(n,k,j,i);
        }
      }}}
    } else if (dim == 2) {
      pnew->data.NewAthenaArray(nx4,nx3,1,nx1);
      for (int n=0; n<nx4; ++n){
      for (int k=okl; k<=oku; ++k){
      for (int j=ojl; j<=oju; ++j){
        for (int i=oil; i<=oiu; ++i){
          pnew->data(n,k,0,i) += pdata->data(n,k,j,i);
        }
      }}}
    } else {
      pnew->data.NewAthenaArray(nx4,nx3,nx2,1);
      for (int n=0; n<nx4; ++n){
      for (int k=okl; k<=oku; ++k){
      for (int j=ojl; j<=oju; ++j){
        for (int i=oil; i<=oiu; ++i){
          pnew->data(n,k,j,0) += pdata->data(n,k,j,i);
        }
      }}}
    }

    ReplaceOutputDataNode(pdata,pnew);
    pdata = pdata->pnext;
  }
 
  // modify array indices
  if (dim == 3) {
    okl = 0;
    oku = 0;
  } else if (dim == 2) {
    ojl = 0;
    oju = 0;
  } else {
    oil = 0;
    oiu = 0;
  }

  return;
}
