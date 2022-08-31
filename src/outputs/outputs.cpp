//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file outputs.cpp
//! \brief implements functions for Athena++ outputs
//!
//! The number and types of outputs are all controlled by the number and values of
//! parameters specified in <output[n]> blocks in the input file.  Each output block must
//! be labelled by a unique integer "n".  Following the convention of the parser
//! implemented in the ParameterInput class, a second output block with the same integer
//! "n" of an earlier block will silently overwrite the values read by the first block.
//! The numbering of the output blocks does not need to be consecutive, and blocks may
//! appear in any order in the input file.  Moreover, unlike the C version of Athena,
//! the total number of <output[n]> blocks does not need to be specified --
//! in Athena++ a new output type will be created for each and every <output[n]> block
//! in the input file.
//!
//! Required parameters that must be specified in an <output[n]> block are:
//!   - variable     = cons,prim,D,d,E,e,m,m1,m2,m3,v,v1=vx,v2=vy,v3=vz,p,
//!                    bcc,bcc1,bcc2,bcc3,b,b1,b2,b3,phi,uov
//!   - file_type    = rst,tab,vtk,hst,hdf5
//!   - dt           = problem time between outputs
//!
//! EXAMPLE of an <output[n]> block for a VTK dump:
//!
//!     <output3>
//!     file_type   = tab       # Tabular data dump
//!     variable    = prim      # variables to be output
//!     data_format = %12.5e    # Optional data format string
//!     dt          = 0.01      # time increment between outputs
//!     x2_slice    = 0.0       # slice in x2
//!     x3_slice    = 0.0       # slice in x3
//!
//!
//! Each <output[n]> block will result in a new node being created in a linked list of
//! OutputType stored in the Outputs class.  During a simulation, outputs are made when
//! the simulation time satisfies the criteria implemented in the MakeOutputs() function.
//!
//! \note
//! To implement a new output type, write a new OutputType derived class, and construct
//! an object of this class in the Outputs constructor at the location indicated by the
//! comment text: 'NEW_OUTPUT_TYPES'. Current summary:
//! - outputs.cpp, OutputType:LoadOutputData() (below): conditionally add new OutputData
//!   node to linked list, depending on the user-input 'variable' string.
//!   Provide direction on how to slice a possible 4D source AthenaArray into separate
//!   3D arrays; automatically enrolls quantity in vtk.cpp, formatted_table.cpp outputs.
//! - athena_hdf5.cpp, ATHDF5Output::WriteOutputFile(): need to allocate space for the new
//!   OutputData node as an HDF5 "variable" inside an existing HDF5 "dataset"
//!   (cell-centered vs. face-centered data).
//! - restart.cpp, RestartOutput::WriteOutputFile(): memcpy array of quantity to pdata
//!   pointer and increment the pointer. pdata points to an allocated region of
//!   memory whose "datasize" is inferred from MeshBlock::GetBlockSizeInBytes(), ---->
//! - mesh/meshblock.cpp, MeshBlock::GetBlockSizeInBytes(): increment std::size_t size by
//!   the size of the new quantity's array(s)
//! - mesh/meshblock.cpp, MeshBlock restart constructor: memcpy quantity
//!   (IN THE SAME ORDER AS THE VARIABLES ARE WRITTEN IN restart.cpp)
//!   from the loaded .rst file to the MeshBlock's appropriate physics member object
//!
//! - history.cpp, HistoryOutput::WriteOutputFile() (3x places):
//!   1) modify NHISTORY_VARS macro
//!      so that the size of data_sum[] can accommodate the new physics, when active.
//!   2) Compute volume-weighted data_sum[i] for the new quantity + etc. factors
//!   3) Provide short string to serve as the column header description of new quantity
//!
//! HDF5 note: packing gas velocity into the "prim" HDF5 dataset will cause VisIt to treat
//! the 3x components as independent scalars instead of a physical vector, unlike how it
//! treats .vtk velocity output from Athena++. The workaround is to import the
//! vis/visit/*.xml expressions file, which can pack these HDF5 scalars into a vector.
//
// TODO(felker): Replace MeshBlock::GetBlockSizeInBytes() by 2x RegisterMeshBlockData()
// overloads. Replace "NEW_OUTPUT_TYPES" region of RestartOutput::WriteOutputFile() with
// automatic loops over registered MeshBlock quantities in pvars_cc, pvars_fc vectors.
//========================================================================================

// C headers

// C++ headers
#include <cstdio>
#include <cstdlib>
#include <cstring>    // strcmp
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>   // std::string, to_string()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../gravity/gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"
#include "outputs.hpp"

//----------------------------------------------------------------------------------------
//! OutputType constructor

OutputType::OutputType(OutputParameters oparams) :
    output_params(oparams),
    pnext_type(),  // Terminate this node in singly linked list with nullptr
    num_vars_(),
    // nested doubly linked list of OutputData:
    pfirst_data_(),  // Initialize head node to nullptr
    plast_data_() { // Initialize tail node to nullptr
}

//----------------------------------------------------------------------------------------
//! Outputs constructor

Outputs::Outputs(Mesh *pm, ParameterInput *pin) {
  pfirst_type_ = nullptr;
  std::stringstream msg;
  InputBlock *pib = pin->pfirst_block;
  OutputType *pnew_type;
  OutputType *plast = pfirst_type_;
  int num_hst_outputs=0, num_rst_outputs=0; // number of history and restart outputs

  // loop over input block names.  Find those that start with "output", read parameters,
  // and construct singly linked list of OutputTypes.
  while (pib != nullptr) {
    if (pib->block_name.compare(0, 6, "output") == 0) {
      OutputParameters op;  // define temporary OutputParameters struct

      // extract integer number of output block.  Save name and number
      std::string outn = pib->block_name.substr(6); // 6 because counting starts at 0!
      op.block_number = atoi(outn.c_str());
      op.block_name.assign(pib->block_name);

      // set time of last output, time between outputs
      op.next_time = pin->GetOrAddReal(op.block_name,"next_time", pm->time);
      op.dt = pin->GetOrAddReal(op.block_name,"dt", 0.0);
      op.dcycle = pin->GetOrAddInteger(op.block_name,"dcycle", 0);

      if (op.dt == 0.0 && op.dcycle == 0) {
        msg << "### FATAL ERROR in Outputs constructor" << std::endl
            << "Either dt or dcycle must be specified in " << op.block_name
            << std::endl;
        ATHENA_ERROR(msg);
      }
      if (op.dt > 0.0 && op.dcycle > 0) {
        msg << "### FATAL ERROR in Outputs constructor" << std::endl
            << "dt and dcycle cannot be specified simultaneously in " << op.block_name
            << std::endl;
        ATHENA_ERROR(msg);
      }
      if (op.dt > 0.0 || op.dcycle > 0) {  // only add output if dt > 0 or dycle > 0
        // set file number, basename, id, and format
        op.file_number = pin->GetOrAddInteger(op.block_name,"file_number",0);
        op.file_basename = pin->GetString("job","problem_id");
        char define_id[10];
        std::snprintf(define_id, sizeof(define_id),
                      "out%d", op.block_number);  // default id="outN"
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
            ATHENA_ERROR(msg);
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
            ATHENA_ERROR(msg);
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
            ATHENA_ERROR(msg);
          }
        }

        // read sum options.  Check for conflicts with slicing.
        op.output_sumx1 = pin->GetOrAddBoolean(op.block_name,"x1_sum",false);
        if ((op.output_slicex1) && (op.output_sumx1)) {
          msg << "### FATAL ERROR in Outputs constructor" << std::endl
              << "Cannot request both slice and sum along x1-direction"
              << " in output block '" << op.block_name << "'" << std::endl;
          ATHENA_ERROR(msg);
        }
        op.output_sumx2 = pin->GetOrAddBoolean(op.block_name,"x2_sum",false);
        if ((op.output_slicex2) && (op.output_sumx2)) {
          msg << "### FATAL ERROR in Outputs constructor" << std::endl
              << "Cannot request both slice and sum along x2-direction"
              << " in output block '" << op.block_name << "'" << std::endl;
          ATHENA_ERROR(msg);
        }
        op.output_sumx3 = pin->GetOrAddBoolean(op.block_name,"x3_sum",false);
        if ((op.output_slicex3) && (op.output_sumx3)) {
          msg << "### FATAL ERROR in Outputs constructor" << std::endl
              << "Cannot request both slice and sum along x3-direction"
              << " in output block '" << op.block_name << "'" << std::endl;
          ATHENA_ERROR(msg);
        }

        // read ghost cell option
        op.include_ghost_zones = pin->GetOrAddBoolean(op.block_name, "ghost_zones",
                                                      false);

        // read orbital system output option
        op.orbital_system_output = pin->GetOrAddBoolean(op.block_name, "orbital_system",
                                                        false);

        // read cartesian mapping option
        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0 ||
            std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
          op.cartesian_vector = pin->GetOrAddBoolean(op.block_name, "cartesian_vector",
                                                   false);
        else
          op.cartesian_vector = false;

        // set output variable and optional data format string used in formatted writes
        if (op.file_type.compare("hst") != 0 && op.file_type.compare("rst") != 0) {
          op.variable = pin->GetString(op.block_name, "variable");
        }
        op.data_format = pin->GetOrAddString(op.block_name, "data_format", "%12.5e");
        op.data_format.insert(0, " "); // prepend with blank to separate columns

        // Construct new OutputType according to file format
        // NEW_OUTPUT_TYPES: Add block to construct new types here
        if (op.file_type.compare("hst") == 0) {
          pnew_type = new HistoryOutput(op);
          num_hst_outputs++;
        } else if (op.file_type.compare("tab") == 0) {
          pnew_type = new FormattedTableOutput(op);
        } else if (op.file_type.compare("vtk") == 0) {
          pnew_type = new VTKOutput(op);
        } else if (op.file_type.compare("rst") == 0) {
          pnew_type = new RestartOutput(op);
          num_rst_outputs++;
        } else if (op.file_type.compare("ath5") == 0
                   || op.file_type.compare("hdf5") == 0) {
#ifdef HDF5OUTPUT
          pnew_type = new ATHDF5Output(op);
#else
          msg << "### FATAL ERROR in Outputs constructor" << std::endl
              << "Executable not configured for HDF5 outputs, but HDF5 file format "
              << "is requested in output block '" << op.block_name << "'" << std::endl;
          ATHENA_ERROR(msg);
#endif
        } else {
          msg << "### FATAL ERROR in Outputs constructor" << std::endl
              << "Unrecognized file format = '" << op.file_type
              << "' in output block '" << op.block_name << "'" << std::endl;
          ATHENA_ERROR(msg);
        }

        // Append type as tail node in singly linked list
        if (pfirst_type_ == nullptr) {
          pfirst_type_ = pnew_type;
        } else {
          plast->pnext_type = pnew_type;
        }
        plast = pnew_type;
      }
    }
    pib = pib->pnext;  // move to next input block name
  }

  // check there were no more than one history or restart files requested
  if (num_hst_outputs > 1 || num_rst_outputs > 1) {
    msg << "### FATAL ERROR in Outputs constructor" << std::endl
        << "More than one history or restart output block detected in input file"
        << std::endl;
    ATHENA_ERROR(msg);
  }

  // Move restarts to the tail end of the OutputType list, so file counters for other
  // output types are up-to-date in restart file
  int pos = 0, found = 0;
  OutputType *pot = pfirst_type_;
  OutputType *prst = pot;
  while (pot != nullptr) {
    if (pot->output_params.file_type.compare("rst") == 0) {
      prst = pot;
      found = 1;
      if (pot->pnext_type == nullptr) found = 2;
      break;
    }
    pos++;
    pot = pot->pnext_type;
  }
  if (found == 1) {
    // remove the restarting block
    pot = pfirst_type_;
    if (pos == 0) { // head node/first block
      pfirst_type_ = pfirst_type_->pnext_type;
    } else {
      for (int j=0; j<pos-1; j++) // seek the list
        pot = pot->pnext_type;
      pot->pnext_type = prst->pnext_type; // remove it
    }
    while (pot->pnext_type != nullptr)
      pot = pot->pnext_type; // find the tail node
    prst->pnext_type = nullptr;
    pot->pnext_type = prst;
  }
  // if found == 2, do nothing; it's already at the tail node/end of the list
}

// destructor - iterates through singly linked list of OutputTypes and deletes nodes

Outputs::~Outputs() {
  OutputType *ptype = pfirst_type_;
  while (ptype != nullptr) {
    OutputType *ptype_old = ptype;
    ptype = ptype->pnext_type;
    delete ptype_old;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void OutputType::LoadOutputData(MeshBlock *pmb)
//! \brief Create doubly linked list of OutputData's containing requested variables

void OutputType::LoadOutputData(MeshBlock *pmb) {
  Hydro *phyd = pmb->phydro;
  Field *pfld = pmb->pfield;
  PassiveScalars *psclr = pmb->pscalars;
  Gravity *pgrav = pmb->pgrav;
  OrbitalAdvection *porb = pmb->porb;
  num_vars_ = 0;
  OutputData *pod;

  // NEW_OUTPUT_TYPES:

  // (lab-frame) density
  if (ContainVariable(output_params.variable, "D") ||
      ContainVariable(output_params.variable, "cons")) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "dens";
    pod->data.InitWithShallowSlice(phyd->u, 4, IDN, 1);
    AppendOutputDataNode(pod);
    num_vars_++;
  }

  // (rest-frame) density
  if (ContainVariable(output_params.variable, "d") ||
      ContainVariable(output_params.variable, "prim")) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "rho";
    pod->data.InitWithShallowSlice(phyd->w, 4, IDN, 1);
    AppendOutputDataNode(pod);
    num_vars_++;
  }

  // total energy
  if (NON_BAROTROPIC_EOS) {
    if (ContainVariable(output_params.variable, "E") ||
        ContainVariable(output_params.variable, "cons")) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "Etot";
      if(porb->orbital_advection_defined
         && !output_params.orbital_system_output) {
        porb->ConvertOrbitalSystem(phyd->w, phyd->u, OrbitalTransform::cons);
        pod->data.InitWithShallowSlice(porb->u_orb, 4, IEN, 1);
      } else {
        pod->data.InitWithShallowSlice(phyd->u, 4, IEN, 1);
      }
      AppendOutputDataNode(pod);
      num_vars_++;
    }

    // pressure
    if (ContainVariable(output_params.variable, "p") ||
        ContainVariable(output_params.variable, "prim")) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "press";
      pod->data.InitWithShallowSlice(phyd->w, 4, IPR, 1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
  }

  // momentum vector
  if (ContainVariable(output_params.variable, "m") ||
      ContainVariable(output_params.variable, "cons")) {
    pod = new OutputData;
    pod->type = "VECTORS";
    pod->name = "mom";
    if(porb->orbital_advection_defined
       && !output_params.orbital_system_output) {
      porb->ConvertOrbitalSystem(phyd->w, phyd->u, OrbitalTransform::cons);
      pod->data.InitWithShallowSlice(porb->u_orb, 4, IM1, 3);
    } else {
      pod->data.InitWithShallowSlice(phyd->u, 4, IM1, 3);
    }
    AppendOutputDataNode(pod);
    num_vars_ += 3;
    if (output_params.cartesian_vector) {
      AthenaArray<Real> src;
      if(porb->orbital_advection_defined
         && !output_params.orbital_system_output) {
        src.InitWithShallowSlice(porb->u_orb, 4, IM1, 3);
      } else {
        src.InitWithShallowSlice(phyd->u, 4, IM1, 3);
      }
      pod = new OutputData;
      pod->type = "VECTORS";
      pod->name = "mom_xyz";
      pod->data.NewAthenaArray(3, phyd->u.GetDim3(), phyd->u.GetDim2(),
                               phyd->u.GetDim1());
      CalculateCartesianVector(src,  pod->data,  pmb->pcoord);
      AppendOutputDataNode(pod);
      num_vars_ += 3;
    }
  }

  // each component of momentum
  if (ContainVariable(output_params.variable, "m1")) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "mom1";
    pod->data.InitWithShallowSlice(phyd->u, 4, IM1, 1);
    AppendOutputDataNode(pod);
    num_vars_++;
  }
  if (ContainVariable(output_params.variable, "m2")) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "mom2";
    if(porb->orbital_advection_defined
       && !output_params.orbital_system_output
       && porb->orbital_direction == 1) {
      porb->ConvertOrbitalSystem(phyd->w, phyd->u, OrbitalTransform::cons);
      pod->data.InitWithShallowSlice(porb->u_orb, 4, IM2, 1);
    } else {
      pod->data.InitWithShallowSlice(phyd->u, 4, IM2, 1);
    }
    AppendOutputDataNode(pod);
    num_vars_++;
  }
  if (ContainVariable(output_params.variable, "m3")) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "mom3";
    if(porb->orbital_advection_defined
       && !output_params.orbital_system_output
       && porb->orbital_direction == 2) {
      porb->ConvertOrbitalSystem(phyd->w, phyd->u, OrbitalTransform::cons);
      pod->data.InitWithShallowSlice(porb->u_orb, 4, IM3, 1);
    } else {
      pod->data.InitWithShallowSlice(phyd->u, 4, IM3, 1);
    }
    AppendOutputDataNode(pod);
    num_vars_++;
  }

  // velocity vector
  if (ContainVariable(output_params.variable, "v") ||
      ContainVariable(output_params.variable, "prim")) {
    pod = new OutputData;
    pod->type = "VECTORS";
    pod->name = "vel";
    if(porb->orbital_advection_defined
       && !output_params.orbital_system_output) {
      porb->ConvertOrbitalSystem(phyd->w, phyd->u, OrbitalTransform::prim);
      pod->data.InitWithShallowSlice(porb->w_orb, 4, IVX, 3);
    } else {
      pod->data.InitWithShallowSlice(phyd->w, 4, IVX, 3);
    }
    AppendOutputDataNode(pod);
    num_vars_ += 3;
    if (output_params.cartesian_vector) {
      AthenaArray<Real> src;
      if(porb->orbital_advection_defined
         && !output_params.orbital_system_output) {
        src.InitWithShallowSlice(porb->w_orb, 4, IVX, 3);
      } else {
        src.InitWithShallowSlice(phyd->w, 4, IVX, 3);
      }
      pod = new OutputData;
      pod->type = "VECTORS";
      pod->name = "vel_xyz";
      pod->data.NewAthenaArray(3, phyd->w.GetDim3(), phyd->w.GetDim2(),
                               phyd->w.GetDim1());
      CalculateCartesianVector(src,  pod->data,  pmb->pcoord);
      AppendOutputDataNode(pod);
      num_vars_ += 3;
    }
  }

  // each component of velocity
  if (ContainVariable(output_params.variable, "vx") ||
      ContainVariable(output_params.variable, "v1")) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "vel1";
    pod->data.InitWithShallowSlice(phyd->w, 4, IVX, 1);
    AppendOutputDataNode(pod);
    num_vars_++;
  }
  if (ContainVariable(output_params.variable, "vy") ||
      ContainVariable(output_params.variable, "v2")) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "vel2";
    if(porb->orbital_advection_defined
       && !output_params.orbital_system_output
       && porb->orbital_direction == 1) {
      porb->ConvertOrbitalSystem(phyd->w, phyd->u, OrbitalTransform::prim);
      pod->data.InitWithShallowSlice(porb->w_orb, 4, IVY, 1);
    } else {
      pod->data.InitWithShallowSlice(phyd->w, 4, IVY, 1);
    }
    AppendOutputDataNode(pod);
    num_vars_++;
  }
  if (ContainVariable(output_params.variable, "vz") ||
      ContainVariable(output_params.variable, "v3")) {
    pod = new OutputData;
    pod->type = "SCALARS";
    pod->name = "vel3";
    if(porb->orbital_advection_defined
       && !output_params.orbital_system_output
       && porb->orbital_direction == 2) {
      porb->ConvertOrbitalSystem(phyd->w, phyd->u, OrbitalTransform::prim);
      pod->data.InitWithShallowSlice(porb->w_orb, 4, IVZ, 1);
    } else {
      pod->data.InitWithShallowSlice(phyd->w, 4, IVZ, 1);
    }
    AppendOutputDataNode(pod);
    num_vars_++;
  }

  if (SELF_GRAVITY_ENABLED) {
    if (ContainVariable(output_params.variable, "phi") ||
        ContainVariable(output_params.variable, "prim") ||
        ContainVariable(output_params.variable, "cons")) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "Phi";
      pod->data.InitWithShallowSlice(pgrav->phi, 4, 0, 1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
  } // endif (SELF_GRAVITY_ENABLED)

  if (NSCALARS > 0) {
    std::string root_name_cons = "s";
    std::string root_name_prim = "r";
    for (int n=0; n<NSCALARS; n++) {
      std::string scalar_name_cons = root_name_cons + std::to_string(n);
      std::string scalar_name_prim = root_name_prim + std::to_string(n);
      if (ContainVariable(output_params.variable, scalar_name_cons) ||
          ContainVariable(output_params.variable, "cons")) {
        pod = new OutputData;
        pod->type = "SCALARS";
        pod->name = scalar_name_cons;
        pod->data.InitWithShallowSlice(psclr->s, 4, n, 1);
        AppendOutputDataNode(pod);
        num_vars_++;
      }
      if (ContainVariable(output_params.variable, scalar_name_prim) ||
          ContainVariable(output_params.variable, "prim")) {
        pod = new OutputData;
        pod->type = "SCALARS";
        pod->name = scalar_name_prim;
        pod->data.InitWithShallowSlice(psclr->r, 4, n, 1);
        AppendOutputDataNode(pod);
        num_vars_++;
      }
    }
  }
  // note, the Bcc variables are stored in a separate HDF5 dataset from the above Output
  // nodes, and it must come after those nodes in the linked list
  if (MAGNETIC_FIELDS_ENABLED) {
    // vector of cell-centered magnetic field
    if (ContainVariable(output_params.variable, "bcc") ||
        ContainVariable(output_params.variable, "prim") ||
        ContainVariable(output_params.variable, "cons")) {
      pod = new OutputData;
      pod->type = "VECTORS";
      pod->name = "Bcc";
      pod->data.InitWithShallowSlice(pfld->bcc, 4, IB1, 3);
      AppendOutputDataNode(pod);
      num_vars_ += 3;
      if (output_params.cartesian_vector) {
        AthenaArray<Real> src;
        src.InitWithShallowSlice(pfld->bcc, 4, IB1, 3);
        pod = new OutputData;
        pod->type = "VECTORS";
        pod->name = "Bcc_xyz";
        pod->data.NewAthenaArray(3, pfld->bcc.GetDim3(), pfld->bcc.GetDim2(),
                                 pfld->bcc.GetDim1());
        CalculateCartesianVector(src, pod->data, pmb->pcoord);
        AppendOutputDataNode(pod);
        num_vars_ += 3;
      }
    }

    // each component of cell-centered magnetic field
    if (ContainVariable(output_params.variable, "bcc1")) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "Bcc1";
      pod->data.InitWithShallowSlice(pfld->bcc, 4, IB1, 1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
    if (ContainVariable(output_params.variable, "bcc2")) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "Bcc2";
      pod->data.InitWithShallowSlice(pfld->bcc, 4, IB2, 1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
    if (ContainVariable(output_params.variable, "bcc3")) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "Bcc3";
      pod->data.InitWithShallowSlice(pfld->bcc, 4, IB3, 1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
    // each component of face-centered magnetic field
    if (ContainVariable(output_params.variable, "b1")
        || ContainVariable(output_params.variable, "b")) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "B1";
      pod->data.InitWithShallowSlice(pfld->b.x1f, 4, 0, 1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
    if (ContainVariable(output_params.variable, "b2")
        || ContainVariable(output_params.variable, "b")) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "B2";
      pod->data.InitWithShallowSlice(pfld->b.x2f, 4, 0, 1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
    if (ContainVariable(output_params.variable, "b3")
        || ContainVariable(output_params.variable, "b")) {
      pod = new OutputData;
      pod->type = "SCALARS";
      pod->name = "B3";
      pod->data.InitWithShallowSlice(pfld->b.x3f, 4, 0, 1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
  } // endif (MAGNETIC_FIELDS_ENABLED)

  bool output_all_uov = ContainVariable(output_params.variable, "uov")
                        || ContainVariable(output_params.variable, "user_out_var");
  for (int n = 0; n < pmb->nuser_out_var; ++n) {
    char abbr_name[16], full_name[32];
    std::snprintf(abbr_name, sizeof(abbr_name), "uov%d", n);
    std::snprintf(full_name, sizeof(full_name), "user_out_var%d", n);
    if (output_all_uov ||
        (pmb->user_out_var_names_[n].length() != 0
         && ContainVariable(output_params.variable, pmb->user_out_var_names_[n]))
        || ContainVariable(output_params.variable, abbr_name)
        || ContainVariable(output_params.variable, full_name)) {
      pod = new OutputData;
      pod->type = "SCALARS";
      if (pmb->user_out_var_names_[n].length() != 0) {
        pod->name = pmb->user_out_var_names_[n];
      } else {
        pod->name = full_name;
      }
      pod->data.InitWithShallowSlice(pmb->user_out_var, 4, n, 1);
      AppendOutputDataNode(pod);
      num_vars_++;
    }
  }

  // throw an error if output variable name not recognized
  if (num_vars_ == 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [OutputType::LoadOutputData]" << std::endl
        << "Output variable '" << output_params.variable << "' not implemented"
        << std::endl;
    ATHENA_ERROR(msg);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutputType::AppendOutputDataNode(OutputData *pnew_data)
//! \brief

void OutputType::AppendOutputDataNode(OutputData *pnew_data) {
  if (pfirst_data_ == nullptr) {
    pfirst_data_ = pnew_data;
  } else {
    pnew_data->pprev = plast_data_;
    plast_data_->pnext = pnew_data;
  }
  // make the input node the new tail node of the doubly linked list
  plast_data_ = pnew_data;
}

//----------------------------------------------------------------------------------------
//! \fn void OutputType::ReplaceOutputDataNode(OutputData *pold, OutputData *pnew)
//! \brief

void OutputType::ReplaceOutputDataNode(OutputData *pold, OutputData *pnew) {
  if (pold == pfirst_data_) {
    pfirst_data_ = pnew;
    if (pold->pnext != nullptr) {    // there is another node in the list
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

//----------------------------------------------------------------------------------------
//! \fn void OutputType::ClearOutputData()
//! \brief

void OutputType::ClearOutputData() {
  OutputData *pdata = pfirst_data_;
  while (pdata != nullptr) {
    OutputData *pdata_old = pdata;
    pdata = pdata->pnext;
    delete pdata_old;
  }
  // reset pointers to head and tail nodes of doubly linked list:
  pfirst_data_ = nullptr;
  plast_data_  = nullptr;
}

//----------------------------------------------------------------------------------------
//! \fn void Outputs::MakeOutputs(Mesh *pm, ParameterInput *pin, bool wtflag)
//! \brief scans through singly linked list of OutputTypes and makes any outputs needed.

void Outputs::MakeOutputs(Mesh *pm, ParameterInput *pin, bool wtflag) {
  // wtflag = only true for making final outputs due to signal or wall-time/cycle/time
  // limit. Used by restart file output to change suffix to .final
  bool first=true;
  OutputType* ptype = pfirst_type_;
  while (ptype != nullptr) {
    if (((pm->time == pm->start_time) // output initial conditions, unless next_time set
         && (ptype->output_params.next_time <= pm->start_time ))
      || (ptype->output_params.dt > 0.0 && pm->time >= ptype->output_params.next_time)
      || (ptype->output_params.dcycle > 0 && pm->ncycle%ptype->output_params.dcycle == 0)
      || (pm->time >= pm->tlim)
      || (wtflag && ptype->output_params.file_type == "rst")) {
      if (first && ptype->output_params.file_type != "hst") {
        pm->ApplyUserWorkBeforeOutput(pin);
        first = false;
      }
      ptype->WriteOutputFile(pm, pin, wtflag);
    }
    ptype = ptype->pnext_type; // move to next OutputType node in singly linked list
  }
}

//----------------------------------------------------------------------------------------
//! \fn void OutputType::TransformOutputData(MeshBlock *pmb)
//! \brief Calls sum and slice functions on each direction in turn, in order to allow
//! mulitple operations performed on the same data set

bool OutputType::TransformOutputData(MeshBlock *pmb) {
  bool flag = true;
  if (output_params.output_slicex3) {
    bool ret = SliceOutputData(pmb,3);
    if (!ret) flag = false;
  }
  if (output_params.output_slicex2) {
    bool ret = SliceOutputData(pmb,2);
    if (!ret) flag = false;
  }
  if (output_params.output_slicex1) {
    bool ret = SliceOutputData(pmb,1);
    if (!ret) flag = false;
  }
  if (output_params.output_sumx3) {
    SumOutputData(pmb,3);
  }
  if (output_params.output_sumx2) {
    SumOutputData(pmb,2);
  }
  if (output_params.output_sumx1) {
    SumOutputData(pmb,1);
  }
  return flag;
}

//----------------------------------------------------------------------------------------
//! \fn bool OutputType::SliceOutputData(MeshBlock *pmb, int dim)
//! \brief perform data slicing and update the data list

bool OutputType::SliceOutputData(MeshBlock *pmb, int dim) {
  int islice(0), jslice(0), kslice(0);

  // Compute i,j,k indices of slice; check if in range of data in this block
  if (dim == 1) {
    if (output_params.x1_slice >= pmb->block_size.x1min &&
        output_params.x1_slice < pmb->block_size.x1max) {
      for (int i=pmb->is+1; i<=pmb->ie+1; ++i) {
        if (pmb->pcoord->x1f(i) > output_params.x1_slice) {
          islice = i - 1;
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
          jslice = j - 1;
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
          kslice = k - 1;
          output_params.kslice = kslice;
          break;
        }
      }
    } else {
      return false;
    }
  }

  // For each node in OutputData doubly linked list, slice arrays containing output data
  OutputData *pdata, *pnew;
  pdata = pfirst_data_;

  while (pdata != nullptr) {
    pnew = new OutputData;
    pnew->type = pdata->type;
    pnew->name = pdata->name;
    int nx4 = pdata->data.GetDim4();
    int nx3 = pdata->data.GetDim3();
    int nx2 = pdata->data.GetDim2();
    int nx1 = pdata->data.GetDim1();

    // Loop over variables and dimensions, extract slice
    if (dim == 3) {
      pnew->data.NewAthenaArray(nx4, 1, nx2, nx1);
      for (int n=0; n<nx4; ++n) {
        for (int j=out_js; j<=out_je; ++j) {
          for (int i=out_is; i<=out_ie; ++i) {
            pnew->data(n,0,j,i) = pdata->data(n,kslice,j,i);
          }
        }
      }
    } else if (dim == 2) {
      pnew->data.NewAthenaArray(nx4,nx3,1,nx1);
      for (int n=0; n<nx4; ++n) {
        for (int k=out_ks; k<=out_ke; ++k) {
          for (int i=out_is; i<=out_ie; ++i) {
            pnew->data(n,k,0,i) = pdata->data(n,k,jslice,i);
          }
        }
      }
    } else {
      pnew->data.NewAthenaArray(nx4,nx3,nx2,1);
      for (int n=0; n<nx4; ++n) {
        for (int k=out_ks; k<=out_ke; ++k) {
          for (int j=out_js; j<=out_je; ++j) {
            pnew->data(n,k,j,0) = pdata->data(n,k,j,islice);
          }
        }
      }
    }

    ReplaceOutputDataNode(pdata,pnew);
    pdata = pnew->pnext;
  }

  // modify array indices
  if (dim == 3) {
    out_ks = 0;
    out_ke = 0;
  } else if (dim == 2) {
    out_js = 0;
    out_je = 0;
  } else {
    out_is = 0;
    out_ie = 0;
  }
  return true;
}

//----------------------------------------------------------------------------------------
//! \fn void OutputType::SumOutputData(OutputData* pod, int dim)
//! \brief perform data summation and update the data list

void OutputType::SumOutputData(MeshBlock* pmb, int dim) {
  // For each node in OutputData doubly linked list, sum arrays containing output data
  OutputData *pdata = pfirst_data_;

  while (pdata != nullptr) {
    OutputData *pnew = new OutputData;
    pnew->type = pdata->type;
    pnew->name = pdata->name;
    int nx4 = pdata->data.GetDim4();
    int nx3 = pdata->data.GetDim3();
    int nx2 = pdata->data.GetDim2();
    int nx1 = pdata->data.GetDim1();

    // Loop over variables and dimensions, sum over specified dimension
    if (dim == 3) {
      pnew->data.NewAthenaArray(nx4, 1, nx2, nx1);
      for (int n=0; n<nx4; ++n) {
        for (int k=out_ks; k<=out_ke; ++k) {
          for (int j=out_js; j<=out_je; ++j) {
            for (int i=out_is; i<=out_ie; ++i) {
              pnew->data(n,0,j,i) += pdata->data(n,k,j,i);
            }
          }
        }
      }
    } else if (dim == 2) {
      pnew->data.NewAthenaArray(nx4, nx3, 1, nx1);
      for (int n=0; n<nx4; ++n) {
        for (int k=out_ks; k<=out_ke; ++k) {
          for (int j=out_js; j<=out_je; ++j) {
            for (int i=out_is; i<=out_ie; ++i) {
              pnew->data(n,k,0,i) += pdata->data(n,k,j,i);
            }
          }
        }
      }
    } else {
      pnew->data.NewAthenaArray(nx4, nx3, nx2, 1);
      for (int n=0; n<nx4; ++n) {
        for (int k=out_ks; k<=out_ke; ++k) {
          for (int j=out_js; j<=out_je; ++j) {
            for (int i=out_is; i<=out_ie; ++i) {
              pnew->data(n,k,j,0) += pdata->data(n,k,j,i);
            }
          }
        }
      }
    }

    ReplaceOutputDataNode(pdata, pnew);
    pdata = pdata->pnext;
  }

  // modify array indices
  if (dim == 3) {
    out_ks = 0;
    out_ke = 0;
  } else if (dim == 2) {
    out_js = 0;
    out_je = 0;
  } else {
    out_is = 0;
    out_ie = 0;
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void OutputType::CalculateCartesianVector(AthenaArray<Real> &src,
//!                               AthenaArray<Real> &dst, Coordinates *pco)
//! \brief Convert vectors in curvilinear coordinates into Cartesian

void OutputType::CalculateCartesianVector(AthenaArray<Real> &src, AthenaArray<Real> &dst,
                                          Coordinates *pco) {
  Real n1x, n1y, n1z, n2x, n2y, n2z, n3x, n3y, n3z;
  if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    if (out_ks == out_ke) { // 2D
      for (int k=out_ks; k<=out_ke; k++) {
        for (int j=out_js; j<=out_je; j++) {
          n1x = std::sin(pco->x2v(j));
          n1z = std::cos(pco->x2v(j));
          n2x = std::cos(pco->x2v(j));
          n2z = -std::sin(pco->x2v(j));
          for (int i=out_is; i<=out_ie; i++) {
            dst(0,k,j,i) = src(0,k,j,i)*n1x + src(1,k,j,i)*n2x;
            dst(1,k,j,i) = src(2,k,j,i);
            dst(2,k,j,i) = src(0,k,j,i)*n1z + src(1,k,j,i)*n2z;
          }
        }
      }
    } else { // 3D
      for (int k=out_ks; k<=out_ke; k++) {
        n3x = -std::sin(pco->x3v(k));
        n3y = std::cos(pco->x3v(k));
        n3z = 0.0;
        for (int j=out_js; j<=out_je; j++) {
          n1x = std::sin(pco->x2v(j))*std::cos(pco->x3v(k));
          n1y = std::sin(pco->x2v(j))*std::sin(pco->x3v(k));
          n1z = std::cos(pco->x2v(j));
          n2x = std::cos(pco->x2v(j))*std::cos(pco->x3v(k));
          n2y = std::cos(pco->x2v(j))*std::sin(pco->x3v(k));
          n2z = -std::sin(pco->x2v(j));
          for (int i=out_is; i<=out_ie; i++) {
            dst(0,k,j,i) = src(0,k,j,i)*n1x + src(1,k,j,i)*n2x + src(2,k,j,i)*n3x;
            dst(1,k,j,i) = src(0,k,j,i)*n1y + src(1,k,j,i)*n2y + src(2,k,j,i)*n3y;
            dst(2,k,j,i) = src(0,k,j,i)*n1z + src(1,k,j,i)*n2z + src(2,k,j,i)*n3z;
          }
        }
      }
    }
  }
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=out_ks; k<=out_ke; k++) {
      for (int j=out_js; j<=out_je; j++) {
        n1x = std::cos(pco->x2v(j));
        n1y = std::sin(pco->x2v(j));
        n2x =-std::sin(pco->x2v(j));
        n2y = std::cos(pco->x2v(j));
        for (int i=out_is; i<=out_ie; i++) {
          dst(0,k,j,i) = src(0,k,j,i)*n1x + src(1,k,j,i)*n2x;
          dst(1,k,j,i) = src(0,k,j,i)*n1y + src(1,k,j,i)*n2y;
          dst(2,k,j,i) = src(2,k,j,i);
        }
      }
    }
  }
  return;
}

bool OutputType::ContainVariable(const std::string &haystack, const std::string &needle) {
  if (haystack.compare(needle) == 0)
    return true;
  if (haystack.find(',' + needle + ',') != std::string::npos)
    return true;
  if (haystack.find(needle + ',') == 0)
    return true;
  if (haystack.find(',' + needle) != std::string::npos
    && haystack.find(',' + needle) == haystack.length() - needle.length() - 1)
    return true;
  return false;
}
