//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh.cpp
//  \brief implementation of functions in MeshBlock class

// C headers
#include <stdlib.h>
#include <string.h>  // memcpy

// C++ headers
#include <algorithm>  // sort
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../bvals/bvals_grav.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../fft/athena_fft.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../gravity/gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "mesh_refinement.hpp"
#include "meshblock_tree.hpp"
#include "mesh.hpp"

//----------------------------------------------------------------------------------------
// MeshBlock constructor: constructs coordinate, boundary condition, hydro, field
//                        and mesh refinement objects.

MeshBlock::MeshBlock(int igid, int ilid, LogicalLocation iloc, RegionSize input_block,
                     enum BoundaryFlag *input_bcs, Mesh *pm, ParameterInput *pin,
                     int igflag, bool ref_flag) {
  std::stringstream msg;
  int root_level;
  pmy_mesh = pm;
  root_level = pm->root_level;
  block_size = input_block;
  prev=NULL;
  next=NULL;
  gid=igid;
  lid=ilid;
  loc=iloc;
  gflag=igflag;
  cost=1.0;

  nuser_out_var = 0;
  nreal_user_meshblock_data_ = 0;
  nint_user_meshblock_data_ = 0;

  // initialize grid indices

  is = NGHOST;
  ie = is + block_size.nx1 - 1;

  if (block_size.nx2 > 1) {
    js = NGHOST;
    je = js + block_size.nx2 - 1;
  } else {
    js = je = 0;
  }

  if (block_size.nx3 > 1) {
    ks = NGHOST;
    ke = ks + block_size.nx3 - 1;
  } else {
    ks = ke = 0;
  }

  if (pm->multilevel==true) {
    cnghost=(NGHOST+1)/2+1;
    cis=cnghost; cie=cis+block_size.nx1/2-1;
    cjs=cje=cks=cke=0;
    if (block_size.nx2>1) // 2D or 3D
      cjs=cnghost, cje=cjs+block_size.nx2/2-1;
    if (block_size.nx3>1) // 3D
      cks=cnghost, cke=cks+block_size.nx3/2-1;
  }

  // construct objects stored in MeshBlock class.  Note in particular that the initial
  // conditions for the simulation are set in problem generator called from main, not
  // in the Hydro constructor

  // mesh-related objects


  // Boundary
  pbval  = new BoundaryValues(this, input_bcs, pin);

  // Coordinates
  if (COORDINATE_SYSTEM == "cartesian") {
    pcoord = new Cartesian(this, pin, false);
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    pcoord = new Cylindrical(this, pin, false);
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    pcoord = new SphericalPolar(this, pin, false);
  } else if (COORDINATE_SYSTEM == "minkowski") {
    pcoord = new Minkowski(this, pin, false);
  } else if (COORDINATE_SYSTEM == "schwarzschild") {
    pcoord = new Schwarzschild(this, pin, false);
  } else if (COORDINATE_SYSTEM == "kerr-schild") {
    pcoord = new KerrSchild(this, pin, false);
  } else if (COORDINATE_SYSTEM == "gr_user") {
    pcoord = new GRUser(this, pin, false);
  }

  if (SELF_GRAVITY_ENABLED) pgrav = new Gravity(this, pin);
  if (SELF_GRAVITY_ENABLED == 1) {
    pgbval = new GravityBoundaryValues(this,input_bcs);
  }

  // Reconstruction (constructor may implicitly depend on Coordinates, and PPM variable
  // floors depend on EOS (but EOS not needed by Reconstruction constructor)
  precon = new Reconstruction(this, pin);

  if (pm->multilevel==true) pmr = new MeshRefinement(this, pin);

  // physics-related objects: may depend on Coordinates for diffusion terms
  phydro = new Hydro(this, pin);
  if (MAGNETIC_FIELDS_ENABLED) pfield = new Field(this, pin);
  peos = new EquationOfState(this, pin);

  // Create user mesh data
  InitUserMeshBlockData(pin);

  return;
}

//----------------------------------------------------------------------------------------
// MeshBlock constructor for restarts

MeshBlock::MeshBlock(int igid, int ilid, Mesh *pm, ParameterInput *pin,
           LogicalLocation iloc, RegionSize input_block, enum BoundaryFlag *input_bcs,
           Real icost, char *mbdata, int igflag) {
  std::stringstream msg;
  pmy_mesh = pm;
  prev=NULL;
  next=NULL;
  gid=igid;
  lid=ilid;
  loc=iloc;
  gflag=igflag;
  cost=icost;
  block_size = input_block;

  nuser_out_var = 0;
  nreal_user_meshblock_data_ = 0;
  nint_user_meshblock_data_ = 0;

  // initialize grid indices
  is = NGHOST;
  ie = is + block_size.nx1 - 1;

  if (block_size.nx2 > 1) {
    js = NGHOST;
    je = js + block_size.nx2 - 1;
  } else {
    js = je = 0;
  }

  if (block_size.nx3 > 1) {
    ks = NGHOST;
    ke = ks + block_size.nx3 - 1;
  } else {
    ks = ke = 0;
  }

  if (pm->multilevel==true) {
    cnghost=(NGHOST+1)/2+1;
    cis=cnghost; cie=cis+block_size.nx1/2-1;
    cjs=cje=cks=cke=0;
    if (block_size.nx2>1) // 2D or 3D
      cjs=cnghost, cje=cjs+block_size.nx2/2-1;
    if (block_size.nx3>1) // 3D
      cks=cnghost, cke=cks+block_size.nx3/2-1;
  }

  // (re-)create mesh-related objects in MeshBlock

  // Boundary
  pbval  = new BoundaryValues(this, input_bcs, pin);

  // (re-)create physics-related objects in MeshBlock
  //phydro = new Hydro(this, pin);
  //if (MAGNETIC_FIELDS_ENABLED) pfield = new Field(this, pin);
  //peos = new EquationOfState(this, pin);

  if (SELF_GRAVITY_ENABLED) pgrav = new Gravity(this, pin);
  if (SELF_GRAVITY_ENABLED == 1) {
    pgbval = new GravityBoundaryValues(this,input_bcs);
  }

  // Coordinates
  if (COORDINATE_SYSTEM == "cartesian") {
    pcoord = new Cartesian(this, pin, false);
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    pcoord = new Cylindrical(this, pin, false);
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    pcoord = new SphericalPolar(this, pin, false);
  } else if (COORDINATE_SYSTEM == "minkowski") {
    pcoord = new Minkowski(this, pin, false);
  } else if (COORDINATE_SYSTEM == "schwarzschild") {
    pcoord = new Schwarzschild(this, pin, false);
  } else if (COORDINATE_SYSTEM == "kerr-schild") {
    pcoord = new KerrSchild(this, pin, false);
  } else if (COORDINATE_SYSTEM == "gr_user") {
    pcoord = new GRUser(this, pin, false);
  }

  // Reconstruction (constructor may implicitly depend on Coordinates)
  precon = new Reconstruction(this, pin);

  if (pm->multilevel==true) pmr = new MeshRefinement(this, pin);

  // (re-)create physics-related objects in MeshBlock
  phydro = new Hydro(this, pin);
  if (MAGNETIC_FIELDS_ENABLED) pfield = new Field(this, pin);
  peos = new EquationOfState(this, pin);
  InitUserMeshBlockData(pin);

  int os=0;
  // load hydro and field data
  memcpy(phydro->u.data(), &(mbdata[os]), phydro->u.GetSizeInBytes());
  // load it into the half-step arrays too
  memcpy(phydro->u1.data(), &(mbdata[os]), phydro->u1.GetSizeInBytes());
  os += phydro->u.GetSizeInBytes();
  if (GENERAL_RELATIVITY) {
    memcpy(phydro->w.data(), &(mbdata[os]), phydro->w.GetSizeInBytes());
    os += phydro->w.GetSizeInBytes();
    memcpy(phydro->w1.data(), &(mbdata[os]), phydro->w1.GetSizeInBytes());
    os += phydro->w1.GetSizeInBytes();
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    memcpy(pfield->b.x1f.data(), &(mbdata[os]), pfield->b.x1f.GetSizeInBytes());
    memcpy(pfield->b1.x1f.data(), &(mbdata[os]), pfield->b1.x1f.GetSizeInBytes());
    os += pfield->b.x1f.GetSizeInBytes();
    memcpy(pfield->b.x2f.data(), &(mbdata[os]), pfield->b.x2f.GetSizeInBytes());
    memcpy(pfield->b1.x2f.data(), &(mbdata[os]), pfield->b1.x2f.GetSizeInBytes());
    os += pfield->b.x2f.GetSizeInBytes();
    memcpy(pfield->b.x3f.data(), &(mbdata[os]), pfield->b.x3f.GetSizeInBytes());
    memcpy(pfield->b1.x3f.data(), &(mbdata[os]), pfield->b1.x3f.GetSizeInBytes());
    os += pfield->b.x3f.GetSizeInBytes();
  }

  // NEW_PHYSICS: add load of new physics from restart file here
  if (SELF_GRAVITY_ENABLED >= 1) {
    memcpy(pgrav->phi.data(), &(mbdata[os]), pgrav->phi.GetSizeInBytes());
    os += pgrav->phi.GetSizeInBytes();
  }

  // load user MeshBlock data
  for (int n=0; n<nint_user_meshblock_data_; n++) {
    memcpy(iuser_meshblock_data[n].data(), &(mbdata[os]),
           iuser_meshblock_data[n].GetSizeInBytes());
    os+=iuser_meshblock_data[n].GetSizeInBytes();
  }
  for (int n=0; n<nreal_user_meshblock_data_; n++) {
    memcpy(ruser_meshblock_data[n].data(), &(mbdata[os]),
           ruser_meshblock_data[n].GetSizeInBytes());
    os+=ruser_meshblock_data[n].GetSizeInBytes();
  }

  return;
}

//----------------------------------------------------------------------------------------
// MeshBlock destructor

MeshBlock::~MeshBlock() {
  if (prev!=NULL) prev->next=next;
  if (next!=NULL) next->prev=prev;

  delete pcoord;
  delete pbval;
  delete precon;
  if (pmy_mesh->multilevel == true) delete pmr;

  delete phydro;
  if (MAGNETIC_FIELDS_ENABLED) delete pfield;
  delete peos;
  if (SELF_GRAVITY_ENABLED) delete pgrav;
  if (SELF_GRAVITY_ENABLED==1) delete pgbval;

  // delete user output variables array
  if (nuser_out_var > 0) {
    user_out_var.DeleteAthenaArray();
    delete [] user_out_var_names_;
  }
  // delete user MeshBlock data
  for (int n=0; n<nreal_user_meshblock_data_; n++)
    ruser_meshblock_data[n].DeleteAthenaArray();
  if (nreal_user_meshblock_data_>0) delete [] ruser_meshblock_data;
  for (int n=0; n<nint_user_meshblock_data_; n++)
    iuser_meshblock_data[n].DeleteAthenaArray();
  if (nint_user_meshblock_data_>0) delete [] iuser_meshblock_data;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::AllocateRealUserMeshBlockDataField(int n)
//  \brief Allocate Real AthenaArrays for user-defned data in MeshBlock

void MeshBlock::AllocateRealUserMeshBlockDataField(int n) {
  if (nreal_user_meshblock_data_!=0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshBlock::AllocateRealUserMeshBlockDataField"
        << std::endl << "User MeshBlock data arrays are already allocated" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  nreal_user_meshblock_data_=n;
  ruser_meshblock_data = new AthenaArray<Real>[n];
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::AllocateIntUserMeshBlockDataField(int n)
//  \brief Allocate integer AthenaArrays for user-defned data in MeshBlock

void MeshBlock::AllocateIntUserMeshBlockDataField(int n) {
  if (nint_user_meshblock_data_!=0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshBlock::AllocateIntusermeshblockDataField"
        << std::endl << "User MeshBlock data arrays are already allocated" << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
  nint_user_meshblock_data_=n;
  iuser_meshblock_data = new AthenaArray<int>[n];
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::AllocateUserOutputVariables(int n)
//  \brief Allocate user-defined output variables

void MeshBlock::AllocateUserOutputVariables(int n) {
  if (n<=0) return;
  if (nuser_out_var!=0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshBlock::AllocateUserOutputVariables"
        << std::endl << "User output variables are already allocated." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
  nuser_out_var=n;
  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = block_size.nx2 + 2*(NGHOST);
  if (block_size.nx3 > 1) ncells3 = block_size.nx3 + 2*(NGHOST);
  user_out_var.NewAthenaArray(nuser_out_var,ncells3,ncells2,ncells1);
  user_out_var_names_ = new std::string[n];
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::SetUserOutputVariableName(int n, const char *name)
//  \brief set the user-defined output variable name

void MeshBlock::SetUserOutputVariableName(int n, const char *name) {
  if (n>=nuser_out_var) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshBlock::SetUserOutputVariableName"
        << std::endl << "User output variable is not allocated." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }
  user_out_var_names_[n]=name;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn size_t MeshBlock::GetBlockSizeInBytes(void)
//  \brief Calculate the block data size required for restart.

size_t MeshBlock::GetBlockSizeInBytes(void) {
  size_t size;

  size=phydro->u.GetSizeInBytes();
  if (GENERAL_RELATIVITY) {
    size+=phydro->w.GetSizeInBytes();
    size+=phydro->w1.GetSizeInBytes();
  }
  if (MAGNETIC_FIELDS_ENABLED)
    size+=(pfield->b.x1f.GetSizeInBytes()+pfield->b.x2f.GetSizeInBytes()
          +pfield->b.x3f.GetSizeInBytes());
  if (SELF_GRAVITY_ENABLED)
    size+=pgrav->phi.GetSizeInBytes();

  // NEW_PHYSICS: modify the size counter here when new physics is introduced

  // calculate user MeshBlock data size
  for (int n=0; n<nint_user_meshblock_data_; n++)
    size+=iuser_meshblock_data[n].GetSizeInBytes();
  for (int n=0; n<nreal_user_meshblock_data_; n++)
    size+=ruser_meshblock_data[n].GetSizeInBytes();

  return size;
}
