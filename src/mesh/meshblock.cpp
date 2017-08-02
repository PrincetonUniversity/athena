//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mesh.cpp
//  \brief implementation of functions in MeshBlock class

// C/C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <algorithm>  // sort
#include <iomanip>
#include <stdlib.h>
#include <string.h>  // memcpy

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
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
           enum BoundaryFlag *input_bcs, Mesh *pm, ParameterInput *pin, bool ref_flag)
{
  std::stringstream msg;
  int root_level;
  pmy_mesh = pm;
  root_level = pm->root_level;
  block_size = input_block;
  for(int i=0; i<6; i++) block_bcs[i] = input_bcs[i];
  prev=NULL;
  next=NULL;
  gid=igid;
  lid=ilid;
  loc=iloc;
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

  if(pm->multilevel==true) {
    cnghost=(NGHOST+1)/2+1;
    cis=cnghost; cie=cis+block_size.nx1/2-1;
    cjs=cje=cks=cke=0;
    if(block_size.nx2>1) // 2D or 3D
      cjs=cnghost, cje=cjs+block_size.nx2/2-1;
    if(block_size.nx3>1) // 3D
      cks=cnghost, cke=cks+block_size.nx3/2-1;
  }

  // construct objects stored in MeshBlock class.  Note in particular that the initial
  // conditions for the simulation are set in problem generator called from main, not
  // in the Hydro constructor

  // mesh-related objects
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

  pbval  = new BoundaryValues(this, pin);
  if (block_bcs[INNER_X2] == POLAR_BNDRY||block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) {
    int level = loc.level - pmy_mesh->root_level;
    int num_north_polar_blocks = pmy_mesh->nrbx3 * (1 << level);
    polar_neighbor_north = new PolarNeighborBlock[num_north_polar_blocks];
  }
  if (block_bcs[OUTER_X2] == POLAR_BNDRY||block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) {
    int level = loc.level - pmy_mesh->root_level;
    int num_south_polar_blocks = pmy_mesh->nrbx3 * (1 << level);
    polar_neighbor_south = new PolarNeighborBlock[num_south_polar_blocks];
  }
  precon = new Reconstruction(this, pin);
  if(pm->multilevel==true) pmr = new MeshRefinement(this, pin);

  // physics-related objects
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
           Real icost, char *mbdata)
{
  std::stringstream msg;
  pmy_mesh = pm;
  prev=NULL;
  next=NULL;
  gid=igid;
  lid=ilid;
  loc=iloc;
  cost=icost;
  block_size = input_block;
  for(int i=0; i<6; i++) block_bcs[i] = input_bcs[i];

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

  if(pm->multilevel==true) {
    cnghost=(NGHOST+1)/2+1;
    cis=cnghost; cie=cis+block_size.nx1/2-1;
    cjs=cje=cks=cke=0;
    if(block_size.nx2>1) // 2D or 3D
      cjs=cnghost, cje=cjs+block_size.nx2/2-1;
    if(block_size.nx3>1) // 3D
      cks=cnghost, cke=cks+block_size.nx3/2-1;
  }

  // (re-)create mesh-related objects in MeshBlock
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

  pbval  = new BoundaryValues(this, pin);
  if (block_bcs[INNER_X2] == POLAR_BNDRY||block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) {
    int level = loc.level - pmy_mesh->root_level;
    int num_north_polar_blocks = pmy_mesh->nrbx3 * (1 << level);
    polar_neighbor_north = new PolarNeighborBlock[num_north_polar_blocks];
  }
  if (block_bcs[OUTER_X2] == POLAR_BNDRY||block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) {
    int level = loc.level - pmy_mesh->root_level;
    int num_south_polar_blocks = pmy_mesh->nrbx3 * (1 << level);
    polar_neighbor_south = new PolarNeighborBlock[num_south_polar_blocks];
  }
  precon = new Reconstruction(this, pin);
  if(pm->multilevel==true) pmr = new MeshRefinement(this, pin);

  // (re-)create physics-related objects in MeshBlock
  phydro = new Hydro(this, pin);
  if (MAGNETIC_FIELDS_ENABLED) pfield = new Field(this, pin);
  peos = new EquationOfState(this, pin);

  InitUserMeshBlockData(pin);

  // load hydro and field data
  int os=0;
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

  // load user MeshBlock data
  for(int n=0; n<nint_user_meshblock_data_; n++) {
    memcpy(iuser_meshblock_data[n].data(), &(mbdata[os]),
           iuser_meshblock_data[n].GetSizeInBytes());
    os+=iuser_meshblock_data[n].GetSizeInBytes();
  }
  for(int n=0; n<nreal_user_meshblock_data_; n++) {
    memcpy(ruser_meshblock_data[n].data(), &(mbdata[os]),
           ruser_meshblock_data[n].GetSizeInBytes());
    os+=ruser_meshblock_data[n].GetSizeInBytes();
  }

  return;
}

//----------------------------------------------------------------------------------------
// MeshBlock destructor

MeshBlock::~MeshBlock()
{
  if(prev!=NULL) prev->next=next;
  if(next!=NULL) next->prev=prev;

  delete pcoord;
  if (block_bcs[INNER_X2] == POLAR_BNDRY||block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) delete[] polar_neighbor_north;
  if (block_bcs[OUTER_X2] == POLAR_BNDRY||block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) delete[] polar_neighbor_south;
  delete pbval;
  delete precon;
  if (pmy_mesh->multilevel == true) delete pmr;

  delete phydro;
  if (MAGNETIC_FIELDS_ENABLED) delete pfield;
  delete peos;

  // delete user output variables array
  if(nuser_out_var > 0) {
    user_out_var.DeleteAthenaArray();
    delete [] user_out_var_names_;
  }
  // delete user MeshBlock data
  for(int n=0; n<nreal_user_meshblock_data_; n++)
    ruser_meshblock_data[n].DeleteAthenaArray();
  if(nreal_user_meshblock_data_>0) delete [] ruser_meshblock_data;
  for(int n=0; n<nint_user_meshblock_data_; n++)
    iuser_meshblock_data[n].DeleteAthenaArray();
  if(nint_user_meshblock_data_>0) delete [] iuser_meshblock_data;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlock::AllocateRealUserMeshBlockDataField(int n)
//  \brief Allocate Real AthenaArrays for user-defned data in MeshBlock

void MeshBlock::AllocateRealUserMeshBlockDataField(int n)
{
  if(nreal_user_meshblock_data_!=0) {
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

void MeshBlock::AllocateIntUserMeshBlockDataField(int n)
{
  if(nint_user_meshblock_data_!=0) {
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

void MeshBlock::AllocateUserOutputVariables(int n)
{
  if(n<=0) return;
  if(nuser_out_var!=0) {
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

void MeshBlock::SetUserOutputVariableName(int n, const char *name)
{
  if(n>=nuser_out_var) {
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

size_t MeshBlock::GetBlockSizeInBytes(void)
{
  size_t size;

  size=phydro->u.GetSizeInBytes();
  if (GENERAL_RELATIVITY) {
    size+=phydro->w.GetSizeInBytes();
    size+=phydro->w1.GetSizeInBytes();
  }
  if (MAGNETIC_FIELDS_ENABLED)
    size+=(pfield->b.x1f.GetSizeInBytes()+pfield->b.x2f.GetSizeInBytes()
          +pfield->b.x3f.GetSizeInBytes());

  // NEW_PHYSICS: modify the size counter here when new physics is introduced

  // calculate user MeshBlock data size
  for(int n=0; n<nint_user_meshblock_data_; n++)
    size+=iuser_meshblock_data[n].GetSizeInBytes();
  for(int n=0; n<nreal_user_meshblock_data_; n++)
    size+=ruser_meshblock_data[n].GetSizeInBytes();

  return size;
}

//----------------------------------------------------------------------------------------
// \!fn void NeighborBlock::SetNeighbor(int irank, int ilevel, int igid, int ilid,
//                          int iox1, int iox2, int iox3, enum NeighborType itype,
//                          int ibid, int itargetid, int ifi1=0, int ifi2=0,
//                          bool ipolar=false)
// \brief Set neighbor information

void NeighborBlock::SetNeighbor(int irank, int ilevel, int igid, int ilid,
  int iox1, int iox2, int iox3, enum NeighborType itype, int ibid, int itargetid,
//[JMSHI
  bool ipolar, bool ishear, int ifi1=0, int ifi2=0)
  //bool ipolar, int ifi1=0, int ifi2=0)
{
  rank=irank; level=ilevel; gid=igid; lid=ilid; ox1=iox1; ox2=iox2; ox3=iox3;
  type=itype; bufid=ibid; targetid=itargetid; polar=ipolar; fi1=ifi1; fi2=ifi2;
  shear=ishear;
//JMSHI]
  if(type==NEIGHBOR_FACE) {
    if(ox1==-1)      fid=INNER_X1;
    else if(ox1==1)  fid=OUTER_X1;
    else if(ox2==-1) fid=INNER_X2;
    else if(ox2==1)  fid=OUTER_X2;
    else if(ox3==-1) fid=INNER_X3;
    else if(ox3==1)  fid=OUTER_X3;
  }
  if(type==NEIGHBOR_EDGE) {
    if(ox3==0)      eid=(   ((ox1+1)>>1) | ((ox2+1)&2));
    else if(ox2==0) eid=(4+(((ox1+1)>>1) | ((ox3+1)&2)));
    else if(ox1==0) eid=(8+(((ox2+1)>>1) | ((ox3+1)&2)));
  }
  return;
}

//----------------------------------------------------------------------------------------
// \!fn void MeshBlock::SearchAndSetNeighbors(MeshBlockTree &tree, int *ranklist, int *nslist)
// \brief Search and set all the neighbor blocks

void MeshBlock::SearchAndSetNeighbors(MeshBlockTree &tree, int *ranklist, int *nslist)
{
  MeshBlockTree* neibt;
  int myox1, myox2=0, myox3=0, myfx1, myfx2, myfx3;
  myfx1=(int)(loc.lx1&1L);
  myfx2=(int)(loc.lx2&1L);
  myfx3=(int)(loc.lx3&1L);
  myox1=((int)(loc.lx1&1L))*2-1;
  if(block_size.nx2>1) myox2=((int)(loc.lx2&1L))*2-1;
  if(block_size.nx3>1) myox3=((int)(loc.lx3&1L))*2-1;
  long int nrbx1=pmy_mesh->nrbx1, nrbx2=pmy_mesh->nrbx2, nrbx3=pmy_mesh->nrbx3;

  int nf1=1, nf2=1;
  if(pmy_mesh->multilevel==true) {
    if(block_size.nx2>1) nf1=2;
    if(block_size.nx3>1) nf2=2;
  }
  int bufid=0;
  nneighbor=0;
  for(int k=0; k<=2; k++) {
    for(int j=0; j<=2; j++) {
      for(int i=0; i<=2; i++)
        nblevel[k][j][i]=-1;
    }
  }
  nblevel[1][1][1]=loc.level;

  // x1 face
  for(int n=-1; n<=1; n+=2) {
    neibt=tree.FindNeighbor(loc,n,0,0,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
    if(neibt==NULL) { bufid+=nf1*nf2; continue;}
    if(neibt->flag==false) { // neighbor at finer level
      int fface=1-(n+1)/2; // 0 for OUTER_X1, 1 for INNER_X1
      nblevel[1][1][n+1]=neibt->loc.level+1;
      for(int f2=0;f2<nf2;f2++) {
        for(int f1=0;f1<nf1;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(fface,f1,f2);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(-n,0,0,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
          //[JMSHI
              fid-nslist[ranklist[fid]], n, 0, 0, NEIGHBOR_FACE, bufid, tbid, false, false, f1,
              //fid-nslist[ranklist[fid]], n, 0, 0, NEIGHBOR_FACE, bufid, tbid, false, f1,
          //JMSHI]
              f2);
          bufid++; nneighbor++;
        }
      }
    }
    else { // neighbor at same or coarser level
      int nlevel=neibt->loc.level;
      int nid=neibt->gid;
      nblevel[1][1][n+1]=nlevel;
      int tbid;
      //[JMSHI
      bool shear=false;
      if(nlevel==loc.level) { // neighbor at same level
        tbid=FindBufferID(-n,0,0,0,0,pmy_mesh->maxneighbor_);
        if ((n == -1 and block_bcs[INNER_X1] == SHEAR_PERIODIC_BNDRY)
            or (n == 1 and block_bcs[OUTER_X1] == SHEAR_PERIODIC_BNDRY)) {
          shear = true; // neighbor is shearing periodic
        }
      //JMSHI]
      }
      else { // neighbor at coarser level
        tbid=FindBufferID(-n,0,0,myfx2,myfx3,pmy_mesh->maxneighbor_);
      }
      neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
      //[JMSHI
          nid-nslist[ranklist[nid]], n, 0, 0, NEIGHBOR_FACE, bufid, tbid, false, shear);
          //nid-nslist[ranklist[nid]], n, 0, 0, NEIGHBOR_FACE, bufid, tbid, false);
      //JMSHI]
      bufid+=nf1*nf2; nneighbor++;
    }
  }
  if(block_size.nx2==1) return;

  // x2 face
  for(int n=-1; n<=1; n+=2) {
    neibt=tree.FindNeighbor(loc,0,n,0,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
    if(neibt==NULL) { bufid+=nf1*nf2; continue;}
    if(neibt->flag==false) { // neighbor at finer level
      int fface=1-(n+1)/2; // 0 for OUTER_X2, 1 for INNER_X2
      nblevel[1][n+1][1]=neibt->loc.level+1;
      for(int f2=0;f2<nf2;f2++) {
        for(int f1=0;f1<nf1;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(f1,fface,f2);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(0,-n,0,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
          //[JMSHI
              fid-nslist[ranklist[fid]], 0, n, 0, NEIGHBOR_FACE, bufid, tbid, false, false, f1,
              //fid-nslist[ranklist[fid]], 0, n, 0, NEIGHBOR_FACE, bufid, tbid, false, f1,
          //JMSHI]
              f2);
          bufid++; nneighbor++;
        }
      }
    }
    else { // neighbor at same or coarser level
      int nlevel=neibt->loc.level;
      int nid=neibt->gid;
      nblevel[1][n+1][1]=nlevel;
      int tbid;
      bool polar=false;
      if(nlevel==loc.level) { // neighbor at same level
        if ((n == -1 and block_bcs[INNER_X2] == POLAR_BNDRY)
            or (n == 1 and block_bcs[OUTER_X2] == POLAR_BNDRY)) {
          polar = true; // neighbor is across top or bottom pole
        }
        tbid=FindBufferID(0,polar?n:-n,0,0,0,pmy_mesh->maxneighbor_);
      }
      else { // neighbor at coarser level
        tbid=FindBufferID(0,-n,0,myfx1,myfx3,pmy_mesh->maxneighbor_);
      }
      neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
      //[JMSHI
          nid-nslist[ranklist[nid]], 0, n, 0, NEIGHBOR_FACE, bufid, tbid, polar, false);
          //nid-nslist[ranklist[nid]], 0, n, 0, NEIGHBOR_FACE, bufid, tbid, polar);
      //JMSHI]
      bufid+=nf1*nf2; nneighbor++;
    }
  }

  // x3 face
  if(block_size.nx3>1) {
    for(int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(loc,0,0,n,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
      if(neibt==NULL) { bufid+=nf1*nf2; continue;}
      if(neibt->flag==false) { // neighbor at finer level
        int fface=1-(n+1)/2; // 0 for OUTER_X3, 1 for INNER_X3
        nblevel[n+1][1][1]=neibt->loc.level+1;
        for(int f2=0;f2<nf2;f2++) {
          for(int f1=0;f1<nf1;f1++) {
            MeshBlockTree* nf=neibt->GetLeaf(f1,f2,fface);
            int fid = nf->gid;
            int nlevel=nf->loc.level;
            int tbid=FindBufferID(0,0,-n,0,0,pmy_mesh->maxneighbor_);
            neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
      //[JMSHI
                fid-nslist[ranklist[fid]], 0, 0, n, NEIGHBOR_FACE, bufid, tbid, false, false,
                //fid-nslist[ranklist[fid]], 0, 0, n, NEIGHBOR_FACE, bufid, tbid, false,
      //JMSHI]
                f1, f2);
            bufid++; nneighbor++;
          }
        }
      }
      else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[n+1][1][1]=nlevel;
        int tbid;
        if(nlevel==loc.level) { // neighbor at same level
          tbid=FindBufferID(0,0,-n,0,0,pmy_mesh->maxneighbor_);
        }
        else { // neighbor at coarser level
          tbid=FindBufferID(0,0,-n,myfx1,myfx2,pmy_mesh->maxneighbor_);
        }
        neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
      //[JMSHI
            nid-nslist[ranklist[nid]], 0, 0, n, NEIGHBOR_FACE, bufid, tbid, false, false);
            //nid-nslist[ranklist[nid]], 0, 0, n, NEIGHBOR_FACE, bufid, tbid, false);
      //JMSHI]
        bufid+=nf1*nf2; nneighbor++;
      }
    }
  }

  // x1x2 edge
  for(int m=-1; m<=1; m+=2) {
    for(int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(loc,n,m,0,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
      if(neibt==NULL) { bufid+=nf2; continue;}
      if(neibt->flag==false) { // neighbor at finer level
        int ff1=1-(n+1)/2; // 0 for OUTER_X1, 1 for INNER_X1
        int ff2=1-(m+1)/2; // 0 for OUTER_X2, 1 for INNER_X2
        nblevel[1][m+1][n+1]=neibt->loc.level+1;
        for(int f1=0;f1<nf2;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(ff1,ff2,f1);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(-n,-m,0,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
          //[JMSHI
              fid-nslist[ranklist[fid]], n, m, 0, NEIGHBOR_EDGE, bufid, tbid, false, false, f1,
              //fid-nslist[ranklist[fid]], n, m, 0, NEIGHBOR_EDGE, bufid, tbid, false, f1,
          //JMSHI]
              0);
          bufid++; nneighbor++;
        }
      }
      else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[1][m+1][n+1]=nlevel;
        int tbid;
        bool polar=false;
        //[JMSHI
        bool shear=false;
        //JMSHI]
        if(nlevel==loc.level) { // neighbor at same level
          if ((m == -1 and block_bcs[INNER_X2] == POLAR_BNDRY)
              or (m == 1 and block_bcs[OUTER_X2] == POLAR_BNDRY)) {
            polar = true; // neighbor is across top or bottom pole
          }
          //[JMSHI
          if ((n == -1 and block_bcs[INNER_X1] == SHEAR_PERIODIC_BNDRY)
              or (n == 1 and block_bcs[OUTER_X1] == SHEAR_PERIODIC_BNDRY)) {
            shear = true; // neighbor is on shearing periodic bcs
          }
          //JMSHI]
          tbid=FindBufferID(-n,polar?m:-m,0,0,0,pmy_mesh->maxneighbor_);
        }
        else { // neighbor at coarser level
          tbid=FindBufferID(-n,polar?m:-m,0,myfx3,0,pmy_mesh->maxneighbor_);
        }
        if(nlevel>=loc.level || (myox1==n && myox2==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
          //[JMSHI
              nid-nslist[ranklist[nid]], n, m, 0, NEIGHBOR_EDGE, bufid, tbid, polar, shear);
              //nid-nslist[ranklist[nid]], n, m, 0, NEIGHBOR_EDGE, bufid, tbid, polar);
          //JMSHI]
          nneighbor++;
        }
        bufid+=nf2;
      }
    }
  }

  // polar neighbors
  if (block_bcs[INNER_X2] == POLAR_BNDRY||block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) {
    int level = loc.level - pmy_mesh->root_level;
    int num_north_polar_blocks = nrbx3 * (1 << level);
    for (int n = 0; n < num_north_polar_blocks; ++n) {
      LogicalLocation neighbor_loc;
      neighbor_loc.lx1 = loc.lx1;
      neighbor_loc.lx2 = loc.lx2;
      neighbor_loc.lx3 = n;
      neighbor_loc.level = loc.level;
      neibt = tree.FindMeshBlock(neighbor_loc);
      int nid = neibt->gid;
      polar_neighbor_north[neibt->loc.lx3].rank = ranklist[nid];
      polar_neighbor_north[neibt->loc.lx3].lid = nid - nslist[ranklist[nid]];
      polar_neighbor_north[neibt->loc.lx3].gid = nid;
      polar_neighbor_north[neibt->loc.lx3].north = true;
    }
  }
  if (block_bcs[OUTER_X2] == POLAR_BNDRY||block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) {
    int level = loc.level - pmy_mesh->root_level;
    int num_south_polar_blocks = nrbx3 * (1 << level);
    for (int n = 0; n < num_south_polar_blocks; ++n) {
      LogicalLocation neighbor_loc;
      neighbor_loc.lx1 = loc.lx1;
      neighbor_loc.lx2 = loc.lx2;
      neighbor_loc.lx3 = n;
      neighbor_loc.level = loc.level;
      neibt = tree.FindMeshBlock(neighbor_loc);
      int nid = neibt->gid;
      polar_neighbor_south[neibt->loc.lx3].rank = ranklist[nid];
      polar_neighbor_south[neibt->loc.lx3].lid = nid - nslist[ranklist[nid]];
      polar_neighbor_south[neibt->loc.lx3].gid = nid;
      polar_neighbor_south[neibt->loc.lx3].north = false;
    }
  }
  if(block_size.nx3==1) return;

  // x1x3 edge
  for(int m=-1; m<=1; m+=2) {
    for(int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(loc,n,0,m,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
      if(neibt==NULL) { bufid+=nf1; continue;}
      if(neibt->flag==false) { // neighbor at finer level
        int ff1=1-(n+1)/2; // 0 for OUTER_X1, 1 for INNER_X1
        int ff2=1-(m+1)/2; // 0 for OUTER_X3, 1 for INNER_X3
        nblevel[m+1][1][n+1]=neibt->loc.level+1;
        for(int f1=0;f1<nf1;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(ff1,f1,ff2);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(-n,0,-m,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
          //[JMSHI
              fid-nslist[ranklist[fid]], n, 0, m, NEIGHBOR_EDGE, bufid, tbid, false, false, f1,
              //fid-nslist[ranklist[fid]], n, 0, m, NEIGHBOR_EDGE, bufid, tbid, false, f1,
          //JMSHI]
              0);
          bufid++; nneighbor++;
        }
      }
      else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[m+1][1][n+1]=nlevel;
        int tbid;
        //[JMSHI
        bool shear=false;
        //JMSHI]
        if(nlevel==loc.level) { // neighbor at same level
          tbid=FindBufferID(-n,0,-m,0,0,pmy_mesh->maxneighbor_);
          //[JMSHI
          if ((n == -1 and block_bcs[INNER_X1] == SHEAR_PERIODIC_BNDRY)
              or (n == 1 and block_bcs[OUTER_X1] == SHEAR_PERIODIC_BNDRY)) {
            shear = true; // neighbor is across top or bottom pole
          }
          //JMSHI]
        }
        else { // neighbor at coarser level
          tbid=FindBufferID(-n,0,-m,myfx2,0,pmy_mesh->maxneighbor_);
        }
        if(nlevel>=loc.level || (myox1==n && myox3==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
          //[JMSHI
              nid-nslist[ranklist[nid]], n, 0, m, NEIGHBOR_EDGE, bufid, tbid, false, shear);
              //nid-nslist[ranklist[nid]], n, 0, m, NEIGHBOR_EDGE, bufid, tbid, false);
          //JMSHI]
          nneighbor++;
        }
        bufid+=nf1;
      }
    }
  }

  // x2x3 edge
  for(int m=-1; m<=1; m+=2) {
    for(int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(loc,0,n,m,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
      if(neibt==NULL) { bufid+=nf1; continue;}
      if(neibt->flag==false) { // neighbor at finer level
        int ff1=1-(n+1)/2; // 0 for OUTER_X2, 1 for INNER_X2
        int ff2=1-(m+1)/2; // 0 for OUTER_X3, 1 for INNER_X3
        nblevel[m+1][n+1][1]=neibt->loc.level+1;
        for(int f1=0;f1<nf1;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(f1,ff1,ff2);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(0,-n,-m,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
          //[JMSHI
              fid-nslist[ranklist[fid]], 0, n, m, NEIGHBOR_EDGE, bufid, tbid, false, false, f1,
              //fid-nslist[ranklist[fid]], 0, n, m, NEIGHBOR_EDGE, bufid, tbid, false, f1,
          //JMSHI]
              0);
          bufid++; nneighbor++;
        }
      }
      else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[m+1][n+1][1]=nlevel;
        int tbid;
        bool polar=false;
        if(nlevel==loc.level) { // neighbor at same level
          if ((n == -1 and block_bcs[INNER_X2] == POLAR_BNDRY)
              or (n == 1 and block_bcs[OUTER_X2] == POLAR_BNDRY)) {
            polar = true; // neighbor is across top or bottom pole
          }
          tbid=FindBufferID(0,polar?n:-n,-m,0,0,pmy_mesh->maxneighbor_);
        }
        else { // neighbor at coarser level
          tbid=FindBufferID(0,-n,-m,myfx1,0,pmy_mesh->maxneighbor_);
        }
        if(nlevel>=loc.level || (myox2==n && myox3==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
          //[JMSHI
              nid-nslist[ranklist[nid]], 0, n, m, NEIGHBOR_EDGE, bufid, tbid, polar, false);
              //nid-nslist[ranklist[nid]], 0, n, m, NEIGHBOR_EDGE, bufid, tbid, polar);
          //JMSHI]
          nneighbor++;
        }
        bufid+=nf1;
      }
    }
  }

  // corners
  for(int l=-1; l<=1; l+=2) {
    for(int m=-1; m<=1; m+=2) {
      for(int n=-1; n<=1; n+=2) {
        neibt=tree.FindNeighbor(loc,n,m,l,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh->root_level);
        if(neibt==NULL) { bufid++; continue;}
        bool polar=false;
        if ((m == -1 and block_bcs[INNER_X2] == POLAR_BNDRY)
            or (m == 1 and block_bcs[OUTER_X2] == POLAR_BNDRY)) {
          polar = true; // neighbor is across top or bottom pole
        }
        if(neibt->flag==false) { // neighbor at finer level
          int ff1=1-(n+1)/2; // 0 for OUTER_X1, 1 for INNER_X1
          int ff2=1-(m+1)/2; // 0 for OUTER_X2, 1 for INNER_X2
          int ff3=1-(l+1)/2; // 0 for OUTER_X3, 1 for INNER_X3
          neibt=neibt->GetLeaf(ff1,ff2,ff3);
        }
        int nlevel=neibt->loc.level;
        nblevel[l+1][m+1][n+1]=nlevel;
        if(nlevel>=loc.level || (myox1==n && myox2==m && myox3==l)) {
          int nid=neibt->gid;
          int tbid=FindBufferID(-n,polar?m:-m,-l,0,0,pmy_mesh->maxneighbor_);
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
          //[JMSHI
              nid-nslist[ranklist[nid]], n, m, l, NEIGHBOR_CORNER, bufid, tbid, polar, false);
              //nid-nslist[ranklist[nid]], n, m, l, NEIGHBOR_CORNER, bufid, tbid, polar);
          //[JMSHI
          nneighbor++;
        }
        bufid++;
      }
    }
  }

  return;
}
