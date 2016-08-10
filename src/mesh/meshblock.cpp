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
//! \file mesh.cpp
//  \brief implementation of functions in MeshBlock class
//======================================================================================

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

//--------------------------------------------------------------------------------------
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

  // allocate user output variables array
  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = block_size.nx2 + 2*(NGHOST);
  if (block_size.nx3 > 1) ncells3 = block_size.nx3 + 2*(NGHOST);
  user_out_var.NewAthenaArray(NUSER_OUT_VAR,ncells3,ncells2,ncells1);

  nreal_user_meshblock_data_ = 0, nint_user_meshblock_data_ = 0; 

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
  pcoord = new Coordinates(this, pin);
  if(ref_flag==false) pcoord->CheckMeshSpacing();
  pbval  = new BoundaryValues(this, pin);
  if (block_bcs[INNER_X2] == POLAR_BNDRY) {
    int level = loc.level - pmy_mesh->root_level;
    int num_north_polar_blocks = pmy_mesh->nrbx3 * (1 << level);
    polar_neighbor_north = new PolarNeighborBlock[num_north_polar_blocks];
  }
  if (block_bcs[OUTER_X2] == POLAR_BNDRY) {
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

//--------------------------------------------------------------------------------------
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

  // allocate user output variables array
  int ncells1 = block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = block_size.nx2 + 2*(NGHOST);
  if (block_size.nx3 > 1) ncells3 = block_size.nx3 + 2*(NGHOST);
  user_out_var.NewAthenaArray(NUSER_OUT_VAR,ncells3,ncells2,ncells1);

  nreal_user_meshblock_data_ = 0, nint_user_meshblock_data_ = 0; 

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
  pcoord = new Coordinates(this, pin);
  pbval  = new BoundaryValues(this, pin);
  if (block_bcs[INNER_X2] == POLAR_BNDRY) {
    int level = loc.level - pmy_mesh->root_level;
    int num_north_polar_blocks = pmy_mesh->nrbx3 * (1 << level);
    polar_neighbor_north = new PolarNeighborBlock[num_north_polar_blocks];
  }
  if (block_bcs[OUTER_X2] == POLAR_BNDRY) {
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
    memcpy(iusermeshblockdata[n].data(), &(mbdata[os]),
           iusermeshblockdata[n].GetSizeInBytes());
    os+=iusermeshblockdata[n].GetSizeInBytes();
  }
  for(int n=0; n<nreal_user_meshblock_data_; n++) {
    memcpy(rusermeshblockdata[n].data(), &(mbdata[os]),
           rusermeshblockdata[n].GetSizeInBytes());
    os+=rusermeshblockdata[n].GetSizeInBytes();
  }

  return;
}

//--------------------------------------------------------------------------------------
// MeshBlock destructor

MeshBlock::~MeshBlock()
{
  if(prev!=NULL) prev->next=next;
  if(next!=NULL) next->prev=prev;

  delete pcoord;
  if (block_bcs[INNER_X2] == POLAR_BNDRY) delete[] polar_neighbor_north;
  if (block_bcs[OUTER_X2] == POLAR_BNDRY) delete[] polar_neighbor_south;
  delete pbval;
  delete precon;
  if (pmy_mesh->multilevel == true) delete pmr;

  delete phydro;
  if (MAGNETIC_FIELDS_ENABLED) delete pfield;
  delete peos;

  // delete user output variables array
  user_out_var.DeleteAthenaArray();
  // delete user MeshBlock data
  for(int n=0; n<nreal_user_meshblock_data_; n++)
    rusermeshblockdata[n].DeleteAthenaArray();
  if(nreal_user_meshblock_data_>0) delete [] rusermeshblockdata;
  for(int n=0; n<nint_user_meshblock_data_; n++)
    iusermeshblockdata[n].DeleteAthenaArray();
  if(nint_user_meshblock_data_>0) delete [] iusermeshblockdata;
}

//--------------------------------------------------------------------------------------
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
  rusermeshblockdata = new AthenaArray<Real>[n];
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void MeshBlock::AllocateIntUserMeshBlockDataField(int n)
//  \brief Allocate integer AthenaArrays for user-defned data in MeshBlock

void MeshBlock::AllocateIntUserMeshBlockDataField(int n)
{
  if(nint_user_meshblock_data_!=0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MeshBlock::AllocateIntusermeshblockDataField"
        << std::endl << "User MeshBlock data arrays are already allocated" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  nint_user_meshblock_data_=n;
  iusermeshblockdata = new AthenaArray<int>[n];
  return;
}

//--------------------------------------------------------------------------------------
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
    size+=iusermeshblockdata[n].GetSizeInBytes();
  for(int n=0; n<nreal_user_meshblock_data_; n++)
    size+=rusermeshblockdata[n].GetSizeInBytes();

  return size;
}

//--------------------------------------------------------------------------------------
// \!fn void NeighborBlock::SetNeighbor(int irank, int ilevel, int igid, int ilid,
//                          int iox1, int iox2, int iox3, enum NeighborType itype,
//                          int ibid, int itargetid, int ifi1=0, int ifi2=0,
//                          bool ipolar=false)
// \brief Set neighbor information

void NeighborBlock::SetNeighbor(int irank, int ilevel, int igid, int ilid,
  int iox1, int iox2, int iox3, enum NeighborType itype, int ibid, int itargetid,
  bool ipolar, int ifi1=0, int ifi2=0)
{
  rank=irank; level=ilevel; gid=igid; lid=ilid; ox1=iox1; ox2=iox2; ox3=iox3;
  type=itype; bufid=ibid; targetid=itargetid; polar=ipolar; fi1=ifi1; fi2=ifi2;
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

//--------------------------------------------------------------------------------------
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
              fid-nslist[ranklist[fid]], n, 0, 0, NEIGHBOR_FACE, bufid, tbid, false, f1,
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
      if(nlevel==loc.level) { // neighbor at same level
        tbid=FindBufferID(-n,0,0,0,0,pmy_mesh->maxneighbor_);
      }
      else { // neighbor at coarser level
        tbid=FindBufferID(-n,0,0,myfx2,myfx3,pmy_mesh->maxneighbor_);
      }
      neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
          nid-nslist[ranklist[nid]], n, 0, 0, NEIGHBOR_FACE, bufid, tbid, false);
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
              fid-nslist[ranklist[fid]], 0, n, 0, NEIGHBOR_FACE, bufid, tbid, false, f1,
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
          nid-nslist[ranklist[nid]], 0, n, 0, NEIGHBOR_FACE, bufid, tbid, polar);
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
                fid-nslist[ranklist[fid]], 0, 0, n, NEIGHBOR_FACE, bufid, tbid, false,
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
            nid-nslist[ranklist[nid]], 0, 0, n, NEIGHBOR_FACE, bufid, tbid, false);
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
              fid-nslist[ranklist[fid]], n, m, 0, NEIGHBOR_EDGE, bufid, tbid, false, f1,
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
        if(nlevel==loc.level) { // neighbor at same level
          if ((m == -1 and block_bcs[INNER_X2] == POLAR_BNDRY)
              or (m == 1 and block_bcs[OUTER_X2] == POLAR_BNDRY)) {
            polar = true; // neighbor is across top or bottom pole
          }
          tbid=FindBufferID(-n,polar?m:-m,0,0,0,pmy_mesh->maxneighbor_);
        }
        else { // neighbor at coarser level
          tbid=FindBufferID(-n,polar?m:-m,0,myfx3,0,pmy_mesh->maxneighbor_);
        }
        if(nlevel>=loc.level || (myox1==n && myox2==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
              nid-nslist[ranklist[nid]], n, m, 0, NEIGHBOR_EDGE, bufid, tbid, polar);
          nneighbor++;
        }
        bufid+=nf2;
      }
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
              fid-nslist[ranklist[fid]], n, 0, m, NEIGHBOR_EDGE, bufid, tbid, false, f1,
              0);
          bufid++; nneighbor++;
        }
      }
      else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[m+1][1][n+1]=nlevel;
        int tbid;
        if(nlevel==loc.level) { // neighbor at same level
          tbid=FindBufferID(-n,0,-m,0,0,pmy_mesh->maxneighbor_);
        }
        else { // neighbor at coarser level
          tbid=FindBufferID(-n,0,-m,myfx2,0,pmy_mesh->maxneighbor_);
        }
        if(nlevel>=loc.level || (myox1==n && myox3==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
              nid-nslist[ranklist[nid]], n, 0, m, NEIGHBOR_EDGE, bufid, tbid, false);
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
              fid-nslist[ranklist[fid]], 0, n, m, NEIGHBOR_EDGE, bufid, tbid, false, f1,
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
              nid-nslist[ranklist[nid]], 0, n, m, NEIGHBOR_EDGE, bufid, tbid, polar);
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
              nid-nslist[ranklist[nid]], n, m, l, NEIGHBOR_CORNER, bufid, tbid, polar);
          nneighbor++;
        }
        bufid++;
      }
    }
  }

  // polar neighbors
  if (block_bcs[INNER_X2] == POLAR_BNDRY) {
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
  if (block_bcs[OUTER_X2] == POLAR_BNDRY) {
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
  return;
}
