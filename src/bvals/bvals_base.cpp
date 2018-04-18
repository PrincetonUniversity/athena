//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_buffer.cpp
//  \brief utility functions for BoundaryValues buffers

// C++ headers
#include <iostream>   // endl
#include <iomanip>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cstring>    // memcpy
#include <cstdlib>
#include <cmath>

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/buffer_utils.hpp"

// this class header
#include "bvals.hpp"

// forward declaration of static members of this class
NeighborIndexes BoundaryBase::ni[56];
int BoundaryBase::bufid[56];
bool BoundaryBase::called_ = false;
int BoundaryBase::maxneighbor_ ;

//----------------------------------------------------------------------------------------
// \!fn void NeighborBlock::SetNeighbor(int irank, int ilevel, int igid, int ilid,
//                          int iox1, int iox2, int iox3, enum NeighborType itype,
//                          int ibid, int itargetid, int ifi1=0, int ifi2=0,
//                          bool ipolar=false)
// \brief Set neighbor information

void NeighborBlock::SetNeighbor(int irank, int ilevel, int igid, int ilid,
  int iox1, int iox2, int iox3, enum NeighborType itype, int ibid, int itargetid,
  bool ipolar, bool ishear, int ifi1=0, int ifi2=0) {
  rank=irank; level=ilevel; gid=igid; lid=ilid; ox1=iox1; ox2=iox2; ox3=iox3;
  type=itype; bufid=ibid; targetid=itargetid; polar=ipolar; shear=ishear;
  fi1=ifi1; fi2=ifi2;
  if (type==NEIGHBOR_FACE) {
    if (ox1==-1)      fid=INNER_X1;
    else if (ox1==1)  fid=OUTER_X1;
    else if (ox2==-1) fid=INNER_X2;
    else if (ox2==1)  fid=OUTER_X2;
    else if (ox3==-1) fid=INNER_X3;
    else if (ox3==1)  fid=OUTER_X3;
  }
  if (type==NEIGHBOR_EDGE) {
    if (ox3==0)      eid=(   ((ox1+1)>>1) | ((ox2+1)&2));
    else if (ox2==0) eid=(4+(((ox1+1)>>1) | ((ox3+1)&2)));
    else if (ox1==0) eid=(8+(((ox2+1)>>1) | ((ox3+1)&2)));
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn BoundaryBase::BoundaryBase(Mesh *pm, LogicalLocation iloc, RegionSize isize,
//                                 enum BoundaryFlag *input_bcs)
//  \brief constructor of BoundaryBase
BoundaryBase::BoundaryBase(Mesh *pm, LogicalLocation iloc, RegionSize isize,
                           enum BoundaryFlag *input_bcs) {
  loc=iloc;
  block_size_=isize;
  pmy_mesh_=pm;
  if (called_==false) {
    int dim=1;
    if (block_size_.nx2>1) dim=2;
    if (block_size_.nx3>1) dim=3;
    maxneighbor_=BufferID(dim, pmy_mesh_->multilevel);
    called_=true;
  }

  for (int i=0; i<6; i++)
    block_bcs[i]=input_bcs[i];
  if (block_bcs[INNER_X2] == POLAR_BNDRY
   || block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) {
    int level = loc.level - pmy_mesh_->root_level;
    // possible loss of precision to 32 bit int, if int64_t nrbx3 is large
    int num_north_polar_blocks = static_cast<int>(pmy_mesh_->nrbx3 * (1 << level));
    polar_neighbor_north = new PolarNeighborBlock[num_north_polar_blocks];
  }
  if (block_bcs[OUTER_X2] == POLAR_BNDRY
   || block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE) {
    int level = loc.level - pmy_mesh_->root_level;
    int num_south_polar_blocks = static_cast<int>(pmy_mesh_->nrbx3 * (1 << level));
    polar_neighbor_south = new PolarNeighborBlock[num_south_polar_blocks];
  }

  if (pmy_mesh_->multilevel==true) { // SMR or AMR
    // allocate surface area array
    int nc1=block_size_.nx1+2*NGHOST;
    sarea_[0].NewAthenaArray(nc1);
    sarea_[1].NewAthenaArray(nc1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn BoundaryBase::~BoundaryBase()
//  \brief destructor of BoundaryBase
BoundaryBase::~BoundaryBase() {
  if (block_bcs[INNER_X2] == POLAR_BNDRY
   || block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE)
    delete [] polar_neighbor_north;
  if (block_bcs[OUTER_X2] == POLAR_BNDRY
   || block_bcs[OUTER_X2] == POLAR_BNDRY_WEDGE)
    delete [] polar_neighbor_south;
  if (pmy_mesh_->multilevel==true) {
    sarea_[0].DeleteAthenaArray();
    sarea_[1].DeleteAthenaArray();
  }
}


//----------------------------------------------------------------------------------------
//! \fn unsigned int BoundaryBase::CreateBufferID(int ox1, int ox2, int ox3,
//                                                       int fi1, int fi2)
//  \brief calculate a buffer identifier

unsigned int BoundaryBase::CreateBufferID(int ox1, int ox2, int ox3,
                                                 int fi1, int fi2) {
  unsigned int ux1=(unsigned)(ox1+1);
  unsigned int ux2=(unsigned)(ox2+1);
  unsigned int ux3=(unsigned)(ox3+1);
  return (ux1<<6) | (ux2<<4) | (ux3<<2) | (fi1<<1) | fi2;
}


//----------------------------------------------------------------------------------------
//! \fn int BoundaryBase::BufferID(int dim, bool multilevel)
//  \brief calculate neighbor indexes and target buffer IDs

int BoundaryBase::BufferID(int dim, bool multilevel) {
  int nf1=1, nf2=1;
  if (multilevel==true) {
    if (dim>=2) nf1=2;
    if (dim>=3) nf2=2;
  }
  int b=0;
  // x1 face
  for (int n=-1; n<=1; n+=2) {
    for (int f2=0;f2<nf2;f2++) {
      for (int f1=0;f1<nf1;f1++) {
        ni[b].ox1=n; ni[b].ox2=0; ni[b].ox3=0;
        ni[b].fi1=f1; ni[b].fi2=f2;
        ni[b].type=NEIGHBOR_FACE;
        b++;
      }
    }
  }
  // x2 face
  if (dim>=2) {
    for (int n=-1; n<=1; n+=2) {
      for (int f2=0;f2<nf2;f2++) {
        for (int f1=0;f1<nf1;f1++) {
          ni[b].ox1=0; ni[b].ox2=n; ni[b].ox3=0;
          ni[b].fi1=f1; ni[b].fi2=f2;
          ni[b].type=NEIGHBOR_FACE;
          b++;
        }
      }
    }
  }
  if (dim==3) {
    // x3 face
    for (int n=-1; n<=1; n+=2) {
      for (int f2=0;f2<nf2;f2++) {
        for (int f1=0;f1<nf1;f1++) {
          ni[b].ox1=0; ni[b].ox2=0; ni[b].ox3=n;
          ni[b].fi1=f1; ni[b].fi2=f2;
          ni[b].type=NEIGHBOR_FACE;
          b++;
        }
      }
    }
  }
  // edges
  // x1x2
  if (dim>=2) {
    for (int m=-1; m<=1; m+=2) {
      for (int n=-1; n<=1; n+=2) {
        for (int f1=0;f1<nf2;f1++) {
          ni[b].ox1=n; ni[b].ox2=m; ni[b].ox3=0;
          ni[b].fi1=f1; ni[b].fi2=0;
          ni[b].type=NEIGHBOR_EDGE;
          b++;
        }
      }
    }
  }
  if (dim==3) {
    // x1x3
    for (int m=-1; m<=1; m+=2) {
      for (int n=-1; n<=1; n+=2) {
        for (int f1=0;f1<nf1;f1++) {
          ni[b].ox1=n; ni[b].ox2=0; ni[b].ox3=m;
          ni[b].fi1=f1; ni[b].fi2=0;
          ni[b].type=NEIGHBOR_EDGE;
          b++;
        }
      }
    }
    // x2x3
    for (int m=-1; m<=1; m+=2) {
      for (int n=-1; n<=1; n+=2) {
        for (int f1=0;f1<nf1;f1++) {
          ni[b].ox1=0; ni[b].ox2=n; ni[b].ox3=m;
          ni[b].fi1=f1; ni[b].fi2=0;
          ni[b].type=NEIGHBOR_EDGE;
          b++;
        }
      }
    }
    // corners
    for (int l=-1; l<=1; l+=2) {
      for (int m=-1; m<=1; m+=2) {
        for (int n=-1; n<=1; n+=2) {
          ni[b].ox1=n; ni[b].ox2=m; ni[b].ox3=l;
          ni[b].fi1=0; ni[b].fi2=0;
          ni[b].type=NEIGHBOR_CORNER;
          b++;
        }
      }
    }
  }

  for (int n=0;n<b;n++)
    bufid[n]=CreateBufferID(ni[n].ox1, ni[n].ox2, ni[n].ox3, ni[n].fi1, ni[n].fi2);

  return b;
}


//----------------------------------------------------------------------------------------
//! \fn int BoundaryBase::FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2)
//  \brief find the boundary buffer ID from the direction

int BoundaryBase::FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2) {
  int bid=CreateBufferID(ox1, ox2, ox3, fi1, fi2);

  for (int i=0;i<maxneighbor_;i++) {
    if (bid==bufid[i]) return i;
  }
  return -1;
}


//----------------------------------------------------------------------------------------
//! \fn unsigned int BoundaryBase::CreateBvalsMPITag(int lid, int phys, int bufid)
//  \brief calculate an MPI tag for Bval communications
// tag = local id of destination (20) + bufid(6) + physics(5)

unsigned int BoundaryBase::CreateBvalsMPITag(int lid, int phys, int bufid) {
  return (lid<<11) | (bufid<<5) | phys;
}


//----------------------------------------------------------------------------------------
// \!fn void BoundaryBase::SearchAndSetNeighbors(MeshBlockTree &tree,
//                                               int *ranklist, int *nslist)
// \brief Search and set all the neighbor blocks

void BoundaryBase::SearchAndSetNeighbors(MeshBlockTree &tree, int *ranklist,
                                         int *nslist) {
  MeshBlockTree* neibt;
  int myox1, myox2=0, myox3=0, myfx1, myfx2, myfx3;
  myfx1=static_cast<int>(loc.lx1&1L);
  myfx2=static_cast<int>(loc.lx2&1L);
  myfx3=static_cast<int>(loc.lx3&1L);
  myox1=(static_cast<int>(loc.lx1&1L))*2-1;
  if (block_size_.nx2>1) myox2=(static_cast<int>(loc.lx2&1L))*2-1;
  if (block_size_.nx3>1) myox3=(static_cast<int>(loc.lx3&1L))*2-1;
  int64_t nrbx1=pmy_mesh_->nrbx1, nrbx2=pmy_mesh_->nrbx2, nrbx3=pmy_mesh_->nrbx3;

  int nf1=1, nf2=1;
  if (pmy_mesh_->multilevel==true) {
    if (block_size_.nx2>1) nf1=2;
    if (block_size_.nx3>1) nf2=2;
  }
  int bufid=0;
  nneighbor=0;
  for (int k=0; k<=2; k++) {
    for (int j=0; j<=2; j++) {
      for (int i=0; i<=2; i++)
        nblevel[k][j][i]=-1;
    }
  }
  nblevel[1][1][1]=loc.level;

  // x1 face
  for (int n=-1; n<=1; n+=2) {
    neibt=tree.FindNeighbor(loc,n,0,0,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh_->root_level);
    if (neibt==NULL) { bufid+=nf1*nf2; continue;}
    if (neibt->flag==false) { // neighbor at finer level
      int fface=1-(n+1)/2; // 0 for OUTER_X1, 1 for INNER_X1
      nblevel[1][1][n+1]=neibt->loc.level+1;
      for (int f2=0;f2<nf2;f2++) {
        for (int f1=0;f1<nf1;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(fface,f1,f2);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(-n,0,0,0,0);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
                                          fid-nslist[ranklist[fid]], n, 0, 0,
                                          NEIGHBOR_FACE, bufid, tbid, false,
                                          false, f1, f2);
          bufid++; nneighbor++;
        }
      }
    } else { // neighbor at same or coarser level
      int nlevel=neibt->loc.level;
      int nid=neibt->gid;
      nblevel[1][1][n+1]=nlevel;
      int tbid;
      bool shear=false;
      if (nlevel==loc.level) { // neighbor at same level
        tbid=FindBufferID(-n,0,0,0,0);
        if ((n == -1 and block_bcs[INNER_X1] == SHEAR_PERIODIC_BNDRY)
            or (n == 1 and block_bcs[OUTER_X1] == SHEAR_PERIODIC_BNDRY)) {
          shear = true; // neighbor is shearing periodic
        }
      } else { // neighbor at coarser level
        tbid=FindBufferID(-n,0,0,myfx2,myfx3);
      }
      neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
          nid-nslist[ranklist[nid]], n, 0, 0, NEIGHBOR_FACE, bufid, tbid, false, shear);
      bufid+=nf1*nf2; nneighbor++;
    }
  }
  if (block_size_.nx2==1) return;

  // x2 face
  for (int n=-1; n<=1; n+=2) {
    neibt=tree.FindNeighbor(loc,0,n,0,block_bcs,nrbx1,nrbx2,nrbx3,pmy_mesh_->root_level);
    if (neibt==NULL) { bufid+=nf1*nf2; continue;}
    if (neibt->flag==false) { // neighbor at finer level
      int fface=1-(n+1)/2; // 0 for OUTER_X2, 1 for INNER_X2
      nblevel[1][n+1][1]=neibt->loc.level+1;
      for (int f2=0;f2<nf2;f2++) {
        for (int f1=0;f1<nf1;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(f1,fface,f2);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(0,-n,0,0,0);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
                                          fid-nslist[ranklist[fid]], 0, n, 0,
                                          NEIGHBOR_FACE, bufid, tbid, false, false,
                                          f1, f2);
          bufid++; nneighbor++;
        }
      }
    } else { // neighbor at same or coarser level
      int nlevel=neibt->loc.level;
      int nid=neibt->gid;
      nblevel[1][n+1][1]=nlevel;
      int tbid;
      bool polar=false;
      if (nlevel==loc.level) { // neighbor at same level
        if ((n == -1 and block_bcs[INNER_X2] == POLAR_BNDRY)
            or (n == 1 and block_bcs[OUTER_X2] == POLAR_BNDRY)) {
          polar = true; // neighbor is across top or bottom pole
        }
        tbid=FindBufferID(0,polar?n:-n,0,0,0);
      } else { // neighbor at coarser level
        tbid=FindBufferID(0,-n,0,myfx1,myfx3);
      }
      neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
          nid-nslist[ranklist[nid]], 0, n, 0, NEIGHBOR_FACE, bufid, tbid, polar, false);
      bufid+=nf1*nf2; nneighbor++;
    }
  }

  // x3 face
  if (block_size_.nx3>1) {
    for (int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(loc, 0, 0, n, block_bcs, nrbx1, nrbx2, nrbx3,
                              pmy_mesh_->root_level);
      if (neibt==NULL) { bufid+=nf1*nf2; continue;}
      if (neibt->flag==false) { // neighbor at finer level
        int fface=1-(n+1)/2; // 0 for OUTER_X3, 1 for INNER_X3
        nblevel[n+1][1][1]=neibt->loc.level+1;
        for (int f2=0;f2<nf2;f2++) {
          for (int f1=0;f1<nf1;f1++) {
            MeshBlockTree* nf=neibt->GetLeaf(f1,f2,fface);
            int fid = nf->gid;
            int nlevel=nf->loc.level;
            int tbid=FindBufferID(0,0,-n,0,0);
            neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
                                            fid-nslist[ranklist[fid]], 0, 0, n,
                                            NEIGHBOR_FACE, bufid, tbid,
                                            false, false, f1, f2);
            bufid++; nneighbor++;
          }
        }
      } else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[n+1][1][1]=nlevel;
        int tbid;
        if (nlevel==loc.level) { // neighbor at same level
          tbid=FindBufferID(0,0,-n,0,0);
        } else { // neighbor at coarser level
          tbid=FindBufferID(0,0,-n,myfx1,myfx2);
        }
        neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
            nid-nslist[ranklist[nid]], 0, 0, n, NEIGHBOR_FACE, bufid, tbid, false, false);
        bufid+=nf1*nf2; nneighbor++;
      }
    }
  }

  // x1x2 edge
  for (int m=-1; m<=1; m+=2) {
    for (int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(loc, n, m, 0, block_bcs, nrbx1, nrbx2, nrbx3,
                              pmy_mesh_->root_level);
      if (neibt==NULL) { bufid+=nf2; continue;}
      bool polar=false;
      if ((m == -1 and block_bcs[INNER_X2] == POLAR_BNDRY)
          or (m == 1 and block_bcs[OUTER_X2] == POLAR_BNDRY)) {
        polar = true; // neighbor is across top or bottom pole
      }
      if (neibt->flag==false) { // neighbor at finer level
        int ff1=1-(n+1)/2; // 0 for OUTER_X1, 1 for INNER_X1
        int ff2=1-(m+1)/2; // 0 for OUTER_X2, 1 for INNER_X2
        if (polar) {
          ff2 = 1 - ff2;
        }
        nblevel[1][m+1][n+1]=neibt->loc.level+1;
        for (int f1=0;f1<nf2;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(ff1,ff2,f1);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(-n,polar?m:-m,0,0,0);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
                                          fid-nslist[ranklist[fid]], n, m, 0,
                                          NEIGHBOR_EDGE, bufid, tbid, polar,
                                          false, f1, 0);
          bufid++; nneighbor++;
        }
      } else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[1][m+1][n+1]=nlevel;
        int tbid;
        bool shear=false;
        if (nlevel==loc.level) { // neighbor at same level
          if ((n == -1 and block_bcs[INNER_X1] == SHEAR_PERIODIC_BNDRY)
              or (n == 1 and block_bcs[OUTER_X1] == SHEAR_PERIODIC_BNDRY)) {
            shear = true; // neighbor is on shearing periodic bcs
          }
          tbid=FindBufferID(-n,polar?m:-m,0,0,0);
        } else { // neighbor at coarser level
          tbid=FindBufferID(-n,polar?m:-m,0,myfx3,0);
        }
        if (nlevel>=loc.level || (myox1==n && myox2==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
                                          nid-nslist[ranklist[nid]], n, m, 0,
                                          NEIGHBOR_EDGE, bufid, tbid, polar, shear);
          nneighbor++;
        }
        bufid+=nf2;
      }
    }
  }

  // polar neighbors
  if (block_bcs[INNER_X2] == POLAR_BNDRY||block_bcs[INNER_X2] == POLAR_BNDRY_WEDGE) {
    int level = loc.level - pmy_mesh_->root_level;
    int num_north_polar_blocks = static_cast<int>(nrbx3 * (1 << level));
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
    int level = loc.level - pmy_mesh_->root_level;
    int num_south_polar_blocks = static_cast<int>(nrbx3 * (1 << level));
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
  if (block_size_.nx3==1) return;

  // x1x3 edge
  for (int m=-1; m<=1; m+=2) {
    for (int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(loc, n, 0, m, block_bcs, nrbx1, nrbx2, nrbx3,
                              pmy_mesh_->root_level);
      if (neibt==NULL) { bufid+=nf1; continue;}
      if (neibt->flag==false) { // neighbor at finer level
        int ff1=1-(n+1)/2; // 0 for OUTER_X1, 1 for INNER_X1
        int ff2=1-(m+1)/2; // 0 for OUTER_X3, 1 for INNER_X3
        nblevel[m+1][1][n+1]=neibt->loc.level+1;
        for (int f1=0;f1<nf1;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(ff1,f1,ff2);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(-n,0,-m,0,0);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
                                          fid-nslist[ranklist[fid]], n, 0, m,
                                          NEIGHBOR_EDGE, bufid, tbid,
                                          false, false, f1, 0);
          bufid++; nneighbor++;
        }
      } else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[m+1][1][n+1]=nlevel;
        int tbid;
        bool shear=false;
        if (nlevel==loc.level) { // neighbor at same level
          tbid=FindBufferID(-n,0,-m,0,0);
          if ((n == -1 and block_bcs[INNER_X1] == SHEAR_PERIODIC_BNDRY)
              or (n == 1 and block_bcs[OUTER_X1] == SHEAR_PERIODIC_BNDRY)) {
            shear = true; //neighbor is on shearing periodic boundary
          }
        } else { // neighbor at coarser level
          tbid=FindBufferID(-n,0,-m,myfx2,0);
        }
        if (nlevel>=loc.level || (myox1==n && myox3==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
                                          nid-nslist[ranklist[nid]], n, 0, m,
                                          NEIGHBOR_EDGE, bufid, tbid, false, shear);
          nneighbor++;
        }
        bufid+=nf1;
      }
    }
  }

  // x2x3 edge
  for (int m=-1; m<=1; m+=2) {
    for (int n=-1; n<=1; n+=2) {
      neibt=tree.FindNeighbor(loc, 0, n, m, block_bcs, nrbx1, nrbx2, nrbx3,
                              pmy_mesh_->root_level);
      if (neibt==NULL) { bufid+=nf1; continue;}
      if (neibt->flag==false) { // neighbor at finer level
        int ff1=1-(n+1)/2; // 0 for OUTER_X2, 1 for INNER_X2
        int ff2=1-(m+1)/2; // 0 for OUTER_X3, 1 for INNER_X3
        nblevel[m+1][n+1][1]=neibt->loc.level+1;
        for (int f1=0;f1<nf1;f1++) {
          MeshBlockTree* nf=neibt->GetLeaf(f1,ff1,ff2);
          int fid = nf->gid;
          int nlevel=nf->loc.level;
          int tbid=FindBufferID(0,-n,-m,0,0);
          neighbor[nneighbor].SetNeighbor(ranklist[fid], nlevel, fid,
                                          fid-nslist[ranklist[fid]], 0, n, m,
                                          NEIGHBOR_EDGE, bufid, tbid,
                                          false, false, f1, 0);
          bufid++; nneighbor++;
        }
      } else { // neighbor at same or coarser level
        int nlevel=neibt->loc.level;
        int nid=neibt->gid;
        nblevel[m+1][n+1][1]=nlevel;
        int tbid;
        bool polar=false;
        if (nlevel==loc.level) { // neighbor at same level
          if ((n == -1 and block_bcs[INNER_X2] == POLAR_BNDRY)
              or (n == 1 and block_bcs[OUTER_X2] == POLAR_BNDRY)) {
            polar = true; // neighbor is across top or bottom pole
          }
          tbid=FindBufferID(0,polar?n:-n,-m,0,0);
        } else { // neighbor at coarser level
          tbid=FindBufferID(0,-n,-m,myfx1,0);
        }
        if (nlevel>=loc.level || (myox2==n && myox3==m)) {
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
                                          nid-nslist[ranklist[nid]], 0, n, m,
                                          NEIGHBOR_EDGE, bufid, tbid, polar, false);
          nneighbor++;
        }
        bufid+=nf1;
      }
    }
  }

  // corners
  for (int l=-1; l<=1; l+=2) {
    for (int m=-1; m<=1; m+=2) {
      for (int n=-1; n<=1; n+=2) {
        neibt=tree.FindNeighbor(loc, n, m, l, block_bcs, nrbx1, nrbx2, nrbx3,
                                pmy_mesh_->root_level);
        if (neibt==NULL) { bufid++; continue;}
        bool polar=false;
        if ((m == -1 and block_bcs[INNER_X2] == POLAR_BNDRY)
            or (m == 1 and block_bcs[OUTER_X2] == POLAR_BNDRY)) {
          polar = true; // neighbor is across top or bottom pole
        }
        if (neibt->flag==false) { // neighbor at finer level
          int ff1=1-(n+1)/2; // 0 for OUTER_X1, 1 for INNER_X1
          int ff2=1-(m+1)/2; // 0 for OUTER_X2, 1 for INNER_X2
          int ff3=1-(l+1)/2; // 0 for OUTER_X3, 1 for INNER_X3
          if (polar) {
            ff2 = 1 - ff2;
          }
          neibt=neibt->GetLeaf(ff1,ff2,ff3);
        }
        int nlevel=neibt->loc.level;
        nblevel[l+1][m+1][n+1]=nlevel;
        if (nlevel>=loc.level || (myox1==n && myox2==m && myox3==l)) {
          int nid=neibt->gid;
          int tbid=FindBufferID(-n,polar?m:-m,-l,0,0);
          neighbor[nneighbor].SetNeighbor(ranklist[nid], nlevel, nid,
                                          nid-nslist[ranklist[nid]], n, m, l,
                                          NEIGHBOR_CORNER, bufid, tbid, polar, false);
          nneighbor++;
        }
        bufid++;
      }
    }
  }

  return;
}
