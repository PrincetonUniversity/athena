//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file meshblock_tree.cpp
//! \brief implementation of functions in the MeshBlockTree class
//!
//! The MeshBlockTree stores the logical grid structure, and is used for neighbor
//! searches, assigning global IDs, etc.  Level is defined as "logical level", where the
//! logical root (single block) level is 0.  Note the logical level of the physical root
//! grid (user-specified root grid) will be greater than zero if it contains more than
//! one MeshBlock

// C headers

// C++ headers
#include <cstdint>    // int64_t
#include <iostream>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "meshblock_tree.hpp"

// Define static member variables
Mesh* MeshBlockTree::pmesh_;
MeshBlockTree* MeshBlockTree::proot_;
int MeshBlockTree::nleaf_;


//----------------------------------------------------------------------------------------
//! \fn bool operator==(const LogicalLocation &l1, const LogicalLocation &l2)
//! \brief overloading the comparison operator for LogicalLocation
bool operator==(const LogicalLocation &l1, const LogicalLocation &l2) {
  return ((l1.level == l2.level) && (l1.lx1 == l2.lx1)
       && (l1.lx2 == l2.lx2) && (l1.lx3 == l2.lx3));
}


//----------------------------------------------------------------------------------------
//! \fn MeshBlockTree::MeshBlockTree(Mesh* pmesh)
//! \brief constructor for the logical root

MeshBlockTree::MeshBlockTree(Mesh* pmesh) : pleaf_(nullptr), gid_(-1) {
  pmesh_ = pmesh;
  proot_ = this;
  loc_.lx1 = 0;
  loc_.lx2 = 0;
  loc_.lx3 = 0;
  loc_.level = 0;
}

//----------------------------------------------------------------------------------------
//! \fn MeshBlockTree::MeshBlockTree(MeshBlockTree *parent, int ox1, int ox2, int ox3)
//! \brief constructor for a leaf

MeshBlockTree::MeshBlockTree(MeshBlockTree *parent, int ox1, int ox2, int ox3)
                           : pleaf_(nullptr), gid_(parent->gid_) {
  loc_.lx1 = (parent->loc_.lx1<<1)+ox1;
  loc_.lx2 = (parent->loc_.lx2<<1)+ox2;
  loc_.lx3 = (parent->loc_.lx3<<1)+ox3;
  loc_.level = parent->loc_.level+1;
}


//----------------------------------------------------------------------------------------
//! \fn MeshBlockTree::~MeshBlockTree()
//! \brief destructor (for both root and leaves)

MeshBlockTree::~MeshBlockTree() {
  if (pleaf_ != nullptr) {
    for (int i=0; i<nleaf_; i++)
      delete pleaf_[i];
    delete [] pleaf_;
  }
}


//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::CreateRootGrid()
//! \brief create the root grid; the root grid can be incomplete (less than 8 leaves)

void MeshBlockTree::CreateRootGrid() {
  if (loc_.level == 0) {
    nleaf_ = 2;
    if (pmesh_->f2) nleaf_ = 4;
    if (pmesh_->f3) nleaf_ = 8;
  }
  if (loc_.level == pmesh_->root_level) return;

  pleaf_ = new MeshBlockTree*[nleaf_];
  for (int n=0; n<nleaf_; n++)
    pleaf_[n] = nullptr;

  std::int64_t levfac = 1LL<<(pmesh_->root_level - loc_.level-1);
  for (int n=0; n<nleaf_; n++) {
    int i = n&1, j = (n>>1)&1, k = (n>>2)&1;
    if ((loc_.lx3*2 + k)*levfac < pmesh_->nrbx3
     && (loc_.lx2*2 + j)*levfac < pmesh_->nrbx2
     && (loc_.lx1*2 + i)*levfac < pmesh_->nrbx1) {
      pleaf_[n] = new MeshBlockTree(this, i, j, k);
      pleaf_[n]->CreateRootGrid();
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::AddMeshBlock(LogicalLocation rloc, int &nnew)
//! \brief add a MeshBlock to the tree, also creates neighboring blocks

void MeshBlockTree::AddMeshBlock(LogicalLocation rloc, int &nnew) {
  if (loc_.level == rloc.level) return; // done

  if (pleaf_ == nullptr) // leaf -> create the finer level
    Refine(nnew);

  // get leaf index
  int sh = rloc.level-loc_.level-1;
  int mx, my, mz;
  mx = ((rloc.lx1>>sh) & 1LL) == 1LL;
  my = ((rloc.lx2>>sh) & 1LL) == 1LL;
  mz = ((rloc.lx3>>sh) & 1LL) == 1LL;
  int n = mx + (my<<1) + (mz<<2);
  pleaf_[n]->AddMeshBlock(rloc, nnew);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::AddMeshBlockWithoutRefine(LogicalLocation rloc)
//! \brief add a MeshBlock to the tree without refinement, used in restarting.
//!        MeshBlockTree::CreateRootGrid must be called before this method
void MeshBlockTree::AddMeshBlockWithoutRefine(LogicalLocation rloc) {
  if (loc_.level == rloc.level) // done
    return;

  if (pleaf_ == nullptr) {
    pleaf_ = new MeshBlockTree*[nleaf_];
    for (int n=0; n<nleaf_; n++)
      pleaf_[n] = nullptr;
  }

  // get leaf index
  int sh = rloc.level-loc_.level-1;
  int mx, my, mz;
  mx = ((rloc.lx1>>sh) & 1LL) == 1LL;
  my = ((rloc.lx2>>sh) & 1LL) == 1LL;
  mz = ((rloc.lx3>>sh) & 1LL) == 1LL;
  int n = mx + (my<<1) + (mz<<2);
  if (pleaf_[n] == nullptr)
    pleaf_[n] = new MeshBlockTree(this, mx, my, mz);
  pleaf_[n]->AddMeshBlockWithoutRefine(rloc);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::Refine(int &nnew)
//! \brief make finer leaves

void MeshBlockTree::Refine(int &nnew) {
  if (pleaf_ != nullptr) return;

  pleaf_ = new MeshBlockTree*[nleaf_];
  for (int n=0; n<nleaf_; n++)
    pleaf_[n] = nullptr;

  for (int n=0; n<nleaf_; n++) {
    int i = n&1, j = (n>>1)&1, k = (n>>2)&1;
    pleaf_[n] = new MeshBlockTree(this, i, j, k);
  }

  std::int64_t nxmax,nymax,nzmax;
  std::int64_t ox, oy, oz, oxmin, oxmax, oymin, oymax, ozmin, ozmax;
  LogicalLocation nloc;
  nloc.level=loc_.level;

  oxmin=-1, oxmax=1, nxmax=(pmesh_->nrbx1<<(loc_.level-pmesh_->root_level));
  if (pmesh_->f2)
    oymin=-1, oymax=1, nymax=(pmesh_->nrbx2<<(loc_.level-pmesh_->root_level));
  else
    oymin=0,  oymax=0, nymax=1;
  if (pmesh_->f3) // 3D
    ozmin=-1, ozmax=1, nzmax=(pmesh_->nrbx3<<(loc_.level-pmesh_->root_level));
  else
    ozmin=0,  ozmax=0, nzmax=1;

  for (oz=ozmin; oz<=ozmax; oz++) {
    nloc.lx3=loc_.lx3+oz;
    if (nloc.lx3<0) {
      if (pmesh_->mesh_bcs[BoundaryFace::inner_x3]!=BoundaryFlag::periodic)
        continue;
      else
        nloc.lx3=nzmax-1;
    }
    if (nloc.lx3>=nzmax) {
      if (pmesh_->mesh_bcs[BoundaryFace::outer_x3]!=BoundaryFlag::periodic)
        continue;
      else
        nloc.lx3=0;
    }
    for (oy=oymin; oy<=oymax; oy++) {
      nloc.lx2=loc_.lx2+oy;
      bool polar=false;
      if (nloc.lx2<0) {
        if (pmesh_->mesh_bcs[BoundaryFace::inner_x2]==BoundaryFlag::periodic) {
          nloc.lx2=nymax-1;
        } else if (pmesh_->mesh_bcs[BoundaryFace::inner_x2]==BoundaryFlag::polar) {
          nloc.lx2=0;
          polar=true;
        } else {
          continue;
        }
      }
      if (nloc.lx2>=nymax) {
        if (pmesh_->mesh_bcs[BoundaryFace::outer_x2]==BoundaryFlag::periodic) {
          nloc.lx2=0;
        } else if (pmesh_->mesh_bcs[BoundaryFace::outer_x2]==BoundaryFlag::polar) {
          nloc.lx2=nymax-1;
          polar=true;
        } else {
          continue;
        }
      }
      if (polar) nloc.lx3=(nloc.lx3+nzmax/2)%nzmax;
      for (ox=oxmin; ox<=oxmax; ox++) {
        if (ox==0 && oy==0 && oz==0) continue;
        nloc.lx1=loc_.lx1+ox;
        if (nloc.lx1<0) {
          if (pmesh_->mesh_bcs[BoundaryFace::inner_x1]!=BoundaryFlag::periodic)
            continue;
          else
            nloc.lx1=nxmax-1;
        }
        if (nloc.lx1>=nxmax) {
          if (pmesh_->mesh_bcs[BoundaryFace::outer_x1]!=BoundaryFlag::periodic)
            continue;
          else
            nloc.lx1=0;
        }
        proot_->AddMeshBlock(nloc, nnew);
      }
    }
  }
  // this block is not a leaf anymore
  gid_=-1;

  nnew+=nleaf_-1;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::Derefine(int &ndel)
//! \brief destroy leaves and make this block a leaf

void MeshBlockTree::Derefine(int &ndel) {
  int s2=0, e2=0, s3=0, e3=0;
  if (pmesh_->f2) s2=-1, e2=1;
  if (pmesh_->f3) s3=-1, e3=1;
  for (int ox3=s3; ox3<=e3; ox3++) {
    for (int ox2=s2; ox2<=e2; ox2++) {
      for (int ox1=-1; ox1<=1; ox1++) {
        MeshBlockTree *bt = proot_->FindNeighbor(loc_, ox1, ox2, ox3,
                                                 pmesh_->mesh_bcs, true);
        if (bt != nullptr) {
          if (bt->pleaf_ != nullptr) {
            int lis, lie, ljs, lje, lks, lke;
            if (ox1==-1)       lis=lie=1;
            else if (ox1==1)   lis=lie=0;
            else              lis=0, lie=1;
            if (pmesh_->f2) {
              if (ox2==-1)     ljs=lje=1;
              else if (ox2==1) ljs=lje=0;
              else            ljs=0, lje=1;
            } else {
              ljs=lje=0;
            }
            if (pmesh_->f3) {
              if (ox3==-1)     lks=lke=1;
              else if (ox3==1) lks=lke=0;
              else            lks=0, lke=1;
            } else {
              lks=lke=0;
            }
            for (int lk=lks; lk<=lke; lk++) {
              for (int lj=ljs; lj<=lje; lj++) {
                for (int li=lis; li<=lie; li++) {
                  int n = li + (lj<<1) + (lk<<2);
                  if (bt->pleaf_[n]->pleaf_ != nullptr) return;
                }
              }
            }
          }
        }
      }
    }
  }

  gid_ = pleaf_[0]->gid_; // now this is a leaf; inherit the first leaf's GID
  for (int n=0; n<nleaf_; n++)
    delete pleaf_[n];
  delete [] pleaf_;
  pleaf_ = nullptr;
  ndel+=nleaf_-1;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::CountMeshBlock(int& count)
//! \brief creates the Location list sorted by Z-ordering

void MeshBlockTree::CountMeshBlock(int& count) {
  if (loc_.level == 0) count=0;

  if (pleaf_ == nullptr) {
    count++;
  } else {
    for (int n=0; n<nleaf_; n++) {
      if (pleaf_[n] != nullptr)
        pleaf_[n]->CountMeshBlock(count);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::GetMeshBlockList(LogicalLocation *list,
//!                                          int *pglist, int& count)
//! \brief creates the Location list sorted by Z-ordering

void MeshBlockTree::GetMeshBlockList(LogicalLocation *list, int *pglist, int& count) {
  if (loc_.level == 0) count=0;

  if (pleaf_ == nullptr) {
    list[count]=loc_;
    if (pglist != nullptr)
      pglist[count]=gid_;
    gid_=count;
    count++;
  } else {
    for (int n=0; n<nleaf_; n++) {
      if (pleaf_[n] != nullptr)
        pleaf_[n]->GetMeshBlockList(list, pglist, count);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn MeshBlockTree* MeshBlockTree::FindNeighbor(LogicalLocation myloc,
//!                    int ox1, int ox2, int ox3, BoundaryFlag *bcs, bool amrflag)
//! \brief find a neighboring block, called from the root of the tree
//!        If it is coarser or same level, return the pointer to that block.
//!        If it is a finer block, return the pointer to its parent.
//!        Note that this function must be called on a completed tree only

MeshBlockTree* MeshBlockTree::FindNeighbor(LogicalLocation myloc,
               int ox1, int ox2, int ox3, BoundaryFlag *bcs, bool amrflag) {
  std::stringstream msg;
  std::int64_t lx, ly, lz;
  int ll;
  int ox, oy, oz;
  MeshBlockTree *bt = proot_;
  lx=myloc.lx1, ly=myloc.lx2, lz=myloc.lx3, ll=myloc.level;

  lx+=ox1; ly+=ox2; lz+=ox3;
  // periodic and polar boundaries
  if (lx<0) {
    if (bcs[BoundaryFace::inner_x1] == BoundaryFlag::periodic
        || bcs[BoundaryFace::inner_x1] == BoundaryFlag::shear_periodic)
      lx=(pmesh_->nrbx1<<(ll-pmesh_->root_level))-1;
    else
      return nullptr;
  }
  if (lx>=pmesh_->nrbx1<<(ll-pmesh_->root_level)) {
    if (bcs[BoundaryFace::outer_x1] == BoundaryFlag::periodic
        || bcs[BoundaryFace::outer_x1] == BoundaryFlag::shear_periodic)
      lx=0;
    else
      return nullptr;
  }
  bool polar = false;
  if (ly<0) {
    if (bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic) {
      ly=(pmesh_->nrbx2<<(ll-pmesh_->root_level))-1;
    } else if (bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) {
      ly=0;
      polar=true;
    } else {
      return nullptr;
    }
  }
  if (ly>=pmesh_->nrbx2<<(ll-pmesh_->root_level)) {
    if (bcs[BoundaryFace::outer_x2] == BoundaryFlag::periodic) {
      ly=0;
    } else if (bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) {
      ly=(pmesh_->nrbx2<<(ll-pmesh_->root_level))-1;
      polar=true;
    } else {
      return nullptr;
    }
  }
  std::int64_t num_x3 = pmesh_->nrbx3<<(ll-pmesh_->root_level);
  if (lz<0) {
    if (bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic)
      lz=num_x3-1;
    else
      return nullptr;
  }
  if (lz>=num_x3) {
    if (bcs[BoundaryFace::outer_x3] == BoundaryFlag::periodic)
      lz=0;
    else
      return nullptr;
  }
  if (ll<1) return proot_; // single grid; return root
  if (polar) lz=(lz+num_x3/2)%num_x3;

  for (int level=0; level<ll; level++) {
    if (bt->pleaf_ == nullptr) { // leaf
      if (level == ll-1) {
        return bt;
      } else {
        msg << "### FATAL ERROR in FindNeighbor" << std::endl
            << "Neighbor search failed. The Block Tree is broken." << std::endl;
        ATHENA_ERROR(msg);
        return nullptr;
      }
    }
    // find a leaf in the next level
    int sh=ll-level-1;
    ox = ((lx>>sh) & 1LL) == 1LL;
    oy = ((ly>>sh) & 1LL) == 1LL;
    oz = ((lz>>sh) & 1LL) == 1LL;
    bt=bt->GetLeaf(ox, oy, oz);
    if (bt == nullptr) {
      msg << "### FATAL ERROR in FindNeighbor" << std::endl
          << "Neighbor search failed. The Block Tree is broken." << std::endl;
      ATHENA_ERROR(msg);
      return nullptr;
    }
  }
  if (bt->pleaf_ == nullptr) // leaf on the same level
    return bt;
  // one level finer: check if it is a leaf
  ox = oy = oz = 0;
  if (ox1 < 0) ox = 1;
  if (ox2 < 0) oy = 1;
  if (ox3 < 0) oz = 1;
  MeshBlockTree *btleaf = bt->GetLeaf(ox, oy, oz);
  if (btleaf->pleaf_ == nullptr)
    return bt;  // return this block
  if (!amrflag) {
    msg << "### FATAL ERROR in FindNeighbor" << std::endl
        << "Neighbor search failed. The Block Tree is broken." << std::endl;
    ATHENA_ERROR(msg);
  }
  return bt;
}

//----------------------------------------------------------------------------------------
//! \fn MeshBlockTree* MeshBlockTree::FindMeshBlock(LogicalLocation tloc)
//! \brief find MeshBlock with LogicalLocation tloc and return a pointer

MeshBlockTree* MeshBlockTree::FindMeshBlock(LogicalLocation tloc) {
  if (tloc.level == loc_.level) return this;
  if (pleaf_ == nullptr) return nullptr;
  // get leaf index
  int sh = tloc.level - loc_.level - 1;
  int mx = (((tloc.lx1>>sh) & 1LL) == 1LL);
  int my = (((tloc.lx2>>sh) & 1LL) == 1LL);
  int mz = (((tloc.lx3>>sh) & 1LL) == 1LL);
  int n = mx + (my<<1) + (mz<<2);
  if (pleaf_[n] == nullptr)
    return nullptr;
  return pleaf_[n]->FindMeshBlock(tloc);
}


//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::CountMGOctets(int *noct)
//! \brief count the number of octets for Multigrid with mesh refinement

void MeshBlockTree::CountMGOctets(int *noct) {
  if (pleaf_ == nullptr) return;

  int lev = loc_.level - pmesh_->root_level;
  if (lev >= 0)
    noct[lev]++;
  for (int n=0; n<nleaf_; n++) {
    if (pleaf_[n] != nullptr)
      pleaf_[n]->CountMGOctets(noct);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::GetMGOctetList(std::vector<MGOctet> *oct,
//!     std::unordered_map<LogicalLocation, int, LogicalLocationHash> *octmap, int *noct)
//! \brief construct lists of octets Multigrid with mesh refinement

void MeshBlockTree::GetMGOctetList(std::vector<MGOctet> *oct,
     std::unordered_map<LogicalLocation, int, LogicalLocationHash> *octmap, int *noct) {
  if (pleaf_ == nullptr) return;

  int lev = loc_.level - pmesh_->root_level;
  int oid = 0;
  if (lev >= 0) {
    oid = noct[lev];
    oct[lev][oid].loc = loc_;
    oct[lev][oid].fleaf = true;
    octmap[lev][loc_] = oid;
    noct[lev]++;
  }
  for (int n=0; n<nleaf_; n++) {
    if (pleaf_[n] != nullptr) {
      if (pleaf_[n]->pleaf_ != nullptr) {
        if (lev >= 0) oct[lev][oid].fleaf = false;
        pleaf_[n]->GetMGOctetList(oct, octmap, noct);
      }
    }
  }

  return;
}
