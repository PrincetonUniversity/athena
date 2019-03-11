//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file meshblocktree.cpp
//  \brief implementation of functions in the MeshBlockTree class
// The MeshBlockTree stores the logical grid structure, and is used for neighbor
// searches, assigning global IDs, etc.  Level is defined as "logical level", where the
// logical root (single block) level is 0.  Note the logical level of the physical root
// grid (user-specified root grid) will be greater than zero if it contains more than
// one MeshBlock

// C headers

// C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "meshblock_tree.hpp"

//----------------------------------------------------------------------------------------
//! \fn MeshBlockTree::MeshBlockTree()
//  \brief constructor for the logical root

MeshBlockTree::MeshBlockTree() {
  flag=true;
  gid=-1;
  loc.lx1=0, loc.lx2=0, loc.lx3=0;
  loc.level=0;
  pparent=nullptr;
  for (int k=0; k<=1; k++) {
    for (int j=0; j<=1; j++) {
      for (int i=0; i<=1; i++) {
        pleaf[k][j][i]=nullptr;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn MeshBlockTree::MeshBlockTree(MeshBlockTree *parent, int ox, int oy, int oz)
//  \brief constructor for a leaf

MeshBlockTree::MeshBlockTree(MeshBlockTree *parent, int ox, int oy, int oz) {
  flag=true;
  pparent=parent;
  gid=pparent->gid;
  loc.lx1=(parent->loc.lx1<<1)+ox;
  loc.lx2=(parent->loc.lx2<<1)+oy;
  loc.lx3=(parent->loc.lx3<<1)+oz;
  loc.level=parent->loc.level+1;
  for (int k=0; k<=1; k++) {
    for (int j=0; j<=1; j++) {
      for (int i=0; i<=1; i++) {
        pleaf[k][j][i]=nullptr;
      }
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn MeshBlockTree::~MeshBlockTree()
//  \brief destructor (for both root and leaves)

MeshBlockTree::~MeshBlockTree() {
  for (int k=0; k<=1; k++) {
    for (int j=0; j<=1; j++) {
      for (int i=0; i<=1; i++) {
        if (pleaf[k][j][i]!=nullptr)
          delete pleaf[k][j][i];
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::CreateRootGrid(std::int64_t nx,std::int64_t ny,
//                                         std::int64_t nz, int nl)
//  \brief create the root grid; the root grid can be incomplete (less than 8 leaves)

void MeshBlockTree::CreateRootGrid(std::int64_t nx, std::int64_t ny, std::int64_t nz,
                                   int nl) {
  if (loc.level == nl) return;

  for (int k=0; k<=1; k++) {
    if ((loc.lx3*2+k)*(1LL<<(nl-loc.level-1)) < nz) {
      for (int j=0; j<=1; j++) {
        if ((loc.lx2*2+j)*(1LL<<(nl-loc.level-1)) < ny) {
          for (int i=0; i<=1; i++) {
            if ((loc.lx1*2+i)*(1LL<<(nl-loc.level-1)) < nx) {
              flag=false; // if there is a leaf, this is a node
              gid=-1;
              pleaf[k][j][i] = new MeshBlockTree(this, i, j, k);
              pleaf[k][j][i]->CreateRootGrid(nx, ny, nz, nl);
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::AddMeshBlock(MeshBlockTree& root, LogicalLocation rloc,
//   int dim, enum BoundaryFlag* mesh_bcs, std::int64_t rbx, std::int64_t rby,
//   std::int64_t rbz, int rl, int &nnew)
//  \brief add a MeshBlock to the tree, also creates neighboring blocks

void MeshBlockTree::AddMeshBlock(MeshBlockTree& root, LogicalLocation rloc, int dim,
                                 enum BoundaryFlag* mesh_bcs, std::int64_t rbx,
                                 std::int64_t rby, std::int64_t rbz,
                                 int rl, int &nnew) {
  int mx, my, mz;
  if (loc.level==rloc.level) return; // done

  if (flag==true) // leaf -> create the finer level
    Refine(root,dim,mesh_bcs,rbx,rby,rbz,rl,nnew);
  // get leaf indexes
  int sh = rloc.level-loc.level-1;
  mx = (((rloc.lx1>>sh) & 1LL) == 1LL);
  my = (((rloc.lx2>>sh) & 1LL) == 1LL);
  mz = (((rloc.lx3>>sh) & 1LL) == 1LL);
  pleaf[mz][my][mx]->AddMeshBlock(root,rloc,dim,mesh_bcs,rbx,rby,rbz,rl,nnew);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::AddMeshBlockWithoutRefine(LogicalLocation rloc,
//                          std::int64_t rbx, std::int64_t rby, std::int64_t rbz, int rl)
//  \brief add a MeshBlock to the tree without refinement, used in restarting

void MeshBlockTree::AddMeshBlockWithoutRefine(LogicalLocation rloc,
                                              std::int64_t rbx, std::int64_t rby,
                                              std::int64_t rbz, int rl) {
  int mx, my, mz;
  if (loc.level==rloc.level) // done
    return;
  if (flag==true) // leaf -> create the finer level
    flag=false;
  // get leaf indexes
  int sh = rloc.level-loc.level-1;
  mx = (((rloc.lx1>>sh) & 1LL) == 1LL);
  my = (((rloc.lx2>>sh) & 1LL) == 1LL);
  mz = (((rloc.lx3>>sh) & 1LL) == 1LL);
  if (pleaf[mz][my][mx]==nullptr)
    pleaf[mz][my][mx] = new MeshBlockTree(this, mx, my, mz);
  pleaf[mz][my][mx]->AddMeshBlockWithoutRefine(rloc,rbx,rby,rbz,rl);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::Refine(MeshBlockTree& root, int dim,
//           enum BoundaryFlag* mesh_bcs, std::int64_t rbx, std::int64_t rby,
//           std::int64_t rbz, int rl, int &nnew)
//  \brief make finer leaves

void MeshBlockTree::Refine(MeshBlockTree& root, int dim, enum BoundaryFlag* mesh_bcs,
                           std::int64_t rbx, std::int64_t rby, std::int64_t rbz, int rl,
                           int &nnew) {
  if (flag==false) return;
  std::int64_t nxmax,nymax,nzmax;
  std::int64_t ox, oy, oz, oxmin, oxmax, oymin, oymax, ozmin, ozmax;
  int xmax,ymax,zmax;
  LogicalLocation nloc;
  nloc.level=loc.level;

  xmax=1, oxmin=-1, oxmax=1, nxmax=(rbx<<(loc.level-rl));
  if (dim>=2) ymax=1, oymin=-1, oymax=1, nymax=(rby<<(loc.level-rl));
  else       ymax=0, oymin=0,  oymax=0, nymax=1;
  if (dim==3) zmax=1, ozmin=-1, ozmax=1, nzmax=(rbz<<(loc.level-rl));
  else       zmax=0, ozmin=0,  ozmax=0, nzmax=1;
  for (int k=0; k<=zmax; k++) {
    for (int j=0; j<=ymax; j++) {
      for (int i=0; i<=xmax; i++) {
        pleaf[k][j][i] = new MeshBlockTree(this, i, j, k);
        nnew++;
      }
    }
  }

  for (oz=ozmin; oz<=ozmax; oz++) {
    nloc.lx3=loc.lx3+oz;
    if (nloc.lx3<0) {
      if (mesh_bcs[BoundaryFace::inner_x3]!=BoundaryFlag::periodic) {
        continue;
      } else {
        nloc.lx3=nzmax-1;
      }
    }
    if (nloc.lx3>=nzmax) {
      if (mesh_bcs[BoundaryFace::outer_x3]!=BoundaryFlag::periodic) {
        continue;
      } else {
        nloc.lx3=0;
      }
    }
    for (oy=oymin; oy<=oymax; oy++) {
      nloc.lx2=loc.lx2+oy;
      bool polar=false;
      if (nloc.lx2<0) {
        if (mesh_bcs[BoundaryFace::inner_x2]==BoundaryFlag::periodic) {
          nloc.lx2=nymax-1;
        } else if (mesh_bcs[BoundaryFace::inner_x2]==BoundaryFlag::polar) {
          nloc.lx2=0;
          polar=true;
        } else {
          continue;
        }
      }
      if (nloc.lx2>=nymax) {
        if (mesh_bcs[BoundaryFace::outer_x2]==BoundaryFlag::periodic) {
          nloc.lx2=0;
        } else if (mesh_bcs[BoundaryFace::outer_x2]==BoundaryFlag::polar) {
          nloc.lx2=nymax-1;
          polar=true;
        } else {
          continue;
        }
      }
      if (polar) nloc.lx3=(nloc.lx3+nzmax/2)%nzmax;
      for (ox=oxmin; ox<=oxmax; ox++) {
        if (ox==0 && oy==0 && oz==0) continue;
        nloc.lx1=loc.lx1+ox;
        if (nloc.lx1<0) {
          if (mesh_bcs[BoundaryFace::inner_x1]!=BoundaryFlag::periodic) {
            continue;
          } else {
            nloc.lx1=nxmax-1;
          }
        }
        if (nloc.lx1>=nxmax) {
          if (mesh_bcs[BoundaryFace::outer_x1]!=BoundaryFlag::periodic) {
            continue;
          } else {
            nloc.lx1=0;
          }
        }
        root.AddMeshBlock(root,nloc,dim,mesh_bcs,rbx,rby,rbz,rl,nnew);
      }
    }
  }
  // this block is not a leaf anymore
  flag=false;
  gid=-1;
  nnew--;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::Derefine(MeshBlockTree& root, int dim,
//                                   enum BoundaryFlag* mesh_bcs, std::int64_t rbx,
//                                   std::int64_t rby, std::int64_t rbz, int rl,
//                                   int &ndel)
//  \brief destroy leaves and make this block a leaf

void MeshBlockTree::Derefine(MeshBlockTree& root, int dim, enum BoundaryFlag* mesh_bcs,
                             std::int64_t rbx, std::int64_t rby, std::int64_t rbz, int rl,
                             int &ndel) {
  int s2=0, e2=0, s3=0, e3=0;
  if (dim>=2) s2=-1, e2=1;
  if (dim==3) s3=-1, e3=1;
  for (int ox3=s3; ox3<=e3; ox3++) {
    for (int ox2=s2; ox2<=e2; ox2++) {
      for (int ox1=-1; ox1<=1; ox1++) {
        MeshBlockTree *bt=
            root.FindNeighbor(loc,ox1,ox2,ox3,mesh_bcs,rbx,rby,rbz,rl,true);
        if (bt!=nullptr) {
          if (bt->flag==false) {
            int lis, lie, ljs, lje, lks, lke;
            if (ox1==-1)       lis=lie=1;
            else if (ox1==1)   lis=lie=0;
            else              lis=0, lie=1;
            if (dim>=2) {
              if (ox2==-1)     ljs=lje=1;
              else if (ox2==1) ljs=lje=0;
              else            ljs=0, lje=1;
            } else {
              ljs=lje=0;
            }
            if (dim==3) {
              if (ox3==-1)     lks=lke=1;
              else if (ox3==1) lks=lke=0;
              else            lks=0, lke=1;
            } else {
              lks=lke=0;
            }
            for (int lk=lks; lk<=lke; lk++) {
              for (int lj=ljs; lj<=lje; lj++) {
                for (int li=lis; li<=lie; li++) {
                  if (bt->pleaf[lk][lj][li]->flag==false) return;
                }
              }
            }
          }
        }
      }
    }
  }

  gid=pleaf[0][0][0]->gid; // inherit the first leaf's GID
  for (int k=0; k<=e3; k++) {
    for (int j=0; j<=e2; j++) {
      for (int i=0; i<=1; i++) {
        delete pleaf[k][j][i];
        pleaf[k][j][i]=nullptr;
        ndel++;
      }
    }
  }
  flag=true; // now this is a leaf
  ndel--;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::CountMeshBlock(int& count)
//  \brief creates the Location list sorted by Z-ordering

void MeshBlockTree::CountMeshBlock(int& count) {
  if (pparent==nullptr) count=0;
  if (flag==true) {
    count++;
  } else {
    for (int k=0; k<=1; k++) {
      for (int j=0; j<=1; j++) {
        for (int i=0; i<=1; i++) {
          if (pleaf[k][j][i]!=nullptr)
            pleaf[k][j][i]->CountMeshBlock(count);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::GetMeshBlockList(LogicalLocation *list,
//                                           int *pglist, int& count)
//  \brief creates the Location list sorted by Z-ordering

void MeshBlockTree::GetMeshBlockList(LogicalLocation *list, int *pglist, int& count) {
  if (pparent==nullptr) count=0;
  if (flag==true) {
    list[count]=loc;
    if (pglist!=nullptr)
      pglist[count]=gid;
    gid=count;
    count++;
  } else {
    for (int k=0; k<=1; k++) {
      for (int j=0; j<=1; j++) {
        for (int i=0; i<=1; i++) {
          if (pleaf[k][j][i]!=nullptr)
            pleaf[k][j][i]->GetMeshBlockList(list, pglist, count);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn MeshBlockTree* MeshBlockTree::FindNeighbor(LogicalLocation myloc, int ox1,
//                                    int ox2, int ox3, enum BoundaryFlag* bcs,
//                                    std::int64_t rbx, std::int64_t rby,std::int64_t rbz,
//                                    int rl)
//  \brief find a neighboring block, called from the root of the tree
//         If it is coarser or same level, return the pointer to that block.
//         If it is a finer block, return the pointer to its parent.
//         Note that this function must be called on a completed tree only

MeshBlockTree* MeshBlockTree::FindNeighbor(LogicalLocation myloc, int ox1, int ox2,
                                           int ox3, enum BoundaryFlag* bcs,
                                           std::int64_t rbx, std::int64_t rby,
                                           std::int64_t rbz, int rl, bool amrflag) {
  std::stringstream msg;
  std::int64_t lx, ly, lz;
  int ll;
  int ox, oy, oz;
  MeshBlockTree *bt = this;
  lx=myloc.lx1, ly=myloc.lx2, lz=myloc.lx3, ll=myloc.level;

  lx+=ox1; ly+=ox2; lz+=ox3;
  // periodic and polar boundaries
  if (lx<0) {
    if (bcs[BoundaryFace::inner_x1]==BoundaryFlag::periodic || bcs[BoundaryFace::inner_x1]==BoundaryFlag::shear_periodic) {
      lx=(rbx<<(ll-rl))-1;
    } else {
      return nullptr;
    }
  }
  if (lx>=rbx<<(ll-rl)) {;
    if (bcs[BoundaryFace::outer_x1]==BoundaryFlag::periodic || bcs[BoundaryFace::outer_x1]==BoundaryFlag::shear_periodic) {
      lx=0;
    } else {
      return nullptr;
    }
  }
  bool polar = false;
  if (ly<0) {
    if (bcs[BoundaryFace::inner_x2]==BoundaryFlag::periodic) {
      ly=(rby<<(ll-rl))-1;
    } else if (bcs[BoundaryFace::inner_x2]==BoundaryFlag::polar) {
      ly=0;
      polar=true;
    } else {
      return nullptr;
    }
  }
  if (ly>=rby<<(ll-rl)) {
    if (bcs[BoundaryFace::outer_x2]==BoundaryFlag::periodic) {
      ly=0;
    } else if (bcs[BoundaryFace::outer_x2]==BoundaryFlag::polar) {
      ly=(rby<<(ll-rl))-1;
      polar=true;
    } else {
      return nullptr;
    }
  }
  std::int64_t num_x3 = rbz<<(ll-rl);
  if (lz<0) {
    if (bcs[BoundaryFace::inner_x3]==BoundaryFlag::periodic) {
      lz=num_x3-1;
    } else {
      return nullptr;
    }
  }
  if (lz>=num_x3) {
    if (bcs[BoundaryFace::outer_x3]==BoundaryFlag::periodic) {
      lz=0;
    } else {
      return nullptr;
    }
  }
  if (polar) lz=(lz+num_x3/2)%num_x3;

  if (ll<1) return this; // single grid; return itself


  for (int level=0; level<ll; level++) {
    if (bt->flag==true) { // leaf
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
    ox = (((lx>>sh) & 1LL) == 1LL);
    oy = (((ly>>sh) & 1LL) == 1LL);
    oz = (((lz>>sh) & 1LL)  == 1LL);
    bt=bt->pleaf[oz][oy][ox];
    if (bt==nullptr) {
      msg << "### FATAL ERROR in FindNeighbor" << std::endl
          << "Neighbor search failed. The Block Tree is broken." << std::endl;
      ATHENA_ERROR(msg);
      return nullptr;
    }
  }
  if (bt->flag==true) // leaf on the same level
    return bt;
  ox=oy=oz=0;
  // one level finer: check if they are leaves
  if (ox1<0) ox=1;
  if (ox2<0) oy=1;
  if (ox3<0) oz=1;
  if (bt->pleaf[oz][oy][ox]->flag==true)
    return bt;  // return this block
  if (amrflag==false) {
    msg << "### FATAL ERROR in FindNeighbor" << std::endl
        << "Neighbor search failed. The Block Tree is broken." << std::endl;
    ATHENA_ERROR(msg);
  }
  return bt;
}

//----------------------------------------------------------------------------------------
//! \fn MeshBlockTree* MeshBlockTree::FindMeshBlock(LogicalLocation tloc)
//  \brief find MeshBlock with LogicalLocation tloc and return a pointer

MeshBlockTree* MeshBlockTree::FindMeshBlock(LogicalLocation tloc) {
  if (tloc.level==loc.level) return this;
  // get leaf indexes
  int sh = tloc.level-loc.level-1;
  int mx = (((tloc.lx1>>sh) & 1LL) == 1LL);
  int my = (((tloc.lx2>>sh) & 1LL) == 1LL);
  int mz = (((tloc.lx3>>sh) & 1LL) == 1LL);
  if (pleaf[mz][my][mx] == nullptr) {
    return nullptr;
  } else {
    return pleaf[mz][my][mx]->FindMeshBlock(tloc);
  }
}
