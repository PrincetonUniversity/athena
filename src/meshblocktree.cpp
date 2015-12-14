//======================================================================================
//! \file MeshBlockTree.cpp
//  \brief implementation of functions in the MeshBlockTree class
//  The MeshBlockTree stores the logical grid structure, and is used for 
//  neighbor search, assigning global IDs, etc.
//  Level is defined as "logical level", where the logical root (single block) is 0,
//  and the physical root (user-specified root level) can be different.
//======================================================================================


#include "athena.hpp"
#include "meshblocktree.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>

//--------------------------------------------------------------------------------------
//! \fn MeshBlockTree::MeshBlockTree()
//  \brief constructor of MeshBlockTree, creates the logical root
MeshBlockTree::MeshBlockTree()
{
  flag=true;
  gid=-1;
  loc.lx1=0, loc.lx2=0, loc.lx3=0;
  loc.level=0;
  pparent=NULL;
  for(int k=0; k<=1; k++) {
    for(int j=0; j<=1; j++) {
      for(int i=0; i<=1; i++) {
        pleaf[k][j][i]=NULL;
      }
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn MeshBlockTree::MeshBlockTree(MeshBlockTree *parent, int ox, int oy, int oz)
//  \brief constructor of MeshBlockTree, creates a leaf
MeshBlockTree::MeshBlockTree(MeshBlockTree *parent, int ox, int oy, int oz)
{
  flag=true;
  gid=-1;
  pparent=parent;
  loc.lx1=(parent->loc.lx1<<1)+ox;
  loc.lx2=(parent->loc.lx2<<1)+oy;
  loc.lx3=(parent->loc.lx3<<1)+oz;
  loc.level=parent->loc.level+1;
  for(int k=0; k<=1; k++) {
    for(int j=0; j<=1; j++) {
      for(int i=0; i<=1; i++) {
        pleaf[k][j][i]=NULL;
      }
    }
  }
}


//--------------------------------------------------------------------------------------
//! \fn MeshBlockTree::~MeshBlockTree()
//  \brief destructor of MeshBlockTree, destroy all the leaves
MeshBlockTree::~MeshBlockTree()
{
  for(int k=0; k<=1; k++) {
    for(int j=0; j<=1; j++) {
      for(int i=0; i<=1; i++) {
        if(pleaf[k][j][i]!=NULL)
          delete pleaf[k][j][i];
      }
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::CreateRootGrid(long int nx, long int ny, long int nz, int nl)
//  \brief create the root grid; the root grid can be incomplete (less than 8 leaves)
void MeshBlockTree::CreateRootGrid(long int nx, long int ny, long int nz, int nl)
{
  long int mx, my, mz;
  if(loc.level == nl) {
    return;
  }
  for(int k=0; k<=1; k++) {
    if((loc.lx3*2+k)*(1L<<(nl-loc.level-1)) < nz) {
      for(int j=0; j<=1; j++) {
        if((loc.lx2*2+j)*(1L<<(nl-loc.level-1)) < ny) {
          for(int i=0; i<=1; i++) {
            if((loc.lx1*2+i)*(1L<<(nl-loc.level-1)) < nx) {
              flag=false; // if there is a leaf, this is not a leaf
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


//--------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::AddMeshBlock(MeshBlockTree& root, LogicalLocation rloc,
//           int dim, int* mesh_bcs, long int rbx, long int rby, long int rbz, int rl)
//  \brief add a MeshBlock to the tree, also creates neighboring blocks
void MeshBlockTree::AddMeshBlock(MeshBlockTree& root, LogicalLocation rloc, int dim,
                    int* mesh_bcs, long int rbx, long int rby, long int rbz, int rl)
{
  int mx, my, mz;
  if(loc.level==rloc.level) // done
    return;
  if(flag==true) // leaf -> create the finer level
    Refine(root,dim,mesh_bcs,rbx,rby,rbz,rl);
  // get leaf indexes
  int sh = rloc.level-loc.level-1;
  mx=(int)((rloc.lx1>>sh)&1L);
  my=(int)((rloc.lx2>>sh)&1L);
  mz=(int)((rloc.lx3>>sh)&1L);
  pleaf[mz][my][mx]->AddMeshBlock(root,rloc,dim,mesh_bcs,rbx,rby,rbz,rl);
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::AddMeshBlockWithoutRefine(MeshBlockTree& root,
//                          LogicalLocation rloc, int dim, int* mesh_bcs,
//                          long int rbx, long int rby, long int rbz, int rl)
//  \brief add a MeshBlock to the tree without refinement, used in restarting
void MeshBlockTree::AddMeshBlockWithoutRefine(MeshBlockTree& root, LogicalLocation rloc,
     int dim, int* mesh_bcs, long int rbx, long int rby, long int rbz, int rl)
{
  int mx, my, mz;
  if(loc.level==rloc.level) // done
    return;
  if(flag==true) // leaf -> create the finer level
    flag=false;
  // get leaf indexes
  int sh = rloc.level-loc.level-1;
  mx=(int)((rloc.lx1>>sh)&1L);
  my=(int)((rloc.lx2>>sh)&1L);
  mz=(int)((rloc.lx3>>sh)&1L);
  if(pleaf[mz][my][mx]==NULL)
    pleaf[mz][my][mx] = new MeshBlockTree(this, mx, my, mz);
  pleaf[mz][my][mx]->AddMeshBlockWithoutRefine(root,rloc,dim,mesh_bcs,rbx,rby,rbz,rl);
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::Refine(MeshBlockTree& root, int dim, int* mesh_bcs,
//                             long int rbx, long int rby, long int rbz, int rl)
//  \brief make finer leaves
void MeshBlockTree::Refine(MeshBlockTree& root, int dim, int* mesh_bcs,
                       long int rbx, long int rby, long int rbz, int rl)
{
  long int nx,ny,nz,nxmax,nymax,nzmax;
  long int ox, oy, oz, oxmin, oxmax, oymin, oymax, ozmin, ozmax;
  int xmax,ymax,zmax;
  LogicalLocation nloc;
  nloc.level=loc.level;

  xmax=1, oxmin=-1, oxmax=1, nxmax=(rbx<<(loc.level-rl));
  if(dim>=2) ymax=1, oymin=-1, oymax=1, nymax=(rby<<(loc.level-rl));
  else       ymax=0, oymin=0,  oymax=0, nymax=1;
  if(dim==3) zmax=1, ozmin=-1, ozmax=1, nzmax=(rbz<<(loc.level-rl));
  else       zmax=0, ozmin=0,  ozmax=0, nzmax=1;

  for(int k=0; k<=zmax; k++) {
    for(int j=0; j<=ymax; j++) {
      for(int i=0; i<=xmax; i++)
        pleaf[k][j][i] = new MeshBlockTree(this, i, j, k);
    }
  }

  for(oz=ozmin;oz<=ozmax;oz++) {
    nloc.lx3=loc.lx3+oz;
    if(nloc.lx3<0) {
      if(mesh_bcs[inner_x3]!=4) continue;
      else nloc.lx3=nzmax-1;
    }
    if(nloc.lx3>=nzmax) {
      if(mesh_bcs[outer_x3]!=4) continue;
      else nloc.lx3=0;
    }
    for(oy=oymin;oy<=oymax;oy++) {
      nloc.lx2=loc.lx2+oy;
      if(nloc.lx2<0) {
        if(mesh_bcs[inner_x2]!=4) continue;
        else nloc.lx2=nymax-1;
      }
      if(nloc.lx2>=nymax) {
        if(mesh_bcs[outer_x2]!=4) continue;
        else nloc.lx2=0;
      }
      for(ox=oxmin;ox<=oxmax;ox++) {
        if(ox==0 && oy==0 && oz==0) continue;
        nloc.lx1=loc.lx1+ox;
        if(nloc.lx1<0) {
          if(mesh_bcs[inner_x1]!=4) continue;
          else nloc.lx1=nxmax-1;
        }
        if(nloc.lx1>=nxmax) {
          if(mesh_bcs[outer_x1]!=4) continue;
          else nloc.lx1=0;
        }
        root.AddMeshBlock(root,nloc,dim,mesh_bcs,rbx,rby,rbz,rl);
      }
    }
  }
  gid=-1;
  flag=false; // this block is not a leaf anymore
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::Derefine(void)
//  \brief delete leaves
void MeshBlockTree::Derefine(void)
{
  for(int k=0; k<=1; k++) {
    for(int j=0; j<=1; j++) {
      for(int i=0; i<=1; i++) {
        if(pleaf[k][j][i]!=NULL)
          delete pleaf[k][j][i];
      }
    }
  }
  flag=true; // this block is now a leaf
  // need to recalculate gid
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::AssignGID(int& id)
//  \brief assign IDs to the leaves and count the total number of the blocks
void MeshBlockTree::AssignGID(int& id)
{
  if(pparent==NULL) // clear id if this is the root of the tree
    id=0;

  if(flag==true) {
    gid=id;
    id++;
  }
  else {
    gid=-1;

    for(int k=0; k<=1; k++) {
      for(int j=0; j<=1; j++) {
        for(int i=0; i<=1; i++) {
          if(pleaf[k][j][i]!=NULL)
            pleaf[k][j][i]->AssignGID(id); // depth-first search
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void MeshBlockTree::GetLocationList(LogicalLocation *list, int& count)
//  \brief creates the Location list sorted by Z-ordering
void MeshBlockTree::GetLocationList(LogicalLocation *list, int& count)
{
  if(pparent==NULL) count=0;
  if(flag==true) {
    list[count]=loc;
    count++;
  }
  else {
    for(int k=0; k<=1; k++) {
      for(int j=0; j<=1; j++) {
        for(int i=0; i<=1; i++) {
          if(pleaf[k][j][i]!=NULL)
            pleaf[k][j][i]->GetLocationList(list, count);
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn MeshBlockTree* MeshBlockTree::FindNeighbor(LogicalLocation myloc, int ox1,
//      int ox2, int ox3, int *bcs, long int rbx, long int rby, long int rbz, int rl)
//  \brief find a neighboring block, called from the root of the tree
//         If it is coarser or same level, return the pointer to that block.
//         If it is a finer block, return the pointer to its parent.
MeshBlockTree* MeshBlockTree::FindNeighbor(LogicalLocation myloc, int ox1, int ox2, int ox3,
                              int *bcs, long int rbx, long int rby, long int rbz, int rl)
{
  std::stringstream msg;
  long int lx, ly, lz;
  int ll;
  int ox, oy, oz;
  MeshBlockTree *bt = this;
  lx=myloc.lx1, ly=myloc.lx2, lz=myloc.lx3, ll=myloc.level;

  lx+=ox1; ly+=ox2; lz+=ox3;
  // periodic and polar boundaries
  if(lx<0) {
    if(bcs[inner_x1]==4) lx=(rbx<<(ll-rl))-1;
    else return NULL;
  }
  if(lx>=rbx<<(ll-rl)) {;
    if(bcs[outer_x1]==4) lx=0;
    else return NULL;
  }
  bool polar = false;
  if(ly<0) {
    if(bcs[inner_x2]==4) ly=(rby<<(ll-rl))-1;
    else if(bcs[inner_x2]==5) {
      ly=0;
      polar=true;
    }
    else return NULL;
  }
  if(ly>=rby<<(ll-rl)) {
    if(bcs[outer_x2]==4) ly=0;
    else if(bcs[outer_x2]==5) {
      ly=(rby<<(ll-rl))-1;
      polar=true;
    }
    else return NULL;
  }
  long int num_x3 = rbz<<(ll-rl);
  if(lz<0) {
    if(bcs[inner_x3]==4) lz=num_x3-1;
    else return NULL;
  }
  if(lz>=num_x3) {
    if(bcs[outer_x3]==4) lz=0;
    else return NULL;
  }
  if(polar) lz=(lz+num_x3/2)%num_x3;

  if(ll<1) return this; // single grid; return itself


  for(int level=0;level<ll;level++) {
    if(bt->flag==true) { // leaf
      if(level == ll-1)
        return bt;
      else {
        msg << "### FATAL ERROR in FindNeighbor" << std::endl
            << "Neighbor search failed. The Block Tree is broken." << std::endl;
        throw std::runtime_error(msg.str().c_str());
        return NULL;
      }
    }
    // find a leaf in the next level
    int sh=ll-level-1;
    ox=(int)((lx>>sh) & 1L);
    oy=(int)((ly>>sh) & 1L);
    oz=(int)((lz>>sh) & 1L);
    bt=bt->pleaf[oz][oy][ox];
    if(bt==NULL) {
      msg << "### FATAL ERROR in FindNeighbor" << std::endl
          << "Neighbor search failed. The Block Tree is broken." << std::endl;
      throw std::runtime_error(msg.str().c_str());
      return NULL;
    }
  }
  if(bt->flag==true) // leaf on the same level
    return bt;
  ox=oy=oz=0;
  // one level finer: check if they are leaves
  if(ox1<0) ox=1;
  if(ox2<0) oy=1;
  if(ox3<0) oz=1;
  if(bt->pleaf[oz][oy][ox]->flag==true)
    return bt;  // return this block
  msg << "### FATAL ERROR in FindNeighbor" << std::endl
      << "Neighbor search failed. The Block Tree is broken." << std::endl;
  throw std::runtime_error(msg.str().c_str());
  return NULL;
}



//--------------------------------------------------------------------------------------
//! \fn MeshBlockTree* MeshBlockTree::GetLeaf(int ox, int oy, int oz)
//  \brief returns the pointer to a leaf
MeshBlockTree* MeshBlockTree::GetLeaf(int ox, int oy, int oz)
{
  return pleaf[oz][oy][ox];
}

