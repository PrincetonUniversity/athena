#ifndef MESH_MESHBLOCK_TREE_HPP_
#define MESH_MESHBLOCK_TREE_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file meshblock_tree.hpp
//  \brief defines the LogicalLocation structure and MeshBlockTree class
//======================================================================================

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../defs.hpp"
#include "../bvals/bvals.hpp"


//--------------------------------------------------------------------------------------
//! \class MeshBlockTree
//  \brief Objects are nodes in an AMR MeshBlock tree structure

class MeshBlockTree {
  friend class Mesh;
  friend class MeshBlock;
  friend class BoundaryBase;
public:
  MeshBlockTree();
  MeshBlockTree(MeshBlockTree *parent, int ox, int oy, int oz);
  ~MeshBlockTree();

  // accessor
  MeshBlockTree* GetLeaf(int ox, int oy, int oz) {return pleaf[oz][oy][ox];}

  // functions
  void CreateRootGrid(int64_t nx, int64_t ny, int64_t nz, int nl);
  void AddMeshBlock(MeshBlockTree& root, LogicalLocation rloc, int dim,
       enum BoundaryFlag* mesh_bcs, int64_t rbx, int64_t rby, int64_t rbz,
       int rl, int &nnew);
  void AddMeshBlockWithoutRefine(LogicalLocation rloc,
                                 int64_t rbx, int64_t rby, int64_t rbz, int rl);
  void Refine(MeshBlockTree& root, int dim, enum BoundaryFlag* mesh_bcs,
              int64_t rbx, int64_t rby, int64_t rbz, int rl, int &nnew);
  void Derefine(MeshBlockTree& root, int dim, enum BoundaryFlag* mesh_bcs,
              int64_t rbx, int64_t rby, int64_t rbz, int rl, int &ndel);
  MeshBlockTree* FindMeshBlock(LogicalLocation tloc);
  void CountMeshBlock(int& count);
  void GetMeshBlockList(LogicalLocation *list, int *pglist, int& count);
  MeshBlockTree* FindNeighbor(LogicalLocation myloc, int ox1, int ox2, int ox3,
                 enum BoundaryFlag* bcs, int64_t rbx, int64_t rby, int64_t rbz,
                 int rl, bool amrflag=false);

private:
  // data
  bool flag; // false: virtual node, has leaves; true: real node, is a leaf
  MeshBlockTree* pparent;
  MeshBlockTree* pleaf[2][2][2];
  LogicalLocation loc;
  int gid;
};

#endif // MESH_MESHBLOCK_TREE_HPP_
