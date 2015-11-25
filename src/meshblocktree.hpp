//======================================================================================
//! \file meshblocktree.hpp
//  \brief defines the LogicalLocation structure and MeshBlockTree class
//======================================================================================
#ifndef MESHBLOCKTREE_HPP
#define MESHBLOCKTREE_HPP
#include "athena.hpp"
#include "athena_arrays.hpp"  // AthenaArray
#include "defs.hpp"

//! \class MeshBlockTree
//  \brief Construct AMR Block tree structure
class MeshBlockTree
{
private:
  bool flag; // false: vitrual, has leaves, true: real, is a leaf
  MeshBlockTree* pparent;
  MeshBlockTree* pleaf[2][2][2];
  LogicalLocation loc;
  int gid;
  friend class Mesh;
  friend class MeshBlock;
public:
  MeshBlockTree();
  MeshBlockTree(MeshBlockTree *parent, int ox, int oy, int oz);
  ~MeshBlockTree();
  void CreateRootGrid(long int nx, long int ny, long int nz, int nl);
  void AddMeshBlock(MeshBlockTree& root, LogicalLocation rloc, int dim, int* mesh_bcs,
                    long int rbx, long int rby, long int rbz, int rl, int &nnew);
  void AddMeshBlockWithoutRefine(LogicalLocation rloc, 
                                 long int rbx, long int rby, long int rbz, int rl);
  void Refine(MeshBlockTree& root, int dim, int* mesh_bcs,
              long int rbx, long int rby, long int rbz, int rl, int &nnew);
  void Derefine(MeshBlockTree& root, int dim, int* mesh_bcs,
              long int rbx, long int rby, long int rbz, int rl, int &ndel);
  MeshBlockTree* FindMeshBlock(LogicalLocation tloc);
  void CountMeshBlock(int& count);
  void GetMeshBlockList(LogicalLocation *list, int *pglist, int& count);
  MeshBlockTree* FindNeighbor(LogicalLocation myloc, int ox1, int ox2, int ox3, int *bcs,
                              long int rbx, long int rby, long int rbz, int rl);
  MeshBlockTree* GetLeaf(int ox, int oy, int oz);
};

#endif

