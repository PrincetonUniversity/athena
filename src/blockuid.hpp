//======================================================================================
//! \file blockuid.hpp
//  \brief defines the unique Block ID class for MeshBlock
//======================================================================================
#ifndef BLOCKUID_HPP
#define BLOCKUID_HPP
#include "athena.hpp"
#include "athena_arrays.hpp"  // AthenaArray
#include "defs.hpp"

typedef unsigned int ID_t;

static const int usize=sizeof(ID_t)*8/3;


//! \class BlockUID
//  \brief the unique ID for MeshBlock and related functions
class BlockUID
{
private:
  int level;
  ID_t uid[IDLENGTH];
public:
  BlockUID() : level(0) {};
  ~BlockUID() {};
  BlockUID(const BlockUID& bid);
  void SetUID(ID_t *suid, int llevel);
  void GetRawUID(ID_t *puid);
  int GetLevel(void);
  BlockUID &operator= (const BlockUID& bid);
  bool operator== (const BlockUID& bid) const;
  bool operator> (const BlockUID& bid) const;
  bool operator< (const BlockUID& bid) const;
  void CreateUIDfromLocation(long int lx, long int ly, long int lz, int llevel);
  void CreateUIDbyRefinement(BlockUID& coarse, int ox, int oy, int oz);
  void CreateUIDbyDerefinement(BlockUID& fine);
  void GetLocation(long int& lx, long int& ly, long int& lz, int& llevel);
};


//--------------------------------------------------------------------------------------
//! \fn inline void BlockUID::GetRawUID(ID_t *puid)
//  \brief returns the unique ID as an integer array, mainly used for output
inline void BlockUID::GetRawUID(ID_t *puid)
{
  for(int i=0;i<IDLENGTH;i++)
    puid[i]=uid[i];
  return;
}

//--------------------------------------------------------------------------------------
//! \fn inline int BlockUID::GetLevel(void)
//  \brief returns the logical level
inline int BlockUID::GetLevel(void)
{
  return level;
}



//! \class BlockTree
//  \brief Construct AMR Block tree structure
class BlockTree
{
private:
  bool flag; // false: vitrual, has leaves, true: real, is a leaf
  BlockTree* pparent;
  BlockTree* pleaf[2][2][2];
  BlockUID uid;
  int gid;
public:
  BlockTree();
  BlockTree(BlockTree *parent, int ox, int oy, int oz);
  ~BlockTree();
  void CreateRootGrid(long int nx, long int ny, long int nz, int nl);
  void Refine(int dim);
  void Derefine(void);
  void AssignGID(int& id);
  void GetIDList(BlockUID *list, int& count);
  BlockTree* FindNeighbor(enum direction dir, BlockUID id,
                          long int rbx, long int rby, long int rbz, int rl);
  BlockTree* GetLeaf(int ox, int oy, int oz);
  NeighborBlock GetNeighbor(void);
};

#endif

