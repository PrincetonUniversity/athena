//======================================================================================
//! \file blockuid.hpp
//  \brief defines the unique Block ID class for MeshBlock
//======================================================================================
#ifndef BLOCKUID_HPP
#define BLOCKUID_HPP
#include "athena_arrays.hpp"  // AthenaArray


typedef unsigned int ID_t;

static const int usize=sizeof(ID_t)*8/3;


//! \class BlockUID
//  \brief the unique ID for MeshBlock and related functions
class BlockUID
{
  private:
  int level;
  AthenaArray<ID_t> uid;
  public:
  BlockUID(int maxlevel);
  ~BlockUID();
  void SetUID(AthenaArray<ID_t> suid, int llevel);
  AthenaArray<ID_t> GetRawUID(void);
  int GetLevel(void);
  BlockUID &operator= (const BlockUID& bid);
  bool operator== (const BlockUID& bid);
  bool operator> (const BlockUID& bid);
  bool operator< (const BlockUID& bid);
  void CreateUIDfromLocation(int lx, int ly, int lz, int llevel);
  void CreateUIDbyRefinement(BlockUID& coarse, int ox, int oy, int oz);
  void CreateUIDbyDerefinement(BlockUID& fine);
  void GetLocation(int& lx, int& ly, int& lz, int& llevel);
};


//--------------------------------------------------------------------------------------
//! \fn inline AthenaArray<int> BlockUID::GetRawUID(void)
//  \brief returns the unique ID as an integer array, mainly used for output
inline AthenaArray<ID_t> BlockUID::GetRawUID(void)
{
  return uid;
}

//--------------------------------------------------------------------------------------
//! \fn inline int BlockUID::GetLevel(void)
//  \brief returns the logical level
inline int BlockUID::GetLevel(void)
{
  return level;
}


#endif

