//======================================================================================
//! \file blockuid.hpp
//  \brief defines the unique Block ID class for MeshBlock
//======================================================================================
#ifndef BLOCKUID_HPP
#define BLOCKUID_HPP
#include "../athena_arrays.hpp"  // AthenaArray
#include "../defs.hpp"

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
  BlockUID() {};
  ~BlockUID() {};
  BlockUID(const BlockUID& bid);
  void SetUID(ID_t *suid, int llevel);
  void GetRawUID(ID_t *puid);
  int GetLevel(void);
  BlockUID &operator= (const BlockUID& bid);
  bool operator== (const BlockUID& bid) const;
  bool operator> (const BlockUID& bid) const;
  bool operator< (const BlockUID& bid) const;
  void CreateUIDfromLocation(int lx, int ly, int lz, int llevel);
  void CreateUIDbyRefinement(BlockUID& coarse, int ox, int oy, int oz);
  void CreateUIDbyDerefinement(BlockUID& fine);
  void GetLocation(int& lx, int& ly, int& lz, int& llevel);
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


#endif

