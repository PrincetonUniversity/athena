//======================================================================================
//! \file blockuid.cpp
//  \brief implementation of functions in the BlockUID class
//  Level is defined as "logical level", where the logical root (single block) is 0,
//  and the physical root (user-specified root level) can be different.
//  Unique ID (+level) gives the absolute location of the block, and can be used to
//  sort the blocks using the Z-ordering.
//======================================================================================

#include "blockuid.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>

//--------------------------------------------------------------------------------------
//! \fn BlockUID::BlockUID(int maxlevel)
//  \brief the constructor of BlockUID, requires the maximum logical level of refinement
BlockUID::BlockUID(int maxlevel)
{
  int length=(maxlevel+usize-1)/usize; // round up
  uid.NewAthenaArray(length);
}

//--------------------------------------------------------------------------------------
//! \fn BlockUID::~BlockUID()
//  \brief the destructor of BlockUID
BlockUID::~BlockUID()
{
  uid.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
//! \fn BlockUID& BlockUID::operator=(const BlockUID& bid)
//  \brief override the assignment operator
BlockUID& BlockUID::operator=(const BlockUID& bid)
{
  int length=uid.GetDim1();
  level=bid.level;
  for(int i=0; i<length; i++)
  { uid(i)=bid.uid(i);}
  return *this;
}


//--------------------------------------------------------------------------------------
//! \fn bool BlockUID::operator==(const BlockUID& bid)
//  \brief override the comparison (equal) operator
bool BlockUID::operator==(const BlockUID& bid)
{
  int length;  
  if(level!=bid.level)
  {return false;}
  length=uid.GetDim1();
  for(int i=0; i<length; i++)
  {
    if(uid(i)!=bid.uid(i))
    { return false;}
  }
  return true;
}

//--------------------------------------------------------------------------------------
//! \fn bool BlockUID::operator<(const BlockUID& bid)
//  \brief override the comparison (less than) operator
bool BlockUID::operator< (const BlockUID& bid)
{
  int length=uid.GetDim1();
  for(int i=0; i<length-1; i++)
  {
    if(uid(i)>bid.uid(i))
    { return false;}
    if(uid(i)<bid.uid(i))
    { return true;}
  }
  return uid(length-1)<bid.uid(length-1);
}

//--------------------------------------------------------------------------------------
//! \fn bool BlockUID::&operator> (const BlockUID& bid)
//  \brief override the comparison (larger than) operator
bool BlockUID::operator> (const BlockUID& bid)
{
  int length=uid.GetDim1();
  for(int i=0; i<length-1; i++)
  {
    if(uid(i)<bid.uid(i))
    { return false;}
    if(uid(i)>bid.uid(i))
    { return true;}
  }
  return uid(length-1)>bid.uid(length-1);
}


//--------------------------------------------------------------------------------------
//! \fn void BlockUID::SetUID(AthenaArray<ID_t> suid, int llevel)
//  \brief set the unique ID directly, mainly used for restarting
void BlockUID::SetUID(AthenaArray<ID_t> suid, int llevel)
{
  uid=suid;
  level=llevel;
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BlockUID::CreateUIDfromLocation(int lx, int ly, int lz, int llevel)
//  \brief create the unique ID from the location, mainly used for initialization
//         lx, ly, lz are the absolute location in llevel
void BlockUID::CreateUIDfromLocation(int lx, int ly, int lz, int llevel)
{
  long int maxid=1<<llevel;
  int length=uid.GetDim1();
  int bx, by, bz, sh;
  ID_t pack;
  std::stringstream msg;

  if(lx<0 || lx>=maxid || ly<0 || ly>=maxid || lz<0 || lz>=maxid)
  {
    msg << "### FATAL ERROR in CreateUIDfromLocation" << std::endl
        << "The block location is beyond the maximum." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  level=llevel;
  for(int i=0; i<length; i++)
  { uid(i)=0;}
  for(int l=1;l<=llevel;l++)
  {
    sh=llevel-l;
    bx=(lx >> sh-2) & 4;
    by=(ly >> sh-1) & 2;
    bz=(lz >> sh) & 1;
    pack=bx | by | bz;
    sh=(usize-(l-1)%usize-1)*3;
    uid((l-1)/usize) |= (pack << sh);
  }
}

//--------------------------------------------------------------------------------------
//! \fn void BlockUID::CreateUIDbyRefinement(BlockUID coarse, int ox, int oy, int oz)
//  \brief create the unique ID from the ID of the finer level, used for refinement
//         ox, oy, oz indicate the location of the finer level (0=left, 1=right, etc.)
void BlockUID::CreateUIDbyRefinement(BlockUID& coarse, int ox, int oy, int oz)
{
  ID_t pack;
  int sh;
  *this=coarse;
  level++;
  pack = (ox << 2) | (oy << 1) | oz;
  sh=(usize-(level-1)%usize-1)*3;
  uid((level-1)/usize) |= (pack << sh);
}

//--------------------------------------------------------------------------------------
//! \fn void BlockUID::CreateUIDbyRefinement(BlockUID& fine)
//  \brief create the unique ID from the ID of the finer block, used for derefinement
void BlockUID::CreateUIDbyDerefinement(BlockUID& fine)
{
  ID_t mask=~7;
  int sh;
  *this=fine;
  level--;
  sh=(usize-level%usize-1)*3;
  uid(level/usize) &= (mask << sh);
}


//--------------------------------------------------------------------------------------
//! \fn void BlockUID::CreateUIDfromLocation(int lx, int ly, int lz, int llevel)
//  \brief get the location from the unique ID
void BlockUID::GetLocation(int& lx, int& ly, int& lz, int& llevel)
{
  int pack, bx, by, bz, l, sh;
  llevel=level;
  lx=ly=lz=0;
  for(l=1;l<=level;l++)
  {
    sh=(usize-(l-1)%usize-1)*3;
    pack = (uid((l-1)/usize) >> sh);
    bx=(pack>>2) & 1;
    by=(pack>>1) & 1;
    bz=pack & 1;
    sh=llevel-l;
    lx|=(bx<<sh);
    ly|=(by<<sh);
    lz|=(bz<<sh);
  }
}

