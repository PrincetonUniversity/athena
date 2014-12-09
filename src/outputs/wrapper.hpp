#ifndef WRAPPER_HPP
#define WRAPPER_HPP
//======================================================================================
//! \file wrapper.hpp
//  \brief small wrapper class for MPI/Serial Output
//======================================================================================

#include <stdio.h>

#ifdef MPI_PARALLEL
typedef MPI_File ResFileHnd;
typedef MPI_Offset ResSize_t;
#else
typedef FILE * ResFileHnd;
typedef size_t ResSize_t;
#endif


class ResFile
{
private:
  ResFileHnd rfh;
public:
  ResFile() {};
  ~ResFile() {};

  // wrapper functions
  int ResFileOpen(const char* fname);
  int ResFileSeek(ResSize_t offset);
  int ResFileSeekToEnd(void);
  int ResFileWrite(void *buf, ResSize_t size, ResSize_t count);
  int ResFileRead(void *buf, ResSize_t size, ResSize_t count);
  int ResFileClose(void);
  ResSize_t ResFileTell(void);
};

#endif
