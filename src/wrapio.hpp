#ifndef WRAPPER_HPP
#define WRAPPER_HPP
//======================================================================================
//! \file wrapper.hpp
//  \brief small wrapper class for MPI/Serial Output
//======================================================================================

#include "athena.hpp"

#include <stdio.h>

#ifdef MPI_PARALLEL
#include <mpi.h>
typedef MPI_File WrapIOFile;
#else
typedef FILE * WrapIOFile;
#endif

typedef long int WrapIOSize_t;

class WrapIO
{
private:
  WrapIOFile fh;
#ifdef MPI_PARALLEL
  MPI_Comm comm;
#endif
public:
#ifdef MPI_PARALLEL
  WrapIO() {comm=MPI_COMM_WORLD;};
  void SetCommunicator(MPI_Comm scomm) { comm=scomm;};
#else
  WrapIO() {};
#endif
  ~WrapIO() {};

  // wrapper functions
  int Open(const char* fname, enum rwmode rw);
  int Seek(WrapIOSize_t offset);
  int Write(const void *buf, WrapIOSize_t size, WrapIOSize_t count);
  int Read(void *buf, WrapIOSize_t size, WrapIOSize_t count);
  int Close(void);
  WrapIOSize_t Tell(void);
};

#endif
