#ifndef WRAPPER_HPP
#define WRAPPER_HPP
//======================================================================================
//! \file wrapper.hpp
//  \brief small wrapper class for MPI/Serial Output
//======================================================================================

#include "../athena.hpp"

#include <stdio.h>

#ifdef MPI_PARALLEL
#include <mpi.h>
typedef MPI_File IOWrapperFile;
#else
typedef FILE * IOWrapperFile;
#endif

typedef long int IOWrapperSize_t;

class IOWrapper
{
private:
  IOWrapperFile fh;
#ifdef MPI_PARALLEL
  MPI_Comm comm;
#endif
public:
#ifdef MPI_PARALLEL
  IOWrapper() {comm=MPI_COMM_WORLD;};
  void SetCommunicator(MPI_Comm scomm) { comm=scomm;};
#else
  IOWrapper() {};
#endif
  ~IOWrapper() {};

  // wrapper functions
  int Open(const char* fname, enum rwmode rw);
  int Seek(IOWrapperSize_t offset);
  int Write(const void *buf, IOWrapperSize_t size, IOWrapperSize_t count);
  int Read(void *buf, IOWrapperSize_t size, IOWrapperSize_t count);
  int Close(void);
  IOWrapperSize_t Tell(void);
};

#endif // WRAPPER_HPP
