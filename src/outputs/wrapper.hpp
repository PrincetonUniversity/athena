#ifndef WRAPPER_HPP
#define WRAPPER_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file wrapper.hpp
//  \brief small wrapper class for MPI/Serial Output.
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
enum rwmode {WRAPPER_READ_MODE, WRAPPER_WRITE_MODE};

class IOWrapper
{
public:
#ifdef MPI_PARALLEL
  IOWrapper() {comm=MPI_COMM_WORLD;};
  void SetCommunicator(MPI_Comm scomm) { comm=scomm;};
#else
  IOWrapper() {};
#endif
  ~IOWrapper() {};

  // wrappers for basic I/O functions
  int Open(const char* fname, enum rwmode rw);
  int Read(void *buf, IOWrapperSize_t size, IOWrapperSize_t count);
  int Read_at(void *buf, IOWrapperSize_t size,
              IOWrapperSize_t count, IOWrapperSize_t offset);
  int Write(const void *buf, IOWrapperSize_t size, IOWrapperSize_t count);
  int Write_at_all(const void *buf, IOWrapperSize_t size,
                   IOWrapperSize_t cnt, IOWrapperSize_t offset);
  int Close(void);
  int Seek(IOWrapperSize_t offset);
  IOWrapperSize_t GetPosition(void);
private:
  IOWrapperFile fh;
#ifdef MPI_PARALLEL
  MPI_Comm comm;
#endif
};
#endif // WRAPPER_HPP
