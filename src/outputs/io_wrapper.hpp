#ifndef IO_WRAPPER_HPP
#define IO_WRAPPER_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file io_wrapper.hpp
//  \brief defines a set of small wrapper functions for MPI versus Serial Output.

// C headers
#include <stdio.h>

// Athena++ headers
#include "../athena.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
typedef MPI_File IOWrapperFile;
#else
typedef FILE * IOWrapperFile;
#endif

typedef uint64_t IOWrapperSize_t;
enum rwmode {IO_WRAPPER_READ_MODE, IO_WRAPPER_WRITE_MODE};

class IOWrapper {
public:
#ifdef MPI_PARALLEL
  IOWrapper() {comm=MPI_COMM_WORLD;}
  void SetCommunicator(MPI_Comm scomm) { comm=scomm;}
#else
  IOWrapper() {}
#endif
  ~IOWrapper() {}

  // wrapper functions for basic I/O tasks
  int Open(const char* fname, enum rwmode rw);
  int Read(void *buf, IOWrapperSize_t size, IOWrapperSize_t count);
  int Read_all(void *buf, IOWrapperSize_t size, IOWrapperSize_t count);
  int Read_at_all(void *buf, IOWrapperSize_t size,
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
#endif // IO_WRAPPER_HPP
