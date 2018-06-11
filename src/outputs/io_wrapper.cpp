//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file io_wrapper.cpp
//  \brief functions that provide wrapper for MPI-IO versus serial input/output

// C headers
#include <stdio.h>
#include <stdlib.h>

// C++ headers
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ classes headers
#include "../athena.hpp"
#include "io_wrapper.hpp"

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Open(const char* fname, enum rwmode rw)
//  \brief wrapper for {MPI_File_open} versus {fopen} including error check

int IOWrapper::Open(const char* fname, enum rwmode rw) {
  std::stringstream msg;

  if (rw==IO_WRAPPER_READ_MODE) {
#ifdef MPI_PARALLEL
    if (MPI_File_open(comm_,const_cast<char*>(fname),MPI_MODE_RDONLY,MPI_INFO_NULL,&fh_)
       !=MPI_SUCCESS) {  // use const_cast to convince the compiler.
#else
    if ((fh_ = fopen(fname,"rb")) == NULL) {
#endif
      msg << "### FATAL ERROR in function [IOWrapper:Open]"
          <<std::endl<< "Input file '" << fname << "' could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
      return false;
    }

  } else if (rw==IO_WRAPPER_WRITE_MODE) {
#ifdef MPI_PARALLEL
    MPI_File_delete(const_cast<char*>(fname), MPI_INFO_NULL); // truncation
    if (MPI_File_open(comm_,const_cast<char*>(fname),MPI_MODE_WRONLY | MPI_MODE_CREATE,
                     MPI_INFO_NULL,&fh_) != MPI_SUCCESS) {
#else
    if ((fh_ = fopen(fname,"wb")) == NULL) {
#endif
      msg << "### FATAL ERROR in function [IOWrapper:Open]"
          <<std::endl<< "Output file '" << fname << "' could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
      return false;
    }
  } else {
    return false;
  }

  return true;
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read(void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
//  \brief wrapper for {MPI_File_read} versus {fread}

size_t IOWrapper::Read(void *buf, IOWrapperSize_t size, IOWrapperSize_t count) {
#ifdef MPI_PARALLEL
  MPI_Status status;
  int nread;
  if (MPI_File_read(fh_,buf,count*size,MPI_BYTE,&status)!=MPI_SUCCESS) return -1;
  if (MPI_Get_count(&status,MPI_BYTE,&nread)==MPI_UNDEFINED) return -1;
  return nread/size;
#else
  return fread(buf,size,count,fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read_all(void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
//  \brief wrapper for {MPI_File_read_all} versus {fread}

size_t IOWrapper::Read_all(void *buf, IOWrapperSize_t size, IOWrapperSize_t count) {
#ifdef MPI_PARALLEL
  MPI_Status status;
  int nread;
  if (MPI_File_read_all(fh_,buf,count*size,MPI_BYTE,&status)!=MPI_SUCCESS) return -1;
  if (MPI_Get_count(&status,MPI_BYTE,&nread)==MPI_UNDEFINED) return -1;
  return nread/size;
#else
  return fread(buf,size,count,fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read_at_all(void *buf, IOWrapperSize_t size,
//                             IOWrapperSize_t count, IOWrapperSize_t offset)
//  \brief wrapper for {MPI_File_read_at_all} versus {fseek+fread}

size_t IOWrapper::Read_at_all(void *buf, IOWrapperSize_t size,
                           IOWrapperSize_t count, IOWrapperSize_t offset) {
#ifdef MPI_PARALLEL
  MPI_Status status;
  int nread;
  if (MPI_File_read_at_all(fh_,offset,buf,count*size,MPI_BYTE,&status)!=MPI_SUCCESS)
    return -1;
  if (MPI_Get_count(&status,MPI_BYTE,&nread)==MPI_UNDEFINED) return -1;
  return nread/size;
#else
  fseek(fh_, offset, SEEK_SET);
  return fread(buf,size,count,fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Write(const void *buf, IOWrapperSize_t size, IOWrapperSize_t cnt)
//  \brief wrapper for {MPI_File_write} versus {fwrite}

size_t IOWrapper::Write(const void *buf, IOWrapperSize_t size, IOWrapperSize_t cnt) {
#ifdef MPI_PARALLEL
  MPI_Status status;
  int nwrite;
  if (MPI_File_write(fh_,const_cast<void*>(buf),cnt*size,MPI_BYTE,&status)!=MPI_SUCCESS)
    return -1;
  if (MPI_Get_count(&status,MPI_BYTE,&nwrite)==MPI_UNDEFINED) return -1;
  return nwrite/size;
#else
  return fwrite(buf,size,cnt,fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Write_at_all(const void *buf, IOWrapperSize_t size,
//                                  IOWrapperSize_t cnt, IOWrapperSize_t offset)
//  \brief wrapper for {MPI_File_write_at_all} versus {fseek+fwrite}.

size_t IOWrapper::Write_at_all(const void *buf, IOWrapperSize_t size,
                            IOWrapperSize_t cnt, IOWrapperSize_t offset) {
#ifdef MPI_PARALLEL
  MPI_Status status;
  int nwrite;
  if (MPI_File_write_at_all(fh_,offset,const_cast<void*>(buf),cnt*size,MPI_BYTE,&status)
     !=MPI_SUCCESS)
    return -1;
  if (MPI_Get_count(&status,MPI_BYTE,&nwrite)==MPI_UNDEFINED) return -1;
  return nwrite/size;
#else
  fseek(fh_, offset, SEEK_SET);
  return fwrite(buf,size,cnt,fh_);
#endif
}


//----------------------------------------------------------------------------------------
//! \fn void IOWrapper::Close(void)
//  \brief wrapper for {MPI_File_close} versus {fclose}

int IOWrapper::Close(void) {
#ifdef MPI_PARALLEL
  return MPI_File_close(&fh_);
#else
  return fclose(fh_);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int IOWrapper::Seek(IOWrapperSize_t offset)
//  \brief wrapper for {MPI_File_seek} versus {fseek}

int IOWrapper::Seek(IOWrapperSize_t offset) {
#ifdef MPI_PARALLEL
  return MPI_File_seek(fh_,offset,MPI_SEEK_SET);
#else
  return fseek(fh_, offset, SEEK_SET);
#endif
}

//----------------------------------------------------------------------------------------
//! \fn IOWrapperSize_t IOWrapper::GetPosition(void)
//  \brief wrapper for {MPI_File_get_position} versus {ftell}

IOWrapperSize_t IOWrapper::GetPosition(void) {
#ifdef MPI_PARALLEL
  MPI_Offset position;
  MPI_File_get_position(fh_,&position);
  return position;
#else
  return ftell(fh_);
#endif
}
