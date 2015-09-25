//======================================================================================
//! \file outputs.hpp
//  \brief provides multiple classes to handle ALL types of data output (fluid, bfield,
//  gravity, radiation, particles, etc.)
//======================================================================================
// Wrapper Functions for MPI/Serial Output

#include "wrapper.hpp"
#include "../athena.hpp"
#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#ifdef MPI_PARALLEL
#include <mpi.h>

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Open(const char* fname, enum rwmode rw)
//  \brief wrap fopen + error check
int IOWrapper::Open(const char* fname, enum rwmode rw)
{
  std::stringstream msg;
  if(rw==readmode) {
    if(MPI_File_open(comm,const_cast<char*>(fname),MPI_MODE_RDONLY,MPI_INFO_NULL,&fh)
       !=MPI_SUCCESS) {  // use const_cast to convince the compiler.
      msg << "### FATAL ERROR in function [IOWrapper:Open]"
          << std::endl << "Input file '" << fname << "' could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
      return false;
    }
  }
  else if(rw==writemode) {
    MPI_File_delete(const_cast<char*>(fname), MPI_INFO_NULL); // truncation
    if(MPI_File_open(comm,const_cast<char*>(fname),MPI_MODE_WRONLY | MPI_MODE_CREATE,
                     MPI_INFO_NULL,&fh)!=MPI_SUCCESS) {
      msg << "### FATAL ERROR in function [IOWrapper:Open]"
          << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
      return false;
    }
  }
  else 
    return false;
  return true;
}

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Seek(IOWrapperSize_t offset)
//  \brief wrap fseek
int IOWrapper::Seek(IOWrapperSize_t offset)
{
  return MPI_File_seek(fh,offset,MPI_SEEK_SET);
}

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Write(const void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
//  \brief wrap fwrite
int IOWrapper::Write(const void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
{
  MPI_Status status;
  int ierr, nwrite;
  if(MPI_File_write(fh,const_cast<void*>(buf),count*size,MPI_BYTE,&status)!=MPI_SUCCESS)
    return -1;
  if(MPI_Get_count(&status,MPI_BYTE,&nwrite)==MPI_UNDEFINED)
    return -1;
  return nwrite/size;
}

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read(void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
//  \brief wrap fread
int IOWrapper::Read(void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
{
  MPI_Status status;
  int ierr, nread;
  if(MPI_File_read(fh,buf,count*size,MPI_BYTE,&status)!=MPI_SUCCESS)
    return -1;
  if(MPI_Get_count(&status,MPI_BYTE,&nread)==MPI_UNDEFINED)
    return -1;
  return nread/size;
}

//--------------------------------------------------------------------------------------
//! \fn void IOWrapper::Close(void)
//  \brief wrap fclose
int IOWrapper::Close(void)
{
  return MPI_File_close(&fh);
}

//--------------------------------------------------------------------------------------
//! \fn IOWrapperSize_t IOWrapper::Tell(void)
//  \brief wrap ftell
IOWrapperSize_t IOWrapper::Tell(void)
{
  MPI_Offset tell;
  MPI_File_get_position(fh,&tell);
  return tell;
}

#else // Serial
//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Open(const char* fname)
//  \brief wrap fopen + error check
int IOWrapper::Open(const char* fname, enum rwmode rw)
{
  std::stringstream msg;
  if(rw==readmode) {
    if ((fh = fopen(fname,"rb")) == NULL) {
      msg << "### FATAL ERROR in function [IOWrapper:Open]"
          << std::endl << "Input file '" << fname << "' could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
      return false;
    }
  }
  else if(rw==writemode) {
    if ((fh = fopen(fname,"wb")) == NULL) {
      msg << "### FATAL ERROR in function [IOWrapper:Open]"
          << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
      return false;
    }
  }
  else
    return false;
  return true;
}

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Seek(IOWrapperSize_t offset)
//  \brief wrap fseek
int IOWrapper::Seek(IOWrapperSize_t offset)
{
  return fseek(fh, offset, SEEK_SET);
}


//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Write(const void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
//  \brief wrap fwrite
int IOWrapper::Write(const void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
{
  return fwrite(buf,size,count,fh);
}

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read(void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
//  \brief wrap fread
int IOWrapper::Read(void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
{
  return fread(buf,size,count,fh);
}

//--------------------------------------------------------------------------------------
//! \fn void IOWrapper::Close(void)
//  \brief wrap fclose
int IOWrapper::Close(void)
{
  return fclose(fh);
}

//--------------------------------------------------------------------------------------
//! \fn IOWrapperSize_t IOWrapper::Tell(void)
//  \brief wrap ftell
IOWrapperSize_t IOWrapper::Tell(void)
{
  return ftell(fh);
}
#endif
