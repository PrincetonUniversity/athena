//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file wrapper.cpp
//  \brief functions that provide wrapper for MPI-IO versus serial input/output
//======================================================================================

// C/C++ headers
#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

// Athena++ classes headers
#include "../athena.hpp"

// this class header
#include "wrapper.hpp"

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Open(const char* fname, enum rwmode rw)
//  \brief wrapper for {MPI_File_open} versus {fopen} including error check

int IOWrapper::Open(const char* fname, enum rwmode rw)
{
  std::stringstream msg;

  if(rw==WRAPPER_READ_MODE) {
#ifdef MPI_PARALLEL
    if(MPI_File_open(comm,const_cast<char*>(fname),MPI_MODE_RDONLY,MPI_INFO_NULL,&fh)
       !=MPI_SUCCESS) {  // use const_cast to convince the compiler.
#else
    if ((fh = fopen(fname,"rb")) == NULL) {
#endif
      msg << "### FATAL ERROR in function [IOWrapper:Open]"
          <<std::endl<< "Input file '" << fname << "' could not be opened" <<std::endl;
      throw std::runtime_error(msg.str().c_str());
      return false;
    }

  } else if(rw==WRAPPER_WRITE_MODE) {
#ifdef MPI_PARALLEL
    MPI_File_delete(const_cast<char*>(fname), MPI_INFO_NULL); // truncation
    if(MPI_File_open(comm,const_cast<char*>(fname),MPI_MODE_WRONLY | MPI_MODE_CREATE,
                     MPI_INFO_NULL,&fh) != MPI_SUCCESS) {
#else
    if ((fh = fopen(fname,"wb")) == NULL) {
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

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read(void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
//  \brief wrapper for {MPI_File_read} versus {fread}

int IOWrapper::Read(void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
{
#ifdef MPI_PARALLEL
  MPI_Status status;
  int ierr, nread;
  if(MPI_File_read(fh,buf,count*size,MPI_BYTE,&status)!=MPI_SUCCESS) return -1;
  if(MPI_Get_count(&status,MPI_BYTE,&nread)==MPI_UNDEFINED) return -1;
  return nread/size;
#else
  return fread(buf,size,count,fh);
#endif
}

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read_all(void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
//  \brief wrapper for {MPI_File_read_all} versus {fread}

int IOWrapper::Read_all(void *buf, IOWrapperSize_t size, IOWrapperSize_t count)
{
#ifdef MPI_PARALLEL
  MPI_Status status;
  int ierr, nread;
  if(MPI_File_read_all(fh,buf,count*size,MPI_BYTE,&status)!=MPI_SUCCESS) return -1;
  if(MPI_Get_count(&status,MPI_BYTE,&nread)==MPI_UNDEFINED) return -1;
  return nread/size;
#else
  return fread(buf,size,count,fh);
#endif
}

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Read_at_all(void *buf, IOWrapperSize_t size,
//                             IOWrapperSize_t count, IOWrapperSize_t offset)
//  \brief wrapper for {MPI_File_read_at_all} versus {fseek+fread}

int IOWrapper::Read_at_all(void *buf, IOWrapperSize_t size,
                           IOWrapperSize_t count, IOWrapperSize_t offset)
{
#ifdef MPI_PARALLEL
  MPI_Status status;
  int ierr, nread;
  if(MPI_File_read_at_all(fh,offset,buf,count*size,MPI_BYTE,&status)!=MPI_SUCCESS)
    return -1;
  if(MPI_Get_count(&status,MPI_BYTE,&nread)==MPI_UNDEFINED) return -1;
  return nread/size;
#else
  fseek(fh, offset, SEEK_SET);
  return fread(buf,size,count,fh);
#endif
}

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Write(const void *buf, IOWrapperSize_t size, IOWrapperSize_t cnt)
//  \brief wrapper for {MPI_File_write} versus {fwrite}

int IOWrapper::Write(const void *buf, IOWrapperSize_t size, IOWrapperSize_t cnt)
{
#ifdef MPI_PARALLEL
  MPI_Status status;
  int ierr, nwrite;
  if(MPI_File_write(fh,const_cast<void*>(buf),cnt*size,MPI_BYTE,&status)!=MPI_SUCCESS)
    return -1;
  if(MPI_Get_count(&status,MPI_BYTE,&nwrite)==MPI_UNDEFINED) return -1;
  return nwrite/size;
#else
  return fwrite(buf,size,cnt,fh);
#endif
}

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Write_at_all(const void *buf, IOWrapperSize_t size,
//                                  IOWrapperSize_t cnt, IOWrapperSize_t offset)
//  \brief wrapper for {MPI_File_write_at_all} versus {fseek+fwrite}.

int IOWrapper::Write_at_all(const void *buf, IOWrapperSize_t size,
                            IOWrapperSize_t cnt, IOWrapperSize_t offset)
{
#ifdef MPI_PARALLEL
  MPI_Status status;
  int ierr, nwrite;
  if(MPI_File_write_at_all(fh,offset,const_cast<void*>(buf),cnt*size,MPI_BYTE,&status)
     !=MPI_SUCCESS)
    return -1;
  if(MPI_Get_count(&status,MPI_BYTE,&nwrite)==MPI_UNDEFINED) return -1;
  return nwrite/size;
#else
  fseek(fh, offset, SEEK_SET);
  return fwrite(buf,size,cnt,fh);
#endif
}


//--------------------------------------------------------------------------------------
//! \fn void IOWrapper::Close(void)
//  \brief wrapper for {MPI_File_close} versus {fclose}

int IOWrapper::Close(void)
{
#ifdef MPI_PARALLEL
  return MPI_File_close(&fh);
#else
  return fclose(fh);
#endif
}

//--------------------------------------------------------------------------------------
//! \fn int IOWrapper::Seek(IOWrapperSize_t offset)
//  \brief wrapper for {MPI_File_seek} versus {fseek}

int IOWrapper::Seek(IOWrapperSize_t offset)
{
#ifdef MPI_PARALLEL
  return MPI_File_seek(fh,offset,MPI_SEEK_SET);
#else
  return fseek(fh, offset, SEEK_SET);
#endif
}

//--------------------------------------------------------------------------------------
//! \fn IOWrapperSize_t IOWrapper::GetPosition(void)
//  \brief wrapper for {MPI_File_get_position} versus {ftell}

IOWrapperSize_t IOWrapper::GetPosition(void)
{
#ifdef MPI_PARALLEL
  MPI_Offset position;
  MPI_File_get_position(fh,&position);
  return position;
#else
  return ftell(fh);
#endif
}
