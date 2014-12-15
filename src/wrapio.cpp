//======================================================================================
//! \file outputs.hpp
//  \brief provides multiple classes to handle ALL types of data output (fluid, bfield,
//  gravity, radiation, particles, etc.)
//======================================================================================
// Wrapper Functions for MPI/Serial Output

#include "wrapio.hpp"
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
//! \fn void WrapIO::Open(char* fname)
//  \brief wrap fopen + error check
int WrapIO::Open(char* fname)
{
  std::stringstream msg;
  if(MPI_File_Open(comm,fname,MPI_MODE_RDWR,&fh)!=MPI_SUCCESS) {
    msg << "### FATAL ERROR in function [WrapIO:WrapIOOpen]"
        << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
    throw std::runtime_error(msg.str().c_str());
    return false;
  }
  return true;
}

//--------------------------------------------------------------------------------------
//! \fn int WrapIO::Seek(WrapIOSize_t offset)
//  \brief wrap fseek
int WrapIO::Seek(WrapIOSize_t offset)
{
  return MPI_File_Seek(fh,offset,MPI_SEEK_SET);
}

//--------------------------------------------------------------------------------------
//! \fn int WrapIO::Write(const void *buf, WrapIOSize_t size, WrapIOSize_t count)
//  \brief wrap fwrite
int WrapIO::Write(const void *buf, WrapIOSize_t size, WrapIOSize_t count)
{
  MPI_Status status;
  int ierr, nwrite;
  if(MPI_File_Write(fh,buf,count*size,MPI_BYTE,&status)!=MPI_SUCCSESS)
    return -1;
  if(MPI_Get_Count(status,MPI_BYTE,&nwrite)==MPI_UNDEFINED)
    return -1;
  return nwrite/size;
}

//--------------------------------------------------------------------------------------
//! \fn int WrapIO::Read(void *buf, WrapIOSize_t size, WrapIOSize_t count)
//  \brief wrap fread
int WrapIO::Read(void *buf, WrapIOSize_t size, WrapIOSize_t count)
{
  MPI_Status status;
  int ierr, nread;
  if(MPI_File_Read(fh,buf,count*size,MPI_BYTE,&status)!=MPI_SUCCSESS)
    return -1;
  if(MPI_Get_Count(status,MPI_BYTE,&nread)==MPI_UNDEFINED)
    return -1;
  return nread/size;
}

//--------------------------------------------------------------------------------------
//! \fn void WrapIO::Close(void)
//  \brief wrap fclose
int WrapIO::Close(void)
{
  return MPI_File_Close(&fh);
}

//--------------------------------------------------------------------------------------
//! \fn WrapIOSize_t WrapIO::Tell(void)
//  \brief wrap ftell
WrapIOSize_t WrapIO::Tell(void)
{
  WrapIOSize_t tell;
  MPI_File_get_position(fh,&tell);
  return tell;
}

#else // Serial
//--------------------------------------------------------------------------------------
//! \fn void WrapIO::Open(const char* fname)
//  \brief wrap fopen + error check
int WrapIO::Open(const char* fname)
{
  std::stringstream msg;
  if ((fh = fopen(fname,"rb+")) == NULL){
    msg << "### FATAL ERROR in function [WrapIO:WrapIOOpen]"
        << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
    throw std::runtime_error(msg.str().c_str());
    return false;
  }
  return true;
}

//--------------------------------------------------------------------------------------
//! \fn int WrapIO::Seek(WrapIOSize_t offset)
//  \brief wrap fseek
int WrapIO::Seek(WrapIOSize_t offset)
{
  return fseek(fh, offset, SEEK_SET);
}


//--------------------------------------------------------------------------------------
//! \fn int WrapIO::Write(const void *buf, WrapIOSize_t size, WrapIOSize_t count)
//  \brief wrap fwrite
int WrapIO::Write(const void *buf, WrapIOSize_t size, WrapIOSize_t count)
{
  return fwrite(buf,size,count,fh);
}

//--------------------------------------------------------------------------------------
//! \fn int WrapIO::Read(void *buf, WrapIOSize_t size, WrapIOSize_t count)
//  \brief wrap fread
int WrapIO::Read(void *buf, WrapIOSize_t size, WrapIOSize_t count)
{
  return fread(buf,size,count,fh);
}

//--------------------------------------------------------------------------------------
//! \fn void WrapIO::Close(void)
//  \brief wrap fclose
int WrapIO::Close(void)
{
  return fclose(fh);
}

//--------------------------------------------------------------------------------------
//! \fn WrapIOSize_t WrapIO::Tell(void)
//  \brief wrap ftell
WrapIOSize_t WrapIO::Tell(void)
{
  return ftell(fh);
}
#endif
