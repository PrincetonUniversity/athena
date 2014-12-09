//======================================================================================
//! \file outputs.hpp
//  \brief provides multiple classes to handle ALL types of data output (fluid, bfield,
//  gravity, radiation, particles, etc.)
//======================================================================================
// Wrapper Functions for MPI/Serial Output

#include "wrapper.hpp"
#include <sstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#ifdef MPI_PARALLEL
//--------------------------------------------------------------------------------------
//! \fn void ResFile::ResFileOpen(char* fname)
//  \brief wrap fopen + error check
int ResFile::ResFileOpen(char* fname)
{
  return false;
}

//--------------------------------------------------------------------------------------
//! \fn int ResFile::ResFileSeek(ResSize_t offset)
//  \brief wrap fseek
int ResFile::ResFileSeek(ResSize_t offset)
{
  return false;
}

//--------------------------------------------------------------------------------------
//! \fn int ResFile::ResFileSeekToEnd(void)
//  \brief wrap fseek, to the EOF
int ResFile::ResFileSeekToEnd(void)
{
  return false;
}

//--------------------------------------------------------------------------------------
//! \fn int ResFile::ResFileWrite(void *buf, ResSize_t size, ResSize_t count)
//  \brief wrap fwrite
int ResFile::ResFileWrite(void *buf, ResSize_t size, ResSize_t count)
{
  return false;
}

//--------------------------------------------------------------------------------------
//! \fn int ResFile::ResFileRead(void *buf, ResSize_t size, ResSize_t count)
//  \brief wrap fread
int ResFile::ResFileRead(void *buf, ResSize_t size, ResSize_t count)
{
  return false;
}

//--------------------------------------------------------------------------------------
//! \fn void ResFile::ResFileClose(void)
//  \brief wrap fclose
int ResFile::ResFileClose(void)
{
  return false;
}

//--------------------------------------------------------------------------------------
//! \fn ResSize_t ResFile::ResFileTell(void)
//  \brief wrap ftell
ResSize_t ResFile::ResFileTell(void)
{
  return false;
}

#else // Serial
//--------------------------------------------------------------------------------------
//! \fn void ResFile::ResFileOpen(const char* fname)
//  \brief wrap fopen + error check
int ResFile::ResFileOpen(const char* fname)
{
  std::stringstream msg;
  if ((rfh = fopen(fname,"rb+")) == NULL){
    msg << "### FATAL ERROR in function [ResFile:ResFileOpen]"
        << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
    throw std::runtime_error(msg.str().c_str());
    return false;
  }
  return true;
}

//--------------------------------------------------------------------------------------
//! \fn int ResFile::ResFileSeek(ResSize_t offset)
//  \brief wrap fseek
int ResFile::ResFileSeek(ResSize_t offset)
{
  return fseek(rfh, offset, SEEK_SET);
}

//--------------------------------------------------------------------------------------
//! \fn int ResFile::ResFileSeekToEnd(void)
//  \brief wrap fseek, to the EOF
int ResFile::ResFileSeekToEnd(void)
{
  return fseek(rfh,0L,SEEK_END);
}

//--------------------------------------------------------------------------------------
//! \fn int ResFile::ResFileWrite(void *buf, ResSize_t size, ResSize_t count)
//  \brief wrap fwrite
int ResFile::ResFileWrite(void *buf, ResSize_t size, ResSize_t count)
{
  return fwrite(buf,size,count,rfh);
}

//--------------------------------------------------------------------------------------
//! \fn int ResFile::ResFileRead(void *buf, ResSize_t size, ResSize_t count)
//  \brief wrap fread
int ResFile::ResFileRead(void *buf, ResSize_t size, ResSize_t count)
{
  return fread(buf,size,count,rfh);
}

//--------------------------------------------------------------------------------------
//! \fn void ResFile::ResFileClose(void)
//  \brief wrap fclose
int ResFile::ResFileClose(void)
{
  return fclose(rfh);
}

//--------------------------------------------------------------------------------------
//! \fn ResSize_t ResFile::ResFileTell(void)
//  \brief wrap ftell
ResSize_t ResFile::ResFileTell(void)
{
  return ftell(rfh);
}
#endif
