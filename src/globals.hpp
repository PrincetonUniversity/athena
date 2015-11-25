#ifndef GLOBALS_HPP
#define GLOBALS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file globals.hpp
//  \brief namespace containing global variables defined in globals.cpp
//======================================================================================

namespace Globals
{
  extern int my_rank, nranks;
}

unsigned int CreateMPITag(int lid, int flag, int phys, int bufid);
unsigned int CreateAMRMPITag(int gid, int ox1, int ox2, int ox3);

#endif // GLOBALS_HPP
