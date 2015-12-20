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

//======================================================================================
//! \file globals.cpp 
//  \brief namespace containing global variables.
//
// Yes, we all know global variables should NEVER be used, but in fact they are ideal
// for, e.g., global constants that are set once and never changed.  To prevent name
// collisions, global variables are wrapped in their own namespace.
//======================================================================================

#include "athena.hpp"

namespace Globals
{
  int my_rank; // MPI rank of this process, set at start of main()
  int nranks;  // total number of MPI ranks, set at start of main()
}

//--------------------------------------------------------------------------------------
//! \fn unsigned int CreateMPITag(int lid, int flag, int phys, int bufid)
//  \brief calculate an MPI tag
unsigned int CreateMPITag(int lid, int flag, int phys, int bufid)
{
// tag = local id of destination (18) + flag (2) + bufid(7) + physics(4)
  return (lid<<13) | (flag<<11) | (bufid<<4) | phys;
}

//--------------------------------------------------------------------------------------
//! \fn unsigned int CreateAMRMPITag(int lid, int ox1, int ox2, int ox3)
//  \brief calculate an MPI tag for AMR block transfer
unsigned int CreateAMRMPITag(int lid, int ox1, int ox2, int ox3)
{
// tag = local id of destination (25) + ox1(1) + ox2(1) + ox3(1) + physics(4)
  return (lid<<7) | (ox1<<6)| (ox2<<5) | (ox3<<4) | tag_amr;
}


