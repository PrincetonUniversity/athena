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

// Athena headers
#include "../athena.hpp"

// C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>  // mkdir()
#include <unistd.h>    // chdir()

//======================================================================================
//! \file change_rundir.cpp 
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn void ChangeRunDir(const char *pdir)
//  \brief change to input run directory; create if it does not exist yet

void ChangeRunDir(const char *pdir)
{
  std::stringstream msg;

  if (pdir == NULL || *pdir == '\0') return;

  mkdir(pdir, 0775);
  if(chdir(pdir)) {
    msg << "### FATAL ERROR in function [ChangeToRunDir]" << std::endl
        << "Cannot cd to directory '" << pdir << "'";
    throw std::runtime_error(msg.str().c_str());
  }

  return;
}
