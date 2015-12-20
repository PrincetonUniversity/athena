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
//! \file get_boundary_flag.cpp 
//======================================================================================

// C++ headers
#include <iostream> 
#include <sstream>
#include <stdexcept>  // runtime_error

// Athena headers
#include "../athena.hpp"
#include "bvals.hpp"

//--------------------------------------------------------------------------------------
//! \fn GetBoundaryFlag(std::string input_string)
//  \brief Parses input string to return integer flag specifying boundary condition

enum BoundaryFlag GetBoundaryFlag(std::string input_string)
{
  if (input_string == "reflecting") {
    return REFLECTING_BNDRY;
  } else if (input_string == "outflow") {
    return OUTFLOW_BNDRY;
  } else if (input_string == "user") {
    return USER_BNDRY;
  } else if (input_string == "periodic") {
    return PERIODIC_BNDRY;
  } else if (input_string == "polar") {
    return POLAR_BNDRY;
  } else if (input_string == "none") {
    return BNDRY_UNDEF;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetBoundaryFlag" << std::endl
        << "Input string=" << input_string << " not valid boundary type" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
}

