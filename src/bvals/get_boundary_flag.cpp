//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file get_boundary_flag.cpp

// C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>

// Athena headers
#include "../athena.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn GetBoundaryFlag(std::string input_string)
//  \brief Parses input string to return integer flag specifying boundary condition

enum BoundaryFlag GetBoundaryFlag(std::string input_string) {
  if (input_string == "reflecting") {
    return REFLECTING_BNDRY;
  } else if (input_string == "outflow") {
    return OUTFLOW_BNDRY;
  } else if (input_string == "user") {
    return USER_BNDRY;
  } else if (input_string == "periodic") {
    return PERIODIC_BNDRY;
  } else if (input_string == "shear_periodic") {
    return SHEAR_PERIODIC_BNDRY;
  } else if (input_string == "polar") {
    return POLAR_BNDRY;
  } else if (input_string == "polar_wedge") {
    return POLAR_BNDRY_WEDGE;
  } else if (input_string == "none") {
    return BNDRY_UNDEF;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetBoundaryFlag" << std::endl
        << "Input string=" << input_string << " not valid boundary type" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
}
