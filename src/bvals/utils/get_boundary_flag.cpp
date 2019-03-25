//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file get_boundary_flag.cpp

// C headers

// C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>

// Athena++ headers
#include "../../athena.hpp"
#include "../bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn GetBoundaryFlag(std::string input_string)
//  \brief Parses input string to return integer flag specifying boundary condition

BoundaryFlag GetBoundaryFlag(const std::string& input_string) {
  if (input_string == "reflecting") {
    return BoundaryFlag::reflect;
  } else if (input_string == "outflow") {
    return BoundaryFlag::outflow;
  } else if (input_string == "user") {
    return BoundaryFlag::user;
  } else if (input_string == "periodic") {
    return BoundaryFlag::periodic;
  } else if (input_string == "shear_periodic") {
    return BoundaryFlag::shear_periodic;
  } else if (input_string == "polar") {
    return BoundaryFlag::polar;
  } else if (input_string == "polar_wedge") {
    return BoundaryFlag::polar_wedge;
  } else if (input_string == "none") {
    return BoundaryFlag::undef;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetBoundaryFlag" << std::endl
        << "Input string=" << input_string << "\n"
        << "is an invalid boundary type" << std::endl;
    ATHENA_ERROR(msg);
  }
}

//----------------------------------------------------------------------------------------
//! \fn GetBoundaryString(BoundaryFlag input_flag)
//  \brief Parses enumerated type BoundaryFlag internal integer representation to return
//  string describing the boundary condition. Used for diagnostics. Inverse of
//  GetBoundaryFlag()

std::string GetBoundaryString(BoundaryFlag input_flag) {
  if (input_flag == BoundaryFlag::reflect) {
    return "reflecting";
  } else if (input_flag == BoundaryFlag::outflow) {
    return "outflow";
  } else if (input_flag == BoundaryFlag::user) {
    return "user";
  } else if (input_flag == BoundaryFlag::periodic) {
    return "periodic";
  } else if (input_flag == BoundaryFlag::shear_periodic) {
    return "shear_periodic";
  } else if (input_flag == BoundaryFlag::polar) {
    return "polar";
  } else if (input_flag == BoundaryFlag::polar_wedge) {
    return "polar_wedge";
  } else if (input_flag == BoundaryFlag::undef) {
    return "none";
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in GetBoundaryString" << std::endl
        << "Input enum class BoundaryFlag=" << static_cast<int>(input_flag) << "\n"
        << "is an invalid boundary type" << std::endl;
    ATHENA_ERROR(msg);
  }
}
