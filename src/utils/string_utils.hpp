#ifndef UTILS_STRING_UTILS_HPP_
#define UTILS_STRING_UTILS_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file string_utils.hpp
//  \brief prototypes of string utility functions
//======================================================================================

// C++ headers
#include <string>     // c_str()
#include <vector>     // vector container

//athena++ header
#include "../athena.hpp"

namespace StringUtils {
  //function to split a string into a vector
  std::vector<std::string> split(std::string str, char delimiter);
  //function to get rid of white space leading/trailing a string
  void trim(std::string &s);
}

#endif //UTILS_STRING_UTILS_HPP_
