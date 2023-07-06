//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file chemistry_utils.cpp
//! \brief namespace containing chemistry utilities.

// this class header
#include "chemistry_utils.hpp"

// C header
#include <math.h>

// C++ header
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // std::runtime_error()
#include <string>     // c_str()
#include <vector>     // vector container

//athena++ header
#include "../../bvals/bvals.hpp"

namespace ChemistryUtility {

  //--------------------------------------------------------------------------------------
  //! \fn int FindStrIndex(const std::string *str_arr, const int len,
  //!                      const std::string name)
  //! \brief find index of string
  int FindStrIndex(const std::string *str_arr, const int len,
      const std::string name) {
    std::vector<int> ifind;
    std::stringstream msg; //error message
    for (int i=0; i<len; i++) {
      if (str_arr[i] == name) {
        ifind.push_back(i);
      }
    }
    if (ifind.size() == 1) {
      return ifind[0];
    } else if (ifind.size() == 0) {
      msg <<  "### FATAL ERROR in ChemNetwork [FindStrIndex]" << std::endl
        << name << " not found." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    } else {
      msg <<  "### FATAL ERROR in ChemNetwork [FindStrIndex]" << std::endl
        << name << " found more than once (" << ifind.size() << ")."  << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }

} // namespace ChemistryUtility
