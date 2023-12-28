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

//this class header
#include "string_utils.hpp"

//athena++ header
#include <sstream>    // stringstream

//======================================================================================
//! \file string_utils.cpp
//  \brief string utility function implementations
//======================================================================================
//
namespace StringUtils {
  //====================================================================================
  //! \fn static std::vector<std::string> split(std::string str, char delimiter)
  //  \brief split a string, and store sub strings in a vector
  //====================================================================================
  std::vector<std::string> split(std::string str, char delimiter) {
    std::vector<std::string> internal;
    std::stringstream ss(str); // Turn the string into a stream.
    std::string tok;

    while(getline(ss, tok, delimiter)) {
      trim(tok);
      if (!tok.empty()) {
        internal.push_back(tok);
      }
    }
    return internal;
  }

  //====================================================================================
  //! \fn static void trim(std::string &s)
  //  \brief get rid of white spaces leading and trailing a string
  //====================================================================================
  void trim(std::string &s) {
    size_t p = s.find_first_not_of(" \t\n");
    s.erase(0, p);

    p = s.find_last_not_of(" \t\n");
    if (p != std::string::npos) {
      s.erase(p+1);
    }
  }
} // namespace StringUtils
