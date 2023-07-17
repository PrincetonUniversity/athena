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
//! \file radiation_implicit.cpp
//  \brief implementation of implicit radiation class
//======================================================================================


// C headers

// C++ headers
#include <algorithm>
#include <cstdio>  // fopen and fwrite
#include <iostream>  // cout
#include <sstream>  // msg
#include <stdexcept> // runtime erro

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../globals.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../integrators/rad_integrators.hpp"
#include "../radiation.hpp"
#include "radiation_implicit.hpp"


// set the default parameters for SRJ
inline void DefaultSRJ(IMRadiation *pimrad) {
  // One example:
  // SRJ_P = 7
  // SRJ_w1 = 2.0, SRJ_w2 = 1.8,SRJ_w3 = 1.6,SRJ_w4 = 1.4,SRJ_w5 = 1.2,
  // SRJ_w6 = 1.0,SRJ_w7 = 0.6
  // SRJ_Q1 = 1,SRJ_Q2 = 1,SRJ_Q3 = 1,SRJ_Q4 = 1,SRJ_Q5 = 1,SRJ_Q6 = 1,SRJ_Q7 = 1

  for (int i=0; i<9; ++i) {
    pimrad->srj_q[i] = 0;
    pimrad->srj_w[i] = 1.0;
  }
  return;
}

IMRadiation::IMRadiation(Mesh *pm, ParameterInput *pin) {
  // read in the parameters
  // maximum number of iterations
  nlimit_ = pin->GetOrAddInteger("radiation","nlimit",100);
  error_limit_ =  pin->GetOrAddReal("radiation","error_limit",1.e-6);
  cfl_rad = pin->GetOrAddReal("radiation","cfl_rad",1.0);
  ite_scheme = pin->GetOrAddInteger("radiation","iteration",1);
  rb_or_not = pin->GetOrAddInteger("radiation","red_or_black",0);

  omega = pin->GetOrAddReal("radiation","omega",1.0);
  srj_p = std::min(9, pin->GetOrAddInteger("radiation","SRJ_P",0));
  srj_level = 0;
  srj_cnt = 0;

  pimraditlist = new IMRadITTaskList(pm);
  pimradhylist = new IMRadHydroTaskList(pm);
  pimradcomptlist = new IMRadComptTaskList(pm);

  // set default parameters
  SetSRJParameters = DefaultSRJ;
}


void IMRadiation::EnrollSRJFunction(SRJFunc MySRJFunction) {
  SetSRJParameters = MySRJFunction;
}
