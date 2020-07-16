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
//! \file network_wrapper.cpp
//  \brief implementation of functions in class NetworkWrapper 
//======================================================================================

// this class header
#include "network/network.hpp"

//athena++ header
#include "../defs.hpp"

//c++ header
#include <iostream>   // endl, ostream
#include <sstream>    // stringstream
#include <stdexcept>    // runtime_error

NetworkWrapper::NetworkWrapper() {}

NetworkWrapper::~NetworkWrapper() {}

int NetworkWrapper::WrapJacobian(const realtype t,
                          const N_Vector y, const N_Vector fy, 
                          SUNMatrix jac, void *user_data,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  Real t1 = t;
  //copy y and fy values
  Real y1[NSCALARS+1];
  Real fy1[NSCALARS+1];
  for (int i=0; i<NSCALARS+1; i++) {
		y1[i] = NV_Ith_S(y, i);
		fy1[i] = NV_Ith_S(fy, i);
  }
  //temporary storage for return
  Real jac1[NSCALARS+1][NSCALARS+1];
  Real tmp11[NSCALARS+1];
  Real tmp21[NSCALARS+1];
  Real tmp31[NSCALARS+1];
  NetworkWrapper *meptr = (NetworkWrapper*) user_data;
  meptr->Jacobian(t1, y1, fy1, jac1, tmp11, tmp21, tmp31);
  //set J, tmp1, tmp2, tmp3 to return
  for (int i=0; i<NSCALARS+1; i++) {
		NV_Ith_S(tmp1, i) = tmp11[i];
		NV_Ith_S(tmp2, i) = tmp21[i];
		NV_Ith_S(tmp3, i) = tmp31[i];
  }
  for (int i=0; i<NSCALARS+1; i++) {
    for (int j=0; j<NSCALARS+1; j++) {
      SM_ELEMENT_D(jac, i, j) = jac1[i][j];
    }
  }
  return 0;
}

int NetworkWrapper::WrapRHS(const realtype t, const N_Vector y,
                     N_Vector ydot, void *user_data) {
  Real t1 = t;
  Real y1[NSCALARS];
  Real E1;
  //set y1 to y
  for (int i=0; i<NSCALARS; i++) {
		y1[i] = NV_Ith_S(y, i);
  }
  if (NON_BAROTROPIC_EOS) {
    //set energy
    E1 = NV_Ith_S(y, NSCALARS);
  }
  //temporary storage for return
  Real ydot1[NSCALARS];
  NetworkWrapper *meptr = (NetworkWrapper*) user_data;
  meptr->RHS(t1, y1, E1, ydot1);
  //set ydot to return
  for (int i=0; i<NSCALARS; i++) {
		NV_Ith_S(ydot, i) = ydot1[i];
  }
  if (NON_BAROTROPIC_EOS) {
    //set dEdt
    NV_Ith_S(ydot, NSCALARS) = meptr->Edot(t1, y1, E1);
  }
  return 0;
}

void NetworkWrapper::Jacobian(const Real t,
               const Real y[NSCALARS+1], const Real fy[NSCALARS+1], 
               Real jac[NSCALARS+1][NSCALARS+1],
               Real tmp1[NSCALARS+1], Real tmp2[NSCALARS+1], Real tmp3[NSCALARS+1]) {
  std::stringstream msg;
  msg << "### FATAL ERROR in function NetworkWrapper::Jacobian: "
      << "Jacobian not specified but used." << std::endl
      << "specify Jacobian or set <chemistry> user_jac=0." << std::endl;
  ATHENA_ERROR(msg);
}
