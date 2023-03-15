//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file network_wrapper.cpp
//! \brief implementation of functions in class NetworkWrapper

// this class header
#include "network/network.hpp"

//c++ header
#include <iostream>   // endl, ostream
#include <sstream>    // stringstream
#include <stdexcept>    // runtime_error

//athena++ header
#include "../defs.hpp"

//----------------------------------------------------------------------------------------
//! \brief NetworkWrapper constructor
NetworkWrapper::NetworkWrapper() {}

//----------------------------------------------------------------------------------------
//! \brief NetworkWrapper destructor
NetworkWrapper::~NetworkWrapper() {}

#ifdef CVODE
//----------------------------------------------------------------------------------------
//! \fn int NetworkWrapper::WrapJacobian(const realtype t, const N_Vector y,
//!            const N_Vector ydot, SUNMatrix jac, void *user_data, N_Vector tmp1,
//!            N_Vector tmp2, N_Vector tmp3)
//! \brief Wrapper for the Jacobian function in CVODE
int NetworkWrapper::WrapJacobian(const realtype t,
                          const N_Vector y, const N_Vector ydot,
                          SUNMatrix jac, void *user_data,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  Real t1 = t;
  //copy y and ydot values
  if (NON_BAROTROPIC_EOS) {//non-isothermal case
    Real y1[NSCALARS+1];
    Real ydot1[NSCALARS+1];
    for (int i=0; i<NSCALARS+1; i++) {
      y1[i] = NV_Ith_S(y, i);
      ydot1[i] = NV_Ith_S(ydot, i);
    }
    //temporary storage for return
    AthenaArray<Real> jac1;
    jac1.NewAthenaArray(NSCALARS+1,NSCALARS+1);
    NetworkWrapper *meptr = static_cast<NetworkWrapper*>(user_data);
    meptr->Jacobian(t1, y1, ydot1, jac1);
    //set J to return
    for (int i=0; i<NSCALARS+1; i++) {
      for (int j=0; j<NSCALARS+1; j++) {
        SM_ELEMENT_D(jac, i, j) = jac1(i,j);
      }
    }
  } else {//isothermal case
    Real y1[NSCALARS];
    Real ydot1[NSCALARS];
    for (int i=0; i<NSCALARS; i++) {
      y1[i] = NV_Ith_S(y, i);
      ydot1[i] = NV_Ith_S(ydot, i);
    }
    //temporary storage for return
    AthenaArray<Real> jac1;
    jac1.NewAthenaArray(NSCALARS, NSCALARS);
    NetworkWrapper *meptr = static_cast<NetworkWrapper*>(user_data);
    meptr->Jacobian_isothermal(t1, y1, ydot1, jac1);
    //set J to return
    for (int i=0; i<NSCALARS; i++) {
      for (int j=0; j<NSCALARS; j++) {
        SM_ELEMENT_D(jac, i, j) = jac1(i,j);
      }
    }
  }
  return 0;
}

//----------------------------------------------------------------------------------------
//! \fn int NetworkWrapper::WrapRHS(const realtype t, const N_Vector y,
//!                                 N_Vector ydot, void *user_data)
//! \brief Wrapper for the RHS function in CVODE
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
  NetworkWrapper *meptr = static_cast<NetworkWrapper*>(user_data);
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
#endif //CVODE

//----------------------------------------------------------------------------------------
//! \fn void NetworkWrapper::Jacobian(const Real t,
//!               const Real y[NSCALARS+1], const Real ydot[NSCALARS+1],
//!               AthenaArray<Real> &jac)
//! \brief Default Jacobian
//TODO (Munan Gong) the default jacobian is not implemented
void __attribute__((weak)) NetworkWrapper::Jacobian(const Real t,
               const Real y[NSCALARS+1], const Real ydot[NSCALARS+1],
               AthenaArray<Real> &jac) {
  std::stringstream msg;
  msg << "### FATAL ERROR in function NetworkWrapper::Jacobian: "
      << "Jacobian not implemented but used." << std::endl
      << "Please implement the Jacobian function or set <chemistry> user_jac=0."
      << std::endl;
  ATHENA_ERROR(msg);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void NetworkWrapper::Jacobian_isothermal(const Real t,
//!               const Real y[NSCALARS], const Real ydot[NSCALARS],
//!              AthenaArray<Real> &jac)
//! \brief Default Jacobian for isothermal EOS
//TODO (Munan Gong) the default jacobian for isothermal EOS is not implemented
void __attribute__((weak)) NetworkWrapper::Jacobian_isothermal(const Real t,
               const Real y[NSCALARS], const Real ydot[NSCALARS],
               AthenaArray<Real> &jac) {
  std::stringstream msg;
  msg << "### FATAL ERROR in function NetworkWrapper::Jacobian_isothermal: "
      << "Jacobian for isothermal EOS not implemented but used." << std::endl
      << "Please implement the Jacobian function or set <chemistry> user_jac=0."
      << std::endl;
  ATHENA_ERROR(msg);
  return;
}
