#ifndef NETWORK_HPP
#define NETWORK_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file network.hpp
//  \brief definitions for chemical  network.
//======================================================================================

// Athena++ classes headers
#include "../../athena.hpp"

//c++ headers
#include <string>

//CVODE headers. 
#include <sundials/sundials_types.h> // realtype type
#include <nvector/nvector_serial.h> // N_Vector type
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix


class PassiveScalars;
class ParameterInput;
class NetworkWrapper;

//! \class NetworkWrapper
//  \brief Wrapper of the chemical network class.
class NetworkWrapper {
public:
  NetworkWrapper();
  virtual ~NetworkWrapper();

  static int WrapJacobian(const realtype t,
                          const N_Vector y, const N_Vector fy, 
                          SUNMatrix jac, void *user_data,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int WrapRHS(const realtype t, const N_Vector y,
                     N_Vector ydot, void *user_data);

  //Jacobian, only necessary when the input parameter user_jac=1 in <chemistry>
  //if user_jac=0 (default), then numerical jacobian is used.
  virtual void Jacobian(const Real t,
               const Real y[NSCALARS], const Real fy[NSCALARS], 
               Real jac[NSCALARS][NSCALARS], Real tmp1[NSCALARS],
               Real tmp2[NSCALARS], Real tmp3[NSCALARS]);

  //------------All functions below has to be overloaded------------
  // Note that the RHS and Jac does NOT have user_data. All parameters should
  // be passed to the class as private variables.
  // initialize the parameters for chemical network. Called before the ODE
  // solver in ODEWrapper::Integrate()
  virtual void InitializeNextStep(const int k, const int j, const int i) = 0;
  // right hand side of ode
  virtual void RHS(const Real t, const Real y[NSCALARS], Real ydot[NSCALARS]) = 0;
};

#endif // NETWORK_HPP
