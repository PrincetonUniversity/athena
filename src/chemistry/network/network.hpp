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
  //the dimension is NSCALSRS+1, because the last equation is energy equation
  //(Edot).
  virtual void Jacobian(const Real t,
               const Real y[NSCALARS+1], const Real fy[NSCALARS+1], 
               Real jac[NSCALARS+1][NSCALARS+1], Real tmp1[NSCALARS+1],
               Real tmp2[NSCALARS+1], Real tmp3[NSCALARS+1]);

  //------------All functions below has to be overloaded------------
  // Note that the RHS and Jac does NOT have user_data. All parameters should
  // be passed to the class as private variables.

  // initialize the parameters for chemical network in the cell,
  // such as reading density and temperature from hydro variables and
  // reading the radiation field. The parameters are set as private variables
  // of the ChemNetwork class.
  // Called before the ODE solver in ODEWrapper::Integrate()
  // input: k, j, i: index of the cell in the MeshBlock.
  virtual void InitializeNextStep(const int k, const int j, const int i) = 0;
  // right hand side of ode.
  // input: 
  //   t: time in code unites
  //   y: chemical abundances. These are dimensionless concentrations: the
  //   variable r in PassiveScalar class.
  //   E: internal energy
  // output:
  //   ydot: time-derivative of abundance y.
  virtual void RHS(const Real t, const Real y[NSCALARS], const Real E, 
                   Real ydot[NSCALARS]) = 0;
  // Energy equation. Currently solved with chemistry as coupled ODE.
  // input:
  //    t: time in code unites
  //    y: chemical abundances, same as in RHS.
  //    E: internal energy
  // return: rate of energy change dE/dt.
  virtual Real Edot(const Real t, const Real y[NSCALARS], const Real E) = 0;
};

#endif // NETWORK_HPP
