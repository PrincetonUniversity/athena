#ifndef CHEMISTRY_NETWORK_NETWORK_HPP_
#define CHEMISTRY_NETWORK_NETWORK_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file network.hpp
//! \brief definitions for chemical network.

// C headers

// C++ headers
#include <string>

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../defs.hpp"

// CVODE headers
#ifdef CVODE
#include <nvector/nvector_serial.h> // N_Vector type
#include <sundials/sundials_types.h> // realtype type
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix
#endif // CVODE

class PassiveScalars;
class ParameterInput;
class NetworkWrapper;

//! \class NetworkWrapper
//! \brief Wrapper of the chemical network class.
class NetworkWrapper {
 public:
  NetworkWrapper();
  virtual ~NetworkWrapper();

#ifdef CVODE
  static int WrapJacobian(const realtype t,
                          const N_Vector y, const N_Vector fy,
                          SUNMatrix jac, void *user_data,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int WrapRHS(const realtype t, const N_Vector y,
                     N_Vector ydot, void *user_data);
#endif // CVODE

  // Jacobian, only necessary when the input parameter user_jac=1 in <chemistry>
  // if user_jac=0 (default), then numerical jacobian is used.
  // the dimension is NSCALSRS+1, because the last equation is energy equation
  // (Edot).
  virtual void Jacobian(const Real t, const Real *y, const Real *ydot,
                        AthenaArray<Real> &jac);

  // Jacobian for isothermal EOS. The dimentions are NSPECIES because the lack
  // of energy equation
  virtual void Jacobian_isothermal(const Real t, const Real *y, const Real *ydot,
                                   AthenaArray<Real> &jac);

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
  virtual void RHS(const Real t, const Real *y, const Real E,
                   Real *ydot) = 0;
  // Energy equation. Currently solved with chemistry as coupled ODE.
  // input:
  //    t: time in code unites
  //    y: chemical abundances, same as in RHS.
  //    E: internal energy
  // return: rate of energy change dE/dt.
  virtual Real Edot(const Real t, const Real *y, const Real E) = 0;
};

#endif // CHEMISTRY_NETWORK_NETWORK_HPP_
