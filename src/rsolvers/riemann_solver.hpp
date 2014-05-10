#ifndef RIEMANN_SOLVER_HPP
#define RIEMANN_SOLVER_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file riemann_solver.hpp
 *  \brief defines class RiemannSolvers
 *  implements various Riemann solvers
 *====================================================================================*/

class Fluid;

//! \class ConvertVariables
//  \brief variable conversion data and functions

class RiemannSolver {
public:
  RiemannSolver(Fluid *pf);
  ~RiemannSolver();

  void HLLC(const int il, const int iu, AthenaArray<Real> &wl, AthenaArray<Real> &wr,
            AthenaArray<Real> &flx);

private:
  Fluid *pmy_fluid_;  // pointer to parent Fluid object

};
#endif
