#ifndef FLUID_INTEGRATOR_HPP
#define FLUID_INTEGRATOR_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file fluid_integrator.hpp
 *  \brief defines class FluidIntegrator
 * implements data and functions to integrate fluid
 *====================================================================================*/

class Fluid;
class RiemannSolver;
class Reconstruction;

//! \class FluidIntegrator
//  \brief member functions implement various integration algorithms for the fluid

class FluidIntegrator {
public:
  FluidIntegrator(Fluid *pf);
  ~FluidIntegrator();

  void Predict(Block *pb);
  Real cfl_number;

private:
  Fluid *pmy_fluid_;  // pointer to parent Fluid object
  RiemannSolver *flux_;
  Reconstruction *lr_states_;

};
#endif
