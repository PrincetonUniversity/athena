#ifndef FLUID_HPP
#define FLUID_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file fluid.hpp
 *  \brief defines derived class Fluid, based on base class Mesh
 *  contains data structures and functions for a Fluid stored on the Mesh
 *====================================================================================*/

// Athena headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray

class MeshBlock;
class ParameterInput;
class FluidIntegrator;
class FluidBoundaryConditions;
class FluidEqnOfState;

//! \class Fluid
//  \brief fluid data and functions

class Fluid {
friend class FluidIntegrator;
public:
  Fluid(MeshBlock *pmb, ParameterInput *pin);
  ~Fluid();
  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Fluid

  AthenaArray<Real> u,w;   // conserved and primitive variables
  AthenaArray<Real> u1,w1; // conserved and primitive variables at the half-time step

  AthenaArray<Real> g, g_inv;  // metric and its inverse

  FluidIntegrator *pf_integrator;   // integration algorithm
  FluidBoundaryConditions *pf_bcs;  // boundary conditions
  FluidEqnOfState *pf_eos;

// conserved to primitive functions implemented in files in /convert_var
  void ConservedToPrimitive(AthenaArray<Real> &c, AthenaArray<Real> &p_old,
      AthenaArray<Real> &p);

  void NewTimeStep(MeshBlock *pmb);    // computes new timestep on a MeshBlock
  void InitFluid(ParameterInput *pin); // problem generator function (files in /pgen)
  Real GetGamma() const {return gamma_;}

private:
  Real gamma_;                       // ratio of specific heats
  AthenaArray<Real> dt1_,dt2_,dt3_;  // scratch arrays used in NewTimeStep
};
#endif
