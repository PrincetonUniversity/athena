#ifndef FLUID_HPP
#define FLUID_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file fluid.hpp
//  \brief defines Fluid class which implements data and functions for thermal fluid
//======================================================================================

// Athena headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray

class MeshBlock;
class ParameterInput;
class FluidIntegrator;
class FluidBCs;
class FluidEqnOfState;
class FluidSourceTerms;

//! \class Fluid
//  \brief fluid data and functions

class Fluid {
//friend class FluidIntegrator;
friend class Field;
public:
  Fluid(MeshBlock *pmb, ParameterInput *pin);
  ~Fluid();
  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Fluid

  AthenaArray<Real> u,w;   // conserved and primitive variables
  AthenaArray<Real> u1,w1; // conserved and primitive variables at intermediate step

  AthenaArray<Real> g, g_inv;  // metric and its inverse

  AthenaArray<Real> ifov;  // internal fluid output variables for analysis

  FluidIntegrator *pf_integrator;  // integration algorithm
  FluidEqnOfState *pf_eos;         // equation of state (including cons->prim func)
  FluidSourceTerms *pf_srcterms;   // physical source terms

  void NewBlockTimeStep(MeshBlock *pmb);    // computes new timestep on a MeshBlock
//  void InitFluid(ParameterInput *pin); // problem generator function (files in /pgen)

private:
  AthenaArray<Real> dt1_,dt2_,dt3_;  // scratch arrays used in NewTimeStep
};
#endif
