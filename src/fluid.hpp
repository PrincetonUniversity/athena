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

class ParameterInput;
class FluidBoundaryConditions;
class ConvertVariables;
class FluidIntegrator;

//! \class Fluid
//  \brief fluid data and functions

class Fluid {
friend class FluidIntegrator;
public:
  Fluid(Block *pb);
  ~Fluid();

  Block* pparent_block;    // ptr to parent Block

  AthenaArray<Real> u,w;   // conserved and primitive variables
  AthenaArray<Real> u1,w1; // conserved and primitive variables at the half-time step

  FluidBoundaryConditions *pf_bcs;       // boundary conditions for fluid
  ConvertVariables *pcons_to_prim;       // convert conserved-to-primitive
  FluidIntegrator *pf_integrator;        // integration algorithm

  void NewTimeStep(Block *pb);           // computes new timestep on a Block
  void InitProblem(ParameterInput *pin); // problem generator function (files in /pgen)
  Real GetGamma() const {return gamma_;}

private:
  Real gamma_;                       // ratio of specific heats
  AthenaArray<Real> dt1_,dt2_,dt3_;  // scratch arrays used in NewTimeStep
};
#endif
