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

  FluidIntegrator *pf_integrator;   // integration algorithm

// conserved to primitive functions implemented in files in /convert_var
  void ConservedToPrimitive(AthenaArray<Real> &c, AthenaArray<Real> &p);

  void NewTimeStep(Block *pb);           // computes new timestep on a Block
  void InitProblem(ParameterInput *pin); // problem generator function (files in /pgen)
  Real GetGamma() const {return gamma_;}

private:
  Real gamma_;                       // ratio of specific heats
  AthenaArray<Real> dt1_,dt2_,dt3_;  // scratch arrays used in NewTimeStep
};
#endif
