#ifndef FLUID_HPP
#define FLUID_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file fluid.hpp
 *  \brief defines class Fluid
 *  Contains data structures and functions for a Fluid stored on the Mesh
 *====================================================================================*/

class ParameterInput;
class Mesh;

//! \struct FluidData
//  \brief conserved/primitive variables, etc. stored in blocks in the Mesh

typedef struct FluidData {
  AthenaArray<Real> u,w;     // conserved and primitive variables
  Real time, dt;
} FluidData;

//! \class Fluid
//  \brief fluid data and functions

class Fluid {
public:
  Fluid(ParameterInput *pin, Mesh *pm);
  ~Fluid();

  FluidData *proot;

  Real GetGamma() const { return gamma_; }
  Real GetGamma_m1() const { return gamma_m1_; }

// public functions implemented in fluid.cpp
// ProblemGenerator is implemented in one of the files in the pgen directory

  void ProblemGenerator(ParameterInput *pin, Domain *pd);
  void Predict(Mesh *pm);

  void PredictVL2(Mesh *pm);
  void hllc(const int il, const int iu,
     AthenaArray<Real> &wl, AthenaArray<Real> &wr, AthenaArray<Real> &flx);
  void PLM(const int k, const int j, const int il, const int iu, const int dir,
           AthenaArray<Real> &w, AthenaArray<Real> &wl, AthenaArray<Real> &wr);
  void ConservedToPrimitive(Domain *pd, AthenaArray<Real> &u, AthenaArray<Real> &w);

private:
  AthenaArray<Real> u1_,w1_; // conserved and primitive variables at the half-time step
  Real gamma_,gamma_m1_;     // ratio of specific heats

  AthenaArray<Real> wl_,wr_,flx_;

// private functions implemented in fluid.cpp

};
#endif
