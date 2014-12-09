#ifndef FLUID_INTEGRATOR_HPP
#define FLUID_INTEGRATOR_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file fluid_integrator.hpp
//  \brief defines class FluidIntegrator, data and functions to integrate fluid
//======================================================================================

// Athena headers
#include "../../athena.hpp"         // Real
#include "../../athena_arrays.hpp"  // AthenaArray

// Forward declarations
class MeshBlock;
class Fluid;
class ParameterInput;

//! \class FluidIntegrator
//  \brief member functions implement various integration algorithms for the fluid

class FluidIntegrator {
public:
  FluidIntegrator(Fluid *pf, ParameterInput *pin);
  ~FluidIntegrator();

  Fluid *pmy_fluid;  // ptr to Fluid containing this FluidIntegrator

  void Predict(MeshBlock *pmb);
  void Correct(MeshBlock *pmb);

  void RiemannSolver(const int k, const int j, const int il, const int iu,
    const int ivx, const AthenaArray<Real> &bx, AthenaArray<Real> &wl,
    AthenaArray<Real> &wr, AthenaArray<Real> &flx);

  void ReconstructionFuncX1(const int n, const int m, const int k, const int j,
    const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void ReconstructionFuncX2(const int n, const int m, const int k, const int j,
    const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

  void ReconstructionFuncX3(const int n, const int m, const int k, const int j, 
    const AthenaArray<Real> &q, AthenaArray<Real> &ql, AthenaArray<Real> &qr);

private:
// scratch space used in integrator
  AthenaArray<Real> wl_, wr_, flx_, src_; 
  AthenaArray<Real> face_area_, cell_volume_;
};
#endif
