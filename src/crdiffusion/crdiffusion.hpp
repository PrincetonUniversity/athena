#ifndef CRDIFFUSION_CRDIFFUSION_HPP_
#define CRDIFFUSION_CRDIFFUSION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file crdiffusion.hpp
//! \brief defines CRDiffusion class which implements data and functions for the implicit
//!        Cosmic-Ray diffusion solver 

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../bvals/cc/bvals_cc.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;
class CRDiffusionBoundaryValues;
class MGCRDiffusion;
class MGCRDiffusionDriver;

constexpr int NCOEFF = 7;
enum CoeffIndex {DXX=0, DXY=1, DYX=1, DXZ=2, DZX=2, DYY=3, DYZ=4, DZZ=5, NLAMBDA=6};

//! \class CRDiffusion
//! \brief gravitational potential data and functions

class CRDiffusion {
 public:
  CRDiffusion(MeshBlock *pmb, ParameterInput *pin);
  ~CRDiffusion();
  
  MeshBlock* pmy_block;
  MGCRDiffusion *pmg;
  AthenaArray<Real> ecr, zeta, coeff;
  AthenaArray<Real> coarse_ecr, empty_flux[3];
  AthenaArray<Real> def;   // defect from the Multigrid solver
  bool output_defect;

  CellCenteredBoundaryVariable crbvar;

  void CalculateCoefficients(const AthenaArray<Real> &w,
                             const AthenaArray<Real> &bcc);
  void CalculateIonizationRate(const AthenaArray<Real> &w);

  friend class MGCRDiffusuionDriver;

 private:
  int refinement_idx_;
  Real Dpara_, Dperp_, Lambda_;
};

#endif // CRDIFFUSION_CRDIFFUSION_HPP_
