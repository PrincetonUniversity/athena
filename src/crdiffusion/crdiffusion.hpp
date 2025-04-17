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

constexpr int NCOEFF = 7, NMATRIX = 19;
enum CoeffIndex {DXX=0, DXY=1, DYX=1, DXZ=2, DZX=2, DYY=3, DYZ=4, DZY=4, DZZ=5,
                 NLAMBDA=6};
enum MatrixIndex {CCC=0, CCM=1, CCP=2, CMC=3, CPC=4, MCC=5, PCC=6, CMM=7, CMP=8, CPM=9,
                  CPP=10, MCM=11, MCP=12, PCM=13, PCP=14, MMC=15, MPC=16, PMC=17, PPC=18};

//! \class CRDiffusion
//! \brief gravitational potential data and functions

class CRDiffusion {
 public:
  CRDiffusion(MeshBlock *pmb, ParameterInput *pin);
  ~CRDiffusion();

  MeshBlock* pmy_block;
  MGCRDiffusion *pmg;
  AthenaArray<Real> ecr, source, zeta, coeff;
  AthenaArray<Real> coarse_ecr, empty_flux[3];
  AthenaArray<Real> def;   // defect from the Multigrid solver
  bool output_defect;

  CellCenteredBoundaryVariable crbvar;

  void CalculateCoefficients(const AthenaArray<Real> &w,
                             const AthenaArray<Real> &bcc);

  friend class MGCRDiffusuionDriver;

 private:
  int refinement_idx_;
  Real Dpara_, Dperp_, Lambda_;
};

#endif // CRDIFFUSION_CRDIFFUSION_HPP_
