#ifndef CR_CR_HPP_
#define CR_CR_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file cr.hpp
//  \brief definitions for cosmic ray class
//======================================================================================

// C headers

// C++ headers
#include <string>

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/cc/bvals_cc.hpp"

class MeshBlock;
class ParameterInput;
class CRIntegrator;

// Array indices for  moments
enum {CRE=0, CRF1=1, CRF2=2, CRF3=3};

class CosmicRay {
  friend class CRIntegrator;
  friend class BoundaryValues;
 public:
  CosmicRay(MeshBlock *pmb, ParameterInput *pin);

  MeshBlock* pmy_block;
  AthenaArray<Real> u_cr, u_cr1, u_cr2;  // cosmic ray energy density and flux
  AthenaArray<Real> coarse_cr_;

  // diffusion coefficients for both normal diffusion term, and advection term
  AthenaArray<Real> sigma_diff, sigma_adv;

  AthenaArray<Real> v_adv;  // streaming velocity
  AthenaArray<Real> v_diff; // the diffuion velocity, need to calculate the flux

  int refinement_idx{-1};

  AthenaArray<Real> flux[3]; // store transport flux, also need for refinement

  Real vmax; // the maximum velocity (effective speed of light)
  Real vlim;
  Real max_opacity;

  CellCenteredBoundaryVariable cr_bvar;
  CRIntegrator *pcrintegrator;

  // Function in problem generators to update opacity, streaming, source terms
  void EnrollOpacityFunction(CROpacityFunc MyOpacityFunction);
  void EnrollStreamingFunction(CRStreamingFunc MyStreamingFunction);
  void EnrollUserCRSource(CRSrcTermFunc my_func);
  bool cr_source_defined;

  // The function pointer for the diffusion coefficient
  CROpacityFunc UpdateOpacity;
  CRStreamingFunc UpdateStreaming;

  AthenaArray<Real> cwidth;
  AthenaArray<Real> cwidth1;
  AthenaArray<Real> cwidth2;
  AthenaArray<Real> b_grad_pc; // array to store B\dot Grad Pc
  AthenaArray<Real> b_angle; //sin\theta,cos\theta,sin\phi,cos\phi of B direction

  int stream_flag; // flag to include streaming or not
  int src_flag;    // flag to include CR source term or not

 private:
  CRSrcTermFunc UserSourceTerm_;
};

#endif // CR_CR_HPP_
