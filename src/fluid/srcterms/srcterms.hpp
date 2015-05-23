#ifndef SRC_TERMS_HPP
#define SRC_TERMS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file srcterms.hpp
//  \brief defines class FluidSourceTerms
//  Contains data and functions that implement physical (not coordinate) source terms
//======================================================================================

// Athena headers
#include "../../athena.hpp"         // Real
#include "../../athena_arrays.hpp"  // AthenaArray

// Declarations
class Fluid;
class ParameterInput;

typedef void (*SrcTermFunc_t)(const Real time, const Real dt,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons);

//! \class FluidSourceTerms
//  \brief data and functions for physical source terms in the fluid

class FluidSourceTerms {
public:
  FluidSourceTerms(Fluid *pf, ParameterInput *pin);
  ~FluidSourceTerms();

  Real GetGM() const {return gm_;}
  Real GetG1() const {return g1_;}
  Real GetG2() const {return g2_;}
  Real GetG3() const {return g3_;}

  void PhysicalSourceTermsX1(int k, int j, const Real dt, 
    const AthenaArray<Real> &flx, const AthenaArray<Real> &p, AthenaArray<Real> &c);
  void PhysicalSourceTermsX2(int k, int j, const Real dt, 
    const AthenaArray<Real> &flx,
    const AthenaArray<Real> &flx_p1, const AthenaArray<Real> &p, AthenaArray<Real> &c);
  void PhysicalSourceTermsX3(int k, int j, const Real dt,
    const AthenaArray<Real> &flx,
    const AthenaArray<Real> &flx_p1, const AthenaArray<Real> &p, AthenaArray<Real> &c);
  void EnrollSrcTermFunction(SrcTermFunc_t my_func);
  void (*UserSourceTerm)(const Real time, const Real dt, const AthenaArray<Real> &prim,
    AthenaArray<Real> &cons);

private:
  Fluid *pmy_fluid_;  // ptr to Fluid containing this FluidSourceTerms
  Real gm_;           // GM for point mass located at origin
  Real g1_, g2_, g3_; // constant acc'n in each direction
};
#endif
