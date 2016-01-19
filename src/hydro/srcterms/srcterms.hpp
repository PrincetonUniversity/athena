#ifndef SRC_TERMS_HPP
#define SRC_TERMS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file srcterms.hpp
//  \brief defines class HydroSourceTerms
//  Contains data and functions that implement physical (not coordinate) source terms
//======================================================================================

// Athena headers
#include "../../athena.hpp"         // Real
#include "../../athena_arrays.hpp"  // AthenaArray

// Declarations
class Hydro;
class ParameterInput;

//! \class HydroSourceTerms
//  \brief data and functions for physical source terms in the hydro

class HydroSourceTerms {
public:
  HydroSourceTerms(Hydro *phyd, ParameterInput *pin);
  ~HydroSourceTerms();

  Real GetGM() const {return gm_;}
  Real GetG1() const {return g1_;}
  Real GetG2() const {return g2_;}
  Real GetG3() const {return g3_;}

  void PhysicalSourceTerms(const Real dt, const AthenaArray<Real> *flx,
                           const AthenaArray<Real> &p, AthenaArray<Real> &c);
  void EnrollSrcTermFunction(SrcTermFunc_t my_func);
  SrcTermFunc_t UserSourceTerm;

private:
  Hydro *pmy_hydro_;  // ptr to Hydro containing this HydroSourceTerms
  Real gm_;           // GM for point mass located at origin
  Real g1_, g2_, g3_; // constant acc'n in each direction
};
#endif
