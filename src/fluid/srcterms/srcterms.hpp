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

//! \struct PointMass
//  \brief node in a linked list contained data about point masses

typedef struct PointMass {
  Real gm;
  ThreeVector position;
  ThreeVector velocity;
  struct PointMass *pnext;   // pointer to next node
} PointMass;

//! \class FluidSourceTerms
//  \brief data and functions for physical source terms in the fluid

class FluidSourceTerms {
public:
  FluidSourceTerms(Fluid *pf, ParameterInput *pin);
  ~FluidSourceTerms();

  void PhysicalSourceTerms(const Real dt, const AthenaArray<Real> &p,
    AthenaArray<Real> &c);

private:
  Fluid *pmy_fluid_;       // ptr to Fluid containing this FluidSourceTerms
  PointMass *pfirst_mass;  // ptr to first PointMass in linked list
  AthenaArray<Real> src_terms_i_, src_terms_j_;
  AthenaArray<Real> volume_i_,    volume_j_;
};
#endif
