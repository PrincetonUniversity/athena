#ifndef RADIATION_RADIATION_HPP_
#define RADIATION_RADIATION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.hpp
//  \brief definitions for Radiation class

// Athena++ headers
#include "../athena.hpp"         // Real
#include "../athena_arrays.hpp"  // AthenaArray

// Forward declarations
class MeshBlock;
class ParameterInput;

//----------------------------------------------------------------------------------------
// Radiation class
// Notes:
//   designed for general relativity

class Radiation {

public:

  // Constructor and destructor
  Radiation(MeshBlock *pmb, ParameterInput *pin);
  ~Radiation();

  // Object pointers
  MeshBlock* pmy_block;  // pointer to containing MeshBlock

  // Parameters
  int nzeta;   // number of polar radiation angles
  int npsi;    // number of azimuthal radiation angles
  int zs, ze;  // start and end zeta-indices
  int ps, pe;  // start and end psi-indices
  int is, ie;  // start and end x1-indices
  int js, je;  // start and end x2-indices
  int ks, ke;  // start and end x3-indices

  // Data arrays
  AthenaArray<Real> zetaf;   // face-centered polar radiation angles
  AthenaArray<Real> zetav;   // volume-centered polar radiation angles
  AthenaArray<Real> dzetaf   // face-to-face polar radiation angle differences
  AthenaArray<Real> psif;    // face-centered azimuthal radiation angles
  AthenaArray<Real> psiv;    // volume-centered azimuthal radiation angles
  AthenaArray<Real> dpsif    // face-to-face azimuthal radiation angle differences
};

#endif // RADIATION_RADIATION_HPP_
