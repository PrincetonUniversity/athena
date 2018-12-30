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

  // Data
  MeshBlock* pmy_block;            // pointer to containing MeshBlock
  int nzeta;                       // number of polar radiation angles
  int npsi;                        // number of azimuthal radiation angles
  AthenaArray<Real> zetaf, zetav;  // face- and volume-centered polar radiation angles
  AthenaArray<Real> psif, psiv;    // face- and volume-centered azimuthal radiation angles
};

#endif // RADIATION_RADIATION_HPP_
