#ifndef CHEM_RAD_CHEM_RAD_HPP_
#define CHEM_RAD_CHEM_RAD_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file chem_rad.hpp
//! \brief definitions for ChemRadiation class

// C headers

// C++ headers
#include <string>

// Athena++ class headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../defs.hpp"
#include "../parameter_input.hpp"

class ChemRadIntegrator;

//! \class Radiation
//! \brief Radiation data and functions

class ChemRadiation {
 public:
  ChemRadiation(MeshBlock *pmb, ParameterInput *pin);
  ~ChemRadiation();

  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Fluid
  ChemRadIntegrator* pchemradintegrator; // ptr to radiation integrator
  std::string integrator;

  AthenaArray<Real> ir; // radiation specific intensity
  // average radiation intensity over all angles, for output
  AthenaArray<Real> ir_avg;

  int nang, nfreq, n_fre_ang; // n_fre_ang=nang*nfreq
  bool output_zone_sec; // output zone/sec for radiation tasklist
};

#endif // CHEM_RAD_CHEM_RAD_HPP_
