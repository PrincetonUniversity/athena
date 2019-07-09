#ifndef RADINTEGRATORS_HPP
#define RADINTEGRATORS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file rad_integrators.hpp
//  \brief definitions for RadIntegrator class
//======================================================================================

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../radiation.hpp" // radiation
#include "../../scalars/scalars.hpp"

class MeshBlock;
class ParameterInput;
class Radiation;
class NeighborBlock;

//! \class RadIntegrator
//  \brief integrate algorithm for radiative transfer


class RadIntegrator {
  friend class Radiation;
  friend class BoundaryValues;
public:
  RadIntegrator(Radiation *prad, ParameterInput *pin);
  ~RadIntegrator();
  
  Radiation *pmy_rad;
  MeshBlock *pmy_mb;

  //average radiation field in all directions to output the radiation field
  //strengths.
  void CopyToOutput();

  //calcuate total column and update radiation
  void UpdateRadiation(int direction);

#ifdef INCLUDE_CHEMISTRY
  int ncol;
  AthenaArray<Real> col, col_avg;
  ChemNetwork* pmy_chemnet;
#endif
};

#endif // RADINTEGRATORS_HPP
