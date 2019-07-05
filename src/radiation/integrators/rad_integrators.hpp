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
  friend class RadiationIntegratorTaskList;
public:
  RadIntegrator(Radiation *prad, ParameterInput *pin);
  ~RadIntegrator();
  
  Radiation *pmy_rad;
  MeshBlock *pmy_mb;

#ifdef INCLUDE_CHEMISTRY
  int ncol;
  AthenaArray<Real> col;
  AthenaArray<Real> col_avg, col_Htot, col_CO, col_H2,  col_C;//TODO:for output
  ChemNetwork* pmy_chemnet;
  //calcuate column within each meshblock
  void GetColMB(int direction);
  //update column density after receive boundary.
  void UpdateCol(int direction);
  //calcuate total column and update radiation
  void UpdateRadiation(int direction);
  //copy column density and average radiation field to output
  void CopyToOutput();
  void SetSixRayNeighbors();
#endif
private:
  Real rad_G0_; //unshielded radiation field strengh, uniform.
  Real unit_length_in_cm_;
  //six_ray: factors for column in current and previous cell
  Real f_cell;
  Real f_prev; 
  //index for direction of rays in six-ray
  NeighborBlock* pfacenb_[6]; //TODO:uniform mesh for now
};

#endif // RADINTEGRATORS_HPP
