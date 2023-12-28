#ifndef CHEM_RAD_INTEGRATORS_RAD_INTEGRATORS_HPP_
#define CHEM_RAD_INTEGRATORS_RAD_INTEGRATORS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file rad_integrators.hpp
//! \brief definitions for ChemRadIntegrator class

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../bvals/sixray/bvals_sixray.hpp"
#include "../../scalars/scalars.hpp"
#include "../chem_rad.hpp" // radiation

class MeshBlock;
class ChemRadiation;
class ChemNetwork;

//! \class ChemRadIntegrator
//! \brief integrate algorithm for radiative transfer

class ChemRadIntegrator {
  friend class ChemRadiation;
  friend class BoundaryValues;
  friend class ChemRadiationIntegratorTaskList;
 public:
  ChemRadIntegrator(ChemRadiation *pchemrad, ParameterInput *pin);
  ~ChemRadIntegrator();

  ChemRadiation *pmy_rad;
  MeshBlock *pmy_mb;
  ChemNetwork* pmy_chemnet;
  int ncol; // number of column densities needed to track
  AthenaArray<Real> col; // column densitites
  // boundary for column densities
  SixRayBoundaryVariable col_bvar;
  AthenaArray<Real> col_avg, col_Htot, col_CO, col_H2,  col_C; // for debug output

  void CopyToOutput();

  void UpdateRadiation();

 private:
  // calculate column densities within the meshblock, for six_ray
  void GetColMB(BoundaryFace direction);
  // update column density after boundary is received
  void UpdateCol(BoundaryFace direction);
};

#endif // CHEM_RAD_INTEGRATORS_RAD_INTEGRATORS_HPP_
