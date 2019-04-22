//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_hydro.cpp
//  \brief implements boundary functions for Hydro variables and utilities to manage
//  primitive/conservative variable relationship in a derived class of the
//  CellCenteredBoundaryVariable base class.

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../hydro/hydro.hpp"
#include "../../../mesh/mesh.hpp"
#include "bvals_hydro.hpp"

//----------------------------------------------------------------------------------------
//! \class HydroBoundaryFunctions

HydroBoundaryVariable::HydroBoundaryVariable(
    MeshBlock *pmb, AthenaArray<Real> *var_hydro, AthenaArray<Real> *coarse_var,
    AthenaArray<Real> *var_flux,
    HydroBoundaryQuantity hydro_type)
    : CellCenteredBoundaryVariable(pmb, var_hydro, coarse_var, var_flux) {
  hydro_type_ = hydro_type;
  flip_across_pole_ = flip_across_pole_hydro;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroBoundaryVariable::SelectCoarseBuffer(HydroBoundaryQuantity type)
//  \brief

void HydroBoundaryVariable::SelectCoarseBuffer(HydroBoundaryQuantity hydro_type) {
  if (pmy_mesh_->multilevel) {
    switch (hydro_type) {
      case (HydroBoundaryQuantity::cons): {
        coarse_buf = &(pmy_block_->phydro->coarse_cons_);
        break;
      }
      case (HydroBoundaryQuantity::prim): {
        coarse_buf = &(pmy_block_->phydro->coarse_prim_);
        break;
      }
    }
  }
  hydro_type_ = hydro_type;
  return;
}

// TODO(felker): make general (but restricted) setter fns in CellCentered and FaceCentered
void HydroBoundaryVariable::SwapHydroQuantity(AthenaArray<Real> &var_hydro,
                                              HydroBoundaryQuantity hydro_type) {
  var_cc = &var_hydro;
  SelectCoarseBuffer(hydro_type);
  return;
}
