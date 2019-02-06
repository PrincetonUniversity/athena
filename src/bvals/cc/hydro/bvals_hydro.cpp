//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_hydro.cpp
//  \brief implements boundary functions for Hydro variables in a derived class of
//  the CellCenteredBoundaryVariable base class.

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../mesh/mesh.hpp"
#include "../bvals_cc.hpp"
#include "bvals_hydro.hpp"

//----------------------------------------------------------------------------------------
//! \class HydroBoundaryFunctions

HydroBoundaryVariable::HydroBoundaryVariable(
    MeshBlock *pmb, enum BoundaryType type, AthenaArray<Real> &var_hydro,
    AthenaArray<Real> *var_flux,
    enum HydroBoundaryType hydro_type)
    // AthenaArray<Real> &prim)
    : CellCenteredBoundaryVariable(pmb, type, var_hydro, var_flux) {
  hydro_type_=hydro_type;
}

// destructor
HydroBoundaryVariable::~HydroBoundaryVariable() {
}

//----------------------------------------------------------------------------------------
//! \fn void HydroBoundaryVariable::SelectCoarseBuffer(enum HydroBoundaryType type)
//  \brief

// 3x calls to long switch for coarse buffer selection: Send(), ReceiveAndSet(), Set()
// +1x call to shortened "non-switch": Receive()
void HydroBoundaryVariable::SelectCoarseBuffer(enum HydroBoundaryType hydro_type) {
  // KGF: currently always true:
  if (hydro_type==HYDRO_CONS || hydro_type==HYDRO_PRIM) {
    nl_=0, nu_=NHYDRO-1;
    flip_across_pole_=flip_across_pole_hydro;
    if (pmy_mesh_->multilevel) {
      if (hydro_type==HYDRO_CONS)
        coarse_buf.InitWithShallowCopy(pmy_block_->pmr->coarse_cons_);
      if (hydro_type==HYDRO_PRIM)
        coarse_buf.InitWithShallowCopy(pmy_block_->pmr->coarse_prim_);
    }
  }
  // Smaller switch used only in ReceiveBoundaryBuffers(void)
  // if (hydro_type==HYDRO_CONS || hydro_type==HYDRO_PRIM) {
  //   pbd=&bd_cc_;
  // }
  hydro_type_ = hydro_type;
  return;
}

// This should be made to be a completely general function in CellCenteredBoundaryVariable
// class, and there should be a counterpart in the FaceCentereDBoundaryVariable class

// Even without the switching between u and w in Hydro, there may be a future case in
// which function calls want to operate on b1, b2, u1, u2, or some general quantity

// However, allowing such generality may be dangerous, and somewhat unnecessary for any
// physical quantity other than "b" or "u". In both of these special cases, the
// alternative is to pass "int nregister", and have the function directly access Hydro
// "u1, u2" or "w, w1" or Field "b1, b2".
void HydroBoundaryVariable::SwapHydroQuantity(AthenaArray<Real> &var_hydro,
                                              enum HydroBoundaryType hydro_type) {
  var_cc.InitWithShallowCopy(var_hydro);
  // KGF: src and dst are completely useless as-is, since they always mirror var_cc
  src.InitWithShallowCopy(var_cc);
  dst.InitWithShallowCopy(var_cc);

  SelectCoarseBuffer(hydro_type);
  return;
}
