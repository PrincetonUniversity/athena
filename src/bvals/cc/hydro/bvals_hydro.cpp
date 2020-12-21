//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_hydro.cpp
//! \brief implements boundary functions for Hydro variables and utilities to manage
//! primitive/conservative variable relationship in a derived class of the
//! CellCenteredBoundaryVariable base class.

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../hydro/hydro.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../utils/buffer_utils.hpp"
#include "bvals_hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn HydroBoundaryVariable::HydroBoundaryVariable
//! \brief

HydroBoundaryVariable::HydroBoundaryVariable(
    MeshBlock *pmb, AthenaArray<Real> *var_hydro, AthenaArray<Real> *coarse_var,
    AthenaArray<Real> *var_flux,
    HydroBoundaryQuantity hydro_type) :
    CellCenteredBoundaryVariable(pmb, var_hydro, coarse_var, var_flux),
    hydro_type_(hydro_type) {
  flip_across_pole_ = flip_across_pole_hydro;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroBoundaryVariable::SelectCoarseBuffer(HydroBoundaryQuantity type)
//! \brief

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

//----------------------------------------------------------------------------------------
//! \fn void HydroBoundaryVariable::SwapHydroQuantity
//! \brief
//! \todo (felker):
//! * make general (but restricted) setter fns in CellCentered and FaceCentered

void HydroBoundaryVariable::SwapHydroQuantity(AthenaArray<Real> &var_hydro,
                                              HydroBoundaryQuantity hydro_type) {
  var_cc = &var_hydro;
  SelectCoarseBuffer(hydro_type);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void HydroBoundaryVariable::SetBoundarySameLevel(Real *buf,
//!                                                      const NeighborBlock& nb)
//! \brief Set hydro boundary received from a block on the same level

void HydroBoundaryVariable::SetBoundarySameLevel(Real *buf,
                                                 const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  AthenaArray<Real> &var = *var_cc;

  if (nb.ni.ox1 == 0)     si = pmb->is,        ei = pmb->ie;
  else if (nb.ni.ox1 > 0) si = pmb->ie + 1,      ei = pmb->ie + NGHOST;
  else              si = pmb->is - NGHOST, ei = pmb->is - 1;
  if (nb.ni.ox2 == 0)     sj = pmb->js,        ej = pmb->je;
  else if (nb.ni.ox2 > 0) sj = pmb->je + 1,      ej = pmb->je + NGHOST;
  else              sj = pmb->js - NGHOST, ej = pmb->js - 1;
  if (nb.ni.ox3 == 0)     sk = pmb->ks,        ek = pmb->ke;
  else if (nb.ni.ox3 > 0) sk = pmb->ke + 1,      ek = pmb->ke + NGHOST;
  else              sk = pmb->ks - NGHOST, ek = pmb->ks - 1;

  int p = 0;

  if (nb.polar) {
    for (int n=nl_; n<=nu_; ++n) {
      Real sign = 1.0;
      if (flip_across_pole_ != nullptr) sign = flip_across_pole_[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma omp simd linear(p)
          for (int i=si; i<=ei; ++i) {
            var(n,k,j,i) = sign * buf[p++];
          }
        }
      }
    }
  } else {
    BufferUtility::UnpackData(buf, var, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  }

  if (pbval_->shearing_box == 2) {
    // 2D shearing box in x-z plane: additional step to shift azimuthal velocity
    int sign[2]{1, -1};
    Real qomL = pbval_->qomL_;
    for (int upper=0; upper<2; upper++) {
      if ((pmb->loc.lx1 == pbval_->loc_shear[upper]) && (sign[upper]*nb.ni.ox1 < 0)) {
        for (int k=sk; k<=ek; ++k) {
          for (int j=sj; j<=ej; ++j) {
            for (int i=si; i<=ei; ++i) {
              if (NON_BAROTROPIC_EOS)
                var(IEN,k,j,i) += (0.5/var(IDN,k,j,i))*(
                    SQR(var(IM3,k,j,i) + sign[upper]*qomL*var(IDN,k,j,i))
                    - SQR(var(IM3,k,j,i)));
              var(IM3,k,j,i) += sign[upper]*qomL*var(IDN,k,j,i);
            }
          }
        }
      }
    }
  }
  return;
}
