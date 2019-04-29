//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_shear_hydro.cpp
//  \brief functions that apply shearing box BCs for hydro variables
//========================================================================================

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../eos/eos.hpp"
#include "../../../field/field.hpp"
#include "../../../globals.hpp"
#include "../../../hydro/hydro.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../parameter_input.hpp"
#include "../../../utils/buffer_utils.hpp"
#include "../../bvals.hpp"
#include "../../bvals_interfaces.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//--------------------------------------------------------------------------------------
//! \fn void HydroBoundaryVariable::AddHydroShearForInit()
//  \brief Send shearing box boundary buffers for hydro variables

void HydroBoundaryVariable::AddHydroShearForInit() {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  AthenaArray<Real> &var = *var_cc;

  int jl = pmb->js - NGHOST;
  int ju = pmb->je + NGHOST;
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  Real qomL = pbval_->qomL_;

  int sign[2]{1, -1};
  int ib[2]{pmb->is - NGHOST, pmb->ie + 1};

  // could call modified ShearQuantities(src=shear_cc_, dst=var, upper), by first loading
  // shear_cc_=var for IDN, IM2 so that order of IM2, IEN update to var doesn't matter.
  // Would need to reassign src=shear_cc_ to updated dst=var for IM2 after? Is it used?
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) {
      // step 1. -- add shear to the periodic boundary values
      for (int k=kl; k<=ku; k++) {
        for (int j=jl; j<=ju; j++) {
          for (int i=0; i<NGHOST; i++) {
            int ii = ib[upper] + i;
            // add shear to conservative
            shear_cc_[upper](IM2,k,j,i) = var(IM2,k,j,ii)
                                          + sign[upper]*qomL*var(IDN,k,j,ii);
            // Update energy, THEN x2 momentum (careful!)
            if (NON_BAROTROPIC_EOS) {
              var(IEN,k,j,ii) += (0.5/var(IDN,k,j,ii))*(SQR(shear_cc_[upper](IM2,k,j,i))
                                                        - SQR(var(IM2,k,j,ii)));
            }
            var(IM2,k,j,ii) = shear_cc_[upper](IM2,k,j,i);
          }
        }
      }
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return;
}
// --------------------------------------------------------------------------------------
// ! \fn void HydroBoundaryVariable::ShearQuantities(AthenaArray<Real> &shear_cc_,
//                                                   bool upper)
//  \brief Apply shear to Hydro x2 momentum and energy

void HydroBoundaryVariable::ShearQuantities(AthenaArray<Real> &shear_cc_, bool upper) {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  AthenaArray<Real> &var = *var_cc;
  int js = pmb->js;
  int je = pmb->je;

  int jl = js - NGHOST;
  int ju = je + NGHOST;
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  Real qomL = pbval_->qomL_;
  int sign[2]{1, -1};
  int ib[2]{pmb->is - NGHOST, pmb->ie + 1};

  for (int k=kl; k<=ku; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=0; i<NGHOST; i++) {
        int ii = ib[upper] + i;
        shear_cc_(IM2,k,j,i) += + sign[upper]*qomL*var(IDN,k,j,ii);
        if (NON_BAROTROPIC_EOS) {
          shear_cc_(IEN,k,j,i) += (0.5/var(IDN,k,j,ii))
                                  *(SQR(shear_cc_(IM2,k,j,i)) - SQR(var(IM2,k,j,ii)));
        }
      }
    }
  }
  return;
}
