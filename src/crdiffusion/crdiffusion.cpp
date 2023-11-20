//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file crdiffusion.cpp
//! \brief implementation of functions in class CRDiffusion

// C headers

// C++ headers
#include <iostream>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../bvals/bvals_interfaces.hpp"
#include "../bvals/cc/bvals_cc.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"
#include "crdiffusion.hpp"
#include "mg_crdiffusion.hpp"


//----------------------------------------------------------------------------------------
//! \fn CRDiffusion::CRDiffusion(MeshBlock *pmb, ParameterInput *pin)
//! \brief CRDiffusion constructor
CRDiffusion::CRDiffusion(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block(pmb), ecr(pmb->ncells3, pmb->ncells2, pmb->ncells1),
    zeta(pmb->ncells3, pmb->ncells2, pmb->ncells1),
    coarse_ecr(pmb->ncc3, pmb->ncc2, pmb->ncc1,
              (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
               AthenaArray<Real>::DataStatus::empty)),
    empty_flux{AthenaArray<Real>(), AthenaArray<Real>(), AthenaArray<Real>()},
    D(NCOEFF, pmb->ncells3, pmb->ncells2, pmb->ncells1),
    nlambda(pmb->ncells3, pmb->ncells2, pmb->ncells1),
    output_defect(false), crbvar(pmb, &ecr, &coarse_ecr, empty_flux, false),
    refinement_idx_(), Dpara_(), Dperp_(), Lambda_() {
  Dpara_ = pin->GetReal("crdiffusion", "Dpara");
  Dperp_ = pin->GetReal("crdiffusion", "Dperp");
  Lambda_ = pin->GetReal("crdiffusion", "Lambda");

  output_defect = pin->GetOrAddBoolean("crdiffusion", "output_defect", false);
  if (output_defect)
    def.NewAthenaArray(pmb->ncells3, pmb->ncells2, pmb->ncells1);

  pmb->RegisterMeshBlockData(ecr);
  // "Enroll" in S/AMR by adding to vector of tuples of pointers in MeshRefinement class
  if (pmb->pmy_mesh->multilevel)
    refinement_idx_ = pmy_block->pmr->AddToRefinement(&ecr, &coarse_ecr);

//  pmg = new MGCRDiffusion(pmb->pmy_mesh->pmgcrd, pmb);

  // Enroll CellCenteredBoundaryVariable object
  crbvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&crbvar);
  pmb->pbval->pcrbvar = &crbvar;
}


//----------------------------------------------------------------------------------------
//! \fn CRDiffusion::~CRDiffusion()
//! \brief CRDiffusion destructor
CRDiffusion::~CRDiffusion() {
  delete pmg;
}


//----------------------------------------------------------------------------------------
//! \fn void CRDiffusion::CalculateCoefficients()
//! \brief Calculate coefficients required for CR calculation
void CRDiffusion::CalculateCoefficients(const AthenaArray<Real> &w,
                                        const AthenaArray<Real> &bcc) {
  int il = pmy_block->is - NGHOST, iu = pmy_block->ie + NGHOST;
  int jl = pmy_block->js, ju = pmy_block->je;
  int kl = pmy_block->ks, ku = pmy_block->ke;
  if (pmy_block->pmy_mesh->f2)
    jl -= NGHOST, ju += NGHOST;
  if (pmy_block->pmy_mesh->f3)
    kl -= NGHOST, ku += NGHOST;
  constexpr Real Dunit = 1.0, nlunit = 1.0;
  Real Dpara = Dpara_ * Dunit, Dperp = Dperp_ * Dunit, Lambda = Lambda_ * nlunit;

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = il; i <= iu; ++i) {
          const Real &bx = bcc(IB1,k,j,i);
          const Real &by = bcc(IB2,k,j,i);
          const Real &bz = bcc(IB3,k,j,i);
          Real ba = std::sqrt(SQR(bx) + SQR(by) + SQR(bz) + TINY_NUMBER);
          Real nx = bx / ba, ny = by / ba, nz = bz / ba;
          D(XX,k,j,i) = Dperp + (Dpara - Dperp) * nx * nx;
          D(XY,k,j,i) =         (Dpara - Dperp) * nx * ny;
          D(XZ,k,j,i) =         (Dpara - Dperp) * nx * nz;
          D(YY,k,j,i) = Dperp + (Dpara - Dperp) * ny * ny;
          D(YZ,k,j,i) =         (Dpara - Dperp) * ny * nz;
          D(ZZ,k,j,i) = Dperp + (Dpara - Dperp) * nz * nz;
          nlambda(k,j,i) = Lambda * w(IDN,k,j,i);
        }
      }
    }
  } else {
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = il; i <= iu; ++i) {
          D(XX,k,j,i) = D(XY,k,j,i) = D(XZ,k,j,i) = D(YY,k,j,i)
                      = D(YZ,k,j,i) = D(ZZ,k,j,i) = Dpara;
          nlambda(k,j,i) = Lambda * w(IDN,k,j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn CRDiffusion::ExpandPhysicalBoundaries()
//! \brief Expand physical boundary values to NGHOST = 2 and to edges/corners.
void CRDiffusion::ExpandPhysicalBoundaries() {
  int is = pmy_block->is, ie = pmy_block->ie,
      js = pmy_block->js, je = pmy_block->je,
      ks = pmy_block->ks, ke = pmy_block->ke;

  // push face boundary values
  if (pmy_block->pbval->nblevel[1][1][0] < 0) {
    for (int k = ks; k <= ke; k++) {
      for (int j = js; j <= je; j++)
        ecr(k, j, is-2) = ecr(k, j, is-1);
    }
  }
  if (pmy_block->pbval->nblevel[1][1][2] < 0) {
    for (int k = ks; k <= ke; k++) {
      for (int j = js; j <= je; j++)
        ecr(k, j, ie+2) = ecr(k, j, ie+1);
    }
  }
  if (pmy_block->pbval->nblevel[1][0][1] < 0) {
    for (int k = ks; k <= ke; k++) {
      for (int i = is; i <= ie; i++)
        ecr(k, js-2, i) = ecr(k, js-1, i);
    }
  }
  if (pmy_block->pbval->nblevel[1][2][1] < 0) {
    for (int k = ks; k <= ke; k++) {
      for (int i = is; i <= ie; i++)
        ecr(k, je+2, i) = ecr(k, je+1, i);
    }
  }
  if (pmy_block->pbval->nblevel[0][1][1] < 0) {
    for (int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++)
        ecr(ks-2, j, i) = ecr(ks-1, j, i);
    }
  }
  if (pmy_block->pbval->nblevel[2][1][1] < 0) {
    for (int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++)
        ecr(ke+2, j, i) = ecr(ke+1, j, i);
    }
  }

  // fill edges
  if (pmy_block->pbval->nblevel[1][0][0] < 0) {
    for (int k = ks; k <= ke; k++) {
      Real p = 0.5*(ecr(k, js-1, is) + ecr(k, js, is-1));
      ecr(k, js-1, is-1) = p;
      ecr(k, js-1, is-2) = p;
      ecr(k, js-2, is-1) = p;
      ecr(k, js-2, is-2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[1][0][2] < 0) {
    for (int k = ks; k <= ke; k++) {
      Real p = 0.5*(ecr(k, js-1, ie) + ecr(k, js, ie+1));
      ecr(k, js-1, ie+1) = p;
      ecr(k, js-1, ie+2) = p;
      ecr(k, js-2, ie+1) = p;
      ecr(k, js-2, ie+2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[1][2][0] < 0) {
    for (int k = ks; k <= ke; k++) {
      Real p = 0.5*(ecr(k, je+1, is) + ecr(k, je, is-1));
      ecr(k, je+1, is-1) = p;
      ecr(k, je+1, is-2) = p;
      ecr(k, je+2, is-1) = p;
      ecr(k, je+2, is-2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[1][2][2] < 0) {
    for (int k = ks; k <= ke; k++) {
      Real p = 0.5*(ecr(k, je+1, ie) + ecr(k, je, ie+1));
      ecr(k, je+1, ie+1) = p;
      ecr(k, je+1, ie+2) = p;
      ecr(k, je+2, ie+1) = p;
      ecr(k, je+2, ie+2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[0][1][0] < 0) {
    for (int j = js; j <= je; j++) {
      Real p = 0.5*(ecr(ks-1, j, is) + ecr(ks, j, is-1));
      ecr(ks-1, j, is-1) = p;
      ecr(ks-1, j, is-2) = p;
      ecr(ks-2, j, is-1) = p;
      ecr(ks-2, j, is-2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[0][1][2] < 0) {
    for (int j = js; j <= je; j++) {
      Real p = 0.5*(ecr(ks-1, j, ie) + ecr(ks, j, ie+1));
      ecr(ks-1, j, ie+1) = p;
      ecr(ks-1, j, ie+2) = p;
      ecr(ks-2, j, ie+1) = p;
      ecr(ks-2, j, ie+2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[2][1][0] < 0) {
    for (int j = js; j <= je; j++) {
      Real p = 0.5*(ecr(ke+1, j, is) + ecr(ke, j, is-1));
      ecr(ke+1, j, is-1) = p;
      ecr(ke+1, j, is-2) = p;
      ecr(ke+2, j, is-1) = p;
      ecr(ke+2, j, is-2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[2][1][2] < 0) {
    for (int j = js; j <= je; j++) {
      Real p = 0.5*(ecr(ke+1, j, ie) + ecr(ke, j, ie+1));
      ecr(ke+1, j, ie+1) = p;
      ecr(ke+1, j, ie+2) = p;
      ecr(ke+2, j, ie+1) = p;
      ecr(ke+2, j, ie+2) = p;
    }
  }
  if (pmy_block->pbval->nblevel[0][0][1] < 0) {
    for (int i = is; i <= ie; i++) {
      Real p = 0.5*(ecr(ks-1, js, i) + ecr(ks, js-1, i));
      ecr(ks-1, js-1, i) = p;
      ecr(ks-1, js-2, i) = p;
      ecr(ks-2, js-1, i) = p;
      ecr(ks-2, js-2, i) = p;
    }
  }
  if (pmy_block->pbval->nblevel[0][2][1] < 0) {
    for (int i = is; i <= ie; i++) {
      Real p = 0.5*(ecr(ks-1, je, i) + ecr(ks, je+1, i));
      ecr(ks-1, je+1, i) = p;
      ecr(ks-1, je+2, i) = p;
      ecr(ks-2, je+1, i) = p;
      ecr(ks-2, je+2, i) = p;
    }
  }
  if (pmy_block->pbval->nblevel[2][0][1] < 0) {
    for (int i = is; i <= ie; i++) {
      Real p = 0.5*(ecr(ke+1, js, i) + ecr(ke, js-1, i));
      ecr(ke+1, js-1, i) = p;
      ecr(ke+1, js-2, i) = p;
      ecr(ke+2, js-1, i) = p;
      ecr(ke+2, js-2, i) = p;
    }
  }
  if (pmy_block->pbval->nblevel[2][2][1] < 0) {
    for (int i = is; i <= ie; i++) {
      Real p = 0.5*(ecr(ke+1, je, i) + ecr(ke, je+1, i));
      ecr(ke+1, je+1, i) = p;
      ecr(ke+1, je+2, i) = p;
      ecr(ke+2, je+1, i) = p;
      ecr(ke+2, je+2, i) = p;
    }
  }

  // fill corners
  if (pmy_block->pbval->nblevel[0][0][0] < 0) {
    Real p = (ecr(ks, js-1, is-1) + ecr(ks-1, js, is-1) + ecr(ks-1, js-1, is))/3.0;
    ecr(ks-1, js-1, is-1) = p;
    ecr(ks-1, js-1, is-2) = p;
    ecr(ks-1, js-2, is-1) = p;
    ecr(ks-1, js-2, is-2) = p;
    ecr(ks-2, js-1, is-1) = p;
    ecr(ks-2, js-1, is-2) = p;
    ecr(ks-2, js-2, is-1) = p;
    ecr(ks-2, js-2, is-2) = p;
  }
  if (pmy_block->pbval->nblevel[0][0][2] < 0) {
    Real p = (ecr(ks, js-1, ie+1) + ecr(ks-1, js, ie+1) + ecr(ks-1, js-1, ie))/3.0;
    ecr(ks-1, js-1, ie+1) = p;
    ecr(ks-1, js-1, ie+2) = p;
    ecr(ks-1, js-2, ie+1) = p;
    ecr(ks-1, js-2, ie+2) = p;
    ecr(ks-2, js-1, ie+1) = p;
    ecr(ks-2, js-1, ie+2) = p;
    ecr(ks-2, js-2, ie+1) = p;
    ecr(ks-2, js-2, ie+2) = p;
  }
  if (pmy_block->pbval->nblevel[0][2][0] < 0) {
    Real p = (ecr(ks, je+1, is-1) + ecr(ks-1, je, is-1) + ecr(ks-1, je+1, is))/3.0;
    ecr(ks-1, je+1, is-1) = p;
    ecr(ks-1, je+1, is-2) = p;
    ecr(ks-1, je+2, is-1) = p;
    ecr(ks-1, je+2, is-2) = p;
    ecr(ks-2, je+1, is-1) = p;
    ecr(ks-2, je+1, is-2) = p;
    ecr(ks-2, je+2, is-1) = p;
    ecr(ks-2, je+2, is-2) = p;
  }
  if (pmy_block->pbval->nblevel[2][0][0] < 0) {
    Real p = (ecr(ke, js-1, is-1) + ecr(ke+1, js, is-1) + ecr(ke+1, js-1, is))/3.0;
    ecr(ke+1, js-1, is-1) = p;
    ecr(ke+1, js-1, is-2) = p;
    ecr(ke+1, js-2, is-1) = p;
    ecr(ke+1, js-2, is-2) = p;
    ecr(ke+2, js-1, is-1) = p;
    ecr(ke+2, js-1, is-2) = p;
    ecr(ke+2, js-2, is-1) = p;
    ecr(ke+2, js-2, is-2) = p;
  }
  if (pmy_block->pbval->nblevel[0][2][2] < 0) {
    Real p = (ecr(ks, je+1, ie+1) + ecr(ks-1, je, ie+1) + ecr(ks-1, je+1, ie))/3.0;
    ecr(ks-1, je+1, ie+1) = p;
    ecr(ks-1, je+1, ie+2) = p;
    ecr(ks-1, je+2, ie+1) = p;
    ecr(ks-1, je+2, ie+2) = p;
    ecr(ks-2, je+1, ie+1) = p;
    ecr(ks-2, je+1, ie+2) = p;
    ecr(ks-2, je+2, ie+1) = p;
    ecr(ks-2, je+2, ie+2) = p;
  }
  if (pmy_block->pbval->nblevel[2][0][2] < 0) {
    Real p = (ecr(ke, js-1, ie+1) + ecr(ke+1, js, ie+1) + ecr(ke+1, js-1, ie))/3.0;
    ecr(ke+1, js-1, ie+1) = p;
    ecr(ke+1, js-1, ie+2) = p;
    ecr(ke+1, js-2, ie+1) = p;
    ecr(ke+1, js-2, ie+2) = p;
    ecr(ke+2, js-1, ie+1) = p;
    ecr(ke+2, js-1, ie+2) = p;
    ecr(ke+2, js-2, ie+1) = p;
    ecr(ke+2, js-2, ie+2) = p;
  }
  if (pmy_block->pbval->nblevel[2][2][0] < 0) {
    Real p = (ecr(ke, je+1, is-1) + ecr(ke+1, je, is-1) + ecr(ke+1, je+1, is))/3.0;
    ecr(ke+1, je+1, is-1) = p;
    ecr(ke+1, je+1, is-2) = p;
    ecr(ke+1, je+2, is-1) = p;
    ecr(ke+1, je+2, is-2) = p;
    ecr(ke+2, je+1, is-1) = p;
    ecr(ke+2, je+1, is-2) = p;
    ecr(ke+2, je+2, is-1) = p;
    ecr(ke+2, je+2, is-2) = p;
  }
  if (pmy_block->pbval->nblevel[2][2][2] < 0) {
    Real p = (ecr(ke, je+1, ie+1) + ecr(ke+1, je, ie+1) + ecr(ke+1, je+1, ie))/3.0;
    ecr(ke+1, je+1, ie+1) = p;
    ecr(ke+1, je+1, ie+2) = p;
    ecr(ke+1, je+2, ie+1) = p;
    ecr(ke+1, je+2, ie+2) = p;
    ecr(ke+2, je+1, ie+1) = p;
    ecr(ke+2, je+1, ie+2) = p;
    ecr(ke+2, je+2, ie+1) = p;
    ecr(ke+2, je+2, ie+2) = p;
  }

  return;
}

