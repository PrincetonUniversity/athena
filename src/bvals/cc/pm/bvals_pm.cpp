//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_pm.cpp
//! \brief implements boundary functions for Particle Mesh variables
//!   and utilities to add particle densities rather than replacce

// C headers

// C++ headers
#include <sstream>    // stringstream

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../orbital_advection/orbital_advection.hpp"
#include "../../../utils/buffer_utils.hpp"
#include "bvals_pm.hpp"

//----------------------------------------------------------------------------------------
//! \fn ParticleMeshBoundaryVariable::ParticleMeshBoundaryVariable
//! \brief

ParticleMeshBoundaryVariable::ParticleMeshBoundaryVariable(
    MeshBlock *pmb, AthenaArray<Real> *var, AthenaArray<Real> *coarse_var,
    AthenaArray<Real> *var_flux, ParticleMesh *ppm) :
    CellCenteredBoundaryVariable(pmb, var, coarse_var, var_flux, false),
    ppm_(ppm) {
  var_buf.NewAthenaArray(var->GetDim4(),pmb->ncells3,pmb->ncells2,pmb->ncells1);
}

//----------------------------------------------------------------------------------------
//! \fn int ParticleMeshBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
//!                                                             const NeighborBlock& nb)
//! \brief Set particle_mesh boundary buffers for sending to a block on the same level

int ParticleMeshBoundaryVariable::LoadBoundaryBufferSameLevel(Real *buf,
                                                              const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;

  // load ghost zones
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
  AthenaArray<Real> &var = *var_cc;
  BufferUtility::PackData(var, buf, nl_, nu_, si, ei, sj, ej, sk, ek, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticleMeshBoundaryVariable::SetBoundarySameLevel(Real *buf,
//!                                                      const NeighborBlock& nb)
//! \brief Unpack PM boundary received from a block on the same level
//!
//! Unpack received to var_buf
//! Add it to var_cc if non-shear-periodic
//! or pass it to SendShearingBoxBoundaryBuffers for shift

void ParticleMeshBoundaryVariable::SetBoundarySameLevel(Real *buf,
                                                 const NeighborBlock& nb) {
  MeshBlock *pmb = pmy_block_;
  int si, sj, sk, ei, ej, ek;
  si = (nb.ni.ox1 > 0) ? (pmb->ie - NGHOST + 1) : pmb->is;
  ei = (nb.ni.ox1 < 0) ? (pmb->is + NGHOST - 1) : pmb->ie;
  sj = (nb.ni.ox2 > 0) ? (pmb->je - NGHOST + 1) : pmb->js;
  ej = (nb.ni.ox2 < 0) ? (pmb->js + NGHOST - 1) : pmb->je;
  sk = (nb.ni.ox3 > 0) ? (pmb->ke - NGHOST + 1) : pmb->ks;
  ek = (nb.ni.ox3 < 0) ? (pmb->ks + NGHOST - 1) : pmb->ke;


  // unpack and add directly onto active zones
  int p = 0;

  AthenaArray<Real> &var = *var_cc;
  for (int n=nl_; n<=nu_; ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
        for (int i=si; i<=ei; ++i) {
          if (nb.shear) var_buf(n,k,j,i) += buf[p++];
          else
            var(n,k,j,i) += buf[p++];
        }
      }
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn int ParticleMeshBoundaryVariable::LoadShearingBoxBoundarySameLevel(
//!                                       AthenaArray<Real> &src, Real *buf, int nb)
//! \brief Set PM shearing boundary buffers for sending to a block on the same level

void ParticleMeshBoundaryVariable::LoadShearingBoxBoundarySameLevel(
                                 AthenaArray<Real> &src, Real *buf, int nb) {
 
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticleMeshBoundaryVariable::SendShearingBoxBoundaryBuffers()
//! \brief Send PM shearing box boundary buffers received

void ParticleMeshBoundaryVariable::SendShearingBoxBoundaryBuffers() {
 
  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool ParticleMeshBoundaryVariable::SetShearingBoxBoundaryBuffers()
//! \brief Add PM density received

void ParticleMeshBoundaryVariable::SetShearingBoxBoundaryBuffers() {

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticleMeshBoundaryVariable::ShearQuantities(AthenaArray<Real> &shear_cc_,
//!                                                   bool upper)
//! \brief Apply shear to ParticleMesh x2 momentum

void ParticleMeshBoundaryVariable::ShearQuantities(AthenaArray<Real> &shear_cc_,
                                                   bool upper) {
 
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ParticleMeshBoundaryVariable::ReceiveAndSetShearingBoxBoundariesWithWait()
//! \brief receive and set the boundary data for initialization

void ParticleMeshBoundaryVariable::ReceiveAndSetShearingBoxBoundariesWithWait() {

  return;
}
