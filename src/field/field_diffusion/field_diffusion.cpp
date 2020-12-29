//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file field_diffusion.cpp
//! \brief Class to implement diffusion processes in the induction equations

// C headers

// C++ headers
#include <algorithm>  // min()
#include <cmath>      // sqrt(), fabs()
#include <cstring>    // strcmp()
#include <iostream>   // endl
#include <limits>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../field.hpp"
#include "field_diffusion.hpp"

//! FieldDiffusion constructor

FieldDiffusion::FieldDiffusion(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block(pmb), field_diffusion_defined(false) {
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;

  // Check if field diffusion
  eta_ohm = pin->GetOrAddReal("problem", "eta_ohm", 0.0);
  eta_hall = pin->GetOrAddReal("problem", "eta_hall", 0.0);
  eta_ad = pin->GetOrAddReal("problem", "eta_ad", 0.0);

  if ((eta_ohm != 0.0) || (eta_hall != 0.0) || (eta_ad != 0.0)) {
    field_diffusion_defined = true;
    // Allocate memory for scratch vectors
    etaB.NewAthenaArray(3, nc3, nc2, nc1);
    e_oa.x1e.NewAthenaArray(nc3+1, nc2+1, nc1);
    e_oa.x2e.NewAthenaArray(nc3+1, nc2, nc1+1);
    e_oa.x3e.NewAthenaArray(nc3, nc2+1, nc1+1);
    e_h.x1e.NewAthenaArray(nc3+1, nc2+1, nc1);
    e_h.x2e.NewAthenaArray(nc3+1, nc2, nc1+1);
    e_h.x3e.NewAthenaArray(nc3, nc2+1, nc1+1);
    pflux.x1f.NewAthenaArray(nc3, nc2, nc1+1);
    pflux.x2f.NewAthenaArray(nc3, nc2+1, nc1);
    pflux.x3f.NewAthenaArray(nc3+1, nc2, nc1);

    jfx.NewAthenaArray(3, nc3+1, nc2+1, nc1+1);
    jfy.NewAthenaArray(3, nc3+1, nc2+1, nc1+1);
    jfz.NewAthenaArray(3, nc3+1, nc2+1, nc1+1);
    jcc.NewAthenaArray(3, nc3+1, nc2+1, nc1+1);

    eta_tot_.NewAthenaArray(nc1);
    bmag_.NewAthenaArray(nc3, nc2, nc1);

    jedge_.x1e.NewAthenaArray(nc3+1, nc2+1, nc1);
    jedge_.x2e.NewAthenaArray(nc3+1, nc2, nc1+1);
    jedge_.x3e.NewAthenaArray(nc3, nc2+1, nc1+1);

    face_area_.NewAthenaArray(nc1);
    face_area_p1_.NewAthenaArray(nc1);
    edge_length_.NewAthenaArray(nc1);
    edge_length_m1_.NewAthenaArray(nc1);
    cell_volume_.NewAthenaArray(nc1);
    dx1_.NewAthenaArray(nc1);
    dx2_.NewAthenaArray(nc1);
    dx3_.NewAthenaArray(nc1);

    if (pmb->pmy_mesh->FieldDiffusivity_ == nullptr)
      CalcMagDiffCoeff_ = ConstDiffusivity;
    else
      CalcMagDiffCoeff_ = pmb->pmy_mesh->FieldDiffusivity_;
  }

  if (field_diffusion_defined && RELATIVISTIC_DYNAMICS) {
    std::stringstream msg;
    msg << "### FATAL ERROR in FieldDiffusion" << std::endl
        << "Diffusion is incompatibile with relativistic dynamics" << std::endl;
    ATHENA_ERROR(msg);
  }
}


//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::CalcDiffusionEMF
//! \brief Calculate diffusion EMF(Ohmic & Ambipolar for now)

void FieldDiffusion::CalcDiffusionEMF(FaceField &bi, const AthenaArray<Real> &bc,
                                      EdgeField &e) {
  Field *pf = pmy_block->pfield;
  Hydro *ph = pmy_block->phydro;
  // Mesh  *pm = pmy_block->pmy_mesh; // unused variable

  if ((eta_ohm == 0.0) && (eta_ad == 0.0)) return;

  SetDiffusivity(ph->w, pf->bcc);

  CalcCurrent(bi);
  ClearEMF(e_oa);
  if (eta_ohm != 0.0) OhmicEMF(bi, bc, e_oa);
  if (eta_ad != 0.0) AmbipolarEMF(bi, bc, e_oa);

  // calculate the Poynting flux pflux and add to energy flux in Hydro class
  if (NON_BAROTROPIC_EOS) PoyntingFlux(e_oa, bc);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::AddEMF(EdgeField &e_src, EdgeField &e_des)
//! \brief Add source EMF to destination EMF

void FieldDiffusion::AddEMF(const EdgeField &e_src, EdgeField &e_des) {
  int size1 = e_src.x1e.GetSize();
  int size2 = e_src.x2e.GetSize();
  int size3 = e_src.x3e.GetSize();

#pragma omp simd
  for (int i=0; i<size1; ++i)
    e_des.x1e(i) += e_src.x1e(i);

#pragma omp simd
  for (int i=0; i<size2; ++i)
    e_des.x2e(i) += e_src.x2e(i);

#pragma omp simd
  for (int i=0; i<size3; ++i)
    e_des.x3e(i) += e_src.x3e(i);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::ClearEMF(EdgeField &e)
//! \brief Clear EMF

// TODO(felker): move out of FieldDiffusion class. Completely general operation
void FieldDiffusion::ClearEMF(EdgeField &e) {
  e.x1e.ZeroClear();
  e.x2e.ZeroClear();
  e.x3e.ZeroClear();
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::SetDiffusivity
//! \brief Set magnetic diffusion coefficients

void FieldDiffusion::SetDiffusivity(const AthenaArray<Real> &w,
                                    const AthenaArray<Real> &bc) {
  MeshBlock *pmb = pmy_block;
  int il = pmb->is - NGHOST; int jl = pmb->js; int kl = pmb->ks;
  int iu = pmb->ie + NGHOST; int ju = pmb->je; int ku = pmb->ke;
  if (pmb->pmy_mesh->f2) {
    jl -= NGHOST; ju += NGHOST;
  }
  if (pmb->pmy_mesh->f3) {
    kl -= NGHOST; ku += NGHOST;
  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        Real Bsq = SQR(bc(IB1,k,j,i)) + SQR(bc(IB2,k,j,i)) + SQR(bc(IB3,k,j,i));
        bmag_(k,j,i) = std::sqrt(Bsq);
      }
    }
  }
  // set diffusivities
  CalcMagDiffCoeff_(this, pmb, w, bmag_, il, iu, jl, ju, kl, ku);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::AddPoyntingFlux
//! \brief Add Poynting flux to the hydro energy flux

void FieldDiffusion::AddPoyntingFlux(FaceField &p_src) {
  MeshBlock *pmb = pmy_block;
  AthenaArray<Real> &x1flux = pmb->phydro->flux[X1DIR];
  AthenaArray<Real> &x2flux = pmb->phydro->flux[X2DIR];
  AthenaArray<Real> &x3flux = pmb->phydro->flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> &f1 = p_src.x1f, &f2 = p_src.x2f, &f3 = p_src.x3f;

  // 1D update:
  if (pmb->block_size.nx2 == 1) {
#pragma omp simd
    for (int i=is; i<=ie+1; ++i) {
      x1flux(IEN,ks,js,i) += f1(ks,js,i);
    }
    return;
  }

  // 2D update:
  if (pmb->block_size.nx3 == 1) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        x1flux(IEN,ks,j,i) += f1(ks,j,i);
        x2flux(IEN,ks,j,i) += f2(ks,j,i);
      }
    }
    return;
  }

  // 3D update:
  for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        x1flux(IEN,k,j,i) += f1(k,j,i);
        x2flux(IEN,k,j,i) += f2(k,j,i);
        x3flux(IEN,k,j,i) += f3(k,j,i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::NewDiffusionDt
//! \brief Get the non-ideal MHD timestep

void FieldDiffusion::NewDiffusionDt(Real &dt_oa, Real &dt_h) {
  MeshBlock *pmb = pmy_block;
  const bool f2 = pmb->pmy_mesh->f2;
  const bool f3 = pmb->pmy_mesh->f3;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  Real fac_oa, fac_h;
  if (f3 ) { // 3D
    fac_oa = 1.0/6.0;
  } else {
    if (f2) { // 2D
      fac_oa = 0.25;
    } else { // 1D
      fac_oa = 0.5;
    }
  }

  if (f2)
    fac_h = 1.0;
  else
    fac_h = 0.5;
  Real real_max = std::numeric_limits<Real>::max();
  dt_oa = real_max;
  dt_h  = real_max;

  AthenaArray<Real> &eta_t = eta_tot_;
  AthenaArray<Real> &len = dx1_, &dx2 = dx2_, &dx3 = dx3_;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        eta_t(i) = 0.0;
      }
      if (eta_ohm > 0.0) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          eta_t(i) += etaB(ohmic,k,j,i);
        }
      }
      if (eta_ad > 0.0) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          eta_t(i) += etaB(DiffProcess::ambipolar,k,j,i);
        }
      }
      pmb->pcoord->CenterWidth1(k, j, is, ie, len);
      pmb->pcoord->CenterWidth2(k, j, is, ie, dx2);
      pmb->pcoord->CenterWidth3(k, j, is, ie, dx3);
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        len(i) = (f2) ? std::min(len(i), dx2(i)):len(i);
        len(i) = (f3) ? std::min(len(i), dx3(i)):len(i);
      }
      if ((eta_ohm > 0.0) || (eta_ad > 0.0)) {
        for (int i=is; i<=ie; ++i)
          dt_oa = std::min(dt_oa, static_cast<Real>(
              fac_oa*SQR(len(i)) / (eta_t(i) + TINY_NUMBER)));
      }
      if (eta_hall > 0.0) {
        for (int i=is; i<=ie; ++i)
          dt_h = std::min(dt_h, static_cast<Real>(
              fac_h*SQR(len(i)) / (std::abs(etaB(hall,k,j,i)) + TINY_NUMBER)));
      }
    }
  }
  return;
}
