//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro.cpp
//  \brief implementation of functions in class Hydro

// C headers

// C++ headers
#include <algorithm>
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../mesh/mesh.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "hydro.hpp"
#include "hydro_diffusion/hydro_diffusion.hpp"
#include "srcterms/hydro_srcterms.hpp"

// constructor, initializes data structures and parameters

Hydro::Hydro(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block(pmb), u(NHYDRO, pmb->ncells3, pmb->ncells2, pmb->ncells1),
    w(NHYDRO, pmb->ncells3, pmb->ncells2, pmb->ncells1),
    u1(NHYDRO, pmb->ncells3, pmb->ncells2, pmb->ncells1),
    w1(NHYDRO, pmb->ncells3, pmb->ncells2, pmb->ncells1),
    hbvar(pmb, &u, &coarse_cons_, flux, HydroBoundaryQuantity::cons),
    hsrc(this, pin),
    hdif(this, pin) {
  int ncells1 = pmb->ncells1, ncells2 = pmb->ncells2, ncells3 = pmb->ncells3;
  Mesh *pm = pmy_block->pmy_mesh;

  // Allocate optional memory primitive/conserved variable registers for time-integrator
  if (pmb->precon->xorder == 4) {
    // fourth-order hydro cell-centered approximations
    u_cc.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);
    w_cc.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);
  }

  // If user-requested time integrator is type 3S*, allocate additional memory registers
  std::string integrator = pin->GetOrAddString("time", "integrator", "vl2");
  if (integrator == "ssprk5_4" || STS_ENABLED)
    // future extension may add "int nregister" to Hydro class
    u2.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);

  flux[X1DIR].NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1+1);
  if (pm->f2_)
    flux[X2DIR].NewAthenaArray(NHYDRO, ncells3, ncells2+1, ncells1);
  if (pm->f3_)
    flux[X3DIR].NewAthenaArray(NHYDRO, ncells3+1, ncells2, ncells1);

  // allocate prolongation buffers
  if (pm->multilevel == true) {
    int ncc1 = pmb->block_size.nx1/2 + 2*NGHOST;
    int ncc2 = 1;
    if (pm->f2_) ncc2 = pmb->block_size.nx2/2 + 2*NGHOST;
    int ncc3 = 1;
    if (pm->f3_) ncc3 = pmb->block_size.nx3/2 + 2*NGHOST;
    coarse_cons_.NewAthenaArray(NHYDRO, ncc3, ncc2, ncc1);
    coarse_prim_.NewAthenaArray(NHYDRO, ncc3, ncc2, ncc1);
    // "Enroll" in S/AMR by adding to vector of tuples of pointers in MeshRefinement class
    pmy_block->pmr->AddToRefinement(&u, &coarse_cons_);
  }

  // enroll HydroBoundaryVariable object
  hbvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&hbvar);
  pmb->pbval->bvars_main_int.push_back(&hbvar);

  // Allocate memory for scratch arrays
  dt1_.NewAthenaArray(ncells1);
  dt2_.NewAthenaArray(ncells1);
  dt3_.NewAthenaArray(ncells1);
  dxw_.NewAthenaArray(ncells1);
  wl_.NewAthenaArray(NWAVE, ncells1);
  wr_.NewAthenaArray(NWAVE, ncells1);
  wlb_.NewAthenaArray(NWAVE, ncells1);
  x1face_area_.NewAthenaArray(ncells1+1);
  if (pm->f2_) {
    x2face_area_.NewAthenaArray(ncells1);
    x2face_area_p1_.NewAthenaArray(ncells1);
  }
  if (pm->f3_) {
    x3face_area_.NewAthenaArray(ncells1);
    x3face_area_p1_.NewAthenaArray(ncells1);
  }
  cell_volume_.NewAthenaArray(ncells1);
  dflx_.NewAthenaArray((NHYDRO), ncells1);
  if (MAGNETIC_FIELDS_ENABLED && RELATIVISTIC_DYNAMICS) { // only used in (SR/GR)MHD
    bb_normal_.NewAthenaArray(ncells1);
    lambdas_p_l_.NewAthenaArray(ncells1);
    lambdas_m_l_.NewAthenaArray(ncells1);
    lambdas_p_r_.NewAthenaArray(ncells1);
    lambdas_m_r_.NewAthenaArray(ncells1);
  }
  if (GENERAL_RELATIVITY) { // only used in GR
    g_.NewAthenaArray(NMETRIC, ncells1);
    gi_.NewAthenaArray(NMETRIC, ncells1);
    cons_.NewAthenaArray(NWAVE, ncells1);
  }
  // for one-time potential calcuation and correction (old Athena)
  if (SELF_GRAVITY_ENABLED == 3) {
    gflx[X1DIR].NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1+1);
    if (pm->f2_)
      gflx[X2DIR].NewAthenaArray(NHYDRO, ncells3, ncells2+1, ncells1);
    if (pm->f3_)
      gflx[X3DIR].NewAthenaArray(NHYDRO, ncells3+1, ncells2, ncells1);

    gflx_old[X1DIR].NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1+1);
    if (pm->f2_)
      gflx_old[X2DIR].NewAthenaArray(NHYDRO, ncells3, ncells2+1, ncells1);
    if (pm->f3_)
      gflx_old[X3DIR].NewAthenaArray(NHYDRO, ncells3+1, ncells2, ncells1);
  }

  // fourth-order hydro
  // 4D scratch arrays
  wl3d_.NewAthenaArray(NWAVE, ncells3, ncells2, ncells1);
  wr3d_.NewAthenaArray(NWAVE, ncells3, ncells2, ncells1);
  scr1_nkji_.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);
  scr2_nkji_.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);
  laplacian_l_fc_.NewAthenaArray(ncells1);
  laplacian_r_fc_.NewAthenaArray(ncells1);

  UserTimeStep_ = pmb->pmy_mesh->UserTimeStep_;
}

//----------------------------------------------------------------------------------------
//! \fn Real Hydro::GetWeightForCT(Real dflx, Real rhol, Real rhor, Real dx, Real dt)
//  \brief Calculate the weighting factor for the constrained transport method

Real Hydro::GetWeightForCT(Real dflx, Real rhol, Real rhor, Real dx, Real dt) {
  Real v_over_c = (1024.0)* dt * dflx / (dx * (rhol + rhor));
  Real tmp_min = std::min(static_cast<Real>(0.5),v_over_c);
  return 0.5 + std::max(static_cast<Real>(-0.5),tmp_min);
}
