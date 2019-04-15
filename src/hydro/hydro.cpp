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
#include "hydro.hpp"
#include "hydro_diffusion/hydro_diffusion.hpp"
#include "srcterms/hydro_srcterms.hpp"

// constructor, initializes data structures and parameters

Hydro::Hydro(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block = pmb;

  // Allocate memory for primitive/conserved variables
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  // Allocate memory registers for primitive/conserved variables for time-integrator
  u.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);
  pmb->RegisterMeshBlockData(u);
  w.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);
  u1.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);
  w1.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);

  // fourth-order hydro cell-centered approximations
  u_cc.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);
  w_cc.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);

  // If user-requested time integrator is type 3S*, allocate additional memory registers
  std::string integrator = pin->GetOrAddString("time", "integrator", "vl2");
  if (integrator == "ssprk5_4" || STS_ENABLED)
    // future extension may add "int nregister" to Hydro class
    u2.NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1);

  flux[X1DIR].NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1+1);
  if (pmy_block->block_size.nx2 > 1)
    flux[X2DIR].NewAthenaArray(NHYDRO, ncells3, ncells2+1, ncells1);
  if (pmy_block->block_size.nx3 > 1)
    flux[X3DIR].NewAthenaArray(NHYDRO, ncells3+1, ncells2, ncells1);

  // allocate prolongation buffers
  if (pmy_block->pmy_mesh->multilevel == true) {
    int ncc1 = pmb->block_size.nx1/2 + 2*NGHOST;
    int ncc2 = 1;
    if (pmb->block_size.nx2 > 1) ncc2 = pmb->block_size.nx2/2 + 2*NGHOST;
    int ncc3 = 1;
    if (pmb->block_size.nx3 > 1) ncc3 = pmb->block_size.nx3/2 + 2*NGHOST;
    coarse_cons_.NewAthenaArray(NHYDRO, ncc3, ncc2, ncc1);
    coarse_prim_.NewAthenaArray(NHYDRO, ncc3, ncc2, ncc1);
    // "Enroll" in S/AMR by adding to vector of tuples of pointers in MeshRefinement class
    pmy_block->pmr->AddToRefinement(&u, &coarse_cons_);
  }
  // create object to interface with BoundaryValues
  phbval  = new HydroBoundaryVariable(pmy_block, &u, &coarse_cons_, flux,
                                      HydroBoundaryQuantity::cons);
  phbval->bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(phbval);
  pmb->pbval->bvars_main_int.push_back(phbval);

  // Allocate memory for scratch arrays
  dt1_.NewAthenaArray(ncells1);
  dt2_.NewAthenaArray(ncells1);
  dt3_.NewAthenaArray(ncells1);
  dxw_.NewAthenaArray(ncells1);
  wl_.NewAthenaArray(NWAVE, ncells1);
  wr_.NewAthenaArray(NWAVE, ncells1);
  wlb_.NewAthenaArray(NWAVE, ncells1);
  x1face_area_.NewAthenaArray(ncells1+1);
  if (pmy_block->block_size.nx2 > 1) {
    x2face_area_.NewAthenaArray(ncells1);
    x2face_area_p1_.NewAthenaArray(ncells1);
  }
  if (pmy_block->block_size.nx3 > 1) {
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
    if (pmy_block->block_size.nx2 > 1)
      gflx[X2DIR].NewAthenaArray(NHYDRO, ncells3, ncells2+1, ncells1);
    if (pmy_block->block_size.nx3 > 1)
      gflx[X3DIR].NewAthenaArray(NHYDRO, ncells3+1, ncells2, ncells1);

    gflx_old[X1DIR].NewAthenaArray(NHYDRO, ncells3, ncells2, ncells1+1);
    if (pmy_block->block_size.nx2 > 1)
      gflx_old[X2DIR].NewAthenaArray(NHYDRO, ncells3, ncells2+1, ncells1);
    if (pmy_block->block_size.nx3 > 1)
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

  // Construct ptrs to objects of various classes needed to integrate hydro/MHD eqns
  psrc  = new HydroSourceTerms(this, pin);

  // ptr to diffusion object
  phdif = new HydroDiffusion(this,pin);
}

// destructor

Hydro::~Hydro() {
  delete psrc;
  delete phdif;
  delete phbval;
}


//----------------------------------------------------------------------------------------
//! \fn Real Hydro::GetWeightForCT(Real dflx, Real rhol, Real rhor, Real dx, Real dt)
//  \brief Calculate the weighting factor for the constrained transport method

Real Hydro::GetWeightForCT(Real dflx, Real rhol, Real rhor, Real dx, Real dt) {
  Real v_over_c = (1024.0)* dt * dflx / (dx * (rhol + rhor));
  Real tmp_min = std::min(static_cast<Real>(0.5),v_over_c);
  return 0.5 + std::max(static_cast<Real>(-0.5),tmp_min);
}
