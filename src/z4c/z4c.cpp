//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file z4c.cpp
//  \brief implementation of functions in the Z4c class

#include <iostream>
#include <fstream>

// Athena++ headers
#include "z4c.hpp"
#include "z4c_macro.hpp"
#include "../athena.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"

#define SQ(X) ((X)*(X))

// constructor, initializes data structures and parameters

char const * const Z4c::Z4c_names[Z4c::N_Z4c] = {
  "z4c.chi",
  "z4c.gxx", "z4c.gxy", "z4c.gxz", "z4c.gyy", "z4c.gyz", "z4c.gzz",
  "z4c.Khat",
  "z4c.Axx", "z4c.Axy", "z4c.Axz", "z4c.Ayy", "z4c.Ayz", "z4c.Azz",
  "z4c.Gamx", "z4c.Gamy", "z4c.Gamz",
  "z4c.Theta",
  "z4c.alpha",
  "z4c.betax", "z4c.betay", "z4c.betaz",
};

char const * const Z4c::ADM_names[Z4c::N_ADM] = {
  "adm.gxx", "adm.gxy", "adm.gxz", "adm.gyy", "adm.gyz", "adm.gzz",
  "adm.Kxx", "adm.Kxy", "adm.Kxz", "adm.Kyy", "adm.Kyz", "adm.Kzz",
  "adm.psi4",
};

char const * const Z4c::Constraint_names[Z4c::N_CON] = {
  "con.C",
  "con.H",
  "con.M",
  "con.Z",
  "con.Mx", "con.My", "con.Mz",
};

char const * const Z4c::Matter_names[Z4c::N_MAT] = {
  "mat.rho",
  "mat.Sx", "mat.Sy", "mat.Sz",
  "mat.Sxx", "mat.Sxy", "mat.Sxz", "mat.Syy", "mat.Syz", "mat.Szz",
};

Z4c::Z4c(MeshBlock *pmb, ParameterInput *pin) :
  pmy_block(pmb),
  empty_flux{AthenaArray<Real>(), AthenaArray<Real>(), AthenaArray<Real>()},
  coarse_u_(N_Z4c, pmb->ncc3, pmb->ncc2, pmb->ncc1,
            (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
             AthenaArray<Real>::DataStatus::empty)),
  //u(N_Z4c, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  storage{{N_Z4c, pmb->ncells3, pmb->ncells2, pmb->ncells1}, // u
          {N_Z4c, pmb->ncells3, pmb->ncells2, pmb->ncells1}, // u1
          {N_Z4c, pmb->ncells3, pmb->ncells2, pmb->ncells1}, // u2
          {N_Z4c, pmb->ncells3, pmb->ncells2, pmb->ncells1}, // rhs
          {N_Z4c, pmb->ncells3, pmb->ncells2, pmb->ncells1}, // rhs1
          {N_ADM, pmb->ncells3, pmb->ncells2, pmb->ncells1}, // adm
          {N_CON, pmb->ncells3, pmb->ncells2, pmb->ncells1}, // con
          {N_MAT, pmb->ncells3, pmb->ncells2, pmb->ncells1}, // mat
  },
  ubvar(pmb, &storage.u, &coarse_u_, empty_flux)
  /*
  u(N_Z4c, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  ubvar(pmb, &storage->u, &coarse_u_, empty_flux)
  */
{
  // pmy_block = pmb;
  Coordinates * pco = pmb->pcoord;

  // Allocate memory for the solution and its time derivative
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  Mesh *pm = pmy_block->pmy_mesh;

  if(pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if(pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  // // inform MeshBlock that this array is the "primary" representation
  // // Used for:
  // // (1) load-balancing
  // // (2) (future) dumping to restart file
  // u = storage.u;
  pmb->RegisterMeshBlockData(storage.u);

  // "Enroll" in SMR/AMR by adding to vector of pointers in MeshRefinement class
  if (pm->multilevel) {
    refinement_idx = pmy_block->pmr->AddToRefinement(&storage.u, &coarse_u_);
  }


  /*
  // Allocate memory for the solution
  storage.u.NewAthenaArray(N_Z4c, ncells3, ncells2, ncells1);
  storage.u1.NewAthenaArray(N_Z4c, ncells3, ncells2, ncells1);
  // If user-requested time integrator is type 3S*, allocate additional memory registers
  std::string integrator = pin->GetOrAddString("time","integrator","vl2");
  if (integrator == "ssprk5_4") storage.u2.NewAthenaArray(N_Z4c, ncells3, ncells2, ncells1);
  storage.rhs.NewAthenaArray(N_Z4c, ncells3, ncells2, ncells1);
  //DEBUG Necesary for Traditional RK4
  storage.rhs1.NewAthenaArray(N_Z4c, ncells3, ncells2, ncells1);
  //ENDDEBUG
  storage.adm.NewAthenaArray(N_ADM, ncells3, ncells2, ncells1);
  storage.con.NewAthenaArray(N_CON, ncells3, ncells2, ncells1);
  storage.mat.NewAthenaArray(N_MAT, ncells3, ncells2, ncells1);
  */
  // enroll CellCenteredBoundaryVariable object
  ubvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&ubvar);
  pmb->pbval->bvars_main_int.push_back(&ubvar);

  dt1_.NewAthenaArray(ncells1);
  dt2_.NewAthenaArray(ncells1);
  dt3_.NewAthenaArray(ncells1);

  // Parameters
  opt.chi_psi_power = pin->GetOrAddReal("z4c", "chi_psi_power", -4.0);
  opt.chi_div_floor = pin->GetOrAddReal("z4c", "chi_div_floor", -1000.0);
  opt.diss = pin->GetOrAddReal("z4c", "diss", 0.0);
  opt.eps_floor = pin->GetOrAddReal("z4c", "eps_floor", 1e-12);
  opt.damp_kappa1 = pin->GetOrAddReal("z4c", "damp_kappa1", 0.0);
  opt.damp_kappa2 = pin->GetOrAddReal("z4c", "damp_kappa2", 0.0);
  // Gauge conditions (default to moving puncture gauge)
  opt.lapse_oplog = pin->GetOrAddReal("z4c", "lapse_oplog", 1.0);
  opt.lapse_harmonicf = pin->GetOrAddReal("z4c", "lapse_harmonicf", 2.0);
  opt.lapse_harmonic = pin->GetOrAddReal("z4c", "lapse_harmonic", 0.0);
  opt.lapse_advect = pin->GetOrAddReal("z4c", "lapse_advect", 1.0);
  opt.shift_advect = pin->GetOrAddReal("z4c", "shift_advect", 1.0);
  opt.shift_eta = pin->GetOrAddReal("z4c", "shift_eta", 0.0);
  // Single puncture parameters
  opt.punc_ADM_mass = pin->GetOrAddReal("z4c", "punc_ADM_mass", 1.0);
  // AwA parameters (default to linear wave test)
  opt.AwA_amplitude = pin->GetOrAddReal("z4c", "AwA_amplitude", 1e-10);
  opt.AwA_d_x = pin->GetOrAddReal("z4c", "AwA_d_x", 1.0);
  opt.AwA_d_y = pin->GetOrAddReal("z4c", "AwA_d_y", 1.0);

  // Set aliases
  SetADMAliases(storage.adm, adm);
  SetConstraintAliases(storage.con, con);
  SetMatterAliases(storage.mat, mat);
  SetZ4cAliases(storage.rhs, rhs);
  SetZ4cAliases(storage.u, z4c);

  // Allocate memory for aux 1D vars
  r.NewAthenaTensor(ncells1);
  detg.NewAthenaTensor(ncells1);
  chi_guarded.NewAthenaTensor(ncells1);
  oopsi4.NewAthenaTensor(ncells1);
  A.NewAthenaTensor(ncells1);
  AA.NewAthenaTensor(ncells1);
  R.NewAthenaTensor(ncells1);
  Ht.NewAthenaTensor(ncells1);
  K.NewAthenaTensor(ncells1);
  KK.NewAthenaTensor(ncells1);
  Ddalpha.NewAthenaTensor(ncells1);
  S.NewAthenaTensor(ncells1);
  M_u.NewAthenaTensor(ncells1);
  Gamma_u.NewAthenaTensor(ncells1);
  DA_u.NewAthenaTensor(ncells1);
  s_u.NewAthenaTensor(ncells1);
  g_uu.NewAthenaTensor(ncells1);
  A_uu.NewAthenaTensor(ncells1);
  AA_dd.NewAthenaTensor(ncells1);
  R_dd.NewAthenaTensor(ncells1);
  Rphi_dd.NewAthenaTensor(ncells1);
  Kt_dd.NewAthenaTensor(ncells1);
  K_ud.NewAthenaTensor(ncells1);
  Ddalpha_dd.NewAthenaTensor(ncells1);
  Ddphi_dd.NewAthenaTensor(ncells1);
  Gamma_ddd.NewAthenaTensor(ncells1);
  Gamma_udd.NewAthenaTensor(ncells1);
  DK_ddd.NewAthenaTensor(ncells1);
  DK_udd.NewAthenaTensor(ncells1);

  dbeta.NewAthenaTensor(ncells1);
  dalpha_d.NewAthenaTensor(ncells1);
  ddbeta_d.NewAthenaTensor(ncells1);
  dchi_d.NewAthenaTensor(ncells1);
  dphi_d.NewAthenaTensor(ncells1);
  dK_d.NewAthenaTensor(ncells1);
  dKhat_d.NewAthenaTensor(ncells1);
  dTheta_d.NewAthenaTensor(ncells1);
  ddalpha_dd.NewAthenaTensor(ncells1);
  dbeta_du.NewAthenaTensor(ncells1);
  ddchi_dd.NewAthenaTensor(ncells1);
  dGam_du.NewAthenaTensor(ncells1);
  dg_ddd.NewAthenaTensor(ncells1);
  dg_duu.NewAthenaTensor(ncells1);
  dK_ddd.NewAthenaTensor(ncells1);
  dA_ddd.NewAthenaTensor(ncells1);
  ddbeta_ddu.NewAthenaTensor(ncells1);
  ddg_dddd.NewAthenaTensor(ncells1);

  Lchi.NewAthenaTensor(ncells1);
  LKhat.NewAthenaTensor(ncells1);
  LTheta.NewAthenaTensor(ncells1);
  Lalpha.NewAthenaTensor(ncells1);
  LGam_u.NewAthenaTensor(ncells1);
  Lbeta_u.NewAthenaTensor(ncells1);
  Lg_dd.NewAthenaTensor(ncells1);
  LA_dd.NewAthenaTensor(ncells1);

  // Setup finite differencing kernel
  // NOTE: this will need to be changed if the Z4c variables become vertex center
  FD.stride[0] = 1;
  FD.stride[1] = 0;
  FD.stride[2] = 0;
  FD.idx[0] = 1.0/pco->dx1v(0);
  FD.idx[1] = 0.;
  FD.idx[2] = 0.;
  if(ncells2 > 1) {
    FD.stride[1] = ncells1;
    FD.idx[1] = 1.0/pco->dx2v(0);
  }
  if(ncells3 > 1) {
    FD.stride[2] = ncells2*ncells1;
    FD.idx[2] = 1.0/pco->dx3v(0);
  }
  FD.diss = opt.diss*pow(2, -2*NGHOST)*(NGHOST % 2 == 0 ? -1 : 1);
}

// destructor

Z4c::~Z4c()
{
  storage.u.DeleteAthenaArray();
  storage.u1.DeleteAthenaArray();
  storage.u2.DeleteAthenaArray();
  storage.rhs.DeleteAthenaArray();
  //DEBUG
  storage.rhs1.DeleteAthenaArray();
  //ENDDEBUG
  storage.adm.DeleteAthenaArray();
  storage.con.DeleteAthenaArray();
  storage.mat.DeleteAthenaArray();

  dt1_.DeleteAthenaArray();
  dt2_.DeleteAthenaArray();
  dt3_.DeleteAthenaArray();

  r.DeleteAthenaTensor();
  detg.DeleteAthenaTensor();
  chi_guarded.DeleteAthenaTensor();
  oopsi4.DeleteAthenaTensor();
  A.DeleteAthenaTensor();
  AA.DeleteAthenaTensor();
  R.DeleteAthenaTensor();
  Ht.DeleteAthenaTensor();
  K.DeleteAthenaTensor();
  KK.DeleteAthenaTensor();
  Ddalpha.DeleteAthenaTensor();
  S.DeleteAthenaTensor();
  M_u.DeleteAthenaTensor();
  Gamma_u.DeleteAthenaTensor();
  DA_u.DeleteAthenaTensor();
  s_u.DeleteAthenaTensor();
  g_uu.DeleteAthenaTensor();
  A_uu.DeleteAthenaTensor();
  AA_dd.DeleteAthenaTensor();
  R_dd.DeleteAthenaTensor();
  Rphi_dd.DeleteAthenaTensor();
  Kt_dd.DeleteAthenaTensor();
  K_ud.DeleteAthenaTensor();
  Ddalpha_dd.DeleteAthenaTensor();
  Ddphi_dd.DeleteAthenaTensor();
  Gamma_ddd.DeleteAthenaTensor();
  Gamma_udd.DeleteAthenaTensor();
  DK_ddd.DeleteAthenaTensor();
  DK_udd.DeleteAthenaTensor();

  dbeta.DeleteAthenaTensor();
  dalpha_d.DeleteAthenaTensor();
  ddbeta_d.DeleteAthenaTensor();
  dchi_d.DeleteAthenaTensor();
  dphi_d.DeleteAthenaTensor();
  dK_d.DeleteAthenaTensor();
  dKhat_d.DeleteAthenaTensor();
  dTheta_d.DeleteAthenaTensor();
  ddalpha_dd.DeleteAthenaTensor();
  dbeta_du.DeleteAthenaTensor();
  ddchi_dd.DeleteAthenaTensor();
  dGam_du.DeleteAthenaTensor();
  dg_ddd.DeleteAthenaTensor();
  dg_duu.DeleteAthenaTensor();
  dK_ddd.DeleteAthenaTensor();
  dA_ddd.DeleteAthenaTensor();
  ddbeta_ddu.DeleteAthenaTensor();
  ddg_dddd.DeleteAthenaTensor();

  Lchi.DeleteAthenaTensor();
  LKhat.DeleteAthenaTensor();
  LTheta.DeleteAthenaTensor();
  Lalpha.DeleteAthenaTensor();
  LGam_u.DeleteAthenaTensor();
  Lbeta_u.DeleteAthenaTensor();
  Lg_dd.DeleteAthenaTensor();
  LA_dd.DeleteAthenaTensor();
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::SetADMAliases(AthenaArray<Real> & u, ADM_vars & z4c)
// \brief Set ADM aliases

void Z4c::SetADMAliases(AthenaArray<Real> & u_adm, Z4c::ADM_vars & adm)
{
  adm.psi4.InitWithShallowSlice(u_adm, I_ADM_psi4);
  adm.g_dd.InitWithShallowSlice(u_adm, I_ADM_gxx);
  adm.K_dd.InitWithShallowSlice(u_adm, I_ADM_Kxx);
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::SetConstraintAliases(AthenaArray<Real> & u, ADM_vars & z4c)
// \brief Set ADM aliases

void Z4c::SetConstraintAliases(AthenaArray<Real> & u_con, Z4c::Constraint_vars & con)
{
  con.C.InitWithShallowSlice(u_con, I_CON_C);
  con.H.InitWithShallowSlice(u_con, I_CON_H);
  con.M.InitWithShallowSlice(u_con, I_CON_M);
  con.Z.InitWithShallowSlice(u_con, I_CON_Z);
  con.M_d.InitWithShallowSlice(u_con, I_CON_Mx);
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::SetMatterAliases(AthenaArray<Real> & u_mat, Matter_vars & mat)
// \brief Set matter aliases

void Z4c::SetMatterAliases(AthenaArray<Real> & u_mat, Z4c::Matter_vars & mat)
{
  mat.rho.InitWithShallowSlice(u_mat, I_MAT_rho);
  mat.S_d.InitWithShallowSlice(u_mat, I_MAT_Sx);
  mat.S_dd.InitWithShallowSlice(u_mat, I_MAT_Sxx);
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::SetZ4cAliases(AthenaArray<Real> & u, Z4c_vars & z4c)
// \brief Set Z4c aliases

void Z4c::SetZ4cAliases(AthenaArray<Real> & u, Z4c::Z4c_vars & z4c)
{
  z4c.chi.InitWithShallowSlice(u, I_Z4c_chi);
  z4c.Khat.InitWithShallowSlice(u, I_Z4c_Khat);
  z4c.Theta.InitWithShallowSlice(u, I_Z4c_Theta);
  z4c.alpha.InitWithShallowSlice(u, I_Z4c_alpha);
  z4c.Gam_u.InitWithShallowSlice(u, I_Z4c_Gamx);
  z4c.beta_u.InitWithShallowSlice(u, I_Z4c_betax);
  z4c.g_dd.InitWithShallowSlice(u, I_Z4c_gxx);
  z4c.A_dd.InitWithShallowSlice(u, I_Z4c_Axx);
}

//----------------------------------------------------------------------------------------
// \!fn Real Z4c::SpatialDet(Real gxx, ... , Real gzz)
// \brief returns determinant of 3-metric

Real Z4c::SpatialDet(Real const gxx, Real const gxy, Real const gxz,
                     Real const gyy, Real const gyz, Real const gzz)
{
  return - SQ(gxz)*gyy + 2*gxy*gxz*gyz - gxx*SQ(gyz) - SQ(gxy)*gzz + gxx*gyy*gzz;
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::SpatialInv(Real const detginv,
//           Real const gxx, Real const gxy, Real const gxz,
//           Real const gyy, Real const gyz, Real const gzz,
//           Real * uxx, Real * uxy, Real * uxz,
//           Real * uyy, Real * uyz, Real * uzz)
// \brief returns inverse of 3-metric

void Z4c::SpatialInv(Real const detginv,
                     Real const gxx, Real const gxy, Real const gxz,
                     Real const gyy, Real const gyz, Real const gzz,
                     Real * uxx, Real * uxy, Real * uxz,
                     Real * uyy, Real * uyz, Real * uzz)
{
  *uxx = (-SQ(gyz) + gyy*gzz)*detginv;
  *uxy = (gxz*gyz  - gxy*gzz)*detginv;
  *uyy = (-SQ(gxz) + gxx*gzz)*detginv;
  *uxz = (-gxz*gyy + gxy*gyz)*detginv;
  *uyz = (gxy*gxz  - gxx*gyz)*detginv;
  *uzz = (-SQ(gxy) + gxx*gyy)*detginv;
  return;
}

//----------------------------------------------------------------------------------------
// \!fn Real Z4c::Trace(Real detginv, Real gxx, ... , Real gzz, Real Axx, ..., Real Azz)
// \brief returns Trace of extrinsic curvature

Real Z4c::Trace(Real const detginv,
                Real const gxx, Real const gxy, Real const gxz,
                Real const gyy, Real const gyz, Real const gzz,
                Real const Axx, Real const Axy, Real const Axz,
                Real const Ayy, Real const Ayz, Real const Azz)
{
  return (detginv*(
       - 2.*Ayz*gxx*gyz + Axx*gyy*gzz +  gxx*(Azz*gyy + Ayy*gzz)
       + 2.*(gxz*(Ayz*gxy - Axz*gyy + Axy*gyz) + gxy*(Axz*gyz - Axy*gzz))
       - Azz*SQ(gxy) - Ayy*SQ(gxz) - Axx*SQ(gyz)
       ));
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::AlgConstr(AthenaArray<Real> & u)
// \brief algebraic constraints projection
//
// This function operates on all grid points of the MeshBlock

void Z4c::AlgConstr(AthenaArray<Real> & u)
{
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);

  GLOOP2(k,j) {

    // compute determinant and "conformal conformal factor"
    GLOOP1(i) {
      detg(i) = SpatialDet(z4c.g_dd,k,j,i);
      detg(i) = detg(i) > 0. ? detg(i) : 1.;
      Real eps = detg(i) - 1.;
      oopsi4(i) = (eps < opt.eps_floor) ? (1. - opt.eps_floor/3.) : (pow(1./detg(i), 1./3.));
    }
    // enforce unitary determinant for conformal metric
    for(int a = 0; a < NDIM; ++a)
    for(int b = a; b < NDIM; ++b) {
      GLOOP1(i) {
        z4c.g_dd(a,b,k,j,i) *= oopsi4(i);
      }
    }

    // compute trace of A
    GLOOP1(i) {
      // note: here we are assuming that det g = 1, which we enforced above
      A(i) = Trace(1.0,
          z4c.g_dd(0,0,k,j,i), z4c.g_dd(0,1,k,j,i), z4c.g_dd(0,2,k,j,i),
          z4c.g_dd(1,1,k,j,i), z4c.g_dd(1,2,k,j,i), z4c.g_dd(2,2,k,j,i),
          z4c.A_dd(0,0,k,j,i), z4c.A_dd(0,1,k,j,i), z4c.A_dd(0,2,k,j,i),
          z4c.A_dd(1,1,k,j,i), z4c.A_dd(1,2,k,j,i), z4c.A_dd(2,2,k,j,i));
    }
    // enforce trace of A to be zero
    for(int a = 0; a < NDIM; ++a)
    for(int b = a; b < NDIM; ++b) {
      GLOOP1(i) {
        z4c.A_dd(a,b,k,j,i) -= (1.0/3.0) * A(i) * z4c.g_dd(a,b,k,j,i);
      }
    }
  }
}
