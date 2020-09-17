//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of core functions in class Radiation

// C++ headers
#include <cmath>      // acos, cos, NAN, sin, sqrt
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str, string

// Athena++ headers
#include "radiation.hpp"
#include "../athena.hpp"                   // Real, indices
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../mesh/mesh.hpp"                // MeshBlock

// Declarations
void DefaultOpacity(MeshBlock *pmb, const AthenaArray<Real> &prim_hydro);

//----------------------------------------------------------------------------------------
// Radiation constructor
// Inputs:
//   pmb: pointer to containing MeshBlock
//   pin: pointer to runtime parameters

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block(pmb),
    coupled_to_matter(pin->GetBoolean("radiation", "coupled")),
    affect_fluid(pin->GetOrAddBoolean("radiation", "affect_fluid", coupled_to_matter)),
    nzeta(pin->GetInteger("radiation", "n_polar")),
    npsi(pin->GetInteger("radiation", "n_azimuthal")),
    nang((nzeta + 2*NGHOST) * (npsi + 2*NGHOST)),
    prim(nang, pmb->ncells3, pmb->ncells2, pmb->ncells1),
    prim1(nang, pmb->ncells3, pmb->ncells2, pmb->ncells1),
    cons(nang, pmb->ncells3, pmb->ncells2, pmb->ncells1),
    cons1(nang, pmb->ncells3, pmb->ncells2, pmb->ncells1),
    cons2(nang, pmb->ncells3, pmb->ncells2, pmb->ncells1),
    flux_x{
      {nang, pmb->ncells3, pmb->ncells2, pmb->ncells1+1},
      {nang, pmb->ncells3, pmb->ncells2+1, pmb->ncells1, (pmb->pmy_mesh->f2 ?
          AthenaArray<Real>::DataStatus::allocated :
          AthenaArray<Real>::DataStatus::empty)},
      {nang, pmb->ncells3+1, pmb->ncells2, pmb->ncells1, (pmb->pmy_mesh->f3 ?
          AthenaArray<Real>::DataStatus::allocated :
          AthenaArray<Real>::DataStatus::empty)}
    },
    coarse_cons(nang, pmb->ncc3, pmb->ncc2, pmb->ncc1, (pmb->pmy_mesh->multilevel ?
        AthenaArray<Real>::DataStatus::allocated : AthenaArray<Real>::DataStatus::empty)),
    coarse_prim(nang, pmb->ncc3, pmb->ncc2, pmb->ncc1, (pmb->pmy_mesh->multilevel ?
        AthenaArray<Real>::DataStatus::allocated : AthenaArray<Real>::DataStatus::empty)),
    rbvar(pmb, &cons, &coarse_cons, flux_x, nzeta, npsi) {

  // Verify temporal and special orders
  std::stringstream msg;
  if (pin->GetString("time", "integrator") != "vl2") {
    msg << "### FATAL ERROR in Radiation\n";
    msg << "only VL2 integration supported\n";
    throw std::runtime_error(msg.str().c_str());
  }
  if (pin->GetInteger("time", "xorder") > 2) {
    msg << "### FATAL ERROR in Radiation\n";
    msg << "only first and second order reconstruction supported\n";
    throw std::runtime_error(msg.str().c_str());
  }

  // Verify numbers of angles
  if (nzeta < 4) {
    msg << "### FATAL ERROR in Radiation constructor\n";
    msg << "too few polar angles\n";
    throw std::runtime_error(msg.str().c_str());
  }
  if (npsi < 4) {
    msg << "### FATAL ERROR in Radiation constructor\n";
    msg << "too few azimuthal angles\n";
    throw std::runtime_error(msg.str().c_str());
  }
  if (npsi%2 != 0) {
    msg << "### FATAL ERROR in Radiation constructor\n";
    msg << "must have even number of azimuthal angles\n";
    throw std::runtime_error(msg.str().c_str());
  }

  // Set object and function pointers
  UserSourceTerm = pmb->pmy_mesh->UserRadSourceTerm_;
  UpdateOpacity = DefaultOpacity;

  // Enroll refinement communication
  if (pmb->pmy_mesh->multilevel) {
    refinement_idx = pmb->pmr->AddToRefinement(&cons, &coarse_cons);
  }

  // Construct objects
  rbvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&rbvar);
  pmb->pbval->bvars_main_int.push_back(&rbvar);

  // Set flags
  if (not coupled_to_matter and UserSourceTerm == nullptr) {
    source_terms_defined = false;
  } else {
    source_terms_defined = true;
  }
  if (coupled_to_matter) {
    moment_fix = pin->GetBoolean("radiation", "moment_fix");
  } else {
    moment_fix = false;
  }

  // Set parameters
  nang_zf = (nzeta + 2*NGHOST + 1) * (npsi + 2*NGHOST);
  nang_pf = (nzeta + 2*NGHOST) * (npsi + 2*NGHOST + 1);
  nang_zpf = (nzeta + 2*NGHOST + 1) * (npsi + 2*NGHOST + 1);
  zs = NGHOST;
  ze = nzeta + NGHOST - 1;
  ps = NGHOST;
  pe = npsi + NGHOST - 1;
  is = pmy_block->is;
  ie = pmy_block->ie;
  js = pmy_block->js;
  je = pmy_block->je;
  ks = pmy_block->ks;
  ke = pmy_block->ke;

  // Set and calculate units
  density_cgs = pin->GetOrAddReal("radiation", "density_cgs", NAN);
  mol_weight = pin->GetOrAddReal("radiation", "mol_weight", NAN);
  arad = arad_cgs * SQR(c_cgs) * SQR(SQR(mol_weight * m_p_cgs * c_cgs / k_b_cgs))
      / density_cgs;

  // Allocate memory for angles
  zetaf.NewAthenaArray(nzeta + 2*NGHOST + 1);
  zetav.NewAthenaArray(nzeta + 2*NGHOST);
  dzetaf.NewAthenaArray(nzeta + 2*NGHOST);
  psif.NewAthenaArray(npsi + 2*NGHOST + 1);
  psiv.NewAthenaArray(npsi + 2*NGHOST);
  dpsif.NewAthenaArray(npsi + 2*NGHOST);
  zeta_length.NewAthenaArray(nzeta + 2*NGHOST, npsi + 2*NGHOST + 1);
  psi_length.NewAthenaArray(nzeta + 2*NGHOST + 1, npsi + 2*NGHOST);
  solid_angle.NewAthenaArray(nzeta + 2*NGHOST, npsi + 2*NGHOST);

  // Construct polar angles, equally spaced in cosine
  Real dczeta = -2.0 / nzeta;
  zetaf(zs) = 0.0;             // set north pole exactly
  zetaf(ze+1) = PI;            // set south pole exactly
  for (int l = zs+1; l <= (nzeta-1)/2+NGHOST; ++l) {
    Real czeta = 1.0 + (l - NGHOST) * dczeta;
    Real zeta = std::acos(czeta);
    zetaf(l) = zeta;                           // set northern active faces
    zetaf(ze+NGHOST+1-l) = PI - zeta;          // set southern active faces
  }
  if (nzeta%2 == 0) {
    zetaf(nzeta/2+NGHOST) = PI/2.0;  // set equator exactly if present
  }
  for (int l = zs-NGHOST; l <= zs-1; ++l) {
    zetaf(l) = -zetaf(2*NGHOST - l);                 // set northern ghost faces
    zetaf(ze+NGHOST+1-l) = 2.0*PI - zetaf(nzeta+l);  // set southern ghost faces
  }
  for (int l = zs-NGHOST; l <= ze+NGHOST; ++l) {
    zetav(l) = (zetaf(l+1) * std::cos(zetaf(l+1)) - std::sin(zetaf(l+1))
        - zetaf(l) * std::cos(zetaf(l)) + std::sin(zetaf(l))) / (std::cos(zetaf(l+1))
        - std::cos(zetaf(l)));
    dzetaf(l) = zetaf(l+1) - zetaf(l);
  }

  // Construct azimuthal angles, equally spaced
  Real dpsi = 2.0*PI / npsi;
  psif(ps) = 0.0;             // set origin exactly
  psif(pe+1) = 2.0*PI;        // set origin exactly
  for (int m = ps+1; m <= pe; ++m) {
    psif(m) = (m - NGHOST) * dpsi;  // set active faces
  }
  for (int m = ps-NGHOST; m <= ps-1; ++m) {
    psif(m) = psif(npsi+m) - 2.0*PI;                  // set beginning ghost faces
    psif(pe+NGHOST+1-m) = psif(2*NGHOST-m) + 2.0*PI;  // set end ghost faces
  }
  for (int m = ps-NGHOST; m <= pe+NGHOST; ++m) {
    psiv(m) = 0.5 * (psif(m) + psif(m+1));
    dpsif(m) = psif(m+1) - psif(m);
  }

  // Calculate angular lengths and areas
  for (int l = zs-NGHOST; l <= ze+NGHOST; ++l) {
    for (int m = ps-NGHOST; m <= pe+NGHOST+1; ++m) {
      zeta_length(l,m) = std::cos(zetaf(l)) - std::cos(zetaf(l+1));
    }
  }
  for (int l = zs-NGHOST; l <= ze+NGHOST+1; ++l) {
    for (int m = ps-NGHOST; m <= pe+NGHOST; ++m) {
      psi_length(l,m) = std::sin(zetaf(l)) * dpsif(m);
    }
  }
  for (int l = zs-NGHOST; l <= ze+NGHOST; ++l) {
    for (int m = ps-NGHOST; m <= pe+NGHOST; ++m) {
      solid_angle(l,m) = (std::cos(zetaf(l)) - std::cos(zetaf(l+1))) * dpsif(m);
    }
  }

  // Allocate memory for angle fluxes
  flux_a[ZETADIR].NewAthenaArray(nang_zf, pmb->ncells3, pmb->ncells2, pmb->ncells1);
  flux_a[PSIDIR].NewAthenaArray(nang_pf, pmb->ncells3, pmb->ncells2, pmb->ncells1);

  // Allocate memory for opacity array
  opacity.NewAthenaArray(NOPA, pmb->ncells3, pmb->ncells2, pmb->ncells1);

  // Allocate memory for moments
  moments_coord.NewAthenaArray(10, pmb->ncells3, pmb->ncells2, pmb->ncells1);
  moments_tetrad.NewAthenaArray(10, pmb->ncells3, pmb->ncells2, pmb->ncells1);
  if (coupled_to_matter) {
    moments_fluid.NewAthenaArray(10, pmb->ncells3, pmb->ncells2, pmb->ncells1);
  }

  // Allocate memory for unit normal components in orthonormal frame
  int num_cells_zeta = ze + NGHOST;
  int num_cells_psi = pe + NGHOST;
  nh_cc_.NewAthenaArray(4, num_cells_zeta, num_cells_psi);
  nh_fc_.NewAthenaArray(4, num_cells_zeta + 1, num_cells_psi);
  nh_cf_.NewAthenaArray(4, num_cells_zeta, num_cells_psi + 1);

  // Calculate unit normal components in orthonormal frame at angle centers
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      nh_cc_(0,l,m) = 1.0;
      nh_cc_(1,l,m) = std::sin(zetav(l)) * std::cos(psiv(m));
      nh_cc_(2,l,m) = std::sin(zetav(l)) * std::sin(psiv(m));
      nh_cc_(3,l,m) = std::cos(zetav(l));
    }
  }

  // Calculate unit normal components in orthonormal frame at zeta-faces
  for (int l = zs; l <= ze+1; ++l) {
    for (int m = ps; m <= pe; ++m) {
      nh_fc_(0,l,m) = 1.0;
      nh_fc_(1,l,m) = std::sin(zetaf(l)) * std::cos(psiv(m));
      nh_fc_(2,l,m) = std::sin(zetaf(l)) * std::sin(psiv(m));
      nh_fc_(3,l,m) = std::cos(zetaf(l));
    }
  }

  // Calculate unit normal components in orthonormal frame at psi-faces
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe+1; ++m) {
      nh_cf_(0,l,m) = 1.0;
      nh_cf_(1,l,m) = std::sin(zetav(l)) * std::cos(psif(m));
      nh_cf_(2,l,m) = std::sin(zetav(l)) * std::sin(psif(m));
      nh_cf_(3,l,m) = std::cos(zetav(l));
    }
  }

  // Allocate memory for unit normal and related components in coordinate frame
  nmu_.NewAthenaArray(4, num_cells_zeta, num_cells_psi, pmb->ncells3, pmb->ncells2,
      pmb->ncells1);
  n0_n_mu_.NewAthenaArray(4, num_cells_zeta, num_cells_psi, pmb->ncells3, pmb->ncells2,
      pmb->ncells1);
  n1_n_0_.NewAthenaArray(num_cells_zeta, num_cells_psi, pmb->ncells3, pmb->ncells2,
      pmb->ncells1 + 1);
  n2_n_0_.NewAthenaArray(num_cells_zeta, num_cells_psi, pmb->ncells3, pmb->ncells2 + 1,
      pmb->ncells1);
  n3_n_0_.NewAthenaArray(num_cells_zeta, num_cells_psi, pmb->ncells3 + 1, pmb->ncells2,
      pmb->ncells1);
  na1_n_0_.NewAthenaArray(num_cells_zeta + 1, num_cells_psi, pmb->ncells3, pmb->ncells2,
      pmb->ncells1);
  na2_n_0_.NewAthenaArray(num_cells_zeta, num_cells_psi + 1, pmb->ncells3, pmb->ncells2,
      pmb->ncells1);

  // Allocate memory for temporary geometric quantities
  AthenaArray<Real> e, e_cov, omega;
  e.NewAthenaArray(4, 4);
  e_cov.NewAthenaArray(4, 4);
  omega.NewAthenaArray(4, 4, 4);

  // Calculate n^mu and n^0 n_mu
  int kl = ks - (pmb->ncells3 > 1 ? NGHOST : 0);
  int ku = ke + (pmb->ncells3 > 1 ? NGHOST : 0);
  int jl = js - (pmb->ncells2 > 1 ? NGHOST : 0);
  int ju = je + (pmb->ncells2 > 1 ? NGHOST : 0);
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            Real n0 = 0.0;
            Real n1 = 0.0;
            Real n2 = 0.0;
            Real n3 = 0.0;
            Real n_0 = 0.0;
            Real n_1 = 0.0;
            Real n_2 = 0.0;
            Real n_3 = 0.0;
            for (int n = 0; n < 4; ++n) {
              n0 += e(n,0) * nh_cc_(n,l,m);
              n1 += e(n,1) * nh_cc_(n,l,m);
              n2 += e(n,2) * nh_cc_(n,l,m);
              n3 += e(n,3) * nh_cc_(n,l,m);
              n_0 += e_cov(n,0) * nh_cc_(n,l,m);
              n_1 += e_cov(n,1) * nh_cc_(n,l,m);
              n_2 += e_cov(n,2) * nh_cc_(n,l,m);
              n_3 += e_cov(n,3) * nh_cc_(n,l,m);
            }
            nmu_(0,l,m,k,j,i) = n0;
            nmu_(1,l,m,k,j,i) = n1;
            nmu_(2,l,m,k,j,i) = n2;
            nmu_(3,l,m,k,j,i) = n3;
            n0_n_mu_(0,l,m,k,j,i) = n0 * n_0;
            n0_n_mu_(1,l,m,k,j,i) = n0 * n_1;
            n0_n_mu_(2,l,m,k,j,i) = n0 * n_2;
            n0_n_mu_(3,l,m,k,j,i) = n0 * n_3;
          }
        }
      }
    }
  }

  // Calculate n^1 n_0
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie+1; ++i) {
        Real x1 = pmb->pcoord->x1f(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            Real n1 = 0.0;
            Real n_0 = 0.0;
            for (int n = 0; n < 4; ++n) {
              n1 += e(n,1) * nh_cc_(n,l,m);
              n_0 += e_cov(n,0) * nh_cc_(n,l,m);
            }
            n1_n_0_(l,m,k,j,i) = n1 * n_0;
          }
        }
      }
    }
  }

  // Calculate n^2 n_0
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je+1; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2f(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            Real n2 = 0.0;
            Real n_0 = 0.0;
            for (int n = 0; n < 4; ++n) {
              n2 += e(n,2) * nh_cc_(n,l,m);
              n_0 += e_cov(n,0) * nh_cc_(n,l,m);
            }
            n2_n_0_(l,m,k,j,i) = n2 * n_0;
          }
        }
      }
    }
  }

  // Calculate n^3 n_0
  for (int k = ks; k <= ke+1; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3f(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            Real n3 = 0.0;
            Real n_0 = 0.0;
            for (int n = 0; n < 4; ++n) {
              n3 += e(n,3) * nh_cc_(n,l,m);
              n_0 += e_cov(n,0) * nh_cc_(n,l,m);
            }
            n3_n_0_(l,m,k,j,i) = n3 * n_0;
          }
        }
      }
    }
  }

  // Calculate n^zeta n_0
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);
        for (int l = zs+1; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            Real na1 = 0.0;
            for (int n = 0; n < 4; ++n) {
              for (int p = 0; p < 4; ++p) {
                na1 += 1.0 / std::sin(zetaf(l)) * nh_fc_(n,l,m) * nh_fc_(p,l,m)
                    * (nh_fc_(0,l,m) * omega(3,n,p) - nh_fc_(3,l,m) * omega(0,n,p));
              }
            }
            Real n_0 = 0.0;
            for (int n = 0; n < 4; ++n) {
              n_0 += e_cov(n,0) * nh_fc_(n,l,m);
            }
            na1_n_0_(l,m,k,j,i) = na1 * n_0;
          }
        }
      }
    }
  }

  // Calculate n^psi n_0
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe+1; ++m) {
            Real na2 = 0.0;
            for (int n = 0; n < 4; ++n) {
              for (int p = 0; p < 4; ++p) {
                na2 += 1.0 / SQR(std::sin(zetav(l))) * nh_cf_(n,l,m) * nh_cf_(p,l,m)
                    * (nh_cf_(2,l,m) * omega(1,n,p) - nh_cf_(1,l,m) * omega(2,n,p));
              }
            }
            Real n_0 = 0.0;
            for (int n = 0; n < 4; ++n) {
              n_0 += e_cov(n,0) * nh_cf_(n,l,m);
            }
            na2_n_0_(l,m,k,j,i) = na2 * n_0;
          }
        }
      }
    }
  }

  // Allocate memory for metric
  g_.NewAthenaArray(NMETRIC, pmb->ncells1);
  gi_.NewAthenaArray(NMETRIC, pmb->ncells1);

  // Allocate memory for reconstruction
  ii_l_.NewAthenaArray(nang_zpf, pmb->ncells1);
  ii_r_.NewAthenaArray(nang_zpf, pmb->ncells1);

  // Allocate memory for flux calculation
  if (coupled_to_matter) {
    norm_to_tet_1_.NewAthenaArray(4, 4, pmb->ncells3, pmb->ncells2, pmb->ncells1);
    if (js != je) {
      norm_to_tet_2_.NewAthenaArray(4, 4, pmb->ncells3, pmb->ncells2, pmb->ncells1);
    }
    if (ks != ke) {
      norm_to_tet_3_.NewAthenaArray(4, 4, pmb->ncells3, pmb->ncells2, pmb->ncells1);
    }
    ii_lr_.NewAthenaArray(nang_zpf, pmb->ncells1);
    jj_f_.NewAthenaArray(pmb->ncells1);
    k_tot_.NewAthenaArray(pmb->ncells1);
    bb_jj_f_.NewAthenaArray(pmb->ncells1);
    neg_u_n_.NewAthenaArray(nang_zpf, pmb->ncells1);
  }

  // Allocate memory for flux divergence calculation
  area_l_.NewAthenaArray(pmb->ncells1 + 1);
  area_r_.NewAthenaArray(pmb->ncells1 + 1);
  vol_.NewAthenaArray(pmb->ncells1 + 1);
  flux_div_.NewAthenaArray(nang, pmb->ncells1 + 1);

  // Allocate memory for source term calculation
  norm_to_tet_.NewAthenaArray(4, 4, pmb->ncells3, pmb->ncells2, pmb->ncells1);
  if (coupled_to_matter) {
    if (affect_fluid) {
      moments_old_.NewAthenaArray(4, pmb->ncells1);
      moments_new_.NewAthenaArray(4, pmb->ncells1);
    }
    u_tet_.NewAthenaArray(4, pmb->ncells1);
    coefficients_.NewAthenaArray(2, pmb->ncells1);
    bad_cell_.NewAthenaArray(pmb->ncells1);
    tt_plus_.NewAthenaArray(pmb->ncells1);
    ee_f_minus_.NewAthenaArray(pmb->ncells1);
    ee_f_plus_.NewAthenaArray(pmb->ncells1);
  }

  // Calculate transformation from normal frame to tetrad frame (cell)
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pmy_block->pcoord->CellMetric(k, j, is, ie, g_, gi_);
      for (int i = is; i <= ie; ++i) {

        // Set Minkowski metric
        Real eta[4][4] = {};
        eta[0][0] = -1.0;
        eta[1][1] = 1.0;
        eta[2][2] = 1.0;
        eta[3][3] = 1.0;

        // Calculate coordinate-to-tetrad transformation
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);

        // Calculate normal-to-coordinate transformation
        Real norm_to_coord[4][4] = {};
        Real alpha = 1.0 / std::sqrt(-gi_(I00,i));
        norm_to_coord[0][0] = 1.0 / alpha;
        norm_to_coord[1][0] = -alpha * gi_(I01,i);
        norm_to_coord[2][0] = -alpha * gi_(I02,i);
        norm_to_coord[3][0] = -alpha * gi_(I03,i);
        norm_to_coord[1][1] = 1.0;
        norm_to_coord[2][2] = 1.0;
        norm_to_coord[3][3] = 1.0;

        // Concatenate transformations
        for (int m = 0; m < 4; ++m) {
          for (int n = 0; n < 4; ++n) {
            norm_to_tet_(m,n,k,j,i) = 0.0;
            for (int p = 0; p < 4; ++p) {
              for (int q = 0; q < 4; ++q) {
                norm_to_tet_(m,n,k,j,i) += eta[m][p] * e_cov(p,q) * norm_to_coord[q][n];
              }
            }
          }
        }
      }
    }
  }

  // Calculate transformation from normal frame to tetrad frame (x1-face)
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pmy_block->pcoord->Face1Metric(k, j, is, ie+1, g_, gi_);
      for (int i = is; i <= ie+1; ++i) {

        // Set Minkowski metric
        Real eta[4][4] = {};
        eta[0][0] = -1.0;
        eta[1][1] = 1.0;
        eta[2][2] = 1.0;
        eta[3][3] = 1.0;

        // Calculate coordinate-to-tetrad transformation
        Real x1 = pmb->pcoord->x1f(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);

        // Calculate normal-to-coordinate transformation
        Real norm_to_coord[4][4] = {};
        Real alpha = 1.0 / std::sqrt(-gi_(I00,i));
        norm_to_coord[0][0] = 1.0 / alpha;
        norm_to_coord[1][0] = -alpha * gi_(I01,i);
        norm_to_coord[2][0] = -alpha * gi_(I02,i);
        norm_to_coord[3][0] = -alpha * gi_(I03,i);
        norm_to_coord[1][1] = 1.0;
        norm_to_coord[2][2] = 1.0;
        norm_to_coord[3][3] = 1.0;

        // Concatenate transformations
        for (int m = 0; m < 4; ++m) {
          for (int n = 0; n < 4; ++n) {
            norm_to_tet_1_(m,n,k,j,i) = 0.0;
            for (int p = 0; p < 4; ++p) {
              for (int q = 0; q < 4; ++q) {
                norm_to_tet_1_(m,n,k,j,i) += eta[m][p] * e_cov(p,q) * norm_to_coord[q][n];
              }
            }
          }
        }
      }
    }
  }

  // Calculate transformation from normal frame to tetrad frame (x2-face)
  if (js != je) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je+1; ++j) {
        pmy_block->pcoord->Face2Metric(k, j, is, ie, g_, gi_);
        for (int i = is; i <= ie; ++i) {

          // Set Minkowski metric
          Real eta[4][4] = {};
          eta[0][0] = -1.0;
          eta[1][1] = 1.0;
          eta[2][2] = 1.0;
          eta[3][3] = 1.0;

          // Calculate coordinate-to-tetrad transformation
          Real x1 = pmb->pcoord->x1v(i);
          Real x2 = pmb->pcoord->x2f(j);
          Real x3 = pmb->pcoord->x3v(k);
          pmb->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);

          // Calculate normal-to-coordinate transformation
          Real norm_to_coord[4][4] = {};
          Real alpha = 1.0 / std::sqrt(-gi_(I00,i));
          norm_to_coord[0][0] = 1.0 / alpha;
          norm_to_coord[1][0] = -alpha * gi_(I01,i);
          norm_to_coord[2][0] = -alpha * gi_(I02,i);
          norm_to_coord[3][0] = -alpha * gi_(I03,i);
          norm_to_coord[1][1] = 1.0;
          norm_to_coord[2][2] = 1.0;
          norm_to_coord[3][3] = 1.0;

          // Concatenate transformations
          for (int m = 0; m < 4; ++m) {
            for (int n = 0; n < 4; ++n) {
              norm_to_tet_2_(m,n,k,j,i) = 0.0;
              for (int p = 0; p < 4; ++p) {
                for (int q = 0; q < 4; ++q) {
                  norm_to_tet_2_(m,n,k,j,i) +=
                      eta[m][p] * e_cov(p,q) * norm_to_coord[q][n];
                }
              }
            }
          }
        }
      }
    }
  }

  // Calculate transformation from normal frame to tetrad frame (x3-face)
  if (ks != ke) {
    for (int k = ks; k <= ke+1; ++k) {
      for (int j = js; j <= je; ++j) {
        pmy_block->pcoord->Face3Metric(k, j, is, ie, g_, gi_);
        for (int i = is; i <= ie; ++i) {

          // Set Minkowski metric
          Real eta[4][4] = {};
          eta[0][0] = -1.0;
          eta[1][1] = 1.0;
          eta[2][2] = 1.0;
          eta[3][3] = 1.0;

          // Calculate coordinate-to-tetrad transformation
          Real x1 = pmb->pcoord->x1v(i);
          Real x2 = pmb->pcoord->x2v(j);
          Real x3 = pmb->pcoord->x3f(k);
          pmb->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);

          // Calculate normal-to-coordinate transformation
          Real norm_to_coord[4][4] = {};
          Real alpha = 1.0 / std::sqrt(-gi_(I00,i));
          norm_to_coord[0][0] = 1.0 / alpha;
          norm_to_coord[1][0] = -alpha * gi_(I01,i);
          norm_to_coord[2][0] = -alpha * gi_(I02,i);
          norm_to_coord[3][0] = -alpha * gi_(I03,i);
          norm_to_coord[1][1] = 1.0;
          norm_to_coord[2][2] = 1.0;
          norm_to_coord[3][3] = 1.0;

          // Concatenate transformations
          for (int m = 0; m < 4; ++m) {
            for (int n = 0; n < 4; ++n) {
              norm_to_tet_3_(m,n,k,j,i) = 0.0;
              for (int p = 0; p < 4; ++p) {
                for (int q = 0; q < 4; ++q) {
                  norm_to_tet_3_(m,n,k,j,i) +=
                      eta[m][p] * e_cov(p,q) * norm_to_coord[q][n];
                }
              }
            }
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
// Radiation destructor

Radiation::~Radiation() {}

//----------------------------------------------------------------------------------------
// Default opacity function
// Inputs:
//   pmb: pointer to MeshBlock
//   prim_hydro: primitive variables
// Outputs: (none)
// Notes:
//   Does nothing; keeps opacities (not absorption coefficients) at their initial values.

void DefaultOpacity(MeshBlock *pmb, const AthenaArray<Real> &prim_hydro) {
  return;
}
