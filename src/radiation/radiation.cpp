//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of functions in class Radiation

// C++ headers
#include <algorithm>  // max
#include <cmath>      // acos, cos, NAN, sin, sqrt
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str, string

// Athena++ headers
#include "radiation.hpp"
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

  // Verify numbers of angles
  std::stringstream msg;
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
  AthenaArray<Real> nh_fc, nh_cf;
  nh_cc_.NewAthenaArray(4, num_cells_zeta, num_cells_psi);
  nh_fc.NewAthenaArray(4, num_cells_zeta + 1, num_cells_psi);
  nh_cf.NewAthenaArray(4, num_cells_zeta, num_cells_psi + 1);

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
      nh_fc(0,l,m) = 1.0;
      nh_fc(1,l,m) = std::sin(zetaf(l)) * std::cos(psiv(m));
      nh_fc(2,l,m) = std::sin(zetaf(l)) * std::sin(psiv(m));
      nh_fc(3,l,m) = std::cos(zetaf(l));
    }
  }

  // Calculate unit normal components in orthonormal frame at psi-faces
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe+1; ++m) {
      nh_cf(0,l,m) = 1.0;
      nh_cf(1,l,m) = std::sin(zetav(l)) * std::cos(psif(m));
      nh_cf(2,l,m) = std::sin(zetav(l)) * std::sin(psif(m));
      nh_cf(3,l,m) = std::cos(zetav(l));
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
                na1 += 1.0 / std::sin(zetaf(l)) * nh_fc(n,l,m) * nh_fc(p,l,m)
                    * (nh_fc(0,l,m) * omega(3,n,p) - nh_fc(3,l,m) * omega(0,n,p));
              }
            }
            Real n_0 = 0.0;
            for (int n = 0; n < 4; ++n) {
              n_0 += e_cov(n,0) * nh_fc(n,l,m);
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
                na2 += 1.0 / SQR(std::sin(zetav(l))) * nh_cf(n,l,m) * nh_cf(p,l,m)
                    * (nh_cf(2,l,m) * omega(1,n,p) - nh_cf(1,l,m) * omega(2,n,p));
              }
            }
            Real n_0 = 0.0;
            for (int n = 0; n < 4; ++n) {
              n_0 += e_cov(n,0) * nh_cf(n,l,m);
            }
            na2_n_0_(l,m,k,j,i) = na2 * n_0;
          }
        }
      }
    }
  }

  // Allocate memory for left and right reconstructed states
  prim_l_.NewAthenaArray(nang_zpf, pmb->ncells1 + 1);
  prim_r_.NewAthenaArray(nang_zpf, pmb->ncells1 + 1);

  // Allocate memory for flux divergence calculation
  area_l_.NewAthenaArray(pmb->ncells1 + 1);
  area_r_.NewAthenaArray(pmb->ncells1 + 1);
  vol_.NewAthenaArray(pmb->ncells1 + 1);
  flux_div_.NewAthenaArray(nang, pmb->ncells1 + 1);

  // Allocate memory for source term calculation in a single cell
  intensity_scr_.NewAthenaArray(nzeta * npsi);
  tran_coef_.NewAthenaArray(nzeta * npsi);
  weight_.NewAthenaArray(nzeta * npsi);
  vncsigma2_.NewAthenaArray(nzeta * npsi);

  // Allocate memory for source term frame transformations
  g_.NewAthenaArray(NMETRIC, pmb->ncells1);
  gi_.NewAthenaArray(NMETRIC, pmb->ncells1);
  norm_to_tet_.NewAthenaArray(4, 4, pmb->ncells3, pmb->ncells2, pmb->ncells1);
  u_tet_.NewAthenaArray(4, pmb->ncells1);
  dt_.NewAthenaArray(pmb->ncells1);
  dtau_.NewAthenaArray(pmb->ncells1);
  weight_sum_.NewAthenaArray(pmb->ncells1);
  n_cm_.NewAthenaArray(4, nzeta * npsi, pmb->ncells1);
  n0_.NewAthenaArray(nzeta * npsi, pmb->ncells3, pmb->ncells2, pmb->ncells1);
  omega_cm_.NewAthenaArray(nzeta * npsi, pmb->ncells1);
  intensity_cm_.NewAthenaArray(nzeta * npsi, pmb->ncells1);
  moments_old_.NewAthenaArray(4, pmb->ncells1);
  moments_new_.NewAthenaArray(4, pmb->ncells1);

  // Calculate transformation from normal frame to tetrad frame
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

  // Calculate n^0
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            int lm_alt = (l - zs) * (pe - ps + 1) + m - ps;
            Real n0 = 0.0;
            for (int n = 0; n < 4; ++n) {
              n0 += e(n,0) * nh_cc_(n,l,m);
            }
            n0_(lm_alt,k,j,i) = n0;
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
// Function for averaging intensities according to integrator weights
// Inputs:
//   cons_out, cons_in_1, cons_in_2: conserved intensity arrays, possibly uninitialized
//   weights: integrator weights
// Outputs:
//   cons_out: weighted intensity
// Notes:
//   Same procedure as in Hydro::WeightedAveU().

void Radiation::WeightedAve(AthenaArray<Real> &cons_out, AthenaArray<Real> &cons_in_1,
    AthenaArray<Real> &cons_in_2, const Real weights[3]) {

  // Apply averaging based on which weights are 0
  if (weights[2] != 0.0) {
    for (int n = 0; n < nang; ++n) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          #pragma omp simd
          for (int i = is; i <= ie; ++i) {
            cons_out(n,k,j,i) = weights[0] * cons_out(n,k,j,i)
                + weights[1] * cons_in_1(n,k,j,i) + weights[2] * cons_in_2(n,k,j,i);
          }
        }
      }
    }
  } else if (weights[1] != 0.0) {
    for (int n = 0; n < nang; ++n) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          #pragma omp simd
          for (int i = is; i <= ie; ++i) {
            cons_out(n,k,j,i) =
                weights[0] * cons_out(n,k,j,i) + weights[1] * cons_in_1(n,k,j,i);
          }
        }
      }
    }
  } else if (weights[0] != 0.0) {
    for (int n = 0; n < nang; ++n) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          #pragma omp simd
          for (int i = is; i <= ie; ++i) {
            cons_out(n,k,j,i) = weights[0] * cons_out(n,k,j,i);
          }
        }
      }
    }
  } else {
    for (int n = 0; n < nang; ++n) {
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          #pragma omp simd
          for (int i = is; i <= ie; ++i) {
            cons_out(n,k,j,i) = 0.0;
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating spatial and angular fluxes
// Inputs:
//   prim_in: primitive intensity
//   order: reconstruction order
// Outputs:
//   this->flux_x, this->flux_a: fluxes set

void Radiation::CalculateFluxes(AthenaArray<Real> &prim_in, int order) {

  // Check order
  if (order != 1 and order != 2) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Radiation\n";
    msg << "only first and second order reconstruction supported\n";
    throw std::runtime_error(msg.str().c_str());
  }

  // Calculate x1-fluxes
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {

          // Reconstruction
          if (order == 1) {
            for (int i = is; i <= ie+1; ++i) {
              prim_l_(lm,i) = prim_in(lm,k,j,i-1);
              prim_r_(lm,i) = prim_in(lm,k,j,i);
            }
          } else {
            for (int i = is; i <= ie+1; ++i) {
              Real x_l = pmy_block->pcoord->x1v(i-1);
              Real x_c = pmy_block->pcoord->x1f(i);
              Real x_r = pmy_block->pcoord->x1v(i);
              Real dx_l = x_c - x_l;
              Real dx_r = x_r - x_c;
              Real q_ll = prim_in(lm,k,j,i-2);
              Real q_l = prim_in(lm,k,j,i-1);
              Real q_r = prim_in(lm,k,j,i);
              Real q_rr = prim_in(lm,k,j,i+1);
              Real dq_l = q_l - q_ll;
              Real dq_c = q_r - q_l;
              Real dq_r = q_rr - q_r;
              Real dq_2_l = dq_l * dq_c;
              Real dq_2_r = dq_c * dq_r;
              Real dq_m_l = (dq_2_l > 0.0) ? 2.0 * dq_2_l / (dq_l + dq_c) : 0.0;
              Real dq_m_r = (dq_2_r > 0.0) ? 2.0 * dq_2_r / (dq_c + dq_r) : 0.0;
              prim_l_(lm,i) = q_l + dx_l * dq_m_l;
              prim_r_(lm,i) = q_r - dx_r * dq_m_r;
            }
          }

          // Upwind flux calculation
          for (int i = is; i <= ie+1; ++i) {
            Real n1_n_0 = n1_n_0_(l,m,k,j,i);
            if (n1_n_0 < 0.0) {
              flux_x[X1DIR](lm,k,j,i) = n1_n_0 * prim_l_(lm,i);
            } else {
              flux_x[X1DIR](lm,k,j,i) = n1_n_0 * prim_r_(lm,i);
            }
          }
        }
      }
    }
  }

  // Calculate x2-fluxes
  if (js != je) {
    for (int l = zs; l <= ze; ++l) {
      for (int m = ps; m <= pe; ++m) {
        int lm = AngleInd(l, m);
        for (int k = ks; k <= ke; ++k) {
          for (int j = js; j <= je+1; ++j) {

            // Reconstruction
            if (order == 1) {
              for (int i = is; i <= ie; ++i) {
                prim_l_(lm,i) = prim_in(lm,k,j-1,i);
                prim_r_(lm,i) = prim_in(lm,k,j,i);
              }
            } else {
              Real x_l = pmy_block->pcoord->x2v(j-1);
              Real x_c = pmy_block->pcoord->x2f(j);
              Real x_r = pmy_block->pcoord->x2v(j);
              Real dx_l = x_c - x_l;
              Real dx_r = x_r - x_c;
              for (int i = is; i <= ie; ++i) {
                Real q_ll = prim_in(lm,k,j-2,i);
                Real q_l = prim_in(lm,k,j-1,i);
                Real q_r = prim_in(lm,k,j,i);
                Real q_rr = prim_in(lm,k,j+1,i);
                Real dq_l = q_l - q_ll;
                Real dq_c = q_r - q_l;
                Real dq_r = q_rr - q_r;
                Real dq_2_l = dq_l * dq_c;
                Real dq_2_r = dq_c * dq_r;
                Real dq_m_l = (dq_2_l > 0.0) ? 2.0 * dq_2_l / (dq_l + dq_c) : 0.0;
                Real dq_m_r = (dq_2_r > 0.0) ? 2.0 * dq_2_r / (dq_c + dq_r) : 0.0;
                prim_l_(lm,i) = q_l + dx_l * dq_m_l;
                prim_r_(lm,i) = q_r - dx_r * dq_m_r;
              }
            }

            // Upwind flux calculation
            for (int i = is; i <= ie; ++i) {
              Real n2_n_0 = n2_n_0_(l,m,k,j,i);
              if (n2_n_0 < 0.0) {
                flux_x[X2DIR](lm,k,j,i) = n2_n_0 * prim_l_(lm,i);
              } else {
                flux_x[X2DIR](lm,k,j,i) = n2_n_0 * prim_r_(lm,i);
              }
            }
          }
        }
      }
    }
  }

  // Calculate x3-fluxes
  if (ks != ke) {
    for (int l = zs; l <= ze; ++l) {
      for (int m = ps; m <= pe; ++m) {
        int lm = AngleInd(l, m);
        for (int k = ks; k <= ke+1; ++k) {
          for (int j = js; j <= je; ++j) {

            // Reconstruction
            if (order == 1) {
              for (int i = is; i <= ie; ++i) {
                prim_l_(lm,i) = prim_in(lm,k-1,j,i);
                prim_r_(lm,i) = prim_in(lm,k,j,i);
              }
            } else {
              Real x_l = pmy_block->pcoord->x3v(k-1);
              Real x_c = pmy_block->pcoord->x3f(k);
              Real x_r = pmy_block->pcoord->x3v(k);
              Real dx_l = x_c - x_l;
              Real dx_r = x_r - x_c;
              for (int i = is; i <= ie; ++i) {
                Real q_ll = prim_in(lm,k-2,j,i);
                Real q_l = prim_in(lm,k-1,j,i);
                Real q_r = prim_in(lm,k,j,i);
                Real q_rr = prim_in(lm,k+1,j,i);
                Real dq_l = q_l - q_ll;
                Real dq_c = q_r - q_l;
                Real dq_r = q_rr - q_r;
                Real dq_2_l = dq_l * dq_c;
                Real dq_2_r = dq_c * dq_r;
                Real dq_m_l = (dq_2_l > 0.0) ? 2.0 * dq_2_l / (dq_l + dq_c) : 0.0;
                Real dq_m_r = (dq_2_r > 0.0) ? 2.0 * dq_2_r / (dq_c + dq_r) : 0.0;
                prim_l_(lm,i) = q_l + dx_l * dq_m_l;
                prim_r_(lm,i) = q_r - dx_r * dq_m_r;
              }
            }

            // Upwind flux calculation
            for (int i = is; i <= ie; ++i) {
              Real n3_n_0 = n3_n_0_(l,m,k,j,i);
              if (n3_n_0 < 0.0) {
                flux_x[X3DIR](lm,k,j,i) = n3_n_0 * prim_l_(lm,i);
              } else {
                flux_x[X3DIR](lm,k,j,i) = n3_n_0 * prim_r_(lm,i);
              }
            }
          }
        }
      }
    }
  }

  // Calculate zeta-fluxes
  for (int l = zs; l <= ze+1; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm_ll = AngleInd(l - 2, m, false, false);
      int lm_l = AngleInd(l - 1, m, false, false);
      int lm_c = AngleInd(l, m, true, false);
      int lm_r = AngleInd(l, m, false, false);
      int lm_rr = AngleInd(l + 1, m, false, false);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {

          // Reconstruction
          if (order == 1) {
            for (int i = is; i <= ie; ++i) {
              prim_l_(lm_c,i) = prim_in(lm_l,k,j,i);
              prim_r_(lm_c,i) = prim_in(lm_r,k,j,i);
            }
          } else {
            Real x_l = zetav(l-1);
            Real x_c = zetaf(l);
            Real x_r = zetav(l);
            Real dx_l = x_c - x_l;
            Real dx_r = x_r - x_c;
            for (int i = is; i <= ie; ++i) {
              Real q_ll = prim_in(lm_ll,k,j,i);
              Real q_l = prim_in(lm_l,k,j,i);
              Real q_r = prim_in(lm_r,k,j,i);
              Real q_rr = prim_in(lm_rr,k,j,i);
              Real dq_l = q_l - q_ll;
              Real dq_c = q_r - q_l;
              Real dq_r = q_rr - q_r;
              Real dq_2_l = dq_l * dq_c;
              Real dq_2_r = dq_c * dq_r;
              Real dq_m_l = (dq_2_l > 0.0) ? 2.0 * dq_2_l / (dq_l + dq_c) : 0.0;
              Real dq_m_r = (dq_2_r > 0.0) ? 2.0 * dq_2_r / (dq_c + dq_r) : 0.0;
              prim_l_(lm_c,i) = q_l + dx_l * dq_m_l;
              prim_r_(lm_c,i) = q_r - dx_r * dq_m_r;
            }
          }

          // Upwind flux calculation
          for (int i = is; i <= ie; ++i) {
            Real na1_n_0 = na1_n_0_(l,m,k,j,i);
            if (na1_n_0 < 0.0) {
              flux_a[ZETADIR](lm_c,k,j,i) = na1_n_0 * prim_l_(lm_c,i);
            } else {
              flux_a[ZETADIR](lm_c,k,j,i) = na1_n_0 * prim_r_(lm_c,i);
            }
          }
        }
      }
    }
  }

  // Calculate psi-fluxes
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe+1; ++m) {
      int lm_ll = AngleInd(l, m - 2, false, false);
      int lm_l = AngleInd(l, m - 1, false, false);
      int lm_c = AngleInd(l, m, false, true);
      int lm_r = AngleInd(l, m, false, false);
      int lm_rr = AngleInd(l, m + 1, false, false);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {

          // Reconstruction
          if (order == 1) {
            for (int i = is; i <= ie; ++i) {
              prim_l_(lm_c,i) = prim_in(lm_l,k,j,i);
              prim_r_(lm_c,i) = prim_in(lm_r,k,j,i);
            }
          } else {
            Real x_l = psiv(m-1);
            Real x_c = psif(m);
            Real x_r = psiv(m);
            Real dx_l = x_c - x_l;
            Real dx_r = x_r - x_c;
            for (int i = is; i <= ie; ++i) {
              Real q_ll = prim_in(lm_ll,k,j,i);
              Real q_l = prim_in(lm_l,k,j,i);
              Real q_r = prim_in(lm_r,k,j,i);
              Real q_rr = prim_in(lm_rr,k,j,i);
              Real dq_l = q_l - q_ll;
              Real dq_c = q_r - q_l;
              Real dq_r = q_rr - q_r;
              Real dq_2_l = dq_l * dq_c;
              Real dq_2_r = dq_c * dq_r;
              Real dq_m_l = (dq_2_l > 0.0) ? 2.0 * dq_2_l / (dq_l + dq_c) : 0.0;
              Real dq_m_r = (dq_2_r > 0.0) ? 2.0 * dq_2_r / (dq_c + dq_r) : 0.0;
              prim_l_(lm_c,i) = q_l + dx_l * dq_m_l;
              prim_r_(lm_c,i) = q_r - dx_r * dq_m_r;
            }
          }

          // Upwind flux calculation
          for (int i = is; i <= ie; ++i) {
            Real na2_n_0 = na2_n_0_(l,m,k,j,i);
            if (na2_n_0 < 0.0) {
              flux_a[PSIDIR](lm_c,k,j,i) = na2_n_0 * prim_l_(lm_c,i);
            } else {
              flux_a[PSIDIR](lm_c,k,j,i) = na2_n_0 * prim_r_(lm_c,i);
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for updating conserved quantities
// Inputs:
//   prim_in: primitive intensity
//   dt: timestep for stage of integration
// Outputs:
//   cons_out: conserved values updated

void Radiation::AddFluxDivergenceToAverage(AthenaArray<Real> &prim_in, const Real dt,
    AthenaArray<Real> &cons_out) {

  // Extract Coordinates and timestep
  Coordinates *pcoord = pmy_block->pcoord;

  // Calculate angle index range (including some but not all unnecessary ghost zones)
  int lms = AngleInd(zs, ps);
  int lme = AngleInd(ze, pe);

  // Go through all cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {

      // Determine poles
      bool left_pole = pcoord->IsPole(j);
      bool right_pole = pcoord->IsPole(j+1);

      // Calculate x1-divergence
      pcoord->Face1Area(k, j, is, ie+1, area_l_);
      for (int lm = lms; lm <= lme; ++lm) {
        for (int i = is; i <= ie; ++i) {
          flux_div_(lm,i) = area_l_(i+1) * flux_x[X1DIR](lm,k,j,i+1)
              - area_l_(i) * flux_x[X1DIR](lm,k,j,i);
        }
      }

      // Add x2-divergence
      if (js != je) {
        pcoord->Face2Area(k, j, is, ie, area_l_);
        pcoord->Face2Area(k, j+1, is, ie, area_r_);
        for (int lm = lms; lm <= lme; ++lm) {
          for (int i = is; i <= ie; ++i) {
            Real left_flux = left_pole ? 0.0 : -area_l_(i) * flux_x[X2DIR](lm,k,j,i);
            Real right_flux = right_pole ? 0.0 : area_r_(i) * flux_x[X2DIR](lm,k,j+1,i);
            flux_div_(lm,i) += left_flux + right_flux;
          }
        }
      }

      // Add x3-divergence
      if (ks != ke) {
        pcoord->Face3Area(k, j, is, ie, area_l_);
        pcoord->Face3Area(k+1, j, is, ie, area_r_);
        for (int lm = lms; lm <= lme; ++lm) {
          for (int i = is; i <= ie; ++i) {
            flux_div_(lm,i) += area_r_(i) * flux_x[X3DIR](lm,k+1,j,i)
                - area_l_(i) * flux_x[X3DIR](lm,k,j,i);
          }
        }
      }

      // Update conserved variables
      pcoord->CellVolume(k, j, is, ie, vol_);
      for (int lm = lms; lm <= lme; ++lm) {
        for (int i = is; i <= ie; ++i) {
          cons_out(lm,k,j,i) -= dt * flux_div_(lm,i) / vol_(i);
        }
      }
    }
  }

  // Go through all angles
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {

      // Determine poles
      bool left_pole = l == zs;
      bool right_pole = l == ze;

      // Calculate angle lengths and solid angles
      int lm = AngleInd(l, m);
      int lm_lc = AngleInd(l, m, true, false);
      int lm_rc = AngleInd(l + 1, m, true, false);
      int lm_cl = AngleInd(l, m, false, true);
      int lm_cr = AngleInd(l, m + 1, false, true);
      Real zeta_length_m = zeta_length(l,m);
      Real zeta_length_p = zeta_length(l,m+1);
      Real psi_length_m = psi_length(l,m);
      Real psi_length_p = psi_length(l+1,m);
      Real omega = solid_angle(l,m);

      // Go through all cells
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i) {

            // Calculate zeta-divergence
            Real left_flux =
                left_pole ? 0.0 : -psi_length_m * flux_a[ZETADIR](lm_lc,k,j,i);
            Real right_flux =
                right_pole ? 0.0 : psi_length_p * flux_a[ZETADIR](lm_rc,k,j,i);
            Real flux_div = left_flux + right_flux;

            // Add psi-divergence
            flux_div += zeta_length_p * flux_a[PSIDIR](lm_cr,k,j,i)
                - zeta_length_m * flux_a[PSIDIR](lm_cl,k,j,i);

            // Update conserved variables
            cons_out(lm,k,j,i) -= dt * flux_div / omega;
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Radiation conversion from primitive to conserved variables
// Inputs:
//   prim_in: primitives
//   pcoord: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   cons_out: conserved quantities

void Radiation::PrimitiveToConserved(const AthenaArray<Real> &prim_in,
    AthenaArray<Real> &cons_out, Coordinates *pcoord, int il, int iu, int jl, int ju,
    int kl, int ku) {
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            cons_out(lm,k,j,i) = n0_n_mu_(0,l,m,k,j,i) * prim_in(lm,k,j,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Radiation inversion from conserved to primitive variables
// Inputs:
//   cons_in: conserved quantities
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   prim_out: primitives
// Notes:
//   Primitives are floored at 0. Conserved quantities are adjusted to match.
//   This should be the only place where angular ghost zones need to be set.

void Radiation::ConservedToPrimitive(AthenaArray<Real> &cons_in,
    AthenaArray<Real> &prim_out, int il, int iu, int jl, int ju, int kl, int ku) {

  // Calculate primitive intensities
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim_out(lm,k,j,i) = cons_in(lm,k,j,i) / n0_n_mu_(0,l,m,k,j,i);
            if (prim_out(lm,k,j,i) < 0.0) {
              prim_out(lm,k,j,i) = 0.0;
              cons_in(lm,k,j,i) = 0.0;
            }
          }
        }
      }
    }
  }

  // Populate angular ghost zones in azimuthal angle
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps-NGHOST; m <= ps-1; ++m) {
      int m_src = pe - ps + 1 + m;
      int lm = AngleInd(l, m);
      int lm_src = AngleInd(l, m_src);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim_out(lm,k,j,i) = prim_out(lm_src,k,j,i);
          }
        }
      }
    }
    for (int m = pe+1; m <= pe+NGHOST; ++m) {
      int m_src = ps - pe - 1 + m;
      int lm = AngleInd(l, m);
      int lm_src = AngleInd(l, m_src);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim_out(lm,k,j,i) = prim_out(lm_src,k,j,i);
          }
        }
      }
    }
  }

  // Populate angular ghost zones in polar angle
  for (int l = zs-NGHOST; l <= zs-1; ++l) {
    for (int m = ps-NGHOST; m <= pe+NGHOST; ++m) {
      int l_src = 2*zs - 1 - l;
      int m_src = (m + npsi/2) % (npsi + 2*NGHOST);
      int lm = AngleInd(l, m);
      int lm_src = AngleInd(l_src, m_src);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim_out(lm,k,j,i) = prim_out(lm_src,k,j,i);
          }
        }
      }
    }
  }
  for (int l = ze+1; l <= ze+NGHOST; ++l) {
    for (int m = ps-NGHOST; m <= pe+NGHOST; ++m) {
      int l_src = 2*ze + 1 - l;
      int m_src = (m + npsi/2) % (npsi + 2*NGHOST);
      int lm = AngleInd(l, m);
      int lm_src = AngleInd(l_src, m_src);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim_out(lm,k,j,i) = prim_out(lm_src,k,j,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Radiation inversion from conserved to primitive variables, including setting moments
// Inputs:
//   cons_in: conserved quantities
//   prim_hydro: up-to-date primitive hydro quantities
//   pcoord: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   prim_out: primitives
// Notes:
//   Defers to ConservedToPrimitive() for actual inversion.
//   Updates all relevant moments to match updated radiation variables.

void Radiation::ConservedToPrimitiveWithMoments(AthenaArray<Real> &cons_in,
    AthenaArray<Real> &prim_out, const AthenaArray<Real> &prim_hydro, Coordinates *pcoord,
    int il, int iu, int jl, int ju, int kl, int ku) {
  ConservedToPrimitive(cons_in, prim_out, il, iu, jl, ju, kl, ku);
  SetMoments(prim_hydro, pcoord, il, iu, jl, ju, kl, ku);
  return;
}

//----------------------------------------------------------------------------------------
// Function for adding all source terms beyond those induced by coordinates
// Inputs:
//   time: time of simulation
//   dt: simulation timestep
//   prim_rad: primitive intensity at beginning of stage
//   prim_hydro: primitive hydro variables at beginning of stage
//   cons_rad: conserved intensity after stage integration
//   cons_hydro: conserved hydro variables after stage integration
// Outputs:
//   cons: conserved intensity updated
//   cons_hydro: conserved hydro variables updated

void Radiation::AddSourceTerms(const Real time, const Real dt,
    const AthenaArray<Real> &prim_rad, const AthenaArray<Real> &prim_hydro,
    AthenaArray<Real> &cons_rad, AthenaArray<Real> &cons_hydro) {

  // Go through outer loops of cells
  if (coupled_to_matter) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        pmy_block->pcoord->CellMetric(k, j, is, ie, g_, gi_);

        // Calculate zeroth and first moments of radiation before coupling
        for (int n = 0; n < 4; ++n) {
          for (int i = is; i <= ie; ++i) {
            moments_old_(n,i) = 0.0;
          }
        }
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            int lm = AngleInd(l, m);
            for (int n = 0; n < 4; ++n) {
              for (int i = is; i <= ie; ++i) {
                moments_old_(n,i) += n0_n_mu_(n,l,m,k,j,i) * cons_rad(lm,k,j,i)
                    / n0_n_mu_(0,l,m,k,j,i) * solid_angle(l,m);
              }
            }
          }
        }

        // Calculate fluid velocity in tetrad frame, as well as time steps
        for (int i = is; i <= ie; ++i) {
          Real uu1 = prim_hydro(IVX,k,j,i);
          Real uu2 = prim_hydro(IVY,k,j,i);
          Real uu3 = prim_hydro(IVZ,k,j,i);
          Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
              + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22) * SQR(uu2)
              + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
          Real uu0 = std::sqrt(1.0 + temp_var);
          u_tet_(0,i) = norm_to_tet_(0,0,k,j,i) * uu0 + norm_to_tet_(0,1,k,j,i) * uu1
              + norm_to_tet_(0,2,k,j,i) * uu2 + norm_to_tet_(0,3,k,j,i) * uu3;
          u_tet_(1,i) = norm_to_tet_(1,0,k,j,i) * uu0 + norm_to_tet_(1,1,k,j,i) * uu1
              + norm_to_tet_(1,2,k,j,i) * uu2 + norm_to_tet_(1,3,k,j,i) * uu3;
          u_tet_(2,i) = norm_to_tet_(2,0,k,j,i) * uu0 + norm_to_tet_(2,1,k,j,i) * uu1
              + norm_to_tet_(2,2,k,j,i) * uu2 + norm_to_tet_(2,3,k,j,i) * uu3;
          u_tet_(3,i) = norm_to_tet_(3,0,k,j,i) * uu0 + norm_to_tet_(3,1,k,j,i) * uu1
              + norm_to_tet_(3,2,k,j,i) * uu2 + norm_to_tet_(3,3,k,j,i) * uu3;
          Real u0 = uu0 * std::sqrt(-gi_(I00,i));
          dt_(i) = dt;
          dtau_(i) = dt / u0;
        }

        // Transform radiation from tetrad to fluid frame
        for (int i = is; i <= ie; ++i) {
          weight_sum_(i) = 0.0;
        }
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            int lm = AngleInd(l, m);
            int lm_alt = (l - zs) * (pe - ps + 1) + m - ps;
            for (int i = is; i <= ie; ++i) {
              Real un_tet = u_tet_(1,i) * nh_cc_(1,l,m) + u_tet_(2,i) * nh_cc_(2,l,m)
                  + u_tet_(3,i) * nh_cc_(3,l,m);
              n_cm_(0,lm_alt,i) = u_tet_(0,i) * nh_cc_(0,l,m) - un_tet;
              n_cm_(1,lm_alt,i) = -u_tet_(1,i) * nh_cc_(0,l,m)
                  + u_tet_(1,i) / (u_tet_(0,i) + 1.0) * un_tet + nh_cc_(1,l,m);
              n_cm_(2,lm_alt,i) = -u_tet_(2,i) * nh_cc_(0,l,m)
                  + u_tet_(2,i) / (u_tet_(0,i) + 1.0) * un_tet + nh_cc_(2,l,m);
              n_cm_(3,lm_alt,i) = -u_tet_(3,i) * nh_cc_(0,l,m)
                  + u_tet_(3,i) / (u_tet_(0,i) + 1.0) * un_tet + nh_cc_(3,l,m);
              omega_cm_(lm_alt,i) = solid_angle(l,m) / SQR(n_cm_(0,lm_alt,i));
              intensity_cm_(lm_alt,i) = prim_rad(lm,k,j,i) * SQR(SQR(n_cm_(0,lm_alt,i)));
              weight_sum_(i) += omega_cm_(lm_alt,i);
            }
          }
        }

        // Calculate radiation-fluid coupling in fluid frame
        for (int n = 0; n < nzeta * npsi; ++n) {
          for (int i = is; i <= ie; ++i) {
            omega_cm_(n,i) /= weight_sum_(i);
            intensity_cm_(n,i) *= 4.0*PI;
          }
        }
        Coupling(prim_hydro, n_cm_, n0_, omega_cm_, dt_, dtau_, k, j, intensity_cm_);
        for (int n = 0; n < nzeta * npsi; ++n) {
          for (int i = is; i <= ie; ++i) {
            intensity_cm_(n,i) /= 4.0*PI;
          }
        }

        // Apply radiation-fluid coupling to radiation in coordinate frame
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            int lm = AngleInd(l, m);
            int lm_alt = (l - zs) * (pe - ps + 1) + m - ps;
            for (int i = is; i <= ie; ++i) {
              Real intensity_coord =
                  intensity_cm_(lm_alt,i) / SQR(SQR(n_cm_(0,lm_alt,i)));
              cons_rad(lm,k,j,i) = intensity_coord * n0_n_mu_(0,l,m,k,j,i);
            }
          }
        }

        // Calculate zeroth and first moments of radiation after coupling
        for (int n = 0; n < 4; ++n) {
          for (int i = is; i <= ie; ++i) {
            moments_new_(n,i) = 0.0;
          }
        }
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            int lm = AngleInd(l, m);
            for (int n = 0; n < 4; ++n) {
              for (int i = is; i <= ie; ++i) {
                moments_new_(n,i) += n0_n_mu_(n,l,m,k,j,i) * cons_rad(lm,k,j,i)
                    / n0_n_mu_(0,l,m,k,j,i) * solid_angle(l,m);
              }
            }
          }
        }

        // Apply radiation-fluid coupling to fluid
        for (int i = is; i <= ie; ++i) {
          cons_hydro(IEN,k,j,i) += moments_old_(0,i) - moments_new_(0,i);
          cons_hydro(IM1,k,j,i) += moments_old_(1,i) - moments_new_(1,i);
          cons_hydro(IM2,k,j,i) += moments_old_(2,i) - moments_new_(2,i);
          cons_hydro(IM3,k,j,i) += moments_old_(3,i) - moments_new_(3,i);
        }
      }
    }
  }

  // Apply user source terms
  if (UserSourceTerm != nullptr) {
    UserSourceTerm(pmy_block, time, dt, prim_rad, cons_rad);
  }
  return;
}

//----------------------------------------------------------------------------------------
// Opacity enrollment
// Inputs:
//   MyOpacityFunction: user-defined function from problem generator
// Outputs: (none)
// Notes:
//   If nothing else enrolled, default function keeps opacities (not absorption
//     coefficients) at their initial values.

void Radiation::EnrollOpacityFunction(OpacityFunc MyOpacityFunction)
{
  UpdateOpacity = MyOpacityFunction;
  return;
}

//----------------------------------------------------------------------------------------
// Indexing function for angles
// Inputs:
//   l: zeta-index
//   m: psi-index
//   zeta_face: flag indicating zeta-index is on faces
//   psi_face: flag indicating psi-index is on faces
// Outputs:
//   returned value: 1D index for both zeta and psi
// Notes:
//   More general version of RadBoundaryVariable::AngleInd().

int Radiation::AngleInd(int l, int m, bool zeta_face, bool psi_face) {
  if (psi_face) {
    return l * (npsi + 2*NGHOST + 1) + m;
  }
  return l * (npsi + 2*NGHOST) + m;
}

//----------------------------------------------------------------------------------------
// Function for calculating conserved intensity corresponding to beam source
// Inputs:
//   pos_1, pos_2, pos_3: coordinates of beam origin
//   width: full proper diameter of beam
//   dir_1, dir_2, dir_3: relative directions of beam center
//   spread: full spread of beam in direction, in degrees
//   dii_dt: injected I per unit time
//   cylindrical: flag indicating coordinates are cylindrical
//   spherical: flag indicating coordinates are spherical
// Outputs:
//   dcons_dt: conserved values (n^0 n_0 I) per unit time set
// Notes:
//   Arrays should be 4D, with first index holding both zeta and psi.
//   Cylindrical coordinates:
//     phi (x2) will be mapped to [-pi, pi]:
//       Beams near phi = 0 should work fine.
//       Beams near phi = pi will likely not be initialized correctly.
//   Spherical coordinates:
//     theta (x2) will be mapped to [0, pi], adjusting phi (x3) as necessary.
//     phi (x3) will be mapped to [-pi, pi]:
//       Beams near phi = 0 should work fine.
//       Beams near phi = pi will likely not be initialized correctly.

void Radiation::CalculateBeamSource(Real pos_1, Real pos_2, Real pos_3, Real width,
    Real dir_1, Real dir_2, Real dir_3, Real spread, Real dii_dt,
    AthenaArray<Real> &dcons_dt, bool cylindrical, bool spherical) {

  // Account for cylindrical/spherical coordinates in beam origin
  if (cylindrical) {
    if (pos_2 > PI) {
      pos_2 -= 2.0*PI;
    }
  }
  if (spherical) {
    if (pos_2 < 0.0) {
      pos_2 = -pos_2;
      pos_3 -= PI;
    }
    if (pos_2 > PI) {
      pos_2 = 2.0*PI - pos_2;
      pos_3 -= PI;
    }
    if (pos_3 > PI) {
      pos_3 -= 2.0*PI;
    }
  }

  // Allocate scratch arrays
  AthenaArray<Real> g, gi, e, e_cov, omega, nh;
  g.NewAthenaArray(NMETRIC, ie + 1);
  gi.NewAthenaArray(NMETRIC, ie + 1);
  e.NewAthenaArray(4, 4);
  e_cov.NewAthenaArray(4, 4);
  omega.NewAthenaArray(4, 4, 4);
  nh.NewAthenaArray(ze + 1, pe + 1, 4);

  // Calculate unit normal components in orthonormal frame
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      nh(l,m,0) = 1.0;
      nh(l,m,1) = std::sin(zetav(l)) * std::cos(psiv(m));
      nh(l,m,2) = std::sin(zetav(l)) * std::sin(psiv(m));
      nh(l,m,3) = std::cos(zetav(l));
    }
  }

  // Calculate minimum angle between directions
  Real mu_min = std::cos(spread/2.0 * PI/180.0);

  // Go through cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pmy_block->pcoord->CellMetric(k, j, is, ie, g, gi);
      for (int i = is; i <= ie; ++i) {

        // Extract position, accounting for cylindrical/spherical coordinates
        Real x1 = pmy_block->pcoord->x1v(i);
        Real x2 = pmy_block->pcoord->x2v(j);
        Real x3 = pmy_block->pcoord->x3v(k);
        if (cylindrical) {
          if (x2 > PI) {
            x2 -= 2.0*PI;
          }
        }
        if (spherical) {
          if (x2 < 0.0) {
            x2 = -x2;
            x3 -= PI;
          }
          if (x2 > PI) {
            x2 = 2.0*PI - x2;
            x3 -= PI;
          }
          if (x3 > PI) {
            x3 -= 2.0*PI;
          }
        }

        // Calculate proper distance to beam origin
        Real dx1 = x1 - pos_1;
        Real dx2 = x2 - pos_2;
        Real dx3 = x3 - pos_3;
        Real dx_sq = g(I11,i) * SQR(dx1) + 2.0 * g(I12,i) * dx1 * dx2
            + 2.0 * g(I13,i) * dx1 * dx3 + g(I22,i) * SQR(dx2)
            + 2.0 * g(I23,i) * dx2 * dx3 + g(I33,i) * SQR(dx3);

        // Set to 0 if too far from beam in space
        if (dx_sq >= SQR(width/2.0)) {
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              dcons_dt(lm,k,j,i) = 0.0;
            }
          }
          continue;
        }

        // Calculate tetrad
        pmy_block->pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);

        // Calculate contravariant time component of direction
        Real temp_a = g(I00,i);
        Real temp_b = 2.0 * (g(I01,i) * dir_1 + g(I02,i) * dir_2 + g(I03,i) * dir_3);
        Real temp_c = g(I11,i) * SQR(dir_1) + 2.0 * g(I12,i) * dir_1 * dir_2
            + 2.0 * g(I13,i) * dir_1 * dir_3 + g(I22,i) * SQR(dir_2)
            + 2.0 * g(I23,i) * dir_2 * dir_3 + g(I33,i) * SQR(dir_3);
        Real dir_0 =
            (-temp_b - std::sqrt(SQR(temp_b) - 4.0 * temp_a * temp_c)) / (2.0 * temp_a);

        // Calculate covariant direction
        Real dir_cov_0, dir_cov_1, dir_cov_2, dir_cov_3;
        pmy_block->pcoord->LowerVectorCell(dir_0, dir_1, dir_2, dir_3, k, j, i,
            &dir_cov_0, &dir_cov_1, &dir_cov_2, &dir_cov_3);

        // Calculate covariant direction in tetrad frame
        Real dir_tet_cov_0 = e(0,0) * dir_cov_0 + e(0,1) * dir_cov_1 + e(0,2) * dir_cov_2
            + e(0,3) * dir_cov_3;
        Real dir_tet_cov_1 = e(1,0) * dir_cov_0 + e(1,1) * dir_cov_1 + e(1,2) * dir_cov_2
            + e(1,3) * dir_cov_3;
        Real dir_tet_cov_2 = e(2,0) * dir_cov_0 + e(2,1) * dir_cov_1 + e(2,2) * dir_cov_2
            + e(2,3) * dir_cov_3;
        Real dir_tet_cov_3 = e(3,0) * dir_cov_0 + e(3,1) * dir_cov_1 + e(3,2) * dir_cov_2
            + e(3,3) * dir_cov_3;

        // Normalize covariant spatial direction in tetrad frame
        dir_tet_cov_1 /= -dir_tet_cov_0;
        dir_tet_cov_2 /= -dir_tet_cov_0;
        dir_tet_cov_3 /= -dir_tet_cov_0;

        // Go through angles
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {

            // Calculate combined angle index
            int lm = AngleInd(l, m);

            // Calculate angle to beam direction
            Real mu = nh(l,m,1) * dir_tet_cov_1 + nh(l,m,2) * dir_tet_cov_2
                + nh(l,m,3) * dir_tet_cov_3;

            // Set to 0 if too far from beam in angle
            if (mu <= mu_min) {
              dcons_dt(lm,k,j,i) = 0.0;
              continue;
            }

            // Set to nonzero value
            dcons_dt(lm,k,j,i) = dii_dt * n0_n_mu_(0,l,m,k,j,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating conserved intensity corresponding to spatially constant source
// Inputs:
//   energy: coordinate-frame energy density
//   u1, u2, u3: contravariant 4-velocity components of isotropic radiation frame
// Outputs:
//   cons_out: conserved values (n^0 n_0 I) set

void Radiation::CalculateConstantRadiation(Real energy, Real u1, Real u2, Real u3,
    AthenaArray<Real> &cons_out) {

  // Allocate scratch arrays
  AthenaArray<Real> g, gi;
  g.NewAthenaArray(NMETRIC, ie + 1);
  gi.NewAthenaArray(NMETRIC, ie + 1);

  // Go through cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      pmy_block->pcoord->CellMetric(k, j, is, ie, g, gi);
      for (int i = is; i <= ie; ++i) {
        CalculateRadiationInCell(energy, u1, u2, u3, k, j, i, g, cons_out);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating conserved intensity corresponding energy density and velocity
// Inputs:
//   energy: coordinate-frame energy density
//   u1, u2, u3: contravariant 4-velocity components of isotropic radiation frame
//   k, j, i: indices for cell to set
//   g: covariant metric
// Outputs:
//   cons_out: conserved values (n^0 n_0 I) set in given cell

void Radiation::CalculateRadiationInCell(Real energy, Real u1, Real u2, Real u3, int k,
    int j, int i, const AthenaArray<Real> &g, AthenaArray<Real> &cons_out) {

  // Calculate contravariant time component of isotropic radiation frame velocity
  Real temp_a = g(I00,i);
  Real temp_b = 2.0 * (g(I01,i) * u1 + g(I02,i) * u2 + g(I03,i) * u3);
  Real temp_c = g(I11,i) * SQR(u1) + 2.0 * g(I12,i) * u1 * u2
    + 2.0 * g(I13,i) * u1 * u3 + g(I22,i) * SQR(u2) + 2.0 * g(I23,i) * u2 * u3
    + g(I33,i) * SQR(u3) + 1.0;
  Real temp_d = std::max(SQR(temp_b) - 4.0 * temp_a * temp_c, 0.0);
  Real u0 =
      (-temp_b - std::sqrt(temp_d)) / (2.0 * temp_a);

  // Calculate covariant radiation velocity
  Real u_0, u_1, u_2, u_3;
  pmy_block->pcoord->LowerVectorCell(u0, u1, u2, u3, k, j, i, &u_0, &u_1, &u_2,
      &u_3);

  // Set conserved quantity at each angle, tracking energy density
  Real energy_sum = 0.0;
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      Real u_n = u_0 * nmu_(0,l,m,k,j,i) + u_1 * nmu_(1,l,m,k,j,i)
          + u_2 * nmu_(2,l,m,k,j,i) + u_3 * nmu_(3,l,m,k,j,i);
      Real ii = 1.0 / SQR(SQR(-u_n));
      cons_out(lm,k,j,i) = n0_n_mu_(0,l,m,k,j,i) * ii;
      energy_sum += SQR(nmu_(0,l,m,k,j,i)) * ii * solid_angle(l,m);
    }
  }

  // Normalize conserved quantities
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      cons_out(lm,k,j,i) *= energy / energy_sum;
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating moments of radiation field
// Inputs:
//   prim_hydro: up-to-date primitive hydro quantities
//   pcoord: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region for which moments should be calculated
// Outputs: (none)
// Notes:
//   Populates moments_coord array with 10 components of radiation stress tensor in
//       coordinate frame.
//   Populates moments_tetrad array with 10 components of radiation stress tensor in
//       orthonormal tetrad frame.
//   Populates moments_fluid array with 10 components of radiation stress tensor in
//       orthonormal fluid frame if radiation is coupled to matter.

void Radiation::SetMoments(const AthenaArray<Real> &prim_hydro, Coordinates *pcoord,
    int il, int iu, int jl, int ju, int kl, int ku) {

  // Zero moments
  for (int n = 0; n < 10; ++n) {
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = il; i <= iu; ++i) {
          moments_coord(n,k,j,i) = 0.0;
          moments_tetrad(n,k,j,i) = 0.0;
          if (coupled_to_matter) {
            moments_fluid(n,k,j,i) = 0.0;
          }
        }
      }
    }
  }

  // Set coordinate-frame components
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int n1 = 0, n12 = 0; n1 < 4; ++n1) {
        for (int n2 = n1; n2 < 4; ++n2, ++n12) {
          for (int k = kl; k <= ku; ++k) {
            for (int j = jl; j <= ju; ++j) {
              for (int i = il; i <= iu; ++i) {
                moments_coord(n12,k,j,i) += nmu_(n1,l,m,k,j,i) * nmu_(n2,l,m,k,j,i)
                    * prim(lm,k,j,i) * solid_angle(l,m);
              }
            }
          }
        }
      }
    }
  }

  // Set tetrad-frame components
  AthenaArray<Real> e, e_cov, omega;
  e.NewAthenaArray(4, 4);
  e_cov.NewAthenaArray(4, 4);
  omega.NewAthenaArray(4, 4, 4);
  Real moments_coord_full[4][4];
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        pcoord->Tetrad(x1, x2, x3, e, e_cov, omega);
        moments_coord_full[0][0] = moments_coord(0,k,j,i);
        moments_coord_full[0][1] = moments_coord_full[1][0] = moments_coord(1,k,j,i);
        moments_coord_full[0][2] = moments_coord_full[2][0] = moments_coord(2,k,j,i);
        moments_coord_full[0][3] = moments_coord_full[3][0] = moments_coord(3,k,j,i);
        moments_coord_full[1][1] = moments_coord(4,k,j,i);
        moments_coord_full[1][2] = moments_coord_full[2][1] = moments_coord(5,k,j,i);
        moments_coord_full[1][3] = moments_coord_full[3][1] = moments_coord(6,k,j,i);
        moments_coord_full[2][2] = moments_coord(7,k,j,i);
        moments_coord_full[2][3] = moments_coord_full[3][2] = moments_coord(8,k,j,i);
        moments_coord_full[3][3] = moments_coord(9,k,j,i);
        for (int n1 = 0, n12 = 0; n1 < 4; ++n1) {
          for (int n2 = n1; n2 < 4; ++n2, ++n12) {
            for (int m1 = 0; m1 < 4; ++m1) {
              for (int m2 = 0; m2 < 4; ++m2) {
                moments_tetrad(n12,k,j,i) +=
                    e_cov(n1,m1) * e_cov(n2,m2) * moments_coord_full[m1][m2];
              }
            }
          }
        }
        moments_tetrad(1,k,j,i) *= -1.0;
        moments_tetrad(2,k,j,i) *= -1.0;
        moments_tetrad(3,k,j,i) *= -1.0;
      }
    }
  }

  // Set fluid-frame components
  if (coupled_to_matter) {
    Real tet_to_fluid[4][4];
    Real moments_tetrad_full[4][4];
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        pmy_block->pcoord->CellMetric(k, j, il, iu, g_, gi_);
        for (int i = il; i <= iu; ++i) {

          // Calculate fluid velocity in tetrad frame
          Real uu1 = prim_hydro(IVX,k,j,i);
          Real uu2 = prim_hydro(IVY,k,j,i);
          Real uu3 = prim_hydro(IVZ,k,j,i);
          Real temp_var = g_(I11,i) * SQR(uu1) + 2.0 * g_(I12,i) * uu1 * uu2
              + 2.0 * g_(I13,i) * uu1 * uu3 + g_(I22) * SQR(uu2)
              + 2.0 * g_(I23,i) * uu2 * uu3 + g_(I33,i) * SQR(uu3);
          Real uu0 = std::sqrt(1.0 + temp_var);
          Real utet0 = norm_to_tet_(0,0,k,j,i) * uu0 + norm_to_tet_(0,1,k,j,i) * uu1
              + norm_to_tet_(0,2,k,j,i) * uu2 + norm_to_tet_(0,3,k,j,i) * uu3;
          Real utet1 = norm_to_tet_(1,0,k,j,i) * uu0 + norm_to_tet_(1,1,k,j,i) * uu1
              + norm_to_tet_(1,2,k,j,i) * uu2 + norm_to_tet_(1,3,k,j,i) * uu3;
          Real utet2 = norm_to_tet_(2,0,k,j,i) * uu0 + norm_to_tet_(2,1,k,j,i) * uu1
              + norm_to_tet_(2,2,k,j,i) * uu2 + norm_to_tet_(2,3,k,j,i) * uu3;
          Real utet3 = norm_to_tet_(3,0,k,j,i) * uu0 + norm_to_tet_(3,1,k,j,i) * uu1
              + norm_to_tet_(3,2,k,j,i) * uu2 + norm_to_tet_(3,3,k,j,i) * uu3;

          // Construct Lorentz boost from tetrad frame to orthonormal fluid frame
          tet_to_fluid[0][0] = utet0;
          tet_to_fluid[0][1] = tet_to_fluid[1][0] = -utet1;
          tet_to_fluid[0][2] = tet_to_fluid[2][0] = -utet2;
          tet_to_fluid[0][3] = tet_to_fluid[3][0] = -utet3;
          tet_to_fluid[1][1] = SQR(utet1) / (1.0 + utet0) + 1.0;
          tet_to_fluid[1][2] = tet_to_fluid[2][1] = utet1 * utet2 / (1.0 + utet0);
          tet_to_fluid[1][3] = tet_to_fluid[3][1] = utet1 * utet3 / (1.0 + utet0);
          tet_to_fluid[2][2] = SQR(utet2) / (1.0 + utet0) + 1.0;
          tet_to_fluid[2][3] = tet_to_fluid[3][2] = utet2 * utet3 / (1.0 + utet0);
          tet_to_fluid[3][3] = SQR(utet3) / (1.0 + utet0) + 1.0;

          // Transform moments
          moments_tetrad_full[0][0] = moments_tetrad(0,k,j,i);
          moments_tetrad_full[0][1] = moments_tetrad_full[1][0] = moments_tetrad(1,k,j,i);
          moments_tetrad_full[0][2] = moments_tetrad_full[2][0] = moments_tetrad(2,k,j,i);
          moments_tetrad_full[0][3] = moments_tetrad_full[3][0] = moments_tetrad(3,k,j,i);
          moments_tetrad_full[1][1] = moments_tetrad(4,k,j,i);
          moments_tetrad_full[1][2] = moments_tetrad_full[2][1] = moments_tetrad(5,k,j,i);
          moments_tetrad_full[1][3] = moments_tetrad_full[3][1] = moments_tetrad(6,k,j,i);
          moments_tetrad_full[2][2] = moments_tetrad(7,k,j,i);
          moments_tetrad_full[2][3] = moments_tetrad_full[3][2] = moments_tetrad(8,k,j,i);
          moments_tetrad_full[3][3] = moments_tetrad(9,k,j,i);
          for (int n1 = 0, n12 = 0; n1 < 4; ++n1) {
            for (int n2 = n1; n2 < 4; ++n2, ++n12) {
              for (int m1 = 0; m1 < 4; ++m1) {
                for (int m2 = 0; m2 < 4; ++m2) {
                  moments_fluid(n12,k,j,i) += tet_to_fluid[n1][m1] * tet_to_fluid[n2][m2]
                      * moments_tetrad_full[m1][m2];
                }
              }
            }
          }
        }
      }
    }
  }
  return;
}

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
