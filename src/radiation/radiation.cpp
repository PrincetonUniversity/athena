//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of core functions in class Radiation

// C++ headers
#include <algorithm>  // max, min
#include <cstdlib>    // abs
#include <cmath>      // acos, cos, isnan, NAN, sin, sqrt
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str, string

// Athena++ headers
#include "radiation.hpp"
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../mesh/mesh.hpp"                // MeshBlock

// Declarations
void DefaultOpacity(MeshBlock *pmb, const AthenaArray<Real> &prim_hydro);
bool FourthPolyRoot(const Real coef4, const Real tconst, Real *root);

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

  // Allocate memory for flux calculation
  rad_l_.NewAthenaArray(nang_zpf, pmb->ncells1 + 1);
  rad_r_.NewAthenaArray(nang_zpf, pmb->ncells1 + 1);

  // Allocate memory for flux divergence calculation
  area_l_.NewAthenaArray(pmb->ncells1 + 1);
  area_r_.NewAthenaArray(pmb->ncells1 + 1);
  vol_.NewAthenaArray(pmb->ncells1 + 1);
  flux_div_.NewAthenaArray(nang, pmb->ncells1 + 1);

  // Allocate memory for source term calculation
  g_.NewAthenaArray(NMETRIC, pmb->ncells1);
  gi_.NewAthenaArray(NMETRIC, pmb->ncells1);
  norm_to_tet_.NewAthenaArray(4, 4, pmb->ncells3, pmb->ncells2, pmb->ncells1);
  moments_old_.NewAthenaArray(4, pmb->ncells1);
  moments_new_.NewAthenaArray(4, pmb->ncells1);
  u_tet_.NewAthenaArray(4, pmb->ncells1);
  coefficients_.NewAthenaArray(2, pmb->ncells1);
  bad_cell_.NewAthenaArray(pmb->ncells1);
  tt_plus_.NewAthenaArray(pmb->ncells1);
  ee_f_minus_.NewAthenaArray(pmb->ncells1);
  ee_f_plus_.NewAthenaArray(pmb->ncells1);

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

  // Get adiabatic index
  Real gamma_adi = pmy_block->peos->GetGamma();

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

        // Calculate fluid velocity in tetrad frame
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
        }

        // Calculate quartic coefficients
        for (int i = is; i <= ie; ++i) {
          Real rho = prim_hydro(IDN,k,j,i);
          Real tt_minus = prim_hydro(IPR,k,j,i) / rho;
          Real k_a = rho * opacity(OPAA,k,j,i);
          Real k_s = rho * opacity(OPAS,k,j,i);
          Real k_tot = k_a + k_s;
          Real var_a = 0.0;
          Real var_b = 0.0;
          ee_f_minus_(i) = 0.0;
          for (int l = zs; l <= ze; ++l) {
            for (int m = ps; m <= pe; ++m) {
              int lm = AngleInd(l, m);
              Real ii_minus = prim_rad(lm,k,j,i);
              Real u_n = -u_tet_(0,i) * nh_cc_(0,l,m) + u_tet_(1,i) * nh_cc_(1,l,m)
                  + u_tet_(2,i) * nh_cc_(2,l,m) + u_tet_(3,i) * nh_cc_(3,l,m);
              Real denominator = 1.0 - dt * k_tot * u_n;
              var_a += ii_minus * SQR(u_n) / denominator * solid_angle(l,m);
              var_b += 1.0 / u_n / denominator * solid_angle(l,m);
              ee_f_minus_(i) += ii_minus * SQR(u_n) * solid_angle(l,m);
            }
          }
          var_b *= dt / (4.0 * PI);
          coefficients_(0,i) =
              -(gamma_adi - 1.0) / rho * var_b * k_a * arad / (1.0 + var_b * k_s);
          coefficients_(1,i) = -tt_minus - (gamma_adi - 1.0) / rho * ee_f_minus_(i)
              + (gamma_adi - 1.0) / rho * var_a / (1.0 + var_b * k_s);
        }

        // Calculate new gas temperature
        for (int i = is; i <= ie; ++i) {
          bad_cell_(i) = false;
          if (std::abs(coefficients_(0,i)) > TINY_NUMBER) {
            bool quartic_flag =
                FourthPolyRoot(coefficients_(0,i), coefficients_(1,i), &tt_plus_(i));
            if (not quartic_flag or std::isnan(tt_plus_(i))) {
              bad_cell_(i) = true;
              tt_plus_(i) = prim_hydro(IPR,k,j,i) / prim_hydro(IDN,k,j,i);
            }
          } else {
            tt_plus_(i) = -coefficients_(1,i);
          }
        }

        // Calculate new radiation energy density
        for (int i = is; i <= ie; ++i) {
          if (not bad_cell_(i)) {
            Real rho = prim_hydro(IDN,k,j,i);
            Real tt_minus = prim_hydro(IPR,k,j,i) / rho;
            ee_f_plus_(i) =
                ee_f_minus_(i) + rho / (gamma_adi - 1.0) * (tt_minus - tt_plus_(i));
            ee_f_plus_(i) = std::max(ee_f_plus_(i), 0.0);
          }
        }

        // Calculate new intensity
        for (int i = is; i <= ie; ++i) {
          if (not bad_cell_(i)) {
            Real rho = prim_hydro(IDN,k,j,i);
            Real k_a = rho * opacity(OPAA,k,j,i);
            Real k_s = rho * opacity(OPAS,k,j,i);
            Real k_tot = k_a + k_s;
            for (int l = zs; l <= ze; ++l) {
              for (int m = ps; m <= pe; ++m) {
                int lm = AngleInd(l, m);
                Real ii_minus = prim_rad(lm,k,j,i);
                Real u_n = -u_tet_(0,i) * nh_cc_(0,l,m) + u_tet_(1,i) * nh_cc_(1,l,m)
                    + u_tet_(2,i) * nh_cc_(2,l,m) + u_tet_(3,i) * nh_cc_(3,l,m);
                Real ii_plus = (ii_minus - dt / (4.0 * PI) / u_n / SQR(u_n)
                    * (k_a * arad * SQR(SQR(tt_plus_(i))) + k_s * ee_f_plus_(i)))
                    / (1.0 - dt * k_tot * u_n);
                cons_rad(lm,k,j,i) += (ii_plus - ii_minus) * n0_n_mu_(0,l,m,k,j,i);
                cons_rad(lm,k,j,i) = std::min(cons_rad(lm,k,j,i), 0.0);
              }
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
//   TODO: move to radiation.hpp?

int Radiation::AngleInd(int l, int m, bool zeta_face, bool psi_face) {
  if (psi_face) {
    return l * (npsi + 2*NGHOST + 1) + m;
  }
  return l * (npsi + 2*NGHOST) + m;
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

//----------------------------------------------------------------------------------------
// Exact solution for fourth order polynomial
// Inputs:
//   coef4: quartic coefficient
//   tconst: constant coefficient
// Outputs:
//   root: solution to equation
//   returned value: flag indicating success
// Notes:
//   Polynomial has the form coef4 * x^4 + x + tconst = 0.

bool FourthPolyRoot(const Real coef4, const Real tconst, Real *root) {

  // Calculate real root of z^3 - 4*tconst/coef4 * z - 1/coef4^2 = 0
  Real asquar = coef4 * coef4;
  Real acubic = coef4 * asquar;
  Real ccubic = tconst * tconst * tconst;
  Real delta1 = 0.25 - 64.0 * ccubic * coef4 / 27.0;
  if (delta1 < 0.0) {
    return false;
  }
  delta1 = std::sqrt(delta1);
  if (delta1 < 0.5) {
    return false;
  }
  Real zroot;
  if (delta1 > 1.0e11) {  // to avoid small number cancellation
    zroot = std::pow(delta1, -2.0/3.0) / 3.0;
  } else {
    zroot = std::pow(0.5 + delta1, 1.0/3.0) - std::pow(-0.5 + delta1, 1.0/3.0);
  }
  if (zroot < 0.0) {
    return false;
  }
  zroot *= std::pow(coef4, -2.0/3.0);

  // Calculate quartic root using cubic root
  Real rcoef = std::sqrt(zroot);
  Real delta2 = -zroot + 2.0 / (coef4 * rcoef);
  if (delta2 < 0.0) {
    return false;
  }
  delta2 = std::sqrt(delta2);
  *root = 0.5 * (delta2 - rcoef);
  if (*root < 0.0) {
    return false;
  }
  return true;
}
