//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of functions in class Radiation

// C++ headers
#include <cmath>      // acos, cos, sin, sqrt
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str, string

// Athena++ headers
#include "radiation.hpp"
#include "../athena_arrays.hpp"            // AthenaArray
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../mesh/mesh.hpp"                // MeshBlock

//----------------------------------------------------------------------------------------
// Radiation constructor
// Inputs:
//   pmb: pointer to containing MeshBlock
//   pin: pointer to runtime parameters

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin) :
    pmy_block(pmb),
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
    rbvar(pmb, &cons, &coarse_cons, flux_x) {

  // Set object and function pointers
  UserSourceTerm = pmb->pmy_mesh->UserRadSourceTerm_;

  // Enroll refinement communication
  if (pmb->pmy_mesh->multilevel) {
    refinement_idx = pmb->pmr->AddToRefinement(&cons, &coarse_cons);
  }

  // Construct objects
  rbvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&rbvar);
  pmb->pbval->bvars_main_int.push_back(&rbvar);

  // Set flags
  if (UserSourceTerm == nullptr) {
    source_terms_defined = false;
  } else {
    source_terms_defined = true;
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

  // Allocate memory for unit normal components in orthonormal frame
  int num_cells_zeta = ze + NGHOST;
  int num_cells_psi = pe + NGHOST;
  AthenaArray<Real> nh_cc, nh_fc, nh_cf;
  nh_cc.NewAthenaArray(4, num_cells_zeta, num_cells_psi);
  nh_fc.NewAthenaArray(4, num_cells_zeta + 1, num_cells_psi);
  nh_cf.NewAthenaArray(4, num_cells_zeta, num_cells_psi + 1);

  // Calculate unit normal components in orthonormal frame at angle centers
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      nh_cc(0,l,m) = 1.0;
      nh_cc(1,l,m) = std::sin(zetav(l)) * std::cos(psiv(m));
      nh_cc(2,l,m) = std::sin(zetav(l)) * std::sin(psiv(m));
      nh_cc(3,l,m) = std::cos(zetav(l));
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
  n0_.NewAthenaArray(num_cells_zeta, num_cells_psi, pmb->ncells3, pmb->ncells2,
      pmb->ncells1);
  n1_.NewAthenaArray(num_cells_zeta, num_cells_psi, pmb->ncells3, pmb->ncells2,
      pmb->ncells1 + 1);
  n2_.NewAthenaArray(num_cells_zeta, num_cells_psi, pmb->ncells3, pmb->ncells2 + 1,
      pmb->ncells1);
  n3_.NewAthenaArray(num_cells_zeta, num_cells_psi, pmb->ncells3 + 1, pmb->ncells2,
      pmb->ncells1);
  na0_.NewAthenaArray(num_cells_zeta, num_cells_psi, pmb->ncells3, pmb->ncells2,
      pmb->ncells1);
  na1_.NewAthenaArray(num_cells_zeta + 1, num_cells_psi, pmb->ncells3, pmb->ncells2,
      pmb->ncells1);
  na2_.NewAthenaArray(num_cells_zeta, num_cells_psi + 1, pmb->ncells3, pmb->ncells2,
      pmb->ncells1);

  // Allocate memory for tetrad and rotation coefficients
  AthenaArray<Real> e, omega;
  e.NewAthenaArray(4, 4);
  omega.NewAthenaArray(4, 4, 4);

  // Calculate n^0
  for (int k = ks-NGHOST; k <= ke+NGHOST; ++k) {
    for (int j = js-NGHOST; j <= je+NGHOST; ++j) {
      for (int i = is-NGHOST; i <= ie+NGHOST; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, omega);
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            n0_(l,m,k,j,i) = 0.0;
            for (int n = 0; n < 4; ++n) {
              n0_(l,m,k,j,i) += e(n,0) * nh_cc(n,l,m);
            }
          }
        }
      }
    }
  }

  // Calculate n^1
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie+1; ++i) {
        Real x1 = pmb->pcoord->x1f(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, omega);
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            n1_(l,m,k,j,i) = 0.0;
            for (int n = 0; n < 4; ++n) {
              n1_(l,m,k,j,i) += e(n,1) * nh_cc(n,l,m);
            }
          }
        }
      }
    }
  }

  // Calculate n^2
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je+1; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2f(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, omega);
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            n2_(l,m,k,j,i) = 0.0;
            for (int n = 0; n < 4; ++n) {
              n2_(l,m,k,j,i) += e(n,2) * nh_cc(n,l,m);
            }
          }
        }
      }
    }
  }

  // Calculate n^3
  for (int k = ks; k <= ke+1; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3f(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, omega);
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            n3_(l,m,k,j,i) = 0.0;
            for (int n = 0; n < 4; ++n) {
              n3_(l,m,k,j,i) += e(n,3) * nh_cc(n,l,m);
            }
          }
        }
      }
    }
  }

  // Calculate -n^ah n^bh omega^0h_{ah,bh}
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, omega);
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {
            na0_(l,m,k,j,i) = 0.0;
            for (int n = 0; n < 4; ++n) {
              for (int p = 0; p < 4; ++p) {
                na0_(l,m,k,j,i) += -nh_cc(n,l,m) * nh_cc(p,l,m) * omega(0,n,p);
              }
            }
          }
        }
      }
    }
  }

  // Calculate n^zeta
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, omega);
        for (int l = zs; l <= ze+1; ++l) {
          for (int m = ps; m <= pe; ++m) {
            na1_(l,m,k,j,i) = 0.0;
            for (int n = 0; n < 4; ++n) {
              for (int p = 0; p < 4; ++p) {
                na1_(l,m,k,j,i) += 1.0 / std::sin(zetaf(l)) * nh_fc(n,l,m)
                    * nh_fc(p,l,m) * (nh_fc(0,l,m) * omega(3,n,p) - nh_fc(3,l,m)
                    * omega(0,n,p));
              }
            }
          }
        }
      }
    }
  }

  // Calculate n^psi
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real x3 = pmb->pcoord->x3v(k);
        pmb->pcoord->Tetrad(x1, x2, x3, e, omega);
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe+1; ++m) {
            na2_(l,m,k,j,i) = 0.0;
            for (int n = 0; n < 4; ++n) {
              for (int p = 0; p < 4; ++p) {
                na2_(l,m,k,j,i) += 1.0 / SQR(std::sin(zetav(l))) * nh_cf(n,l,m)
                    * nh_cf(p,l,m) * (nh_cf(2,l,m) * omega(1,n,p) - nh_cf(1,l,m)
                    * omega(2,n,p));
              }
            }
          }
        }
      }
    }
  }

  // Deallocate memory for tetrad and rotation coefficients
  e.DeleteAthenaArray();
  omega.DeleteAthenaArray();

  // Deallocate unit normal components in orthonormal frame
  nh_cc.DeleteAthenaArray();
  nh_fc.DeleteAthenaArray();
  nh_cf.DeleteAthenaArray();

  // Allocate memory for left and right reconstructed states
  prim_l_.NewAthenaArray(nang_zpf, pmb->ncells1 + 1);
  prim_r_.NewAthenaArray(nang_zpf, pmb->ncells1 + 1);

  // Allocate memory for flux divergence calculation
  area_l_.NewAthenaArray(pmb->ncells1 + 1);
  area_r_.NewAthenaArray(pmb->ncells1 + 1);
  vol_.NewAthenaArray(pmb->ncells1 + 1);
  flux_div_.NewAthenaArray(nang, pmb->ncells1 + 1);
}

//----------------------------------------------------------------------------------------
// Radiation destructor

Radiation::~Radiation() {
  zetaf.DeleteAthenaArray();
  zetav.DeleteAthenaArray();
  dzetaf.DeleteAthenaArray();
  psif.DeleteAthenaArray();
  psiv.DeleteAthenaArray();
  dpsif.DeleteAthenaArray();
  zeta_length.DeleteAthenaArray();
  psi_length.DeleteAthenaArray();
  solid_angle.DeleteAthenaArray();
  flux_a[ZETADIR].DeleteAthenaArray();
  flux_a[PSIDIR].DeleteAthenaArray();
  n0_.DeleteAthenaArray();
  n1_.DeleteAthenaArray();
  n2_.DeleteAthenaArray();
  n3_.DeleteAthenaArray();
  na0_.DeleteAthenaArray();
  na1_.DeleteAthenaArray();
  na2_.DeleteAthenaArray();
  prim_l_.DeleteAthenaArray();
  prim_r_.DeleteAthenaArray();
  area_l_.DeleteAthenaArray();
  area_r_.DeleteAthenaArray();
  vol_.DeleteAthenaArray();
  flux_div_.DeleteAthenaArray();
}

//----------------------------------------------------------------------------------------
// Function for averaging intensities according to integrator weights
// Inputs:
//   cons_out, cons_in_1, cons_in_2: conserved intensity arrays, possibly uninitialized
//   weights: integrator weights
// Outputs:
//   cons_out: weighted intensity
// Notes:
//   same procedure as in Hydro::WeightedAveU()

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
            Real n1 = n1_(l,m,k,j,i);
            if (n1 > 0.0) {
              flux_x[X1DIR](lm,k,j,i) = n1 * prim_l_(lm,i);
            } else {
              flux_x[X1DIR](lm,k,j,i) = n1 * prim_r_(lm,i);
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
              Real n2 = n2_(l,m,k,j,i);
              if (n2 > 0.0) {
                flux_x[X2DIR](lm,k,j,i) = n2 * prim_l_(lm,i);
              } else {
                flux_x[X2DIR](lm,k,j,i) = n2 * prim_r_(lm,i);
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
              Real n3 = n3_(l,m,k,j,i);
              if (n3 > 0.0) {
                flux_x[X3DIR](lm,k,j,i) = n3 * prim_l_(lm,i);
              } else {
                flux_x[X3DIR](lm,k,j,i) = n3 * prim_r_(lm,i);
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
            Real na1 = na1_(l,m,k,j,i);
            if (na1 > 0.0) {
              flux_a[ZETADIR](lm_c,k,j,i) = na1 * prim_l_(lm_c,i);
            } else {
              flux_a[ZETADIR](lm_c,k,j,i) = na1 * prim_r_(lm_c,i);
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
            Real na2 = na2_(l,m,k,j,i);
            if (na2 > 0.0) {
              flux_a[PSIDIR](lm_c,k,j,i) = na2 * prim_l_(lm_c,i);
            } else {
              flux_a[PSIDIR](lm_c,k,j,i) = na2 * prim_r_(lm_c,i);
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
//   weight: fraction of full timestep
// Outputs:
//   cons_out: conserved values updated

void Radiation::AddFluxDivergenceToAverage(AthenaArray<Real> &prim_in, const Real weight,
    AthenaArray<Real> &cons_out) {

  // Extract Coordinates and timestep
  Coordinates *pcoord = pmy_block->pcoord;
  Real dt = pmy_block->pmy_mesh->dt;

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
          cons_out(lm,k,j,i) -= weight * dt * flux_div_(lm,i) / vol_(i);
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
            cons_out(lm,k,j,i) -= weight * dt * flux_div / omega;
          }
        }
      }
    }
  }

  // Add coordinate source term
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i) {
            cons_out(lm,k,j,i) += weight * dt * na0_(l,m,k,j,i) * prim_in(lm,k,j,i);
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
            cons_out(lm,k,j,i) = n0_(l,m,k,j,i) * prim_in(lm,k,j,i);
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
//   pcoord: pointer to Coordinates
//   il,iu,jl,ju,kl,ku: index bounds of region to be updated
// Outputs:
//   prim_out: primitives
// Notes:
//   This should be the only place where angular ghost zones need to be set.

void Radiation::ConservedToPrimitive(AthenaArray<Real> &cons_in,
    AthenaArray<Real> &prim_out, Coordinates *pcoord, int il, int iu, int jl, int ju,
    int kl, int ku) {

  // Calculate primitive intensities
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int k = kl; k <= ku; ++k) {
        for (int j = jl; j <= ju; ++j) {
          for (int i = il; i <= iu; ++i) {
            prim_out(lm,k,j,i) = cons_in(lm,k,j,i) / n0_(l,m,k,j,i);
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
// Function for adding all source terms beyond those induced by coordinates
// Inputs:
//   time: time of simulation
//   dt: simulation timestep
//   prim_in: primitive intensity
//   cons_out: conserved intensity
// Outputs:
//   cons_out: conserved intensity set
void Radiation::AddSourceTerms(const Real time, const Real dt,
    const AthenaArray<Real> &prim_in, AthenaArray<Real> &cons_out) {
  if (UserSourceTerm != NULL) {
    UserSourceTerm(pmy_block, time, dt, prim_in, cons_out);
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
//   spherical: flag indicating coordinates are spherical
// Outputs:
//   prim_vals: primitive values (I) set
//   cons_vals: conserved values (n^0 I) set
// Notes:
//   arrays should be 4D, with first index holding both zeta and psi
//   spherical coordinates:
//     theta (x2) will be mapped to [0, pi], adjusting phi (x3) as necessary
//     phi (x3) will be mapped to [-pi, pi]:
//       beams near phi = 0 should work fine
//       beams near phi = pi will likely not be initialized correctly

void Radiation::CalculateBeamSource(Real pos_1, Real pos_2, Real pos_3, Real width,
    Real dir_1, Real dir_2, Real dir_3, Real spread, AthenaArray<Real> &prim_vals,
    AthenaArray<Real> &cons_vals, bool spherical) {

  // Account for spherical coordinates in beam origin
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
  AthenaArray<Real> g, gi, e, omega, nh;
  g.NewAthenaArray(NMETRIC, ie + 1);
  gi.NewAthenaArray(NMETRIC, ie + 1);
  e.NewAthenaArray(4, 4);
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

        // Extract position, accounting for spherical coordinates
        Real x1 = pmy_block->pcoord->x1v(i);
        Real x2 = pmy_block->pcoord->x2v(j);
        Real x3 = pmy_block->pcoord->x3v(k);
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
              prim_vals(lm,k,j,i) = 0.0;
              cons_vals(lm,k,j,i) = 0.0;
            }
          }
          continue;
        }

        // Normalize direction
        Real dir_mag_sq = g(I11,i) * SQR(dir_1) + 2.0 * g(I12,i) * dir_1 * dir_2
            + 2.0 * g(I13,i) * dir_1 * dir_3 + g(I22,i) * SQR(dir_2)
            + 2.0 * g(I23,i) * dir_2 * dir_3 + g(I33,i) * SQR(dir_3);
        Real dir_norm = 1.0 / std::sqrt(dir_mag_sq);
        Real v1 = dir_1 * dir_norm;
        Real v2 = dir_2 * dir_norm;
        Real v3 = dir_3 * dir_norm;

        // Calculate tetrad
        pmy_block->pcoord->Tetrad(x1, x2, x3, e, omega);

        // Go through angles
        for (int l = zs; l <= ze; ++l) {
          for (int m = ps; m <= pe; ++m) {

            // Calculate combined angle index
            int lm = AngleInd(l, m);

            // Calculate unit normal components in coordinate frame
            Real n1 = 0.0;
            Real n2 = 0.0;
            Real n3 = 0.0;
            for (int n = 0; n < 4; ++n) {
              n1 += e(n,1) * nh(l,m,n);
              n2 += e(n,2) * nh(l,m,n);
              n3 += e(n,3) * nh(l,m,n);
            }

            // Calculate angle to beam direction
            Real mu = g(I11,i) * v1 * n1 + g(I12,i) * (v1 * n2 + v2 * n1)
                + g(I13,i) * (v1 * n3 + v3 * n1) + g(I22,i) * v2 * n2
                + g(I23,i) * (v2 * n3 + v3 * n2) + g(I33,i) * v3 * n3;

            // Set to 0 if too far from beam in angle
            if (mu <= mu_min) {
              prim_vals(lm,k,j,i) = 0.0;
              cons_vals(lm,k,j,i) = 0.0;
              continue;
            }

            // Set to nonzero value
            prim_vals(lm,k,j,i) = 1.0;
            cons_vals(lm,k,j,i) = n0_(l,m,k,j,i);
          }
        }
      }
    }
  }

  // Deallocate scratch arrays
  g.DeleteAthenaArray();
  gi.DeleteAthenaArray();
  e.DeleteAthenaArray();
  omega.DeleteAthenaArray();
  nh.DeleteAthenaArray();
  return;
}

//----------------------------------------------------------------------------------------
// Function for calculating moments of radiation field
// Inputs: (none)
// Outputs:
//   moments: 4D array:
//     index 0:
//       0: energy density (\int I n^0 n^0 d\Omega)
//     index 1: k
//     index 2: j
//     index 3: i
// Notes:
//   intensity defined by prim values when function is called

void Radiation::SetMoments(AthenaArray<Real> &moments) {
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        moments(0,k,j,i) = 0.0;
      }
    }
  }
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd(l, m);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i) {
            moments(0,k,j,i) += prim(lm,k,j,i) * SQR(n0_(l,m,k,j,i)) * solid_angle(l,m);
          }
        }
      }
    }
  }
  return;
}
