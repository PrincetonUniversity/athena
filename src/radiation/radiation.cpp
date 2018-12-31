//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file radiation.cpp
//  \brief implementation of functions in class Radiation

// C++ headers
#include <cmath>      // acos(), cos(), sin()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str(), string

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

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin) {

  // Set object pointers
  pmy_block = pmb;

  // Set user source term
  UserSourceTerm = pmb->pmy_mesh->UserRadSourceTerm_;

  // Set flags
  if (UserSourceTerm == NULL) {
    source_terms_defined = false;
  } else {
    source_terms_defined = true;
  }

  // Set parameters
  nzeta = pin->GetInteger("radiation", "n_polar");
  npsi = pin->GetInteger("radiation", "n_azimuthal");
  nang = (nzeta + 2*NGHOST) * (npsi + 2*NGHOST);
  nang_zf = (nzeta + 2*NGHOST + 1) * (npsi + 2*NGHOST);
  nang_pf = (nzeta + 2*NGHOST) * (npsi + 2*NGHOST + 1);
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

  // Allocate memory for intensities
  int num_cells_1 = ie + NGHOST + 1;
  int num_cells_2 = 1;
  if (js != je) {
    num_cells_2 = je + NGHOST + 1;
  }
  int num_cells_3 = 1;
  if (ks != ke) {
    num_cells_3 = ke + NGHOST + 1;
  }
  prim.NewAthenaArray(nang, num_cells_3, num_cells_2, num_cells_1);
  prim1.NewAthenaArray(nang, num_cells_3, num_cells_2, num_cells_1);
  cons.NewAthenaArray(nang, num_cells_3, num_cells_2, num_cells_1);
  cons1.NewAthenaArray(nang, num_cells_3, num_cells_2, num_cells_1);
  cons2.NewAthenaArray(nang, num_cells_3, num_cells_2, num_cells_1);

  // Allocate memory for fluxes
  flux_x[X1DIR].NewAthenaArray(nang, num_cells_3, num_cells_2, num_cells_1 + 1);
  if (js != je) {
    flux_x[X2DIR].NewAthenaArray(nang, num_cells_3, num_cells_2 + 1, num_cells_1);
  }
  if (ks != ke) {
    flux_x[X3DIR].NewAthenaArray(nang, num_cells_3 + 1, num_cells_2, num_cells_1);
  }
  flux_a[ZETADIR].NewAthenaArray(nang_zf, num_cells_3, num_cells_2, num_cells_1);
  flux_a[PSIDIR].NewAthenaArray(nang_pf, num_cells_3, num_cells_2, num_cells_1);

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
  n0_.NewAthenaArray(num_cells_zeta, num_cells_psi, num_cells_3, num_cells_2,
      num_cells_1);
  n1_.NewAthenaArray(num_cells_zeta, num_cells_psi, num_cells_3, num_cells_2,
      num_cells_1 + 1);
  n2_.NewAthenaArray(num_cells_zeta, num_cells_psi, num_cells_3, num_cells_2 + 1,
      num_cells_1);
  n3_.NewAthenaArray(num_cells_zeta, num_cells_psi, num_cells_3 + 1, num_cells_2,
      num_cells_1);
  na0_.NewAthenaArray(num_cells_zeta, num_cells_psi, num_cells_3, num_cells_2,
      num_cells_1);
  na1_.NewAthenaArray(num_cells_zeta + 1, num_cells_psi, num_cells_3, num_cells_2,
      num_cells_1);
  na2_.NewAthenaArray(num_cells_zeta, num_cells_psi + 1, num_cells_3, num_cells_2,
      num_cells_1);

  // Allocate memory for tetrad and rotation coefficients
  AthenaArray<Real> e, omega;
  e.NewAthenaArray(4, 4);
  omega.NewAthenaArray(4, 4, 4);

  // Calculate n^0
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
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
                na2_(l,m,k,j,i) += 1.0 / SQR(std::sin(zetaf(l))) * nh_cf(n,l,m)
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
  prim_l_.NewAthenaArray(nang, num_cells_1 + 1);
  prim_r_.NewAthenaArray(nang, num_cells_1 + 1);

  // Allocate memory for flux divergence calculation
  area_l_.NewAthenaArray(num_cells_1 + 1);
  area_r_.NewAthenaArray(num_cells_1 + 1);
  vol_.NewAthenaArray(num_cells_1 + 1);
  flux_div_.NewAthenaArray(nang, num_cells_1 + 1);
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
  prim.DeleteAthenaArray();
  prim1.DeleteAthenaArray();
  cons.DeleteAthenaArray();
  cons1.DeleteAthenaArray();
  cons2.DeleteAthenaArray();
  flux_x[X1DIR].DeleteAthenaArray();
  flux_x[X2DIR].DeleteAthenaArray();
  flux_x[X3DIR].DeleteAthenaArray();
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

void Radiation::WeightedAveCons(AthenaArray<Real> &cons_out, AthenaArray<Real> &cons_in_1,
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
  if (order != 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in Radiation\n";
    msg << "only first order supported\n";
    throw std::runtime_error(msg.str().c_str());
  }

  // Calculate x1-fluxes
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe; ++m) {
      int lm = AngleInd_(l, m);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {

          // Reconstruction
          if (order == 1) {
            for (int i = is; i <= ie+1; ++i) {
              prim_l_(lm,i) = prim_in(lm,k,j,i-1);
              prim_r_(lm,i) = prim_in(lm,k,j,i);
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
        int lm = AngleInd_(l, m);
        for (int k = ks; k <= ke; ++k) {
          for (int j = js; j <= je+1; ++j) {

            // Reconstruction
            if (order == 1) {
              for (int i = is; i <= ie; ++i) {
                prim_l_(lm,i) = prim_in(lm,k,j-1,i);
                prim_r_(lm,i) = prim_in(lm,k,j,i);
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
        int lm = AngleInd_(l, m);
        for (int k = ks; k <= ke+1; ++k) {
          for (int j = js; j <= je; ++j) {

            // Reconstruction
            if (order == 1) {
              for (int i = is; i <= ie; ++i) {
                prim_l_(lm,i) = prim_in(lm,k-1,j,i);
                prim_r_(lm,i) = prim_in(lm,k,j,i);
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
      int lm = AngleInd_(l, m, true, false);
      int lm_l = AngleInd_(l - 1, m, true, false);
      int lm_r = AngleInd_(l, m, true, false);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {

          // Reconstruction
          if (order == 1) {
            for (int i = is; i <= ie; ++i) {
              prim_l_(lm,i) = prim_in(lm_l,k,j,i);
              prim_r_(lm,i) = prim_in(lm_r,k,j,i);
            }
          }

          // Upwind flux calculation
          for (int i = is; i <= ie; ++i) {
            Real na1 = na1_(l,m,k,j,i);
            if (na1 > 0.0) {
              flux_a[ZETADIR](lm,k,j,i) = na1 * prim_l_(lm,i);
            } else {
              flux_a[ZETADIR](lm,k,j,i) = na1 * prim_r_(lm,i);
            }
          }
        }
      }
    }
  }

  // Calculate psi-fluxes
  for (int l = zs; l <= ze; ++l) {
    for (int m = ps; m <= pe+1; ++m) {
      int lm = AngleInd_(l, m, false, true);
      int lm_l = AngleInd_(l, m - 1, false, true);
      int lm_r = AngleInd_(l, m, false, true);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {

          // Reconstruction
          if (order == 1) {
            for (int i = is; i <= ie; ++i) {
              prim_l_(lm,i) = prim_in(lm_l,k,j,i);
              prim_r_(lm,i) = prim_in(lm_r,k,j,i);
            }
          }

          // Upwind flux calculation
          for (int i = is; i <= ie; ++i) {
            Real na2 = na2_(l,m,k,j,i);
            if (na2 > 0.0) {
              flux_a[PSIDIR](lm,k,j,i) = na2 * prim_l_(lm,i);
            } else {
              flux_a[PSIDIR](lm,k,j,i) = na2 * prim_r_(lm,i);
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
  int lms = AngleInd_(zs, ps);
  int lme = AngleInd_(ze, pe);

  // Go through all cells
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {

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
            flux_div_(lm,i) += area_r_(i) * flux_x[X2DIR](lm,k,j+1,i)
                - area_l_(i) * flux_x[X2DIR](lm,k,j,i);
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

      // Calculate angle lengths and solid angles
      int lm = AngleInd_(l, m);
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
            Real flux_div = psi_length_p * flux_a[ZETADIR](l+1,k,j,i)
                - psi_length_m * flux_a[ZETADIR](l,k,j,i);

            // Add psi-divergence
            flux_div += zeta_length_p * flux_a[PSIDIR](m+1,k,j,i)
                - zeta_length_m * flux_a[PSIDIR](m,k,j,i);

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
      int lm = AngleInd_(l, m);
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
      int lm = AngleInd_(l, m);
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
      int lm = AngleInd_(l, m);
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
      int lm = AngleInd_(l, m);
      int lm_src = AngleInd_(l, m_src);
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
      int lm = AngleInd_(l, m);
      int lm_src = AngleInd_(l, m_src);
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
      int lm = AngleInd_(l, m);
      int lm_src = AngleInd_(l_src, m_src);
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
      int lm = AngleInd_(l, m);
      int lm_src = AngleInd_(l_src, m_src);
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
      int lm = AngleInd_(l, m);
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

//----------------------------------------------------------------------------------------
// Indexing function for angles
// Inputs:
//   l: zeta-index
//   m: psi-index
//   zeta_face: flag indicating zeta-index is on faces
//   psi_face: flag indicating psi-index is on faces
// Outputs:
//   returned value: 1D index for both zeta and psi

int Radiation::AngleInd_(int l, int m, bool zeta_face, bool psi_face) {
  if (psi_face) {
    return l * (npsi + 2*NGHOST + 1) + m;
  }
  return l * (npsi + 2*NGHOST) + m;
}
