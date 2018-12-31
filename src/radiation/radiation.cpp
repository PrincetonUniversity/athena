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

  // Set flags
  source_terms_defined = false;

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

  // Allocate memory for angles
  zetaf.NewAthenaArray(nzeta + 2*NGHOST + 1);
  zetav.NewAthenaArray(nzeta + 2*NGHOST);
  dzetaf.NewAthenaArray(nzeta + 2*NGHOST);
  psif.NewAthenaArray(npsi + 2*NGHOST + 1);
  psiv.NewAthenaArray(npsi + 2*NGHOST);
  dpsif.NewAthenaArray(npsi + 2*NGHOST);

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

  // Allocate memory for intensities
  int num_cells_1 = ie + NGHOST;
  int num_cells_2 = 1;
  if (js != je) {
    num_cells_2 = je + NGHOST;
  }
  int num_cells_3 = 1;
  if (ks != ke) {
    num_cells_3 = ke + NGHOST;
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
}

//----------------------------------------------------------------------------------------
// Indexing function for intensities
// Inputs:
//   l: zeta-index
//   m: psi-index
//   zeta_face: flag indicating zeta-index is on faces
//   psi_face: flag indicating psi-index is on faces
// Outputs:
//   returned value: 1D index for both zeta and psi

int Radiation::IntInd(int l, int m, bool zeta_face, bool psi_face) {
  if (psi_face) {
    return l * (npsi + 2*NGHOST + 1) + m;
  }
  return l * (npsi + 2*NGHOST) + m;
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

  // Prepare index bounds
  int is = pmy_block->is;
  int ie = pmy_block->ie;
  int js = pmy_block->js;
  int je = pmy_block->je;
  int ks = pmy_block->ks;
  int ke = pmy_block->ke;

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

void Radiation::ConservedToPrimitive(AthenaArray<Real> &cons_in,
    AthenaArray<Real> &prim_out, Coordinates *pcoord, int il, int iu, int jl, int ju,
    int kl, int ku) {
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
  return;
}
