//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file field_diffusion_fourth.cpp
//! \brief fourth-order accurate Ohmic resistive EMF and Poynting flux
//!
//! Companions to the second-order operators in diffusivity.cpp for the fourth-order
//! finite-volume MHD scheme (time/xorder=4, uniform Cartesian square cells).
//! Strategy (all steps are O(h^4) accurate):
//!  1. the current density J = curl(B) is evaluated at edge-centered points with
//!     fourth-order finite differences of the point-valued face fields b_fc (computed
//!     each stage by Field::FaceAveragedToCellAveragedField);
//!  2. the point-valued resistive EMF E = eta_O J is stored on edges (e_pt_);
//!  3. the edge-AVERAGED EMF required by the fourth-order CT update is obtained with
//!     the one-dimensional Laplacian correction along the edge direction,
//!     <E>_edge = E + h^2/24 d^2E/dl^2;
//!  4. the face-averaged resistive Poynting flux S = E x B is assembled from
//!     fourth-order interpolations of e_pt_ and bcc_center to face-centered points,
//!     followed by the transverse Laplacian face-averaging correction.
//! The diffusivity eta_O is interpolated to edges with the same second-order averages
//! as the second-order operators (exact for constant coefficients, which is the only
//! supported use of Ohmic resistivity at fourth order at present).

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../field.hpp"
#include "field_diffusion.hpp"

namespace {
constexpr Real ONE_24TH = 1.0/24.0;
} // namespace

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::OhmicEMFFourth
//! \brief fourth-order edge-averaged Ohmic EMF, E = eta_O J, from point-valued b_fc

void FieldDiffusion::OhmicEMFFourth(EdgeField &e) {
  MeshBlock *pmb = pmy_block;
  Field *pf = pmb->pfield;
  const FaceField &bf = pf->b_fc;
  const bool f2 = pmb->pmy_mesh->f2;
  const bool f3 = pmb->pmy_mesh->f3;
  const int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je,
            ks = pmb->ks, ke = pmb->ke;
  Coordinates *pco = pmb->pcoord;
  const Real h1 = pco->dx1f(is), h2 = pco->dx2f(js), h3 = pco->dx3f(ks);
  AthenaArray<Real> &e1 = e.x1e, &e2 = e.x2e, &e3 = e.x3e;
  AthenaArray<Real> &E1 = e_pt_.x1e, &E2 = e_pt_.x2e, &E3 = e_pt_.x3e;
  const int oh = DiffProcess::ohmic;

  // 1D update: J2 = -dB3/dx1, J3 = dB2/dx1 at x1-faces; edge averages are trivial
  if (!f2) {
#pragma omp simd
    for (int i=is; i<=ie+1; ++i) {
      Real eta_O = 0.5*(etaB(oh,ks,js,i-1) + etaB(oh,ks,js,i));
      Real dB3dx = (27.0*(bf.x3f(ks,js,i) - bf.x3f(ks,js,i-1))
                    - (bf.x3f(ks,js,i+1) - bf.x3f(ks,js,i-2)))/(24.0*h1);
      Real dB2dx = (27.0*(bf.x2f(ks,js,i) - bf.x2f(ks,js,i-1))
                    - (bf.x2f(ks,js,i+1) - bf.x2f(ks,js,i-2)))/(24.0*h1);
      E2(ks,js,i) = -eta_O*dB3dx;
      E3(ks,js,i) = eta_O*dB2dx;
      e2(ks  ,js  ,i) += E2(ks,js,i);
      e2(ke+1,js  ,i)  = e2(ks,js,i);
      e3(ks  ,js  ,i) += E3(ks,js,i);
      e3(ks  ,je+1,i)  = e3(ks,js,i);
    }
    return;
  }

  // point-valued EMF at edge centers, over ranges wide enough for the edge-averaging
  // corrections below and the face interpolations in PoyntingFluxFourth()
  const int kl1 = f3 ? ks-2 : ks, ku1 = f3 ? ke+3 : ks;  // E1 (x1v, x2f, x3f)
  for (int k=kl1; k<=ku1; ++k) {
    for (int j=js-2; j<=je+3; ++j) {
#pragma omp simd
      for (int i=is-1; i<=ie+2; ++i) {
        Real J1 = (27.0*(bf.x3f(k,j,i) - bf.x3f(k,j-1,i))
                   - (bf.x3f(k,j+1,i) - bf.x3f(k,j-2,i)))/(24.0*h2);
        Real eta_O;
        if (f3) {
          J1 -= (27.0*(bf.x2f(k,j,i) - bf.x2f(k-1,j,i))
                 - (bf.x2f(k+1,j,i) - bf.x2f(k-2,j,i)))/(24.0*h3);
          eta_O = 0.25*(etaB(oh,k  ,j,i) + etaB(oh,k  ,j-1,i)
                        + etaB(oh,k-1,j,i) + etaB(oh,k-1,j-1,i));
        } else {
          eta_O = 0.5*(etaB(oh,k,j,i) + etaB(oh,k,j-1,i));
        }
        E1(k,j,i) = eta_O*J1;
      }
    }
  }

  const int kl2 = f3 ? ks-2 : ks, ku2 = f3 ? ke+3 : ks;  // E2 (x1f, x2v, x3f)
  for (int k=kl2; k<=ku2; ++k) {
    for (int j=js-1; j<=je+2; ++j) {
#pragma omp simd
      for (int i=is-2; i<=ie+3; ++i) {
        Real J2 = -(27.0*(bf.x3f(k,j,i) - bf.x3f(k,j,i-1))
                    - (bf.x3f(k,j,i+1) - bf.x3f(k,j,i-2)))/(24.0*h1);
        Real eta_O;
        if (f3) {
          J2 += (27.0*(bf.x1f(k,j,i) - bf.x1f(k-1,j,i))
                 - (bf.x1f(k+1,j,i) - bf.x1f(k-2,j,i)))/(24.0*h3);
          eta_O = 0.25*(etaB(oh,k  ,j,i) + etaB(oh,k  ,j,i-1)
                        + etaB(oh,k-1,j,i) + etaB(oh,k-1,j,i-1));
        } else {
          eta_O = 0.5*(etaB(oh,k,j,i) + etaB(oh,k,j,i-1));
        }
        E2(k,j,i) = eta_O*J2;
      }
    }
  }

  const int kl3 = f3 ? ks-1 : ks, ku3 = f3 ? ke+2 : ks;  // E3 (x1f, x2f, x3v)
  for (int k=kl3; k<=ku3; ++k) {
    for (int j=js-2; j<=je+3; ++j) {
#pragma omp simd
      for (int i=is-2; i<=ie+3; ++i) {
        Real J3 = (27.0*(bf.x2f(k,j,i) - bf.x2f(k,j,i-1))
                   - (bf.x2f(k,j,i+1) - bf.x2f(k,j,i-2)))/(24.0*h1)
                  - (27.0*(bf.x1f(k,j,i) - bf.x1f(k,j-1,i))
                     - (bf.x1f(k,j+1,i) - bf.x1f(k,j-2,i)))/(24.0*h2);
        Real eta_O = 0.25*(etaB(oh,k,j  ,i) + etaB(oh,k,j  ,i-1)
                           + etaB(oh,k,j-1,i) + etaB(oh,k,j-1,i-1));
        E3(k,j,i) = eta_O*J3;
      }
    }
  }

  // edge-averaged EMF: one-dimensional Laplacian correction along the edge direction
  if (!f3) { // 2D: e3 edges run along x3 (no variation); mirror 2nd-order row copies
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        e1(ks  ,j,i) += E1(ks,j,i)
                        + ONE_24TH*(E1(ks,j,i-1) - 2.0*E1(ks,j,i) + E1(ks,j,i+1));
        e1(ke+1,j,i)  = e1(ks,j,i);
        e2(ks  ,j,i) += E2(ks,j,i)
                        + ONE_24TH*(E2(ks,j-1,i) - 2.0*E2(ks,j,i) + E2(ks,j+1,i));
        e2(ke+1,j,i)  = e2(ks,j,i);
        e3(ks,  j,i) += E3(ks,j,i);
      }
    }
    return;
  }

  // 3D update:
  for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        e1(k,j,i) += E1(k,j,i)
                     + ONE_24TH*(E1(k,j,i-1) - 2.0*E1(k,j,i) + E1(k,j,i+1));
        e2(k,j,i) += E2(k,j,i)
                     + ONE_24TH*(E2(k,j-1,i) - 2.0*E2(k,j,i) + E2(k,j+1,i));
        e3(k,j,i) += E3(k,j,i)
                     + ONE_24TH*(E3(k-1,j,i) - 2.0*E3(k,j,i) + E3(k+1,j,i));
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::PoyntingFluxFourth
//! \brief fourth-order face-averaged resistive Poynting flux S = E x B from the
//! point-valued EMF (e_pt_, filled by OhmicEMFFourth) and bcc_center

void FieldDiffusion::PoyntingFluxFourth() {
  MeshBlock *pmb = pmy_block;
  Field *pf = pmb->pfield;
  const AthenaArray<Real> &bc = pf->bcc_center;
  const bool f2 = pmb->pmy_mesh->f2;
  const bool f3 = pmb->pmy_mesh->f3;
  const int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je,
            ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &E1 = e_pt_.x1e, &E2 = e_pt_.x2e, &E3 = e_pt_.x3e;
  AthenaArray<Real> &f1 = pflux.x1f, &f2fl = pflux.x2f, &f3fl = pflux.x3f;
  AthenaArray<Real> &g1 = pfpt_.x1f, &g2 = pfpt_.x2f, &g3 = pfpt_.x3f;

  pflux.x1f.ZeroClear();
  pflux.x2f.ZeroClear();
  pflux.x3f.ZeroClear();

  // 1D update: no transverse averaging
  if (!f2) {
#pragma omp simd
    for (int i=is; i<=ie+1; ++i) {
      Real B2f = (9.0*(bc(IB2,ks,js,i-1) + bc(IB2,ks,js,i))
                  - (bc(IB2,ks,js,i-2) + bc(IB2,ks,js,i+1)))/16.0;
      Real B3f = (9.0*(bc(IB3,ks,js,i-1) + bc(IB3,ks,js,i))
                  - (bc(IB3,ks,js,i-2) + bc(IB3,ks,js,i+1)))/16.0;
      f1(ks,js,i) = -B2f*E3(ks,js,i) + B3f*E2(ks,js,i);
    }
    return;
  }

  //--- x1-face Poynting flux, S1 = E2 B3 - E3 B2
  {
    const int pkl = f3 ? ks-1 : ks, pku = f3 ? ke+1 : ks;
    for (int k=pkl; k<=pku; ++k) {
      for (int j=js-1; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          Real B2f = (9.0*(bc(IB2,k,j,i-1) + bc(IB2,k,j,i))
                      - (bc(IB2,k,j,i-2) + bc(IB2,k,j,i+1)))/16.0;
          Real B3f = (9.0*(bc(IB3,k,j,i-1) + bc(IB3,k,j,i))
                      - (bc(IB3,k,j,i-2) + bc(IB3,k,j,i+1)))/16.0;
          // E3 at (x1f, x2v): interpolate corner values along x2
          Real E3f = (9.0*(E3(k,j,i) + E3(k,j+1,i))
                      - (E3(k,j-1,i) + E3(k,j+2,i)))/16.0;
          // E2 at (x1f, x3v): 2D value is already there; 3D interpolate along x3
          Real E2f = f3 ? (9.0*(E2(k,j,i) + E2(k+1,j,i))
                           - (E2(k-1,j,i) + E2(k+2,j,i)))/16.0
                        : E2(ks,j,i);
          g1(k,j,i) = -B2f*E3f + B3f*E2f;
        }
      }
    }
    for (int k=ks; k<=(f3 ? ke : ks); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          Real corr = (g1(k,j-1,i) - 2.0*g1(k,j,i) + g1(k,j+1,i));
          if (f3) corr += (g1(k-1,j,i) - 2.0*g1(k,j,i) + g1(k+1,j,i));
          f1(k,j,i) = g1(k,j,i) + ONE_24TH*corr;
        }
      }
    }
  }

  //--- x2-face Poynting flux, S2 = E3 B1 - E1 B3
  {
    const int pkl = f3 ? ks-1 : ks, pku = f3 ? ke+1 : ks;
    for (int k=pkl; k<=pku; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is-1; i<=ie+1; ++i) {
          Real B1f = (9.0*(bc(IB1,k,j-1,i) + bc(IB1,k,j,i))
                      - (bc(IB1,k,j-2,i) + bc(IB1,k,j+1,i)))/16.0;
          Real B3f = (9.0*(bc(IB3,k,j-1,i) + bc(IB3,k,j,i))
                      - (bc(IB3,k,j-2,i) + bc(IB3,k,j+1,i)))/16.0;
          // E3 at (x1v, x2f): interpolate corner values along x1
          Real E3f = (9.0*(E3(k,j,i) + E3(k,j,i+1))
                      - (E3(k,j,i-1) + E3(k,j,i+2)))/16.0;
          // E1 at (x2f, x3v): 2D value is already there; 3D interpolate along x3
          Real E1f = f3 ? (9.0*(E1(k,j,i) + E1(k+1,j,i))
                           - (E1(k-1,j,i) + E1(k+2,j,i)))/16.0
                        : E1(ks,j,i);
          g2(k,j,i) = B1f*E3f - B3f*E1f;
        }
      }
    }
    for (int k=ks; k<=(f3 ? ke : ks); ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          Real corr = (g2(k,j,i-1) - 2.0*g2(k,j,i) + g2(k,j,i+1));
          if (f3) corr += (g2(k-1,j,i) - 2.0*g2(k,j,i) + g2(k+1,j,i));
          f2fl(k,j,i) = g2(k,j,i) + ONE_24TH*corr;
        }
      }
    }
  }

  //--- x3-face Poynting flux, S3 = E1 B2 - E2 B1 (3D only)
  if (f3) {
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js-1; j<=je+1; ++j) {
#pragma omp simd
        for (int i=is-1; i<=ie+1; ++i) {
          Real B1f = (9.0*(bc(IB1,k-1,j,i) + bc(IB1,k,j,i))
                      - (bc(IB1,k-2,j,i) + bc(IB1,k+1,j,i)))/16.0;
          Real B2f = (9.0*(bc(IB2,k-1,j,i) + bc(IB2,k,j,i))
                      - (bc(IB2,k-2,j,i) + bc(IB2,k+1,j,i)))/16.0;
          // E1 at (x2v, x3f): interpolate edge values along x2
          Real E1f = (9.0*(E1(k,j,i) + E1(k,j+1,i))
                      - (E1(k,j-1,i) + E1(k,j+2,i)))/16.0;
          // E2 at (x1v, x3f): interpolate edge values along x1
          Real E2f = (9.0*(E2(k,j,i) + E2(k,j,i+1))
                      - (E2(k,j,i-1) + E2(k,j,i+2)))/16.0;
          g3(k,j,i) = B2f*E1f - B1f*E2f;
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          Real corr = (g3(k,j,i-1) - 2.0*g3(k,j,i) + g3(k,j,i+1))
                      + (g3(k,j-1,i) - 2.0*g3(k,j,i) + g3(k,j+1,i));
          f3fl(k,j,i) = g3(k,j,i) + ONE_24TH*corr;
        }
      }
    }
  }
  return;
}
