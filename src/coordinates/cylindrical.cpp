//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cylindrical.cpp
//! \brief implements functions for cylindrical (r-phi-z) coordinates in a derived class
//! of the Coordinates abstract base class.

// C headers

// C++ headers
#include <cmath>   // pow()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../cr/cr.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "coordinates.hpp"

//----------------------------------------------------------------------------------------
//! Cylindrical coordinates initialization

void Coordinates::Initialize(ParameterInput *pin) {
  // initialize volume-averaged coordinates and spacing
  // x1-direction: x1v = (\int r dV / \int dV) = d(r^3/3)d(r^2/2)
  for (int i=il-ng; i<=iu+ng; ++i) {
    x1v(i) = (TWO_3RD)*(std::pow(x1f(i+1), 3) - std::pow(x1f(i), 3)) /
             (std::pow(x1f(i+1), 2) - std::pow(x1f(i), 2));
    // reduces to eq for centroid: R_i + dR_i^2/(12*R_i)
    // see Mignone (2014) eq 17, e.g.
  }
  for (int i=il-ng; i<=iu+ng-1; ++i) {
    dx1v(i) = x1v(i+1) - x1v(i);
  }

  // x2-direction: x2v = (\int phi dV / \int dV) = dphi/2
  if (pmy_block->block_size.nx2 == 1) {
    x2v(jl) = 0.5*(x2f(jl+1) + x2f(jl));
    dx2v(jl) = dx2f(jl);
  } else {
    for (int j=jl-ng; j<=ju+ng; ++j) {
      x2v(j) = 0.5*(x2f(j+1) + x2f(j));
    }
    for (int j=jl-ng; j<=ju+ng-1; ++j) {
      dx2v(j) = x2v(j+1) - x2v(j);
    }
  }

  // x3-direction: x3v = (\int z dV / \int dV) = dz/2
  if (pmy_block->block_size.nx3 == 1) {
    x3v(kl) = 0.5*(x3f(kl+1) + x3f(kl));
    dx3v(kl) = dx3f(kl);
  } else {
    for (int k=kl-ng; k<=ku+ng; ++k) {
      x3v(k) = 0.5*(x3f(k+1) + x3f(k));
    }
    for (int k=kl-ng; k<=ku+ng-1; ++k) {
      dx3v(k) = x3v(k+1) - x3v(k);
    }
  }

  // initialize geometry coefficients
  // x1-direction
  for (int i=il-ng; i<=iu+ng; ++i) {
    h2v(i) = x1v(i);
    h2f(i) = x1f(i);
    h31v(i) = 1.0;
    h31f(i) = 1.0;
    dh2vd1(i) = 1.0;
    dh2fd1(i) = 1.0;
    dh31vd1(i) = 0.0;
    dh31fd1(i) = 0.0;
  }

  // x2-direction
  if (pmy_block->block_size.nx2 == 1) {
    h32v(jl) = 1.0;
    h32f(jl) = 1.0;
    dh32vd2(jl) = 0.0;
    dh32fd2(jl) = 0.0;
  } else {
    for (int j=jl-ng; j<=ju+ng; ++j) {
      h32v(j) = 1.0;
      h32f(j) = 1.0;
      dh32vd2(j) = 0.0;
      dh32fd2(j) = 0.0;
    }
  }

  // initialize area-averaged coordinates used with MHD AMR
  if ((pmy_block->pmy_mesh->multilevel) && MAGNETIC_FIELDS_ENABLED) {
    for (int i=il-ng; i<=iu+ng; ++i) {
      x1s2(i) = 0.5*(x1f(i+1) + x1f(i));
      x1s3(i) = (2.0/3.0)*(std::pow(x1f(i+1),3) - std::pow(x1f(i),3))
                /(SQR(x1f(i+1)) - SQR(x1f(i)));
    }
    if (pmy_block->block_size.nx2 == 1) {
      x2s1(jl) = x2s3(jl) = x2v(jl);
    } else {
      for (int j=jl-ng; j<=ju+ng; ++j) {
        x2s1(j) = x2s3(j) = 0.5*(x2f(j+1) + x2f(j));
      }
    }
    if (pmy_block->block_size.nx3 == 1) {
      x3s1(kl) = x3s2(kl) = x3v(kl);
    } else {
      for (int k=kl-ng; k<=ku+ng; ++k) {
        x3s1(k) = x3s2(k) = 0.5*(x3f(k+1) + x3f(k));
      }
    }
  }

  // Allocate memory for internal scratch arrays to store partial calculations
  if (coarse_flag) {
    coord_area3_i_.NewAthenaArray(nc1);
    // Compute and store constant coefficients needed for face-areas, cell-volumes, etc.
    // This helps improve performance.
#pragma omp simd
    for (int i=il-ng; i<=iu+ng; ++i) {
      Real rm = x1f(i  );
      Real rp = x1f(i+1);
      // dV = 0.5*(R_{i+1}^2 - R_{i}^2)
      coord_area3_i_(i)= 0.5*(rp*rp - rm*rm);
    }
  } else {
    coord_area3_i_.NewAthenaArray(nc1);
    coord_area3vc_i_.NewAthenaArray(nc1);
    coord_vol_i_.NewAthenaArray(nc1);
    coord_src1_i_.NewAthenaArray(nc1);
    coord_src2_i_.NewAthenaArray(nc1);
    phy_src1_i_.NewAthenaArray(nc1);
    phy_src2_i_.NewAthenaArray(nc1);

    // Compute and store constant coefficients needed for face-areas, cell-volumes, etc.
    // This helps improve performance.
#pragma omp simd
    for (int i=il-ng; i<=iu+ng; ++i) {
      Real rm = x1f(i  );
      Real rp = x1f(i+1);
      // dV = 0.5*(R_{i+1}^2 - R_{i}^2)
      coord_area3_i_(i)= 0.5*(rp*rp - rm*rm);
      // dV = 0.5*(R_{i+1}^2 - R_{i}^2)
      coord_vol_i_(i) = coord_area3_i_(i);
      // (A1^{+} - A1^{-})/dV
      coord_src1_i_(i) = dx1f(i)/coord_vol_i_(i);
      // (dR/2)/(R_c dV)
      coord_src2_i_(i) = dx1f(i)/((rm + rp)*coord_vol_i_(i));
      // Rf_{i}/R_{i}/Rf_{i}^2
      phy_src1_i_(i) = 1.0/(x1v(i)*x1f(i));
    }
#pragma omp simd
    for (int i=il-ng; i<=iu+(ng-1); ++i) {
      // Rf_{i+1}/R_{i}/Rf_{i+1}^2
      phy_src2_i_(i) = 1.0/(x1v(i)*x1f(i+1));
      // dV = 0.5*(R_{i+1}^2 - R_{i}^2)
      coord_area3vc_i_(i)= 0.5*(SQR(x1v(i+1)) - SQR(x1v(i)));
    }
  }
}


//----------------------------------------------------------------------------------------
// EdgeXLength functions: compute physical length at cell edge-X as vector
// Only overwrite functions that differ from implementation in base class

// Edge2(i,j,k) located at (i-1/2,j,k-1/2), i.e. (x1f(i), x2v(j), x3f(k))

void Coordinates::Edge2Length(const int k, const int j, const int il, const int iu,
                              AthenaArray<Real> &len) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    len(i) = x1f(i)*dx2f(j);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetEdgeXLength functions: return length of edge-X length at (i,j,k)

Real Coordinates::GetEdge2Length(const int k, const int j, const int i) {
  return x1f(i)*dx2f(j);
}

//----------------------------------------------------------------------------------------
// VolCenterXLength functions: compute physical length connecting cell centers as vector
// VolCenter2(i,j,k) located at (i,j+1/2,k), i.e. (x1v(i), x2f(j+1), x3v(k))

void Coordinates::VolCenter2Length(const int k, const int j, const int il, const int iu,
                                   AthenaArray<Real> &len) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // length2 = r d(theta)
    len(i) = x1v(i)*dx2v(j);
  }
  return;
}

//----------------------------------------------------------------------------------------
// CenterWidthX functions: return physical width in X-dir at (i,j,k) cell-center

void Coordinates::CenterWidth2(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &dx2) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    dx2(i) = x1v(i)*dx2f(j);
  }
  return;
}

//----------------------------------------------------------------------------------------
// FaceXArea functions: compute area of face with normal in X-dir as vector

void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
                            AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // area1 = r dphi dz
    area(i) = x1f(i)*dx2f(j)*dx3f(k);
  }
  return;
}

void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
                            AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // area3 = dr r dphi = d(r^2/2) dphi
    area(i) = coord_area3_i_(i)*dx2f(j);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetFaceXArea functions: return area of face with normal in X-dir at (i,j,k)

Real Coordinates::GetFace1Area(const int k, const int j, const int i) {
  // area1 = r dphi dz
  return x1f(i)*dx2f(j)*dx3f(k);
}

Real Coordinates::GetFace3Area(const int k, const int j, const int i) {
  // area3 = dr r dphi = d(r^2/2) dphi
  return coord_area3_i_(i)*dx2f(j);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
// VolCenterFaceXArea functions: compute area of face with normal in X-dir as vector
// where the faces are joined by cell centers (for non-ideal MHD)

void Coordinates::VolCenterFace1Area(const int k, const int j, const int il, const int iu,
                                     AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // area1 = r dphi dz
    area(i) = x1v(i)*dx2v(j)*dx3v(k);
  }
  return;
}

void Coordinates::VolCenterFace3Area(const int k, const int j, const int il, const int iu,
                                     AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // area3 = dr r dtheta = d(r^2/2) dtheta
    area(i) = coord_area3vc_i_(i)*dx2v(j);
  }
  return;
}
// Cell Volume function: compute volume of cell as vector

void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
                             AthenaArray<Real> &vol) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // volume = dr dz r dphi = d(r^2/2) dphi dz
    vol(i) = coord_vol_i_(i)*dx2f(j)*dx3f(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetCellVolume: returns cell volume at (i,j,k)

Real Coordinates::GetCellVolume(const int k, const int j, const int i) {
  return coord_vol_i_(i)*dx2f(j)*dx3f(k);
}

//----------------------------------------------------------------------------------------
//! Coordinate (Geometric) source term function

void Coordinates::AddCoordTermsDivergence(
    const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &u) {
  Real iso_cs = pmy_block->peos->GetIsoSoundSpeed();

  HydroDiffusion &hd = pmy_block->phydro->hdif;
  bool do_hydro_diffusion = (hd.hydro_diffusion_defined &&
                             (hd.nu_iso > 0.0 || hd.nu_aniso > 0.0));

  for (int k=pmy_block->ks; k<=pmy_block->ke; ++k) {
    for (int j=pmy_block->js; j<=pmy_block->je; ++j) {
#pragma omp simd
      for (int i=pmy_block->is; i<=pmy_block->ie; ++i) {
        // src_1 = <M_{phi phi}><1/r>
        // Skinner and Ostriker (2010) eq. 11a
        Real m_pp = prim(IDN,k,j,i)*prim(IM2,k,j,i)*prim(IM2,k,j,i);
        if (NON_BAROTROPIC_EOS) {
          m_pp += prim(IEN,k,j,i);
        } else {
          m_pp += (iso_cs*iso_cs)*prim(IDN,k,j,i);
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          m_pp += 0.5*(SQR(bcc(IB1,k,j,i)) - SQR(bcc(IB2,k,j,i)) + SQR(bcc(IB3,k,j,i)) );
        }
        if (!STS_ENABLED) {
          if (do_hydro_diffusion)
            m_pp += 0.5*(hd.visflx[X2DIR](IM2,k,j+1,i) + hd.visflx[X2DIR](IM2,k,j,i));
        }
        u(IM1,k,j,i) += dt*coord_src1_i_(i)*m_pp;

        // src_2 = -< M_{phi r} ><1/r>
        const Real& x_i   = x1f(i);
        const Real& x_ip1 = x1f(i+1);
        // Ju PhD thesis equation 2.14
        u(IM2,k,j,i) -= dt*coord_src2_i_(i)*(x_i*flux[X1DIR](IM2,k,j,i)
                                             + x_ip1*flux[X1DIR](IM2,k,j,i+1));
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! Coordinate (Geometric) source term function for STS

void Coordinates::AddCoordTermsDivergence_STS(
    const Real dt, int stage, const AthenaArray<Real> *flux,
    AthenaArray<Real> &u, AthenaArray<Real> &flux_div) {

  HydroDiffusion &hd = pmy_block->phydro->hdif;
  bool do_hydro_diffusion = (hd.hydro_diffusion_defined &&
                             (hd.nu_iso > 0.0 || hd.nu_aniso > 0.0));

  if (do_hydro_diffusion) {
    for (int k=pmy_block->ks; k<=pmy_block->ke; ++k) {
      for (int j=pmy_block->js; j<=pmy_block->je; ++j) {
#pragma omp simd
        for (int i=pmy_block->is; i<=pmy_block->ie; ++i) {
          // src_1 = <M_{phi phi}><1/r>
          Real m_pp = 0.5*(hd.visflx[X2DIR](IM2,k,j+1,i) + hd.visflx[X2DIR](IM2,k,j,i));
          u(IM1,k,j,i) += dt*coord_src1_i_(i)*m_pp;

          // src_2 = -< M_{phi r} ><1/r>
          const Real& x_i   = x1f(i);
          const Real& x_ip1 = x1f(i+1);
          u(IM2,k,j,i) -= dt*coord_src2_i_(i)*(x_i*flux[X1DIR](IM2,k,j,i)
                                               + x_ip1*flux[X1DIR](IM2,k,j,i+1));

          if (stage == 1 && pmy_block->pmy_mesh->sts_integrator=="rkl2") {
            flux_div(IM1,k,j,i) += 0.5*pmy_block->pmy_mesh->dt*coord_src1_i_(i)*m_pp;
            flux_div(IM2,k,j,i) -= 0.5*pmy_block->pmy_mesh->dt*coord_src2_i_(i)
                                   * (x_i*flux[X1DIR](IM2,k,j,i)
                                      + x_ip1*flux[X1DIR](IM2,k,j,i+1));
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Coordinate (Geometric) source term function for cosmic rays

void Coordinates::AddCRCoordTermsDivergence(
  const AthenaArray<Real> &u_input, AthenaArray<Real> &coord_src) {
  if (CR_ENABLED) {
    CosmicRay *pcr=pmy_block->pcr;
    for (int k=pmy_block->ks; k<=pmy_block->ke; ++k) {
      for (int j=pmy_block->js; j<=pmy_block->je; ++j) {
#pragma omp simd
        for (int i=pmy_block->is; i<=pmy_block->ie; ++i) {
          // src_1 = <M_{phi phi}><1/r>
          Real m_pp =  u_input(CRE,k,j,i)/3.0;
          coord_src(CRF1,k,j,i) =  pcr->vmax * coord_src1_i_(i)*m_pp;
          // 0 for other components
          coord_src(CRE,k,j,i) = 0.0;
          coord_src(CRF2,k,j,i) = 0.0;
          coord_src(CRF3,k,j,i) = 0.0;
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
// Coordinate (Geometric) source term for Grad Pc

void Coordinates::CRGradPcCoordTermsDivergence(
    const AthenaArray<Real> &u_cr, AthenaArray<Real> &grad_pc) {
  // Go through cellscosmicay
  if (CR_ENABLED) {
    //CosmicRay *pcr=pmy_block->pcr;
    for (int k=pmy_block->ks; k<=pmy_block->ke; ++k) {
      for (int j=pmy_block->js; j<=pmy_block->je; ++j) {
        for (int i=pmy_block->is; i<=pmy_block->ie; ++i) {
          // src_1 = <M_{phi phi}><1/r>
          Real m_ii = u_cr(CRE,k,j,i)/3.0;
          grad_pc(0,k,j,i) -= coord_src1_i_(i)*m_ii;
        }
      }// end j
    }// end k
  }// end CR_ENABLED
  return;
}

//--------------------------------------------------------------------------
// For radiation angles
void Coordinates::AxisDirection(int *axisx, int *axisy, int *axisz) {
  *axisx = 0;
  *axisy = 1;
  *axisz = 2;
}

void Coordinates::ConvertAngle(MeshBlock *pmb, const int nang,
                               AthenaArray<Real> &mu) {
  if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
    int n1z = pmy_block->ncells1, n2z = pmy_block->ncells2, n3z = pmy_block->ncells3;
    // int ndim=1;
    // if (n2z > 1) ndim = 2;
    // if (n3z > 1) ndim = 3;

    for (int k=0; k<n3z; ++k) {
      for (int j=0; j<n2z; ++j) {
        const Real& x2 = x2v(j);
        Real cosx2 = cos(x2);
        Real sinx2 = sin(x2);
        if (n2z == 1) {
          cosx2 = 1.0;
          sinx2 = 0.0;
        }
        for (int i=0; i<n1z; ++i) {
          Real *miux = &(mu(0,k,j,i,0));
          Real *miuy = &(mu(1,k,j,i,0));
          for (int mi=0; mi<nang; ++mi) {
            Real miux0 = miux[mi];
            Real miuy0 = miuy[mi];
            miux[mi] = miux0 * cosx2 + miuy0 * sinx2;
            miuy[mi] = miuy0 * cosx2 - miux0 * sinx2;
          }
        }
      }
    }
  }
}
