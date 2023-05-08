// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file spherical_polar.cpp
//! \brief implements functions for spherical polar (r-theta-phi) coordinates in a
//! derived class of the Coordinates abstract base class.

// C headers

// C++ headers
#include <cmath>  // pow(), trig functions
#include <iomanip>
#include <iostream>   // endl
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "coordinates.hpp"

//----------------------------------------------------------------------------------------
// Spherical polar coordinates constructor

SphericalPolar::SphericalPolar(MeshBlock *pmb, ParameterInput *pin, bool flag)
    : Coordinates(pmb, pin, flag) {
  RegionSize& block_size = pmy_block->block_size;
  // check that Mesh's polar coordinate range does not exceed [0, pi], even in 2D
  // (use 2D cylindrical coordinates to create a circular Mesh)
  if (block_size.nx2 > 1
      && (pm->mesh_size.x2min < static_cast<Real>(0.0)
          || pm->mesh_size.x2max > static_cast<Real>(PI))) {
    std::stringstream msg;
    msg << "### FATAL ERROR in SphericalPolar constructor" << std::endl
        << "2D or 3D spherical-polar coordinates requires that\n"
        << "x2 does not exceed the following limits:\n"
        << std::setprecision(std::numeric_limits<Real>::max_digits10 -1)
        << "x2min=" << std::scientific << 0.0 << "\n"
        << "x2max=" << PI << "\n"
        << "Current x2 domain limits are: \n"
        << "x2min=" << pm->mesh_size.x2min << "\n"
        << "x2max=" << pm->mesh_size.x2max << std::endl;
    ATHENA_ERROR(msg);
  }

  // initialize volume-averaged coordinates and spacing
  // x1-direction: x1v = (\int r dV / \int dV) = d(r^4/4)/d(r^3/3)
  for (int i=il-ng; i<=iu+ng; ++i) {
    x1v(i) = 0.75*(std::pow(x1f(i+1), 4) - std::pow(x1f(i), 4)) /
             (std::pow(x1f(i+1), 3) - std::pow(x1f(i), 3));
    // reduces to eq for centroid: R_i + 2*R_i*dR_i^2/(12*R_i^2 + dR_i^2)
    // see Mignone (2014) eq 17, e.g.
  }
  for (int i=il-ng; i<=iu+ng-1; ++i) {
    dx1v(i) = x1v(i+1) - x1v(i);
  }

  // x2-direction: x2v = (\int sin[theta] theta dV / \int dV) =
  //   d(sin[theta] - theta cos[theta])/d(-cos[theta])
  if (pmb->block_size.nx2 == 1) {
    x2v(jl) = 0.5*(x2f(jl+1) + x2f(jl));
    dx2v(jl) = dx2f(jl);
  } else {
    for (int j=jl-ng; j<=ju+ng; ++j) {
      x2v(j) = ((std::sin(x2f(j+1)) - x2f(j+1)*std::cos(x2f(j+1))) -
                (std::sin(x2f(j  )) - x2f(j  )*std::cos(x2f(j  ))))/
               (std::cos(x2f(j  )) - std::cos(x2f(j+1)));
    }
    for (int j=jl-ng; j<=ju+ng-1; ++j) {
      dx2v(j) = x2v(j+1) - x2v(j);
    }
  }

  // x3-direction: x3v = (\int phi dV / \int dV) = dphi/2
  if (pmb->block_size.nx3 == 1) {
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
    h31v(i) = x1v(i);
    h31f(i) = x1f(i);
    dh2vd1(i) = 1.0;
    dh2fd1(i) = 1.0;
    dh31vd1(i) = 1.0;
    dh31fd1(i) = 1.0;
  }

  // x2-direction
  if (pmb->block_size.nx2 == 1) {
    h32v(jl) = std::sin(x2v(jl));
    h32f(jl) = std::sin(x2f(jl));
    dh32vd2(jl) = std::cos(x2v(jl));
    dh32fd2(jl) = std::cos(x2f(jl));
  } else {
    for (int j=jl-ng; j<=ju+ng; ++j) {
      h32v(j) = std::sin(x2v(j));
      h32f(j) = std::sin(x2f(j));
      dh32vd2(j) = std::cos(x2v(j));
      dh32fd2(j) = std::cos(x2f(j));
    }
  }

  // initialize area-averaged coordinates used with MHD AMR
  if ((pmb->pmy_mesh->multilevel) && MAGNETIC_FIELDS_ENABLED) {
    for (int i=il-ng; i<=iu+ng; ++i) {
      x1s2(i) = x1s3(i) = (2.0/3.0)*(std::pow(x1f(i+1),3) - std::pow(x1f(i),3))
                /(SQR(x1f(i+1)) - SQR(x1f(i)));
    }
    if (pmb->block_size.nx2 == 1) {
      x2s1(jl) = x2s3(jl) = x2v(jl);
    } else {
      for (int j=jl-ng; j<=ju+ng; ++j) {
        x2s1(j) = x2s3(j) = x2v(j);
      }
    }
    if (pmb->block_size.nx3 == 1) {
      x3s1(kl) = x3s2(kl) = x3v(kl);
    } else {
      for (int k=kl-ng; k<=ku+ng; ++k) {
        x3s1(k) = x3s2(k) = x3v(k);
      }
    }
  }

  // Allocate memory for internal scratch arrays to store partial calculations
  // (note this is skipped if object is for coarse mesh with AMR)
  if (!coarse_flag) {
    coord_area1_i_.NewAthenaArray(nc1+1);
    coord_area2_i_.NewAthenaArray(nc1);
    coord_area3_i_.NewAthenaArray(nc1);
    coord_vol_i_.NewAthenaArray(nc1);
    coord_src1_i_.NewAthenaArray(nc1);
    coord_src2_i_.NewAthenaArray(nc1);
    phy_src1_i_.NewAthenaArray(nc1);
    phy_src2_i_.NewAthenaArray(nc1);

    coord_area1_j_.NewAthenaArray(nc2);
    coord_area2_j_.NewAthenaArray(nc2+1);
    coord_vol_j_.NewAthenaArray(nc2);
    coord_src1_j_.NewAthenaArray(nc2);
    coord_src2_j_.NewAthenaArray(nc2);
    coord_src3_j_.NewAthenaArray(nc2);
    // non-ideal MHD
    coord_area1vc_i_.NewAthenaArray(nc1);
    coord_area2vc_i_.NewAthenaArray(nc1);
    coord_area3vc_i_.NewAthenaArray(nc1);
    coord_area1vc_j_.NewAthenaArray(nc2);
    coord_area2vc_j_.NewAthenaArray(nc2);
    // Compute and store constant coefficients needed for face-areas, cell-volumes, etc.
    // This helps improve performance.
#pragma omp simd
    for (int i=il-ng; i<=iu+ng; ++i) {
      Real rm = x1f(i  );
      Real rp = x1f(i+1);
      // R^2
      coord_area1_i_(i) = rm*rm;
      // 0.5*(R_{i+1}^2 - R_{i}^2)
      coord_area2_i_(i) = 0.5*(rp*rp - rm*rm);
      // 0.5*(R_{i+1}^2 - R_{i}^2)
      coord_area3_i_(i) = coord_area2_i_(i);
      // dV = (R_{i+1}^3 - R_{i}^3)/3
      coord_vol_i_(i) = (ONE_3RD)*(rp*rp*rp - rm*rm*rm);
      // (A1^{+} - A1^{-})/dV
      coord_src1_i_(i) = coord_area2_i_(i)/coord_vol_i_(i);
      // (dR/2)/(R_c dV)
      coord_src2_i_(i) = dx1f(i)/((rm + rp)*coord_vol_i_(i));
      // Rf_{i}^2/R_{i}^2/Rf_{i}^2
      phy_src1_i_(i) = 1.0/SQR(x1v(i));
      // Rf_{i+1}^2/R_{i}^2/Rf_{i+1}^2
      phy_src2_i_(i) = phy_src1_i_(i);
      // R^2 at the volume center for non-ideal MHD
      coord_area1vc_i_(i) = SQR(x1v(i));
    }
    coord_area1_i_(iu+ng+1) = x1f(iu+ng+1)*x1f(iu+ng+1);
#pragma omp simd
    for (int i=il-ng; i<=iu+ng-1; ++i) {//non-ideal MHD
      // 0.5*(R_{i+1}^2 - R_{i}^2)
      coord_area2vc_i_(i)= 0.5*(SQR(x1v(i+1))-SQR(x1v(i)));
      // 0.5*(R_{i+1}^2 - R_{i}^2)
      coord_area3vc_i_(i)= coord_area2vc_i_(i);
    }
    if (pmb->block_size.nx2 > 1) {
#pragma omp simd
      for (int j=jl-ng; j<=ju+ng; ++j) {
        Real sm = std::abs(std::sin(x2f(j  )));
        Real sp = std::abs(std::sin(x2f(j+1)));
        Real cm = std::cos(x2f(j  ));
        Real cp = std::cos(x2f(j+1));
        // d(sin theta) = d(-cos theta)
        coord_area1_j_(j) = std::abs(cm - cp);
        // sin theta
        coord_area2_j_(j) = sm;
        // d(sin theta) = d(-cos theta)
        coord_vol_j_(j) = coord_area1_j_(j);
        // (A2^{+} - A2^{-})/dV
        coord_src1_j_(j) = (sp - sm)/coord_vol_j_(j);
        // (dS/2)/(S_c dV)
        coord_src2_j_(j) = (sp - sm)/((sm + sp)*coord_vol_j_(j));
        // < cot theta > = (|sin th_p| - |sin th_m|) / |cos th_m - cos th_p|
        coord_src3_j_(j) = (sp - sm)/coord_vol_j_(j);
        // sin theta at the volume center for non-ideal MHD
        coord_area2vc_j_(j)= std::abs(std::sin(x2v(j)));
      }
#pragma omp simd
      for (int j=jl-ng; j<=ju+ng-1; ++j) {
        // d(sin theta) = d(-cos theta) at the volume center for non-ideal MHD
        coord_area1vc_j_(j)= std::abs(cos(x2v(j))-cos(x2v(j+1)));
      }
      coord_area2_j_(ju+ng+1) = std::abs(sin(x2f(ju+ng+1)));
      if (IsPole(jl))   // inner polar boundary
        coord_area1vc_j_(jl-1)= 2.0-std::cos(x2v(jl-1))-std::cos(x2v(jl));
      if (IsPole(ju+1))   // outer polar boundary
        coord_area1vc_j_(ju)  = 2.0+std::cos(x2v(ju))+std::cos(x2v(ju+1));
    } else {
      Real sm = std::abs(std::sin(x2f(jl  )));
      Real sp = std::abs(std::sin(x2f(jl+1)));
      Real cm = std::cos(x2f(jl  ));
      Real cp = std::cos(x2f(jl+1));
      coord_area1_j_(jl) = std::abs(cm - cp);
      coord_area2_j_(jl) = sm;
      coord_area1vc_j_(jl)= coord_area1_j_(jl);
      coord_area2vc_j_(jl)= std::sin(x2v(jl));
      coord_vol_j_(jl) = coord_area1_j_(jl);
      coord_src1_j_(jl) = (sp - sm)/coord_vol_j_(jl);
      coord_src2_j_(jl) = (sp - sm)/((sm + sp)*coord_vol_j_(jl));
      coord_src3_j_(jl) = (sp - sm)/coord_vol_j_(jl);
      coord_area2_j_(jl+1) = sp;
    }
  }
}


//----------------------------------------------------------------------------------------
// EdgeXLength functions: compute physical length at cell edge-X as vector
// Edge1(i,j,k) located at (i,j-1/2,k-1/2), i.e. (x1v(i), x2f(j), x3f(k))
// Edge2(i,j,k) located at (i-1/2,j,k-1/2), i.e. (x1f(i), x2v(j), x3f(k))

void SphericalPolar::Edge2Length(const int k, const int j, const int il, const int iu,
                                 AthenaArray<Real> &len) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // length2 = r d(theta)
    len(i) = x1f(i)*dx2f(j);
  }
  return;
}

// Edge3(i,j,k) located at (i-1/2,j-1/2,k), i.e. (x1f(i), x2f(j), x3v(k))

void SphericalPolar::Edge3Length(const int k, const int j, const int il, const int iu,
                                 AthenaArray<Real> &len) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // length3 = r std::sin(theta) d(phi)
    len(i) = x1f(i)*coord_area2_j_(j)*dx3f(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetEdgeXLength functions: return length of edge-X at (i,j,k)

Real SphericalPolar::GetEdge2Length(const int k, const int j, const int i) {
  return x1f(i)*dx2f(j);
}

Real SphericalPolar::GetEdge3Length(const int k, const int j, const int i) {
  return x1f(i)*coord_area2_j_(j)*dx3f(k);
}

//----------------------------------------------------------------------------------------
// VolCenterXLength functions: compute physical length connecting cell centers as vector
// VolCenter1(i,j,k) located at (i+1/2,j,k), i.e. (x1f(i+1), x2v(j), x3v(k))
// VolCenter2(i,j,k) located at (i,j+1/2,k), i.e. (x1v(i), x2f(j+1), x3v(k))

void SphericalPolar::VolCenter2Length(const int k, const int j, const int il,
                                      const int iu, AthenaArray<Real> &len) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // length2 = r d(theta)
    len(i) = x1v(i)*dx2v(j);
  }
  return;
}

// VolCenter3(i,j,k) located at (i,j,k+1/2), i.e. (x1v(i), x2v(j), x3f(k+1))
void SphericalPolar::VolCenter3Length(const int k, const int j, const int il,
                                      const int iu, AthenaArray<Real> &len) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // length3 = r std::sin(theta) d(phi)
    len(i) = x1v(i)*coord_area2vc_j_(j)*dx3v(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// CenterWidthX functions: return physical width in X-dir at (i,j,k) cell-center

void SphericalPolar::CenterWidth2(const int k, const int j, const int il, const int iu,
                                  AthenaArray<Real> &dx2) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    dx2(i) = x1v(i)*dx2f(j);
  }
  return;
}

void SphericalPolar::CenterWidth3(const int k, const int j, const int il, const int iu,
                                  AthenaArray<Real> &dx3) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    dx3(i) = x1v(i)*std::abs(std::sin(x2v(j)))*dx3f(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// FaceXArea functions: compute area of face with normal in X-dir as vector

void SphericalPolar::Face1Area(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // area1 = r^2 sin[theta] dtheta dphi = r^2 d(-cos[theta]) dphi
    area(i) = coord_area1_i_(i)*coord_area1_j_(j)*dx3f(k);
  }
  return;
}

void SphericalPolar::Face2Area(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // area2 = dr r sin[theta] dphi = d(r^2/2) sin[theta] dphi
    area(i) = coord_area2_i_(i)*coord_area2_j_(j)*dx3f(k);
  }
  return;
}

void SphericalPolar::Face3Area(const int k, const int j, const int il, const int iu,
                               AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // area3 = dr r dtheta = d(r^2/2) dtheta
    area(i) = coord_area3_i_(i)*dx2f(j);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetFaceXArea functions: return area of face with normal in X-dir at (i,j,k)

Real SphericalPolar::GetFace1Area(const int k, const int j, const int i) {
  return (coord_area1_i_(i)*coord_area1_j_(j)*dx3f(k));
}

Real SphericalPolar::GetFace2Area(const int k, const int j, const int i) {
  return (coord_area2_i_(i)*coord_area2_j_(j)*dx3f(k));
}

Real SphericalPolar::GetFace3Area(const int k, const int j, const int i) {
  return (coord_area3_i_(i)*dx2f(j));
}

//----------------------------------------------------------------------------------------
// VolCenterFaceXArea functions: compute area of face with normal in X-dir as vector
// where the faces are joined by cell centers (for non-ideal MHD)

void SphericalPolar::VolCenterFace1Area(const int k, const int j, const int il,
                                        const int iu, AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // area1 = r^2 sin[theta] dtheta dphi = r^2 d(-cos[theta]) dphi
    area(i) = coord_area1vc_i_(i)*coord_area1vc_j_(j)*dx3v(k);
  }
  return;
}
void SphericalPolar::VolCenterFace2Area(const int k, const int j, const int il,
                                        const int iu, AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // area2 = dr r sin[theta] dphi = d(r^2/2) sin[theta] dphi
    area(i) = coord_area2vc_i_(i)*coord_area2vc_j_(j)*dx3v(k);
  }
  return;
}
void SphericalPolar::VolCenterFace3Area(const int k, const int j, const int il,
                                        const int iu, AthenaArray<Real> &area) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // area3 = dr r dtheta = d(r^2/2) dtheta
    area(i) = coord_area3vc_i_(i)*dx2v(j);
  }
  return;
}

//----------------------------------------------------------------------------------------
// Cell Volume function: compute volume of cell as vector

void SphericalPolar::CellVolume(const int k, const int j, const int il, const int iu,
                                AthenaArray<Real> &vol) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    // volume = r^2 std::sin(theta) dr dtheta dphi = d(r^3/3) d(-cos theta) dphi
    vol(i) = coord_vol_i_(i)*coord_vol_j_(j)*dx3f(k);
  }
  return;
}

//----------------------------------------------------------------------------------------
// GetCellVolume: returns cell volume at (i,j,k)

Real SphericalPolar::GetCellVolume(const int k, const int j, const int i) {
  return coord_vol_i_(i)*coord_vol_j_(j)*dx3f(k);
}

//----------------------------------------------------------------------------------------
// Coordinate (Geometric) source term function

void SphericalPolar::AddCoordTermsDivergence(const Real dt, const AthenaArray<Real> *flux,
                                   const AthenaArray<Real> &prim,
                                   const AthenaArray<Real> &bcc, AthenaArray<Real> &u) {
  Real iso_cs = pmy_block->peos->GetIsoSoundSpeed();
  bool use_x2_fluxes = pmy_block->block_size.nx2 > 1;

  HydroDiffusion &hd = pmy_block->phydro->hdif;
  bool do_hydro_diffusion = (hd.hydro_diffusion_defined &&
                             (hd.nu_iso > 0.0 || hd.nu_aniso > 0.0));

  // Go through cells
  for (int k=pmy_block->ks; k<=pmy_block->ke; ++k) {
    for (int j=pmy_block->js; j<=pmy_block->je; ++j) {
#pragma omp simd
      for (int i=pmy_block->is; i<=pmy_block->ie; ++i) {
        // src_1 = < M_{theta theta} + M_{phi phi} ><1/r>
        Real m_ii = prim(IDN,k,j,i)*(SQR(prim(IM2,k,j,i)) + SQR(prim(IM3,k,j,i)));
        if (NON_BAROTROPIC_EOS) {
          m_ii += 2.0*prim(IEN,k,j,i);
        } else {
          m_ii += 2.0*(iso_cs*iso_cs)*prim(IDN,k,j,i);
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          m_ii += SQR(bcc(IB1,k,j,i));
        }
        if (!STS_ENABLED) {
          if (do_hydro_diffusion) {
            m_ii += 0.5*(hd.visflx[X2DIR](IM2,k,j+1,i) + hd.visflx[X2DIR](IM2,k,j,i));
            m_ii += 0.5*(hd.visflx[X3DIR](IM3,k+1,j,i) + hd.visflx[X3DIR](IM3,k,j,i));
          }
        }

        u(IM1,k,j,i) += dt*coord_src1_i_(i)*m_ii;

        // src_2 = -< M_{theta r} ><1/r>
        u(IM2,k,j,i) -= dt*coord_src2_i_(i)*
                        (coord_area1_i_(i)*flux[X1DIR](IM2,k,j,i)
                         + coord_area1_i_(i+1)*flux[X1DIR](IM2,k,j,i+1));

        // src_3 = -< M_{phi r} ><1/r>
        u(IM3,k,j,i) -= dt*coord_src2_i_(i)*
                        (coord_area1_i_(i)*flux[X1DIR](IM3,k,j,i)
                         + coord_area1_i_(i+1)*flux[X1DIR](IM3,k,j,i+1));

        // src_2 = < M_{phi phi} ><cot theta/r>
        Real m_pp = prim(IDN,k,j,i)*SQR(prim(IM3,k,j,i));
        if (NON_BAROTROPIC_EOS) {
          m_pp += prim(IEN,k,j,i);
        } else {
          m_pp += (iso_cs*iso_cs)*prim(IDN,k,j,i);
        }
        if (MAGNETIC_FIELDS_ENABLED) {
          m_pp += 0.5*( SQR(bcc(IB1,k,j,i)) + SQR(bcc(IB2,k,j,i))
                        - SQR(bcc(IB3,k,j,i)) );
        }
        if (!STS_ENABLED) {
          if (do_hydro_diffusion)
            m_pp += 0.5*(hd.visflx[X3DIR](IM3,k+1,j,i) + hd.visflx[X3DIR](IM3,k,j,i));
        }

        u(IM2,k,j,i) += dt*coord_src1_i_(i)*coord_src1_j_(j)*m_pp;

        // src_3 = -< M_{phi theta} ><cot theta/r>
        if (use_x2_fluxes) {
          u(IM3,k,j,i) -= dt*coord_src1_i_(i)*coord_src2_j_(j)*
                          (coord_area2_j_(j)*flux[X2DIR](IM3,k,j,i)
                           + coord_area2_j_(j+1)*flux[X2DIR](IM3,k,j+1,i));
        } else {
          Real m_ph = prim(IDN,k,j,i) * prim(IM3,k,j,i) * prim(IM2,k,j,i);
          if (MAGNETIC_FIELDS_ENABLED) {
            m_ph -= bcc(IB3,k,j,i) * bcc(IB2,k,j,i);
          }
          if (!STS_ENABLED) {
            if (do_hydro_diffusion)
              m_ph += 0.5*(hd.visflx[X2DIR](IM3,k,j+1,i) + hd.visflx[X2DIR](IM3,k,j,i));
          }

          u(IM3,k,j,i) -= dt*coord_src1_i_(i)*coord_src3_j_(j)*m_ph;
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
// Coordinate (Geometric) source term function for STS

void SphericalPolar::AddCoordTermsDivergence_STS(const Real dt, int stage,
                                   const AthenaArray<Real> *flux,
                                   AthenaArray<Real> &u, AthenaArray<Real> &flux_div) {
  Real iso_cs = pmy_block->peos->GetIsoSoundSpeed();
  bool use_x2_fluxes = pmy_block->block_size.nx2 > 1;

  HydroDiffusion &hd = pmy_block->phydro->hdif;
  bool do_hydro_diffusion = (hd.hydro_diffusion_defined &&
                             (hd.nu_iso > 0.0 || hd.nu_aniso > 0.0));

  // Go through cells
  if (do_hydro_diffusion) {
    for (int k=pmy_block->ks; k<=pmy_block->ke; ++k) {
      for (int j=pmy_block->js; j<=pmy_block->je; ++j) {
#pragma omp simd
        for (int i=pmy_block->is; i<=pmy_block->ie; ++i) {
          // src_1 = < M_{theta theta} + M_{phi phi} ><1/r>
          Real m_ii = 0.5*(hd.visflx[X2DIR](IM2,k,j+1,i) + hd.visflx[X2DIR](IM2,k,j,i));
          m_ii += 0.5*(hd.visflx[X3DIR](IM3,k+1,j,i) + hd.visflx[X3DIR](IM3,k,j,i));

          u(IM1,k,j,i) += dt*coord_src1_i_(i)*m_ii;

          // src_2 = -< M_{theta r} ><1/r>
          u(IM2,k,j,i) -= dt*coord_src2_i_(i)*
                          (coord_area1_i_(i)*flux[X1DIR](IM2,k,j,i)
                           + coord_area1_i_(i+1)*flux[X1DIR](IM2,k,j,i+1));

          // src_3 = -< M_{phi r} ><1/r>
          u(IM3,k,j,i) -= dt*coord_src2_i_(i)*
                          (coord_area1_i_(i)*flux[X1DIR](IM3,k,j,i)
                           + coord_area1_i_(i+1)*flux[X1DIR](IM3,k,j,i+1));

          if (stage == 1 && pmy_block->pmy_mesh->sts_integrator == "rkl2") {
            flux_div(IM1,k,j,i) += 0.5*pmy_block->pmy_mesh->dt*coord_src1_i_(i)*m_ii;

            // src_2 = -< M_{theta r} ><1/r>
            flux_div(IM2,k,j,i) -= 0.5*pmy_block->pmy_mesh->dt*coord_src2_i_(i)*
                            (coord_area1_i_(i)*flux[X1DIR](IM2,k,j,i)
                             + coord_area1_i_(i+1)*flux[X1DIR](IM2,k,j,i+1));

            // src_3 = -< M_{phi r} ><1/r>
            flux_div(IM3,k,j,i) -= 0.5*pmy_block->pmy_mesh->dt*coord_src2_i_(i)*
                            (coord_area1_i_(i)*flux[X1DIR](IM3,k,j,i)
                             + coord_area1_i_(i+1)*flux[X1DIR](IM3,k,j,i+1));
          }

          // src_2 = < M_{phi phi} ><cot theta/r>
          Real m_pp = 0.5*(hd.visflx[X3DIR](IM3,k+1,j,i) + hd.visflx[X3DIR](IM3,k,j,i));

          u(IM2,k,j,i) += dt*coord_src1_i_(i)*coord_src1_j_(j)*m_pp;

          if (stage == 1 && pmy_block->pmy_mesh->sts_integrator == "rkl2") {
            // src_2 = -< M_{theta r} ><1/r>
            flux_div(IM2,k,j,i) += 0.5*pmy_block->pmy_mesh->dt
                                   * coord_src1_i_(i)*coord_src1_j_(j)*m_pp;
          }

          // src_3 = -< M_{phi theta} ><cot theta/r>
          if (use_x2_fluxes) {
            u(IM3,k,j,i) -= dt*coord_src1_i_(i)*coord_src2_j_(j)*
                            (coord_area2_j_(j)*flux[X2DIR](IM3,k,j,i)
                             + coord_area2_j_(j+1)*flux[X2DIR](IM3,k,j+1,i));
            if (stage == 1 && pmy_block->pmy_mesh->sts_integrator == "rkl2") {
              // src_2 = -< M_{theta r} ><1/r>
              flux_div(IM3,k,j,i) -= 0.5*pmy_block->pmy_mesh->dt
                                     * coord_src1_i_(i)*coord_src2_j_(j)*
                                       (coord_area2_j_(j)*flux[X2DIR](IM3,k,j,i)
                                        + coord_area2_j_(j+1)*flux[X2DIR](IM3,k,j+1,i));
            }
          } else {
            Real m_ph = 0.5*(hd.visflx[X2DIR](IM3,k,j+1,i) + hd.visflx[X2DIR](IM3,k,j,i));

            u(IM3,k,j,i) -= dt*coord_src1_i_(i)*coord_src3_j_(j)*m_ph;

            if (stage == 1 && pmy_block->pmy_mesh->sts_integrator == "rkl2") {
              // src_2 = -< M_{theta r} ><1/r>
              flux_div(IM3,k,j,i) -= 0.5*pmy_block->pmy_mesh->dt
                                     * coord_src1_i_(i)*coord_src3_j_(j)*m_ph;
            }
          }
        }
      }
    }
  }
  return;
}
