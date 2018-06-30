//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file spherical_polar.cpp
//  \brief implements functions for spherical polar (r-theta-phi) coordinates in a
//  derived class of the Coordinates abstract base class.

// C/C++ headers
#include <math.h>  // pow, trig functions

// Athena++ headers
#include "coordinates.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"

//----------------------------------------------------------------------------------------
// Spherical polar coordinates constructor

SphericalPolar::SphericalPolar(MeshBlock *pmb, ParameterInput *pin, bool flag)
  : Coordinates(pmb, pin, flag) {
  pmy_block = pmb;
  coarse_flag=flag;
  int il, iu, jl, ju, kl, ku, ng;
  if (coarse_flag==true) {
    il = pmb->cis; jl = pmb->cjs; kl = pmb->cks;
    iu = pmb->cie; ju = pmb->cje; ku = pmb->cke;
    ng=pmb->cnghost;
  } else {
    il = pmb->is; jl = pmb->js; kl = pmb->ks;
    iu = pmb->ie; ju = pmb->je; ku = pmb->ke;
    ng=NGHOST;
  }
  Mesh *pm=pmy_block->pmy_mesh;
  RegionSize& block_size = pmy_block->block_size;

  // allocate arrays for volume-centered coordinates and positions of cells
  int ncells1 = (iu-il+1) + 2*ng;
  int ncells2 = 1, ncells3 = 1;
  if (block_size.nx2 > 1) ncells2 = (ju-jl+1) + 2*ng;
  if (block_size.nx3 > 1) ncells3 = (ku-kl+1) + 2*ng;
  dx1v.NewAthenaArray(ncells1);
  dx2v.NewAthenaArray(ncells2);
  dx3v.NewAthenaArray(ncells3);
  x1v.NewAthenaArray(ncells1);
  x2v.NewAthenaArray(ncells2);
  x3v.NewAthenaArray(ncells3);
  // allocate arrays for volume- and face-centered geometry coefficients of cells
  h2f.NewAthenaArray(ncells1);
  dh2fd1.NewAthenaArray(ncells1);
  h31f.NewAthenaArray(ncells1);
  dh31fd1.NewAthenaArray(ncells1);
  h32f.NewAthenaArray(ncells2);
  dh32fd2.NewAthenaArray(ncells2);
  h2v.NewAthenaArray(ncells1);
  dh2vd1.NewAthenaArray(ncells1);
  h31v.NewAthenaArray(ncells1);
  dh31vd1.NewAthenaArray(ncells1);
  h32v.NewAthenaArray(ncells2);
  dh32vd2.NewAthenaArray(ncells2);

  // allocate arrays for area weighted positions for AMR/SMR MHD
  if ((pm->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
    x1s2.NewAthenaArray(ncells1);
    x1s3.NewAthenaArray(ncells1);
    x2s1.NewAthenaArray(ncells2);
    x2s3.NewAthenaArray(ncells2);
    x3s1.NewAthenaArray(ncells3);
    x3s2.NewAthenaArray(ncells3);
  }

  // initialize volume-averaged coordinates and spacing
  // x1-direction: x1v = (\int r dV / \int dV) = d(r^4/4)/d(r^3/3)
  for (int i=il-ng; i<=iu+ng; ++i) {
    x1v(i) = 0.75*(pow(x1f(i+1),4) - pow(x1f(i),4))/(pow(x1f(i+1),3) - pow(x1f(i),3));
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
      x2v(j) = ((sin(x2f(j+1)) - x2f(j+1)*cos(x2f(j+1))) -
                (sin(x2f(j  )) - x2f(j  )*cos(x2f(j  ))))/
                (cos(x2f(j  )) - cos(x2f(j+1)));
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
    h32v(jl) = sin(x2v(jl));
    h32f(jl) = sin(x2f(jl));
    dh32vd2(jl) = cos(x2v(jl));
    dh32fd2(jl) = cos(x2f(jl));
  } else {
    for (int j=jl-ng; j<=ju+ng; ++j) {
      h32v(j) = sin(x2v(j));
      h32f(j) = sin(x2f(j));
      dh32vd2(j) = cos(x2v(j));
      dh32fd2(j) = cos(x2f(j));
    }
  }

  // initialize area-averaged coordinates used with MHD AMR
  if ((pmb->pmy_mesh->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
    for (int i=il-ng; i<=iu+ng; ++i) {
      x1s2(i) = x1s3(i) = (2.0/3.0)*(pow(x1f(i+1),3) - pow(x1f(i),3))
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
  if (coarse_flag==false) {
    coord_area1_i_.NewAthenaArray(ncells1+1);
    coord_area2_i_.NewAthenaArray(ncells1);
    coord_area3_i_.NewAthenaArray(ncells1);
    coord_vol_i_.NewAthenaArray(ncells1);
    coord_src1_i_.NewAthenaArray(ncells1);
    coord_src2_i_.NewAthenaArray(ncells1);
    phy_src1_i_.NewAthenaArray(ncells1);
    phy_src2_i_.NewAthenaArray(ncells1);

    coord_area1_j_.NewAthenaArray(ncells2);
    coord_area2_j_.NewAthenaArray(ncells2+1);
    coord_vol_j_.NewAthenaArray(ncells2);
    coord_src1_j_.NewAthenaArray(ncells2);
    coord_src2_j_.NewAthenaArray(ncells2);
    coord_src3_j_.NewAthenaArray(ncells2);
    // non-ideal MHD
    coord_area1vc_i_.NewAthenaArray(ncells1);
    coord_area2vc_i_.NewAthenaArray(ncells1);
    coord_area3vc_i_.NewAthenaArray(ncells1);
    coord_area1vc_j_.NewAthenaArray(ncells2);
    coord_area2vc_j_.NewAthenaArray(ncells2);
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
        Real sm = fabs(sin(x2f(j  )));
        Real sp = fabs(sin(x2f(j+1)));
        Real cm = cos(x2f(j  ));
        Real cp = cos(x2f(j+1));
        // d(sin theta) = d(-cos theta)
        coord_area1_j_(j) = fabs(cm - cp);
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
        // d(sin theta) = d(-cos theta) at the volume center for non-ideal MHD
        coord_area1vc_j_(j)= fabs(cos(x2v(j))-cos(x2v(j+1)));
        // sin theta at the volume center for non-ideal MHD
        coord_area2vc_j_(j)= fabs(sin(x2v(j)));
      }
      coord_area2_j_(ju+ng+1) = fabs(sin(x2f(ju+ng+1)));
      if (IsPole(jl))   // inner polar boundary
        coord_area1vc_j_(jl-1)= 2.0-cos(x2v(jl-1))-cos(x2v(jl));
      if (IsPole(ju))   // outer polar boundary
        coord_area1vc_j_(ju)  = 2.0+cos(x2v(ju))+cos(x2v(ju+1));
    } else {
      Real sm = fabs(sin(x2f(jl  )));
      Real sp = fabs(sin(x2f(jl+1)));
      Real cm = cos(x2f(jl  ));
      Real cp = cos(x2f(jl+1));
      coord_area1_j_(jl) = fabs(cm - cp);
      coord_area2_j_(jl) = sm;
      coord_area1vc_j_(jl)= coord_area1_j_(jl);
      coord_area2vc_j_(jl)= sin(x2v(jl));
      coord_vol_j_(jl) = coord_area1_j_(jl);
      coord_src1_j_(jl) = (sp - sm)/coord_vol_j_(jl);
      coord_src2_j_(jl) = (sp - sm)/((sm + sp)*coord_vol_j_(jl));
      coord_src3_j_(jl) = (sp - sm)/coord_vol_j_(jl);
      coord_area2_j_(jl+1) = sp;
    }
  }
}

// destructor

SphericalPolar::~SphericalPolar() {
  dx1v.DeleteAthenaArray();
  dx2v.DeleteAthenaArray();
  dx3v.DeleteAthenaArray();
  x1v.DeleteAthenaArray();
  x2v.DeleteAthenaArray();
  x3v.DeleteAthenaArray();
  if ((pmy_block->pmy_mesh->multilevel==true) && MAGNETIC_FIELDS_ENABLED) {
    x1s2.DeleteAthenaArray();
    x1s3.DeleteAthenaArray();
    x2s1.DeleteAthenaArray();
    x2s3.DeleteAthenaArray();
    x3s1.DeleteAthenaArray();
    x3s2.DeleteAthenaArray();
  }
  if (coarse_flag==false) {
    coord_area1_i_.DeleteAthenaArray();
    coord_area2_i_.DeleteAthenaArray();
    coord_area3_i_.DeleteAthenaArray();
    coord_vol_i_.DeleteAthenaArray();
    coord_src1_i_.DeleteAthenaArray();
    coord_src2_i_.DeleteAthenaArray();
    phy_src1_i_.DeleteAthenaArray();
    phy_src2_i_.DeleteAthenaArray();

    coord_area1_j_.DeleteAthenaArray();
    coord_area2_j_.DeleteAthenaArray();
    coord_vol_j_.DeleteAthenaArray();
    coord_src1_j_.DeleteAthenaArray();
    coord_src2_j_.DeleteAthenaArray();
    coord_src3_j_.DeleteAthenaArray();

    // non-ideal MHD
    coord_area1vc_i_.DeleteAthenaArray();
    coord_area2vc_i_.DeleteAthenaArray();
    coord_area3vc_i_.DeleteAthenaArray();
    coord_area1vc_j_.DeleteAthenaArray();
    coord_area2vc_j_.DeleteAthenaArray();
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
    // length3 = r sin(theta) d(phi)
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
    // length3 = r sin(theta) d(phi)
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
    dx3(i) = x1v(i)*fabs(sin(x2v(j)))*dx3f(k);
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
    // volume = r^2 sin(theta) dr dtheta dphi = d(r^3/3) d(-cos theta) dphi
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

void SphericalPolar::CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
                                   const AthenaArray<Real> &prim,
                                   const AthenaArray<Real> &bcc, AthenaArray<Real> &u) {
  Real iso_cs = pmy_block->peos->GetIsoSoundSpeed();
  bool use_x2_fluxes = pmy_block->block_size.nx2 > 1;

  HydroDiffusion *phd = pmy_block->phydro->phdif;
  bool do_hydro_diffusion = (phd->hydro_diffusion_defined &&
                            (phd->nu_iso>0.0 || phd->nu_aniso>0.0));

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
        if (do_hydro_diffusion) {
          m_ii += 0.5*(phd->visflx[X2DIR](IM2,k,j+1,i)+phd->visflx[X2DIR](IM2,k,j,i));
          m_ii += 0.5*(phd->visflx[X3DIR](IM3,k+1,j,i)+phd->visflx[X3DIR](IM3,k,j,i));
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
        if (do_hydro_diffusion)
          m_pp += 0.5*(phd->visflx[X3DIR](IM3,k+1,j,i)+phd->visflx[X3DIR](IM3,k,j,i));

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
          if (do_hydro_diffusion)
              m_ph += 0.5*(phd->visflx[X2DIR](IM3,k,j+1,i)+phd->visflx[X2DIR](IM3,k,j,i));

          u(IM3,k,j,i) -= dt*coord_src1_i_(i)*coord_src3_j_(j)*m_ph;
        }
      }
    }
  }

  return;
}
