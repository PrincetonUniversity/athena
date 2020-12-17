//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file field_loop_poles.cpp
//! \brief Advection of a field loop THROUGH the poles in spherical_polar coordinates.
//!
//! Originally developed by ZZ.  Sets up constant uniform-density flow in x-direction
//! through poles, and follows advection of loop.  Set xz>0 (xz<0) for loop through
//! upper (lower) pole.  Works in 2D and 3D.
//========================================================================================

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <fstream>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

namespace {
void VelProfileCyl(const Real rad, const Real phi, const Real z,
                   Real &v1, Real &v2, Real &v3);
Real A3(const Real x1, const Real x2, const Real x3);
Real A2(const Real x1, const Real x2, const Real x3);
Real A1(const Real x1, const Real x2, const Real x3);

// problem parameters which are global variables with internal linkage
Real vy0, rho0, isocs2, gamma_gas;
Real xc, yc, zc, beta, b0;
} // namespace

// User-defined boundary conditions along inner/outer edges (not poles)
void LoopInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void LoopOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void LoopInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void LoopOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);

int RefinementCondition(MeshBlock *pmb);


//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem","rho0");
  vy0 = pin->GetOrAddReal("problem","vy0",0.0);

  // Get parameters of initial pressure and cooling parameters
  isocs2=SQR(pin->GetReal("hydro","iso_sound_speed"));
  if (NON_BAROTROPIC_EOS) {
    gamma_gas = pin->GetReal("hydro","gamma");
  }

  // Get loop center for field loop tests;
  xc = pin->GetOrAddReal("problem","xc",1.0);
  yc = pin->GetOrAddReal("problem","yc",0.0);
  zc = pin->GetOrAddReal("problem","zc",0.0);
  if (MAGNETIC_FIELDS_ENABLED) {
    beta = pin->GetReal("problem","beta");
    b0=std::sqrt(2.*isocs2*rho0/beta);
  }

  if (adaptive)
    EnrollUserRefinementCondition(RefinementCondition);

  // setup boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, LoopInnerX1);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, LoopOuterX1);
  }
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, LoopInnerX2);
  }
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, LoopOuterX2);
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes field loop advection through pole.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Real rad, phi, z;
  Real v1, v2, v3;
  // Set initial magnetic fields
  if (MAGNETIC_FIELDS_ENABLED) {
    AthenaArray<Real> a1, a2, a3;
    // nxN != ncellsN, in general. Allocate to extend through 2*ghost, regardless # dim
    int nx1 = block_size.nx1 + 2*NGHOST;
    int nx2 = block_size.nx2 + 2*NGHOST;
    int nx3 = block_size.nx3 + 2*NGHOST;
    a1.NewAthenaArray(nx3, nx2, nx1);
    a2.NewAthenaArray(nx3, nx2, nx1);
    a3.NewAthenaArray(nx3, nx2, nx1);

    int level = loc.level;
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          if ((pbval->nblevel[1][0][1]>level && j==js)
              || (pbval->nblevel[1][2][1]>level && j==je+1)
              || (pbval->nblevel[0][1][1]>level && k==ks)
              || (pbval->nblevel[2][1][1]>level && k==ke+1)
              || (pbval->nblevel[0][0][1]>level && j==js   && k==ks)
              || (pbval->nblevel[0][2][1]>level && j==je+1 && k==ks)
              || (pbval->nblevel[2][0][1]>level && j==js   && k==ke+1)
              || (pbval->nblevel[2][2][1]>level && j==je+1 && k==ke+1)) {
            Real x1l = pcoord->x1f(i)+0.25*pcoord->dx1f(i);
            Real x1r = pcoord->x1f(i)+0.75*pcoord->dx1f(i);
            a1(k,j,i) = 0.5*(A1(x1l, pcoord->x2f(j), pcoord->x3f(k)) +
                             A1(x1r, pcoord->x2f(j), pcoord->x3f(k)));
          } else {
            a1(k,j,i) = A1(pcoord->x1v(i), pcoord->x2f(j), pcoord->x3f(k));
          }

          if ((pbval->nblevel[1][1][0]>level && i==is)
              || (pbval->nblevel[1][1][2]>level && i==ie+1)
              || (pbval->nblevel[0][1][1]>level && k==ks)
              || (pbval->nblevel[2][1][1]>level && k==ke+1)
              || (pbval->nblevel[0][1][0]>level && i==is   && k==ks)
              || (pbval->nblevel[0][1][2]>level && i==ie+1 && k==ks)
              || (pbval->nblevel[2][1][0]>level && i==is   && k==ke+1)
              || (pbval->nblevel[2][1][2]>level && i==ie+1 && k==ke+1)) {
            Real x2l = pcoord->x2f(j)+0.25*pcoord->dx2f(j);
            Real x2r = pcoord->x2f(j)+0.75*pcoord->dx2f(j);
            a2(k,j,i) = 0.5*(A2(pcoord->x1f(i), x2l, pcoord->x3f(k)) +
                             A2(pcoord->x1f(i), x2r, pcoord->x3f(k)));
          } else {
            a2(k,j,i) = A2(pcoord->x1f(i), pcoord->x2v(j), pcoord->x3f(k));
          }

          if ((pbval->nblevel[1][1][0]>level && i==is)
              || (pbval->nblevel[1][1][2]>level && i==ie+1)
              || (pbval->nblevel[1][0][1]>level && j==js)
              || (pbval->nblevel[1][2][1]>level && j==je+1)
              || (pbval->nblevel[1][0][0]>level && i==is   && j==js)
              || (pbval->nblevel[1][0][2]>level && i==ie+1 && j==js)
              || (pbval->nblevel[1][2][0]>level && i==is   && j==je+1)
              || (pbval->nblevel[1][2][2]>level && i==ie+1 && j==je+1)) {
            Real x3l = pcoord->x3f(k)+0.25*pcoord->dx3f(k);
            Real x3r = pcoord->x3f(k)+0.75*pcoord->dx3f(k);
            a3(k,j,i) = 0.5*(A3(pcoord->x1f(i), pcoord->x2f(j), x3l) +
                             A3(pcoord->x1f(i), pcoord->x2f(j), x3r));
          } else {
            a3(k,j,i) = A3(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3v(k));
          }
        }
      }
    }

    // Initialize interface fields
    AthenaArray<Real> area, len, len_p1;
    area.NewAthenaArray(ncells1);
    len.NewAthenaArray(ncells1);
    len_p1.NewAthenaArray(ncells1);

    // for 1,2,3-D
    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary
      int jl=js; int ju=je+1;
      if (pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) jl=js+1;
      if (pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) ju=je;
      for (int j=jl; j<=ju; ++j) {
        pcoord->Face2Area(k,j,is,ie,area);
        pcoord->Edge3Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = -1.0*(len(i+1)*a3(k,j,i+1) - len(i)*a3(k,j,i))/area(i);
        }
      }
    }
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        pcoord->Face3Area(k,j,is,ie,area);
        pcoord->Edge2Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = (len(i+1)*a2(k,j,i+1) - len(i)*a2(k,j,i))/area(i);
        }
      }
    }

    // for 2D and 3D
    if (block_size.nx2 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face1Area(k,j,is,ie+1,area);
          pcoord->Edge3Length(k,j  ,is,ie+1,len);
          pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) = (len_p1(i)*a3(k,j+1,i) - len(i)*a3(k,j,i))/area(i);
          }
        }
      }
      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face3Area(k,j,is,ie,area);
          pcoord->Edge1Length(k,j  ,is,ie,len);
          pcoord->Edge1Length(k,j+1,is,ie,len_p1);
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) -= (len_p1(i)*a1(k,j+1,i) - len(i)*a1(k,j,i))/area(i);
          }
        }
      }
    }
    // for 3D only
    if (block_size.nx3 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          pcoord->Face1Area(k,j,is,ie+1,area);
          pcoord->Edge2Length(k  ,j,is,ie+1,len);
          pcoord->Edge2Length(k+1,j,is,ie+1,len_p1);
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) -= (len_p1(i)*a2(k+1,j,i) - len(i)*a2(k,j,i))/area(i);
          }
        }
      }
      for (int k=ks; k<=ke; ++k) {
        // reset loop limits for polar boundary
        int jl=js; int ju=je+1;
        if (pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar) jl=js+1;
        if (pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar) ju=je;
        for (int j=jl; j<=ju; ++j) {
          pcoord->Face2Area(k,j,is,ie,area);
          pcoord->Edge1Length(k  ,j,is,ie,len);
          pcoord->Edge1Length(k+1,j,is,ie,len_p1);
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) += (len_p1(i)*a1(k+1,j,i) - len(i)*a1(k,j,i))/area(i);
          }
        }
      }
    }
  }

  //  Initialize density
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = rho0;
        VelProfileCyl(pcoord->x1v(i),pcoord->x2v(j),pcoord->x3v(k),v1,v2,v3);
        phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*v1;
        phydro->u(IM2,k,j,i) = phydro->u(IDN,k,j,i)*v2;
        phydro->u(IM3,k,j,i) = phydro->u(IDN,k,j,i)*v3;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = isocs2*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))
                                       +SQR(phydro->u(IM2,k,j,i))
                                       +SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
          if (MAGNETIC_FIELDS_ENABLED) {
            phydro->u(IEN,k,j,i) +=
                0.5*(SQR(0.5*(pfield->b.x1f(k,j,i+1) + pfield->b.x1f(k,j,i)))
                     + SQR(0.5*(pfield->b.x2f(k,j+1,i) + pfield->b.x2f(k,j,i)))
                     + SQR(0.5*(pfield->b.x3f(k+1,j,i) + pfield->b.x3f(k,j,i))));
          }
        }
      }
    }
  }

  return;
}

namespace {
//----------------------------------------------------------------------------------------
//! \f transforms uniform velocity in x-directin in spherical polar coords

void VelProfileCyl(const Real x1, const Real x2, const Real x3,
                   Real &v1, Real &v2, Real &v3) {
  v1 = vy0*std::sin(x2)*std::sin(x3);
  v2 = vy0*std::cos(x2)*std::sin(x3);
  v3 = vy0*std::cos(x3);
  return;
}

//----------------------------------------------------------------------------------------
//! \f compute 3-compnent of vector potential

Real A3(const Real x1, const Real x2, const Real x3) {
  Real a3=0.0;
  return a3;
}

//----------------------------------------------------------------------------------------
//! \f compute 2-compnent of vector potential

Real A2(const Real x1, const Real x2, const Real x3) {
  Real a2=0.0;
  Real az=0.0;
  Real x=x1*std::abs(std::sin(x2))*std::cos(x3);
  Real y=x1*std::abs(std::sin(x2))*std::sin(x3);
  if (x2<0.0||x2>PI) {
    x=-x;
    y=-y;
  }
  Real z=x1*std::cos(x2);
  if (std::sqrt(SQR(x-xc)+SQR(y-yc))<=0.5 && std::abs(z-zc)<0.2) {
    az=b0*(0.5-std::sqrt(SQR(x-xc)+SQR(y-yc)));
  }
  a2=-az*std::abs(std::sin(x2));
  return a2;
}

//----------------------------------------------------------------------------------------
//! \f compute 1-compnent of vector potential

Real A1(const Real x1, const Real x2, const Real x3) {
  Real a1=0.0;
  Real az=0.0;
  Real x=x1*std::abs(std::sin(x2))*std::cos(x3);
  Real y=x1*std::abs(std::sin(x2))*std::sin(x3);
  if (x2<0.0||x2>PI) {
    x=-x;
    y=-y;
  }
  Real z=x1*std::cos(x2);
  if (std::sqrt(SQR(x-xc)+SQR(y-yc))<=0.5 && std::abs(z-zc)<0.2) {
    az=b0*(0.5-std::sqrt(SQR(x-xc)+SQR(y-yc)));
  }
  a1=az*std::cos(x2);
  return a1;
}
} // namespace

//----------------------------------------------------------------------------------------
//! \brief: User-defined boundary Conditions: LoopInnerX1

void LoopInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // Real rad, phi, z;
  Real v1, v2, v3;
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        prim(IDN,k,j,il-i) = rho0;
        VelProfileCyl(pco->x1v(il-i),pco->x2v(j),pco->x3v(k),v1,v2,v3);
        prim(IM1,k,j,il-i) = v1;
        prim(IM2,k,j,il-i) = v2;
        prim(IM3,k,j,il-i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,j,il-i) = isocs2*prim(IDN,k,j,il-i);
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(il-i)) = b.x1f(k,j,il);
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(il-i)) = b.x2f(k,j,il);
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(il-i)) = b.x3f(k,j,il);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \brief: User-defined boundary Conditions: LoopOuterX1

void LoopOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        prim(IDN,k,j,iu+i) = rho0;
        VelProfileCyl(pco->x1v(iu+i),pco->x2v(j),pco->x3v(k),v1,v2,v3);
        prim(IM1,k,j,iu+i) = v1;
        prim(IM2,k,j,iu+i) = v2;
        prim(IM3,k,j,iu+i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,j,iu+i) = isocs2*prim(IDN,k,j,iu+i);
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x1f(k,j,(iu+i+1)) = b.x1f(k,j,(iu+1));
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x2f(k,j,(iu+i)) = b.x2f(k,j,iu);
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=1; i<=ngh; ++i) {
          b.x3f(k,j,(iu+i)) = b.x3f(k,j,iu);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \brief: User-defined boundary Conditions: LoopInnerX2

void LoopInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        prim(IDN,k,jl-j,i) = rho0;
        VelProfileCyl(pco->x1v(i),pco->x2v(jl-j),pco->x3v(k),v1,v2,v3);
        prim(IM1,k,jl-j,i) = v1;
        prim(IM2,k,jl-j,i) = v2;
        prim(IM3,k,jl-j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,jl-j,i) = isocs2*prim(IDN,k,jl-j,i);
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          b.x1f(k,(jl-j),i) = b.x1f(k,jl,i);
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(jl-j),i) = b.x2f(k,jl,i);
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(jl-j),i) = b.x3f(k,jl,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \brief: User-defined boundary Conditions: LoopOuterX2

void LoopOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  //  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        prim(IDN,k,ju+j,i) = rho0;
        VelProfileCyl(pco->x1v(i),pco->x2v(ju+j),pco->x3v(k),v1,v2,v3);
        prim(IM1,k,ju+j,i) = v1;
        prim(IM2,k,ju+j,i) = v2;
        prim(IM3,k,ju+j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,ju+j,i) = isocs2*prim(IDN,k,ju+j,i);
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu+1; ++i) {
          b.x1f(k,(ju+j  ),i) = b.x1f(k,(ju  ),i);
        }
      }
    }

    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x2f(k,(ju+j+1),i) = b.x2f(k,(ju+1),i);
        }
      }
    }

    for (int k=kl; k<=ku+1; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          b.x3f(k,(ju+j  ),i) = b.x3f(k,(ju  ),i);
        }
      }
    }
  }
}


// refinement condition: check the field amplitude
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &bc = pmb->pfield->bcc;
  Real maxb = 0.0;
  for (int k=pmb->ks-1; k<=pmb->ke+1; k++) {
    for (int j=pmb->js-1; j<=pmb->je+1; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        Real b = std::sqrt(SQR(bc(IB1,k,j,i)) + SQR(bc(IB2,k,j,i))+SQR(bc(IB3,k,j,i)));
        maxb = std::max(maxb, b);
      }
    }
  }

  if (maxb > 0.2*b0) return 1;
  if (maxb < 0.1*b0) return -1;
  return 0;
}
