//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file field_loop_poles.c
//  \brief Advection of a field loop THROUGH the poles in spherical_polar coordinates.
//
//  Originally developed by ZZ.  Sets up constant uniform-density flow in x-direction
//  through poles, and follows advection of loop.  Set xz>0 (xz<0) for loop through
//  upper (lower) pole.  Works in 2D and 3D.
//========================================================================================

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"

static Real DenProfileCyl(const Real rad, const Real phi, const Real z);
static void VelProfileCyl(const Real rad, const Real phi, const Real z,
                          Real &v1, Real &v2, Real &v3);
static Real A3(const Real x1, const Real x2, const Real x3);
static Real A2(const Real x1, const Real x2, const Real x3);
static Real A1(const Real x1, const Real x2, const Real x3);

// User-defined boundary conditions along inner/outer edges (not poles)
void LoopInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void LoopOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void LoopInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void LoopOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

// problem parameters which are useful to make global to this file
static Real vy0, rho0, isocs2, gamma_gas;
static Real xc, yc, zc, beta, b0;

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
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

  // setup boundary condition
  if (mesh_bcs[INNER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(INNER_X1, LoopInnerX1);
  }
  if (mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(OUTER_X1, LoopOuterX1);
  }
  if (mesh_bcs[INNER_X2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(INNER_X2, LoopInnerX2);
  }
  if (mesh_bcs[OUTER_X2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(OUTER_X2, LoopOuterX2);
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Initializes field loop advection through pole.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad, phi, z;
  Real v1, v2, v3;
  // Set initial magnetic fields
  if (MAGNETIC_FIELDS_ENABLED) {
    AthenaArray<Real> a1,a2,a3;
    int nx1 = (ie-is)+1 + 2*(NGHOST);
    int nx2 = (je-js)+1 + 2*(NGHOST);
    int nx3 = (ke-ks)+1 + 2*(NGHOST);
    a1.NewAthenaArray(nx3,nx2,nx1);
    a2.NewAthenaArray(nx3,nx2,nx1);
    a3.NewAthenaArray(nx3,nx2,nx1);

    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          a1(k,j,i) = A1(pcoord->x1v(i), pcoord->x2f(j), pcoord->x3f(k));
          a2(k,j,i) = A2(pcoord->x1f(i), pcoord->x2v(j), pcoord->x3f(k));
          a3(k,j,i) = A3(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3v(k));
        }
      }
    }

    // Initialize interface fields
    AthenaArray<Real> area,len,len_p1;
    area.NewAthenaArray(nx1);
    len.NewAthenaArray(nx1);
    len_p1.NewAthenaArray(nx1);

    // for 1,2,3-D
    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary
      int jl=js; int ju=je+1;
      if (pbval->block_bcs[INNER_X2] == 5) jl=js+1;
      if (pbval->block_bcs[OUTER_X2] == 5) ju=je;
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
        if (pbval->block_bcs[INNER_X2] == 5) jl=js+1;
        if (pbval->block_bcs[OUTER_X2] == 5) ju=je;
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

    a1.DeleteAthenaArray();
    a2.DeleteAthenaArray();
    a3.DeleteAthenaArray();
    area.DeleteAthenaArray();
    len.DeleteAthenaArray();
    len_p1.DeleteAthenaArray();
  }

  //  Initialize density
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = rho0 ;
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

//----------------------------------------------------------------------------------------
//! \f transforms uniform velocity in x-directin in spherical polar coords

static void VelProfileCyl(const Real x1, const Real x2, const Real x3,
                          Real &v1, Real &v2, Real &v3) {
  v1 = vy0*sin(x2)*sin(x3);
  v2 = vy0*cos(x2)*sin(x3);
  v3 = vy0*cos(x3);
  return;
}

//----------------------------------------------------------------------------------------
//! \f compute 3-compnent of vector potential

static Real A3(const Real x1, const Real x2, const Real x3) {
  Real a3=0.0;
  return a3;
}

//----------------------------------------------------------------------------------------
//! \f compute 2-compnent of vector potential

static Real A2(const Real x1, const Real x2, const Real x3) {
  Real a2=0.0;
  Real az=0.0;
  Real x=x1*fabs(sin(x2))*cos(x3);
  Real y=x1*fabs(sin(x2))*sin(x3);
  if (x2<0.0||x2>PI) {
   x=-x;
   y=-y;
  }
  Real z=x1*cos(x2);
  if (std::sqrt(SQR(x-xc)+SQR(y-yc))<=0.5 && fabs(z-zc)<0.2) {
    az=b0*(0.5-std::sqrt(SQR(x-xc)+SQR(y-yc)));
  }
  a2=-az*fabs(sin(x2));
  return a2;
}

//----------------------------------------------------------------------------------------
//! \f compute 1-compnent of vector potential

static Real A1(const Real x1, const Real x2, const Real x3) {
  Real a1=0.0;
  Real az=0.0;
  Real x=x1*fabs(sin(x2))*cos(x3);
  Real y=x1*fabs(sin(x2))*sin(x3);
  if (x2<0.0||x2>PI) {
   x=-x;
   y=-y;
  }
  Real z=x1*cos(x2);
  if (std::sqrt(SQR(x-xc)+SQR(y-yc))<=0.5 && fabs(z-zc)<0.2) {
    az=b0*(0.5-std::sqrt(SQR(x-xc)+SQR(y-yc)));
  }
  a1=az*cos(x2);
  return a1;
}

//----------------------------------------------------------------------------------------
//!\f: User-defined boundary Conditions: LoopInnerX1

void LoopInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        prim(IDN,k,j,is-i) = rho0;
        VelProfileCyl(pco->x1v(is-i),pco->x2v(j),pco->x3v(k),v1,v2,v3);
        prim(IM1,k,j,is-i) = v1;
        prim(IM2,k,j,is-i) = v2;
        prim(IM3,k,j,is-i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,j,is-i) = isocs2*prim(IDN,k,j,is-i);
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x2f(k,j,(is-i)) = b.x2f(k,j,is);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
      }
    }}
  }
}

//----------------------------------------------------------------------------------------
//!\f: User-defined boundary Conditions: LoopOuterX1

void LoopOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=ngh; ++i) {
        prim(IDN,k,j,ie+i) = rho0;
        VelProfileCyl(pco->x1v(ie+i),pco->x2v(j),pco->x3v(k),v1,v2,v3);
        prim(IM1,k,j,ie+i) = v1;
        prim(IM2,k,j,ie+i) = v2;
        prim(IM3,k,j,ie+i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,j,ie+i) = isocs2*prim(IDN,k,j,ie+i);
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
      }
    }}
  }
}

//----------------------------------------------------------------------------------------
//!\f: User-defined boundary Conditions: LoopInnerX2

void LoopInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=is; i<=ie; ++i) {
        prim(IDN,k,js-j,i) = rho0;
        VelProfileCyl(pco->x1v(i),pco->x2v(js-j),pco->x3v(k),v1,v2,v3);
        prim(IM1,k,js-j,i) = v1;
        prim(IM2,k,js-j,i) = v2;
        prim(IM3,k,js-j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,js-j,i) = isocs2*prim(IDN,k,js-j,i);
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(js-j),i) = b.x1f(k,js,i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(js-j),i) = b.x2f(k,js,i);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(js-j),i) = b.x3f(k,js,i);
      }
    }}
  }
}

//----------------------------------------------------------------------------------------
//!\f: User-defined boundary Conditions: LoopOuterX2

void LoopOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  Real rad,phi,z;
  Real v1, v2, v3;
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=is; i<=ie; ++i) {
        prim(IDN,k,je+j,i) = rho0;
        VelProfileCyl(pco->x1v(i),pco->x2v(je+j),pco->x3v(k),v1,v2,v3);
        prim(IM1,k,je+j,i) = v1;
        prim(IM2,k,je+j,i) = v2;
        prim(IM3,k,je+j,i) = v3;
        if (NON_BAROTROPIC_EOS)
          prim(IEN,k,je+j,i) = isocs2*prim(IDN,k,je+j,i);
      }
    }
  }
  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(je+j  ),i) = b.x1f(k,(je  ),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(je+j+1),i) = b.x2f(k,(je+1),i);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(je+j  ),i) = b.x3f(k,(je  ),i);
      }
    }}
  }
}
