//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file noh.c
//  \brief Spherical Noh implosion problem, from Liska & Wendroff, section 4.5 (fig 4.7)
//
//  Tests code on VERY strong shock, also sensitive to carbuncle instability.
// REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)
//========================================================================================

// C++ headers
#include <cmath>      // sqrt()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// BCs on outer edges of grid in each dimension
void Noh3DOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void Noh3DOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void Noh3DOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

// made global to share with BC functions
static Real gmma, gmma1;

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Enroll boundary value function pointers
  EnrollUserBoundaryFunction(OUTER_X1, Noh3DOuterX1);
  EnrollUserBoundaryFunction(OUTER_X2, Noh3DOuterX2);
  if (mesh_size.nx3 > 1)
    EnrollUserBoundaryFunction(OUTER_X3, Noh3DOuterX3);
  return;
}


//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Noh spherical implosion test
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  gmma  = peos->GetGamma();
  gmma1 = gmma - 1.0;

  // Initialize the grid: d=1, v=-1.0 in radial direction, p=10^-6
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    Real rad;
    if (block_size.nx3 > 1) {
      rad = std::sqrt(SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j)) + SQR(pcoord->x3v(k)));
      phydro->u(IM3,k,j,i) = -pcoord->x3v(k)/rad;
    } else {
      rad = std::sqrt(SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j)));
      phydro->u(IM3,k,j,i) = 0.0;
    }
    phydro->u(IDN,k,j,i) = 1.0;
    phydro->u(IM1,k,j,i) = -pcoord->x1v(i)/rad;
    phydro->u(IM2,k,j,i) = -pcoord->x2v(j)/rad;
    phydro->u(IEN,k,j,i) = 1.0e-6/gmma1 + 0.5;
  }}}
}


//----------------------------------------------------------------------------------------
//! \fn void Noh3DOuterX1()
//  \brief Sets boundary condition on right X1 boundary (oib) for noh3d test
//
// Quantities at this boundary are held fixed at the time-dependent upstream state

void Noh3DOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=1;  i<=ngh; ++i) {
      Real rad,f_t;
      if (pmb->block_size.nx3 > 1) {
        rad = std::sqrt(SQR(pco->x1v(ie+i)) + SQR(pco->x2v(j))
              + SQR(pco->x3v(k)));
        f_t = SQR(1.0 + pmb->pmy_mesh->time/rad);
      } else {
        rad = std::sqrt(SQR(pco->x1v(ie+i)) + SQR(pco->x2v(j)));
        f_t = (1.0 + pmb->pmy_mesh->time/rad);
      }
      Real d0 = 1.0*f_t;

      prim(IDN,k,j,ie+i)  = d0;
      prim(IVX,k,j,ie+i) = -pco->x1v(ie+i)/rad;
      prim(IVY,k,j,ie+i) = -pco->x2v(j   )/rad;
      if (pmb->block_size.nx3 > 1) {
        prim(IVZ,k,j,ie+i) = -pco->x3v(k)/rad;
        prim(IPR,k,j,ie+i) = 1.0e-6*pow(f_t,(1.0+gmma));
      } else {
        prim(IVZ,k,j,ie+i) = 0.0;
        prim(IPR,k,j,ie+i)= 1.0e-6;
      }
    }
  }}
}

//----------------------------------------------------------------------------------------
//! \fn void Noh3DOuterX2()
//  \brief Sets boundary condition on right X2 boundary (ojb) for noh3d test
//
// Quantities at this boundary are held fixed at the time-dependent upstream state

void Noh3DOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=ngh; ++j) {
    for (int i=is; i<=ie; ++i) {
      Real rad,f_t;
      if (pmb->block_size.nx3 > 1) {
        rad = std::sqrt(SQR(pco->x1v(i)) + SQR(pco->x2v(je+j))
              + SQR(pco->x3v(k)));
        f_t = SQR(1.0 + pmb->pmy_mesh->time/rad);
      } else {
        rad = std::sqrt(SQR(pco->x1v(i)) + SQR(pco->x2v(je+j)));
        f_t = (1.0 + pmb->pmy_mesh->time/rad);
      }
      Real d0 = 1.0*f_t;

      prim(IDN,k,je+j,i)  = d0;
      prim(IVX,k,je+j,i) = -pco->x1v(i)/rad;
      prim(IVY,k,je+j,i) = -pco->x2v(je+j)/rad;
      if (pmb->block_size.nx3 > 1) {
        prim(IVZ,k,je+j,i) = -pco->x3v(k)/rad;
        prim(IPR,k,je+j,i) = 1.0e-6*pow(f_t,(1.0+gmma));
      } else {
        prim(IVZ,k,je+j,i) = 0.0;
        prim(IPR,k,je+j,i)= 1.0e-6;
      }
    }
  }}
}

//----------------------------------------------------------------------------------------
//! \fn void Noh3DOuterX3()
//  \brief Sets boundary condition on right X3 boundary (okb) for noh3d test
//
// Quantities at this boundary are held fixed at the time-dependent upstream state

void Noh3DOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
       Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int k=1; k<=ngh; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      Real rad = std::sqrt(SQR(pco->x1v(i)) + SQR(pco->x2v(j))
              + SQR(pco->x3v(ke+k)));
      Real f_t = SQR(1.0 + pmb->pmy_mesh->time/rad);
      Real d0 = 1.0*f_t;

      prim(IDN,ke+k,j,i)  = d0;
      prim(IVX,ke+k,j,i) = -pco->x1v(i)/rad;
      prim(IVY,ke+k,j,i) = -pco->x2v(j)/rad;
      prim(IVZ,ke+k,j,i) = -pco->x3v(ke+k)/rad;
      prim(IPR,ke+k,j,i) = 1.0e-6*pow(f_t,(1.0+gmma));
    }
  }}
}
