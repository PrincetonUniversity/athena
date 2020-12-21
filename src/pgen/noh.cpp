//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file noh.cpp
//! \brief Spherical Noh implosion problem, from Liska & Wendroff, section 4.5 (fig 4.7)
//!
//! Tests code on VERY strong shock, also sensitive to carbuncle instability.
//! REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// BCs on outer edges of grid in each dimension
void Noh3DOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void Noh3DOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void Noh3DOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh);

// made global to share with BC functions
namespace {
Real gmma, gmma1;
} // namespace

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Enroll boundary value function pointers
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Noh3DOuterX1);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x2, Noh3DOuterX2);
  if (mesh_size.nx3 > 1)
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, Noh3DOuterX3);
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
          rad = std::sqrt(SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j))
                          + SQR(pcoord->x3v(k)));
          phydro->u(IM3,k,j,i) = -pcoord->x3v(k)/rad;
        } else {
          rad = std::sqrt(SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j)));
          phydro->u(IM3,k,j,i) = 0.0;
        }
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = -pcoord->x1v(i)/rad;
        phydro->u(IM2,k,j,i) = -pcoord->x2v(j)/rad;
        phydro->u(IEN,k,j,i) = 1.0e-6/gmma1 + 0.5;
      }
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void Noh3DOuterX1()
//  \brief Sets boundary condition on right X1 boundary (oib) for noh3d test
//
// Quantities at this boundary are held fixed at the time-dependent upstream state

void Noh3DOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1;  i<=ngh; ++i) {
        Real rad,f_tm;
        if (pmb->block_size.nx3 > 1) {
          rad = std::sqrt(SQR(pco->x1v(iu+i)) + SQR(pco->x2v(j))
                          + SQR(pco->x3v(k)));
          f_tm = SQR(1.0 + pmb->pmy_mesh->time/rad);
        } else {
          rad = std::sqrt(SQR(pco->x1v(iu+i)) + SQR(pco->x2v(j)));
          f_tm = (1.0 + pmb->pmy_mesh->time/rad);
        }
        Real d0 = 1.0*f_tm;

        prim(IDN,k,j,iu+i)  = d0;
        prim(IVX,k,j,iu+i) = -pco->x1v(iu+i)/rad;
        prim(IVY,k,j,iu+i) = -pco->x2v(j   )/rad;
        if (pmb->block_size.nx3 > 1) {
          prim(IVZ,k,j,iu+i) = -pco->x3v(k)/rad;
          prim(IPR,k,j,iu+i) = 1.0e-6*std::pow(f_tm,(1.0+gmma));
        } else {
          prim(IVZ,k,j,iu+i) = 0.0;
          prim(IPR,k,j,iu+i)= 1.0e-6;
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Noh3DOuterX2()
//  \brief Sets boundary condition on right X2 boundary (ojb) for noh3d test
//
// Quantities at this boundary are held fixed at the time-dependent upstream state

void Noh3DOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        Real rad,f_tm;
        if (pmb->block_size.nx3 > 1) {
          rad = std::sqrt(SQR(pco->x1v(i)) + SQR(pco->x2v(ju+j))
                          + SQR(pco->x3v(k)));
          f_tm = SQR(1.0 + pmb->pmy_mesh->time/rad);
        } else {
          rad = std::sqrt(SQR(pco->x1v(i)) + SQR(pco->x2v(ju+j)));
          f_tm = (1.0 + pmb->pmy_mesh->time/rad);
        }
        Real d0 = 1.0*f_tm;

        prim(IDN,k,ju+j,i)  = d0;
        prim(IVX,k,ju+j,i) = -pco->x1v(i)/rad;
        prim(IVY,k,ju+j,i) = -pco->x2v(ju+j)/rad;
        if (pmb->block_size.nx3 > 1) {
          prim(IVZ,k,ju+j,i) = -pco->x3v(k)/rad;
          prim(IPR,k,ju+j,i) = 1.0e-6*std::pow(f_tm,(1.0+gmma));
        } else {
          prim(IVZ,k,ju+j,i) = 0.0;
          prim(IPR,k,ju+j,i)= 1.0e-6;
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Noh3DOuterX3()
//  \brief Sets boundary condition on right X3 boundary (okb) for noh3d test
//
// Quantities at this boundary are held fixed at the time-dependent upstream state

void Noh3DOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        Real rad = std::sqrt(SQR(pco->x1v(i)) + SQR(pco->x2v(j))
                             + SQR(pco->x3v(ku+k)));
        Real f_tm = SQR(1.0 + pmb->pmy_mesh->time/rad);
        Real d0 = 1.0*f_tm;

        prim(IDN,ku+k,j,i)  = d0;
        prim(IVX,ku+k,j,i) = -pco->x1v(i)/rad;
        prim(IVY,ku+k,j,i) = -pco->x2v(j)/rad;
        prim(IVZ,ku+k,j,i) = -pco->x3v(ku+k)/rad;
        prim(IPR,ku+k,j,i) = 1.0e-6*std::pow(f_tm,(1.0+gmma));
      }
    }
  }
}
