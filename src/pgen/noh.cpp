//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file noh.c
//  \brief Spherical Noh implosion problem, from Liska & Wendroff, section 4.5 (fig 4.7)
//
//  Tests code on VERY strong shock, also sensitive to carbuncle instability.
// REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../hydro/eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// BCs on outer edges of grid in each dimension
void Noh3DOuterX1(MeshBlock *pmb, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke);
void Noh3DOuterX2(MeshBlock *pmb, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke);
void Noh3DOuterX3(MeshBlock *pmb, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke);

// made global to share with BC functions
static Real gmma, gmma1;


//======================================================================================
//! \fn void Mesh::InitUserMeshProperties(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshProperties(ParameterInput *pin)
{
// Enroll boundary value function pointers
  EnrollUserBoundaryFunction(OUTER_X1, Noh3DOuterX1);
  EnrollUserBoundaryFunction(OUTER_X2, Noh3DOuterX2);
  if (mesh_size.nx3 > 1)
    EnrollUserBoundaryFunction(OUTER_X3, Noh3DOuterX3);
  return;
}


//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================

void Mesh::TerminateUserMeshProperties(void)
{
  // nothing to do
  return;
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Noh spherical implosion test
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  gmma  = phydro->peos->GetGamma();
  gmma1 = gmma - 1.0;

// Initialize the grid: d=1, v=-1.0 in radial direction, p=10^-6

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    Real rad;
    if (block_size.nx3 > 1) {
      rad = sqrt(SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j)) + SQR(pcoord->x3v(k)));
      phydro->u(IM3,k,j,i) = -pcoord->x3v(k)/rad;
    } else {
      rad = sqrt(SQR(pcoord->x1v(i)) + SQR(pcoord->x2v(j)));
      phydro->u(IM3,k,j,i) = 0.0;
    }
    phydro->u(IDN,k,j,i) = 1.0;
    phydro->u(IM1,k,j,i) = -pcoord->x1v(i)/rad;
    phydro->u(IM2,k,j,i) = -pcoord->x2v(j)/rad;
    phydro->u(IEN,k,j,i) = 1.0e-6/gmma1 + 0.5;
  }}}
}


//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
//======================================================================================

void MeshBlock::UserWorkInLoop(void)
{
  // nothing to do
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void Noh3DOuterX1()
//  \brief Sets boundary condition on right X1 boundary (oib) for noh3d test
//
// Quantities at this boundary are held fixed at the time-dependent upstream state

void Noh3DOuterX1(MeshBlock *pmb, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=1;  i<=(NGHOST); ++i) {
      Real rad,f_t;
      if (pmb->block_size.nx3 > 1) {
        rad = sqrt(SQR(pmb->pcoord->x1v(ie+i)) + SQR(pmb->pcoord->x2v(j)) 
              + SQR(pmb->pcoord->x3v(k)));
        f_t = SQR(1.0 + pmb->pmy_mesh->time/rad);
      } else {
        rad = sqrt(SQR(pmb->pcoord->x1v(ie+i)) + SQR(pmb->pcoord->x2v(j)));
        f_t = (1.0 + pmb->pmy_mesh->time/rad);
      }
      Real d0 = 1.0*f_t;
   
      a(IDN,k,j,ie+i)  = d0;
      a(IVX,k,j,ie+i) = -pmb->pcoord->x1v(ie+i)/rad;
      a(IVY,k,j,ie+i) = -pmb->pcoord->x2v(j   )/rad;
      if (pmb->block_size.nx3 > 1) {
        a(IVZ,k,j,ie+i) = -pmb->pcoord->x3v(k)/rad;
        a(IEN,k,j,ie+i) = 1.0e-6*pow(f_t,(1.0+gmma));
      } else {
        a(IVZ,k,j,ie+i) = 0.0;
        a(IEN,k,j,ie+i)= 1.0e-6;
      }
    }
  }}
}

//--------------------------------------------------------------------------------------
//! \fn void Noh3DOuterX2()
//  \brief Sets boundary condition on right X2 boundary (ojb) for noh3d test
//
// Quantities at this boundary are held fixed at the time-dependent upstream state

void Noh3DOuterX2(MeshBlock *pmb, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
    for (int i=is; i<=ie; ++i) {
      Real rad,f_t;
      if (pmb->block_size.nx3 > 1) {
        rad = sqrt(SQR(pmb->pcoord->x1v(i)) + SQR(pmb->pcoord->x2v(je+j)) 
              + SQR(pmb->pcoord->x3v(k)));
        f_t = SQR(1.0 + pmb->pmy_mesh->time/rad);
      } else {
        rad = sqrt(SQR(pmb->pcoord->x1v(i)) + SQR(pmb->pcoord->x2v(je+j)));
        f_t = (1.0 + pmb->pmy_mesh->time/rad);
      }
      Real d0 = 1.0*f_t;

      a(IDN,k,je+j,i)  = d0;
      a(IVX,k,je+j,i) = -pmb->pcoord->x1v(i)/rad;
      a(IVY,k,je+j,i) = -pmb->pcoord->x2v(je+j)/rad;
      if (pmb->block_size.nx3 > 1) {
        a(IVZ,k,je+j,i) = -pmb->pcoord->x3v(k)/rad;
        a(IEN,k,je+j,i) = 1.0e-6*pow(f_t,(1.0+gmma));
      } else {
        a(IVZ,k,je+j,i) = 0.0;
        a(IEN,k,je+j,i)= 1.0e-6;
      }
    }
  }}
}

//--------------------------------------------------------------------------------------
//! \fn void Noh3DOuterX3()
//  \brief Sets boundary condition on right X3 boundary (okb) for noh3d test
//
// Quantities at this boundary are held fixed at the time-dependent upstream state

void Noh3DOuterX3(MeshBlock *pmb, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      Real rad = sqrt(SQR(pmb->pcoord->x1v(i)) + SQR(pmb->pcoord->x2v(j)) 
              + SQR(pmb->pcoord->x3v(ke+k)));
      Real f_t = SQR(1.0 + pmb->pmy_mesh->time/rad);
      Real d0 = 1.0*f_t;

      a(IDN,ke+k,j,i)  = d0;
      a(IVX,ke+k,j,i) = -pmb->pcoord->x1v(i)/rad;
      a(IVY,ke+k,j,i) = -pmb->pcoord->x2v(j)/rad;
      a(IVZ,ke+k,j,i) = -pmb->pcoord->x3v(ke+k)/rad;
      a(IEN,ke+k,j,i) = 1.0e-6*pow(f_t,(1.0+gmma));
    }
  }}
}
