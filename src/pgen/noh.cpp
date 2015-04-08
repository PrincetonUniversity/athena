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

// Primary header
#include "../mesh.hpp"

// Athena headers
#include "../athena.hpp"           // enums, Real
#include "../athena_arrays.hpp"    // AthenaArray
#include "../parameter_input.hpp"  // ParameterInput
#include "../fluid/fluid.hpp"      // Fluid
#include "../fluid/eos/eos.hpp"    // GetGamma
#include "../bvals/bvals.hpp"      // Boundary Enroll

// BCs on outer edges of grid in each dimension
void noh3d_oib(MeshBlock *pmb, AthenaArray<Real> &a,
               int is, int ie, int js, int je, int ks, int ke);
void noh3d_ojb(MeshBlock *pmb, AthenaArray<Real> &a,
               int is, int ie, int js, int je, int ks, int ke);
void noh3d_okb(MeshBlock *pmb, AthenaArray<Real> &a,
               int is, int ie, int js, int je, int ks, int ke);

// made global to share with BC functions
static Real gmma, gmma1;

//======================================================================================
//! \file noh.c
//  \brief Spherical Noh implosion problem, from Liska & Wendroff, section 4.5 (fig 4.7)
//
//  Tests code on VERY strong shock, also sensitive to carbuncle instability.
//
// PRIVATE FUNCTION PROTOTYPES:
//   - void noh3d_oib() - sets BCs on R-x1 boundary
//   - void noh3d_ojb() - sets BCs on R-x2 boundary
//   - void noh3d_okb() - sets BCs on R-x3 boundary
//   - void scat_plot() - makes scatter plot of density
//
// REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)
//======================================================================================

void Mesh::ProblemGenerator(Fluid *pfl, Field *pfd, ParameterInput *pin)
{
  MeshBlock *pmb = pfl->pmy_block;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  gmma  = pfl->pf_eos->GetGamma();
  gmma1 = gmma - 1.0;

// Initialize the grid: d=1, v=-1.0 in radial direction, p=10^-6

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
  for (int i=is; i<=ie; i++) {
    Real rad;
    if (pmb->block_size.nx3 > 1) {
      rad = sqrt(SQR(pmb->x1v(i)) + SQR(pmb->x2v(j)) + SQR(pmb->x3v(k)));
      pfl->u(IM3,k,j,i) = -pmb->x3v(k)/rad;
    } else {
      rad = sqrt(SQR(pmb->x1v(i)) + SQR(pmb->x2v(j)));
      pfl->u(IM3,k,j,i) = 0.0;
    }
    pfl->u(IDN,k,j,i) = 1.0;
    pfl->u(IM1,k,j,i) = -pmb->x1v(i)/rad;
    pfl->u(IM2,k,j,i) = -pmb->x2v(j)/rad;
    pfl->u(IEN,k,j,i) = 1.0e-6/gmma1 + 0.5;
  }}}

// Enroll boundary value function pointers
  pmb->pbval->EnrollFluidBoundaryFunction(outer_x1, noh3d_oib);
  pmb->pbval->EnrollFluidBoundaryFunction(outer_x2, noh3d_ojb);
  if (pmb->block_size.nx3 > 1) {
    pmb->pbval->EnrollFluidBoundaryFunction(outer_x3, noh3d_okb);
  }

}

//--------------------------------------------------------------------------------------
//! \fn void noh3d_oib()
//  \brief Sets boundary condition on right X1 boundary (oib) for noh3d test
//
// Quantities at this boundary are held fixed at the time-dependent upstream state

void noh3d_oib(MeshBlock *pmb, AthenaArray<Real> &a,
               int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
    for (int i=1;  i<=(NGHOST); ++i) {
      Real rad,f_t;
      if (pmb->block_size.nx3 > 1) {
        rad = sqrt(SQR(pmb->x1v(ie+i)) + SQR(pmb->x2v(j)) + SQR(pmb->x3v(k)));
        f_t = SQR(1.0 + pmb->pmy_mesh->time/rad);
      } else {
        rad = sqrt(SQR(pmb->x1v(ie+i)) + SQR(pmb->x2v(j)));
        f_t = (1.0 + pmb->pmy_mesh->time/rad);
      }
      Real d0 = 1.0*f_t;
   
      a(IDN,k,j,ie+i)  = d0;
      a(IM1,k,j,ie+i) = -pmb->x1v(ie+i)*d0/rad;
      a(IM2,k,j,ie+i) = -pmb->x2v(j   )*d0/rad;
      if (pmb->block_size.nx3 > 1) {
        a(IM3,k,j,ie+i) = -pmb->x3v(k)*d0/rad;
        a(IEN,k,j,ie+i) = 1.0e-6*pow(f_t,(1.0+gmma))/gmma1 + 0.5*d0;
      } else {
        a(IM3,k,j,ie+i) = 0.0;
        a(IEN,k,j,ie+i)= 1.0e-6/gmma1 + 0.5*d0;
      }
    }
  }}
}

//--------------------------------------------------------------------------------------
//! \fn void noh3d_oib()
//  \brief Sets boundary condition on right X2 boundary (ojb) for noh3d test
//
// Quantities at this boundary are held fixed at the time-dependent upstream state

void noh3d_ojb(MeshBlock *pmb, AthenaArray<Real> &a,
               int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
    for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
      Real rad,f_t;
      if (pmb->block_size.nx3 > 1) {
        rad = sqrt(SQR(pmb->x1v(i)) + SQR(pmb->x2v(je+j)) + SQR(pmb->x3v(k)));
        f_t = SQR(1.0 + pmb->pmy_mesh->time/rad);
      } else {
        rad = sqrt(SQR(pmb->x1v(i)) + SQR(pmb->x2v(je+j)));
        f_t = (1.0 + pmb->pmy_mesh->time/rad);
      }
      Real d0 = 1.0*f_t;

      a(IDN,k,je+j,i)  = d0;
      a(IM1,k,je+j,i) = -pmb->x1v(i)*d0/rad;
      a(IM2,k,je+j,i) = -pmb->x2v(je+j)*d0/rad;
      if (pmb->block_size.nx3 > 1) {
        a(IM3,k,je+j,i) = -pmb->x3v(k)*d0/rad;
        a(IEN,k,je+j,i) = 1.0e-6*pow(f_t,(1.0+gmma))/gmma1 + 0.5*d0;
      } else {
        a(IM3,k,je+j,i) = 0.0;
        a(IEN,k,je+j,i)= 1.0e-6/gmma1 + 0.5*d0;
      }
    }
  }}
}

//--------------------------------------------------------------------------------------
//! \fn void noh3d_oib()
//  \brief Sets boundary condition on right X3 boundary (okb) for noh3d test
//
// Quantities at this boundary are held fixed at the time-dependent upstream state

void noh3d_okb(MeshBlock *pmb, AthenaArray<Real> &a,
               int is, int ie, int js, int je, int ks, int ke)
{
  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js-(NGHOST); j<=je+(NGHOST); ++j) {
    for (int i=is-(NGHOST); i<=ie+(NGHOST); ++i) {
      Real rad = sqrt(SQR(pmb->x1v(i)) + SQR(pmb->x2v(j)) + SQR(pmb->x3v(ke+k)));
      Real f_t = SQR(1.0 + pmb->pmy_mesh->time/rad);
      Real d0 = 1.0*f_t;

      a(IDN,ke+k,j,i)  = d0;
      a(IM1,ke+k,j,i) = -pmb->x1v(i)*d0/rad;
      a(IM2,ke+k,j,i) = -pmb->x2v(j)*d0/rad;
      a(IM3,ke+k,j,i) = -pmb->x3v(ke+k)*d0/rad;
      a(IEN,ke+k,j,i) = 1.0e-6*pow(f_t,(1.0+gmma))/gmma1 + 0.5*d0;
    }
  }}
}
