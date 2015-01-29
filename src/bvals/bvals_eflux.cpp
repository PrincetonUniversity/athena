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
#include "bvals.hpp"

// Athena headers
#include "../athena.hpp"         // macros, Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../mesh.hpp"           // MeshBlock
#include "../field/field.hpp"          // InterfaceField

//======================================================================================
//! \file bvals_eflux.cpp
//  \brief implements defualt boundary conditions for electric fields
//======================================================================================


//--------------------------------------------------------------------------------------
//! \fn void DefaultEFluxInnerX1(MeshBlock *pmb, InterfaceField &ei, InterfaceField &w)
//  \brief  Default boundary conditions for electric fields for inner x1 boundary 
void DefaultEFluxInnerX1(MeshBlock *pmb, InterfaceField &ei, InterfaceField &w)
{
  int is = pmb->is;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  if(pmb->block_size.nx3 > 1) {// 3D
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        ei.x3f(X3E2,k,j,is-1) = ei.x3f(X3E2,k,j,is); 
        w.x3f(k,j,is-1) = w.x3f(k,j,is);
      }
    }
  }
  if(pmb->block_size.nx2 > 1) {// 2D or 3D
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        ei.x2f(X2E3,k,j,is-1) = ei.x2f(X2E3,k,j,is); 
        w.x2f(k,j,is-1) = w.x2f(k,j,is);
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void DefaultEFluxOuterX1(MeshBlock *pmb, InterfaceField &ei, InterfaceField &w)
//  \brief  Default boundary conditions for electric fields for outer x1 boundary 
void DefaultEFluxOuterX1(MeshBlock *pmb, InterfaceField &ei, InterfaceField &w)
{
  int ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  if(pmb->block_size.nx3 > 1) {// 3D
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        ei.x3f(X3E2,k,j,ie+1) = ei.x3f(X3E2,k,j,ie); 
        w.x3f(k,j,ie+1) = w.x3f(k,j,ie);
      }
    }
  }
  if(pmb->block_size.nx2 > 1) {// 2D or 3D
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je+1; ++j) {
        ei.x2f(X2E3,k,j,ie+1) = ei.x2f(X2E3,k,j,ie); 
        w.x2f(k,j,ie+1) = w.x2f(k,j,ie);
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void DefaultEFluxInnerX2(MeshBlock *pmb, InterfaceField &ei, InterfaceField &w)
//  \brief  Default boundary conditions for electric fields for inner x2 boundary 
void DefaultEFluxInnerX2(MeshBlock *pmb, InterfaceField &ei, InterfaceField &w)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js;
  int ks = pmb->ks, ke = pmb->ke;
  if(pmb->block_size.nx3 > 1) {// 3D
    for (int k=ks; k<=ke+1; ++k) {
      for (int i=is; i<=ie; ++i) {
        ei.x3f(X3E1,k,js-1,i) = ei.x3f(X3E1,k,js,i); 
        w.x3f(k,js-1,i) = w.x3f(k,js,i);
      }
    }
  }
  if(pmb->block_size.nx2 > 1) {// 2D or 3D
    for (int k=ks; k<=ke; ++k) {
      for (int i=is; i<=ie+1; ++i) {
        ei.x1f(X1E3,k,js-1,i) = ei.x1f(X1E3,k,js,i); 
        w.x1f(k,js-1,i) = w.x1f(k,js,i);
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void DefaultEFluxOuterX2(MeshBlock *pmb, InterfaceField &ei, InterfaceField &w)
//  \brief  Default boundary conditions for electric fields for outer x2 boundary 
void DefaultEFluxOuterX2(MeshBlock *pmb, InterfaceField &ei, InterfaceField &w)
{
  int is = pmb->is, ie = pmb->ie;
  int je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  if(pmb->block_size.nx3 > 1) {// 3D
    for (int k=ks; k<=ke+1; ++k) {
      for (int i=is; i<=ie; ++i) {
        ei.x3f(X3E1,k,je+1,i) = ei.x3f(X3E1,k,je,i); 
        w.x3f(k,je+1,i) = w.x3f(k,je,i);
      }
    }
  }
  if(pmb->block_size.nx2 > 1) {// 2D or 3D
    for (int k=ks; k<=ke; ++k) {
      for (int i=is; i<=ie+1; ++i) {
        ei.x1f(X1E3,k,je+1,i) = ei.x1f(X1E3,k,je,i); 
        w.x1f(k,je+1,i) = w.x1f(k,je,i);
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void DefaultEFluxInnerX3(MeshBlock *pmb, InterfaceField &ei, InterfaceField &w)
//  \brief  Default boundary conditions for electric fields for inner x3 boundary 
void DefaultEFluxInnerX3(MeshBlock *pmb, InterfaceField &ei, InterfaceField &w)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks;
  if(pmb->block_size.nx3 > 1) {// 3D
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie+1; ++i) {
        ei.x1f(X1E2,ks-1,j,i) = ei.x1f(X1E2,ks,j,i); 
        w.x1f(ks-1,j,i) = w.x1f(ks,j,i);
      }
    }
    for (int j=js; j<=je+1; ++j) {
      for (int i=is; i<=ie; ++i) {
        ei.x2f(X2E1,ks-1,j,i) = ei.x2f(X2E1,ks,j,i); 
        w.x2f(ks-1,j,i) = w.x2f(ks,j,i);
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void DefaultEFluxOuterX3(MeshBlock *pmb, InterfaceField &ei, InterfaceField &w)
//  \brief  Default boundary conditions for electric fields for inner x3 boundary 
void DefaultEFluxOuterX3(MeshBlock *pmb, InterfaceField &ei, InterfaceField &w)
{
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ke = pmb->ke;
  if(pmb->block_size.nx3 > 1) {// 3D
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie+1; ++i) {
        ei.x1f(X1E2,ke+1,j,i) = ei.x1f(X1E2,ke,j,i); 
        w.x1f(ke+1,j,i) = w.x1f(ke,j,i);
      }
    }
    for (int j=js; j<=je+1; ++j) {
      for (int i=is; i<=ie; ++i) {
        ei.x2f(X2E1,ke+1,j,i) = ei.x2f(X2E1,ke,j,i); 
        w.x2f(ke+1,j,i) = w.x2f(ke,j,i);
      }
    }
  }
  return;
}

