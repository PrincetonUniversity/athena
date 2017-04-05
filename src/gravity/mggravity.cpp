//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mggravity.cpp
//  \brief implementation of functions in class MGGravity

// Athena++ headers
#include "mggravity.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../multigrid/multigrid.hpp"

//----------------------------------------------------------------------------------------
//! \fn  void MGGravity::Smooth(int color)
//  \brief Red-Black Gauss-Seidel Smoother
void MGGravity::Smooth(int color)
{
  int ns=ngh, ne=ngh+(1<<lev)-1;
  int c=color;
  AthenaArray<Real> &u=u_[lev];
  AthenaArray<Real> &src=src_[lev];
  int ll=nlev_-1-lev;
  int is, ie, js, je, ks, ke;
  if(ngh_==2 && color==0) {
    is=js=ks=ngh_-1;
    ie=is+(nx_>>ll), je=js+(ny_>>ll), ke=ks+(nz_>>ll);
    c=1;
  }
  else {
    is=js=ks=ngh_;
    ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
  }
  Real dx = rdx_/(Real)(1<<ll);
  Real dx2 = SQR(dx), isix=omega_/6.0;
#pragma ivdep
  for(int k=ks; k<=ke; k++) {
#pragma ivdep
    for(int j=js; j<=je; j++) {
#pragma ivdep
      for(int i=is+c; i<=ie; i+=2)
        u(0,k,j,i)-=((6.0*u(0,k,j,i)-u(0,k+1,j,i)-u(0,k,j+1,i)-u(0,k,j,i+1)
                     -u(0,k-1,j,i)-u(0,k,j-1,i)-u(0,k,j,i-1))+src(0,k,j,i)*dx2)*isix;
      c^=1:
    }
    c^=1:
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MGGravity::CalculateResidual(void)
//  \brief calculate the residual

void MGGravity::CalculateResidual(void)
{
  int ns=ngh, ne=ngh+(1<<current_level_)-1;
  AthenaArray<Real> &u=u_[current_level_];
  AthenaArray<Real> &src=src_[current_level_];
  AthenaArray<Real> &res=res_[current_level_];
  int ll=ncurrent_level__-1-current_level_;
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+(nx_>>ll)-1, je=js+(ny_>>ll)-1, ke=ks+(nz_>>ll)-1;
  Real dx = rdx_/(Real)(1<<ll);
  Real idx2 = 1.0/SQR(dx);
#pragma ivdep
  for(int k=ks; k<=ke; k++) {
#pragma ivdep
    for(int j=js; j<=je; j++) {
#pragma ivdep
      for(int i=is; i<=ie; i++)
        res(0,k,j,i)=(6.0*u(0,k,j,i)-u(0,k+1,j,i)-u(0,k,j+1,i)-u(0,k,j,i+1)
                         -u(0,k-1,j,i)-u(0,k,j-1,i)-u(0,k,j,i-1))*idx2+src(0,k,j,i);
    }
  }
  return;
}

