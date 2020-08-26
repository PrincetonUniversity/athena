//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mgbval_multipole.cpp
//  \brief Multipole expansion (up to 4 and 16) boundary functions for Multigrid

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../defs.hpp"

//----------------------------------------------------------------------------------------
//! \fn MGMultipole4InnerX1(AthenaArray<Real> &dst, Real time, int nvar,
//                          int is, int ie, int js, int je, int ks, int ke, int ngh,
//                          Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief  Multipole4 boundary condition in the inner-X1 direction

void MGMultipole4InnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                         int is, int ie, int js, int je, int ks, int ke, int ngh,
                         Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=0; i<ngh; i++)
          dst(n,k,j,is-i-1) = dst(n,k,j,is);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGMultipole4OuterX1(AthenaArray<Real> &dst, Real time, int nvar,
//                          int is, int ie, int js, int je, int ks, int ke, int ngh,
//                          Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief  Multipole4 boundary condition in the outer-X1 direction

void MGMultipole4OuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=0; i<ngh; i++)
          dst(n,k,j,ie+i+1) = dst(n,k,j,ie);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGMultipole4InnerX2(AthenaArray<Real> &dst, Real time, int nvar,
//                          int is, int ie, int js, int je, int ks, int ke, int ngh,
//                          Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief  Multipole4 boundary condition in the inner-X2 direction

void MGMultipole4InnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                         int is, int ie, int js, int je, int ks, int ke, int ngh,
                         Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=0; j<ngh; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,k,js-j-1,i) = dst(n,k,js,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGMultipole4OuterX2(AthenaArray<Real> &dst, Real time, int nvar,
//                          int is, int ie, int js, int je, int ks, int ke, int ngh,
//                          Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief  Multipole4 boundary condition in the outer-X2 direction

void MGMultipole4OuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                         int is, int ie, int js, int je, int ks, int ke, int ngh,
                         Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=0; j<ngh; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,k,je+j+1,i) = dst(n,k,je,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGMultipole4InnerX3(AthenaArray<Real> &dst, Real time, int nvar,
//                          int is, int ie, int js, int je, int ks, int ke, int ngh,
//                          Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief  Multipole4 boundary condition in the inner-X3 direction

void MGMultipole4InnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                         int is, int ie, int js, int je, int ks, int ke, int ngh,
                         Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=0; k<ngh; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,ks-k-1,j,i) = dst(n,ks,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGMultipole4OuterX3(AthenaArray<Real> &dst, Real time, int nvar,
//                          int is, int ie, int js, int je, int ks, int ke, int ngh,
//                          Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief  Multipole4 boundary condition in the outer-X3 direction

void MGMultipole4OuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                         int is, int ie, int js, int je, int ks, int ke, int ngh,
                         Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=0; k<ngh; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,ke+k+1,j,i) = dst(n,ke,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGMultipole16InnerX1(AthenaArray<Real> &dst, Real time, int nvar,
//                            int is, int ie, int js, int je, int ks, int ke, int ngh,
//                            Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief  Multipole16 boundary condition in the inner-X1 direction

void MGMultipole16InnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                          int is, int ie, int js, int je, int ks, int ke, int ngh,
                          Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=0; i<ngh; i++)
          dst(n,k,j,is-i-1) = dst(n,k,j,is);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGMultipole16OuterX1(AthenaArray<Real> &dst, Real time, int nvar,
//                           int is, int ie, int js, int je, int ks, int ke, int ngh,
//                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief  Multipole16 boundary condition in the outer-X1 direction

void MGMultipole16OuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                          int is, int ie, int js, int je, int ks, int ke, int ngh,
                          Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=0; i<ngh; i++)
          dst(n,k,j,ie+i+1) = dst(n,k,j,ie);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGMultipole16InnerX2(AthenaArray<Real> &dst, Real time, int nvar,
//                           int is, int ie, int js, int je, int ks, int ke, int ngh,
//                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief  Multipole16 boundary condition in the inner-X2 direction

void MGMultipole16InnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                          int is, int ie, int js, int je, int ks, int ke, int ngh,
                          Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=0; j<ngh; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,k,js-j-1,i) = dst(n,k,js,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGMultipole16OuterX2(AthenaArray<Real> &dst, Real time, int nvar,
//                           int is, int ie, int js, int je, int ks, int ke, int ngh,
//                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief  Multipole16 boundary condition in the outer-X2 direction

void MGMultipole16OuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                          int is, int ie, int js, int je, int ks, int ke, int ngh,
                          Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=ks; k<=ke; k++) {
      for (int j=0; j<ngh; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,k,je+j+1,i) = dst(n,k,je,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGMultipole16InnerX3(AthenaArray<Real> &dst, Real time, int nvar,
//                           int is, int ie, int js, int je, int ks, int ke, int ngh,
//                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief  Multipole16 boundary condition in the inner-X3 direction

void MGMultipole16InnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                          int is, int ie, int js, int je, int ks, int ke, int ngh,
                          Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=0; k<ngh; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,ks-k-1,j,i) = dst(n,ks,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn MGMultipole16OuterX3(AthenaArray<Real> &dst, Real time, int nvar,
//                           int is, int ie, int js, int je, int ks, int ke, int ngh,
//                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz)
//  \brief  Multipole16 boundary condition in the outer-X3 direction

void MGMultipole16OuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                          int is, int ie, int js, int je, int ks, int ke, int ngh,
                          Real x0, Real y0, Real z0, Real dx, Real dy, Real dz) {
  for (int n=0; n<nvar; n++) {
    for (int k=0; k<ngh; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++)
          dst(n,ke+k+1,j,i) = dst(n,ke,j,i);
      }
    }
  }
  return;
}



