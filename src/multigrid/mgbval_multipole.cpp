//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mgbval_multipole.cpp
//! \brief Multipole expansion (up to 4 and 16) boundary functions for Multigrid

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../defs.hpp"
#include "multigrid.hpp"


//----------------------------------------------------------------------------------------
//! \fn void MGMultipoleInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
//!            int is, int ie, int js, int je, int ks, int ke, int ngh,
//!            const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
//!            const AthenaArray<Real> &mporigin, int mporder)
//! \brief Multipole boundary condition in the inner-X1 direction
//!  *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MGMultipoleInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
       int is, int ie, int js, int je, int ks, int ke, int ngh,
       const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
       const AthenaArray<Real> &mporigin, int mporder) {
  int i = is;
  Real xorigin = mporigin(0);
  Real yorigin = mporigin(1);
  Real zorigin = mporigin(2);
  Real x = coord.x1f(i) - xorigin, x2 = x*x;
  if (mporder == 2) {
    for (int k = ks; k <= ke; ++k) {
      Real z = coord.x3v(k) - zorigin;
      Real z2 = z*z, zx = z*x;
#pragma ivdep
      for (int j = js; j <= je; ++j) {
        Real y = coord.x2v(j) - yorigin;
        Real y2 = y*y, yz = y*z, xy = x*y;
        Real r2 = x2 + y2 + z2;
        Real ir2 = 1.0/r2, ir1 = std::sqrt(ir2);
        Real ir3 = ir2*ir1, ir5 = ir3*ir2;
        Real phis = ir1*mpcoeff(0)
          + ir3*(mpcoeff(1)*y + mpcoeff(2)*z + mpcoeff(3)*x)
          + ir5*(mpcoeff(4)*xy + mpcoeff(5)*yz + (3.0*z2-r2)*mpcoeff(6)
               + mpcoeff(7)*zx + mpcoeff(8)*0.5*(x2-y2));
        dst(0,k,j,i-1) = 2.0*phis - dst(0,k,j,i);
      }
    }
  } else if (mporder == 4) {
    for (int k = ks; k <= ke; ++k) {
      Real z = coord.x3v(k) - zorigin;
      Real z2 = z*z, zx = z*x;
#pragma ivdep
      for (int j = js; j <= je; ++j) {
        Real y = coord.x2v(j) - yorigin;
        Real y2 = y*y, yz = y*z, xy = x*y;
        Real r2 = x2 + y2 + z2;
        Real ir2 = 1.0/r2, ir1 = std::sqrt(ir2);
        Real ir3 = ir2*ir1, ir5 = ir3*ir2, ir7 = ir5*ir2, ir9 = ir7*ir2;
        Real hx2my2 = 0.5*(x2-y2);
        Real x2mty2 = x2-3.0*y2;
        Real tx2my2 = 3.0*x2-y2;
        Real phis = ir1*mpcoeff(0)
          + ir3*(mpcoeff(1)*y + mpcoeff(2)*z + mpcoeff(3)*x)
          + ir5*(mpcoeff(4)*xy + mpcoeff(5)*yz + (3.0*z2-r2)*mpcoeff(6)
               + mpcoeff(7)*zx + mpcoeff(8)*0.5*(x2-y2))
          + ir7*(y*tx2my2*mpcoeff(9) + x*x2mty2*mpcoeff(15)
               + xy*z*mpcoeff(10) + z*hx2my2*mpcoeff(14)
               + (5.0*z2-r2)*(y*mpcoeff(11) + x*mpcoeff(13))
               + z*(z2-3.0*r2)*mpcoeff(12))
          + ir9*(xy*hx2my2*mpcoeff(16) + 0.125*(x2*x2mty2-y2*tx2my2)*mpcoeff(24)
               + yz*tx2my2*mpcoeff(17) + zx*x2mty2*mpcoeff(23)
               + (7.0*z2-r2)*(xy*mpcoeff(18) + hx2my2*mpcoeff(22))
               + (7.0*z2-3.0*r2)*(yz*mpcoeff(19) + zx*mpcoeff(21))
               + (35.0*z2*z2-30.0*z2*r2+3.0*r2*r2)*mpcoeff(20));
        dst(0,k,j,i-1) = 2.0*phis - dst(0,k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGMultipoleOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
//!            int is, int ie, int js, int je, int ks, int ke, int ngh,
//!            const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
//!            const AthenaArray<Real> &mporigin, int mporder)
//!  *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MGMultipoleOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
       int is, int ie, int js, int je, int ks, int ke, int ngh,
       const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
       const AthenaArray<Real> &mporigin, int mporder) {
  int i = ie;
  Real xorigin = mporigin(0);
  Real yorigin = mporigin(1);
  Real zorigin = mporigin(2);
  Real x = coord.x1f(i+1) - xorigin, x2 = x*x;
  if (mporder == 2) {
    for (int k = ks; k <= ke; ++k) {
      Real z = coord.x3v(k) - zorigin;
      Real z2 = z*z, zx = z*x;
#pragma ivdep
      for (int j = js; j <= je; ++j) {
        Real y = coord.x2v(j) - yorigin;
        Real y2 = y*y, yz = y*z, xy = x*y;
        Real r2 = x2 + y2 + z2;
        Real ir2 = 1.0/r2, ir1 = std::sqrt(ir2);
        Real ir3 = ir2*ir1, ir5 = ir3*ir2;
        Real phis = ir1*mpcoeff(0)
          + ir3*(mpcoeff(1)*y + mpcoeff(2)*z + mpcoeff(3)*x)
          + ir5*(mpcoeff(4)*xy + mpcoeff(5)*yz + (3.0*z2-r2)*mpcoeff(6)
               + mpcoeff(7)*zx + mpcoeff(8)*0.5*(x2-y2));
        dst(0,k,j,i+1) = 2.0*phis - dst(0,k,j,i);
      }
    }
  } else if (mporder == 4) {
    for (int k = ks; k <= ke; ++k) {
      Real z = coord.x3v(k) - zorigin;
      Real z2 = z*z, zx = z*x;
#pragma ivdep
      for (int j = js; j <= je; ++j) {
        Real y = coord.x2v(j) - yorigin;
        Real y2 = y*y, yz = y*z, xy = x*y;
        Real r2 = x2 + y2 + z2;
        Real ir2 = 1.0/r2, ir1 = std::sqrt(ir2);
        Real ir3 = ir2*ir1, ir5 = ir3*ir2, ir7 = ir5*ir2, ir9 = ir7*ir2;
        Real hx2my2 = 0.5*(x2-y2);
        Real x2mty2 = x2-3.0*y2;
        Real tx2my2 = 3.0*x2-y2;
        Real phis = ir1*mpcoeff(0)
          + ir3*(mpcoeff(1)*y + mpcoeff(2)*z + mpcoeff(3)*x)
          + ir5*(mpcoeff(4)*xy + mpcoeff(5)*yz + (3.0*z2-r2)*mpcoeff(6)
               + mpcoeff(7)*zx + mpcoeff(8)*0.5*(x2-y2))
          + ir7*(y*tx2my2*mpcoeff(9) + x*x2mty2*mpcoeff(15)
               + xy*z*mpcoeff(10) + z*hx2my2*mpcoeff(14)
               + (5.0*z2-r2)*(y*mpcoeff(11) + x*mpcoeff(13))
               + z*(z2-3.0*r2)*mpcoeff(12))
          + ir9*(xy*hx2my2*mpcoeff(16) + 0.125*(x2*x2mty2-y2*tx2my2)*mpcoeff(24)
               + yz*tx2my2*mpcoeff(17) + zx*x2mty2*mpcoeff(23)
               + (7.0*z2-r2)*(xy*mpcoeff(18) + hx2my2*mpcoeff(22))
               + (7.0*z2-3.0*r2)*(yz*mpcoeff(19) + zx*mpcoeff(21))
               + (35.0*z2*z2-30.0*z2*r2+3.0*r2*r2)*mpcoeff(20));
        dst(0,k,j,i+1) = 2.0*phis - dst(0,k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGMultipoleInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
//!            int is, int ie, int js, int je, int ks, int ke, int ngh,
//!            const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
//!            const AthenaArray<Real> &mporigin, int mporder)
//! \brief Multipole boundary condition in the inner-X2 direction
//!  *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MGMultipoleInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
       int is, int ie, int js, int je, int ks, int ke, int ngh,
       const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
       const AthenaArray<Real> &mporigin, int mporder) {
  int j = js;
  Real xorigin = mporigin(0);
  Real yorigin = mporigin(1);
  Real zorigin = mporigin(2);
  Real y = coord.x2f(j) - yorigin, y2 = y*y;
  if (mporder == 2) {
    for (int k = ks; k <= ke; ++k) {
      Real z = coord.x3v(k) - zorigin;
      Real z2 = z*z, yz = y*z;
#pragma ivdep
      for (int i = is; i <= ie; ++i) {
        Real x = coord.x1v(i) - xorigin, x2 = x*x, xy = x*y, zx = z*x;
        Real r2 = x2 + y2 + z2;
        Real ir2 = 1.0/r2, ir1 = std::sqrt(ir2);
        Real ir3 = ir2*ir1, ir5 = ir3*ir2;
        Real phis = ir1*mpcoeff(0)
          + ir3*(mpcoeff(1)*y + mpcoeff(2)*z + mpcoeff(3)*x)
          + ir5*(mpcoeff(4)*xy + mpcoeff(5)*yz + (3.0*z2-r2)*mpcoeff(6)
               + mpcoeff(7)*zx + mpcoeff(8)*0.5*(x2-y2));
        dst(0,k,j-1,i) = 2.0*phis - dst(0,k,j,i);
      }
    }
  } else if (mporder == 4) {
    for (int k = ks; k <= ke; ++k) {
      Real z = coord.x3v(k) - zorigin;
      Real z2 = z*z, yz = y*z;
#pragma ivdep
      for (int i = is; i <= ie; ++i) {
        Real x = coord.x1v(i) - xorigin, x2 = x*x, xy = x*y, zx = z*x;
        Real r2 = x2 + y2 + z2;
        Real ir2 = 1.0/r2, ir1 = std::sqrt(ir2);
        Real ir3 = ir2*ir1, ir5 = ir3*ir2, ir7 = ir5*ir2, ir9 = ir7*ir2;
        Real hx2my2 = 0.5*(x2-y2);
        Real x2mty2 = x2-3.0*y2;
        Real tx2my2 = 3.0*x2-y2;
        Real phis = ir1*mpcoeff(0)
          + ir3*(mpcoeff(1)*y + mpcoeff(2)*z + mpcoeff(3)*x)
          + ir5*(mpcoeff(4)*xy + mpcoeff(5)*yz + (3.0*z2-r2)*mpcoeff(6)
               + mpcoeff(7)*zx + mpcoeff(8)*0.5*(x2-y2))
          + ir7*(y*tx2my2*mpcoeff(9) + x*x2mty2*mpcoeff(15)
               + xy*z*mpcoeff(10) + z*hx2my2*mpcoeff(14)
               + (5.0*z2-r2)*(y*mpcoeff(11) + x*mpcoeff(13))
               + z*(z2-3.0*r2)*mpcoeff(12))
          + ir9*(xy*hx2my2*mpcoeff(16) + 0.125*(x2*x2mty2-y2*tx2my2)*mpcoeff(24)
               + yz*tx2my2*mpcoeff(17) + zx*x2mty2*mpcoeff(23)
               + (7.0*z2-r2)*(xy*mpcoeff(18) + hx2my2*mpcoeff(22))
               + (7.0*z2-3.0*r2)*(yz*mpcoeff(19) + zx*mpcoeff(21))
               + (35.0*z2*z2-30.0*z2*r2+3.0*r2*r2)*mpcoeff(20));
        dst(0,k,j-1,i) = 2.0*phis - dst(0,k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGMultipoleOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
//!            int is, int ie, int js, int je, int ks, int ke, int ngh,
//!            const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
//!            const AthenaArray<Real> &mporigin, int mporder)
//! \brief Multipole boundary condition in the outer-X2 direction
//!  *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MGMultipoleOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
       int is, int ie, int js, int je, int ks, int ke, int ngh,
       const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
       const AthenaArray<Real> &mporigin, int mporder) {
  int j = je;
  Real xorigin = mporigin(0);
  Real yorigin = mporigin(1);
  Real zorigin = mporigin(2);
  Real y = coord.x2f(j+1) - yorigin, y2 = y*y;
  if (mporder == 2) {
    for (int k = ks; k <= ke; ++k) {
      Real z = coord.x3v(k) - zorigin;
      Real z2 = z*z, yz = y*z;
#pragma ivdep
      for (int i = is; i <= ie; ++i) {
        Real x = coord.x1v(i) - xorigin, x2 = x*x, xy = x*y, zx = z*x;
        Real r2 = x2 + y2 + z2;
        Real ir2 = 1.0/r2, ir1 = std::sqrt(ir2);
        Real ir3 = ir2*ir1, ir5 = ir3*ir2;
        Real phis = ir1*mpcoeff(0)
          + ir3*(mpcoeff(1)*y + mpcoeff(2)*z + mpcoeff(3)*x)
          + ir5*(mpcoeff(4)*xy + mpcoeff(5)*yz + (3.0*z2-r2)*mpcoeff(6)
               + mpcoeff(7)*zx + mpcoeff(8)*0.5*(x2-y2));
        dst(0,k,j+1,i) = 2.0*phis - dst(0,k,j,i);
      }
    }
  } else if (mporder == 4) {
    for (int k = ks; k <= ke; ++k) {
      Real z = coord.x3v(k) - zorigin;
      Real z2 = z*z, yz = y*z;
#pragma ivdep
      for (int i = is; i <= ie; ++i) {
        Real x = coord.x1v(i) - xorigin, x2 = x*x, xy = x*y, zx = z*x;
        Real r2 = x2 + y2 + z2;
        Real ir2 = 1.0/r2, ir1 = std::sqrt(ir2);
        Real ir3 = ir2*ir1, ir5 = ir3*ir2, ir7 = ir5*ir2, ir9 = ir7*ir2;
        Real hx2my2 = 0.5*(x2-y2);
        Real x2mty2 = x2-3.0*y2;
        Real tx2my2 = 3.0*x2-y2;
        Real phis = ir1*mpcoeff(0)
          + ir3*(mpcoeff(1)*y + mpcoeff(2)*z + mpcoeff(3)*x)
          + ir5*(mpcoeff(4)*xy + mpcoeff(5)*yz + (3.0*z2-r2)*mpcoeff(6)
               + mpcoeff(7)*zx + mpcoeff(8)*0.5*(x2-y2))
          + ir7*(y*tx2my2*mpcoeff(9) + x*x2mty2*mpcoeff(15)
               + xy*z*mpcoeff(10) + z*hx2my2*mpcoeff(14)
               + (5.0*z2-r2)*(y*mpcoeff(11) + x*mpcoeff(13))
               + z*(z2-3.0*r2)*mpcoeff(12))
          + ir9*(xy*hx2my2*mpcoeff(16) + 0.125*(x2*x2mty2-y2*tx2my2)*mpcoeff(24)
               + yz*tx2my2*mpcoeff(17) + zx*x2mty2*mpcoeff(23)
               + (7.0*z2-r2)*(xy*mpcoeff(18) + hx2my2*mpcoeff(22))
               + (7.0*z2-3.0*r2)*(yz*mpcoeff(19) + zx*mpcoeff(21))
               + (35.0*z2*z2-30.0*z2*r2+3.0*r2*r2)*mpcoeff(20));
        dst(0,k,j+1,i) = 2.0*phis - dst(0,k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGMultipoleInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
//!            int is, int ie, int js, int je, int ks, int ke, int ngh,
//!            const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
//!            const AthenaArray<Real> &mporigin, int mporder)
//! \brief Multipole boundary condition in the inner-X3 direction
//!  *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MGMultipoleInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
       int is, int ie, int js, int je, int ks, int ke, int ngh,
       const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
       const AthenaArray<Real> &mporigin, int mporder) {
  int k = ks;
  Real xorigin = mporigin(0);
  Real yorigin = mporigin(1);
  Real zorigin = mporigin(2);
  Real z = coord.x3f(k) - zorigin, z2 = z*z;
  if (mporder == 2) {
    for (int j = js; j <= je; ++j) {
      Real y = coord.x2v(j) - yorigin, y2 = y*y, yz = y*z;
#pragma ivdep
      for (int i = is; i <= ie; ++i) {
        Real x = coord.x1v(i) - xorigin, x2 = x*x, xy = x*y, zx = z*x;
        Real r2 = x2 + y2 + z2;
        Real ir2 = 1.0/r2, ir1 = std::sqrt(ir2);
        Real ir3 = ir2*ir1, ir5 = ir3*ir2;
        Real phis = ir1*mpcoeff(0)
          + ir3*(mpcoeff(1)*y + mpcoeff(2)*z + mpcoeff(3)*x)
          + ir5*(mpcoeff(4)*xy + mpcoeff(5)*yz + (3.0*z2-r2)*mpcoeff(6)
               + mpcoeff(7)*zx + mpcoeff(8)*0.5*(x2-y2));
        dst(0,k-1,j,i) = 2.0*phis - dst(0,k,j,i);
      }
    }
  } else if (mporder == 4) {
    for (int j = js; j <= je; ++j) {
      Real y = coord.x2v(j) - yorigin, y2 = y*y, yz = y*z;
#pragma ivdep
      for (int i = is; i <= ie; ++i) {
        Real x = coord.x1v(i) - xorigin, x2 = x*x, xy = x*y, zx = z*x;
        Real r2 = x2 + y2 + z2;
        Real ir2 = 1.0/r2, ir1 = std::sqrt(ir2);
        Real ir3 = ir2*ir1, ir5 = ir3*ir2, ir7 = ir5*ir2, ir9 = ir7*ir2;
        Real hx2my2 = 0.5*(x2-y2);
        Real x2mty2 = x2-3.0*y2;
        Real tx2my2 = 3.0*x2-y2;
        Real phis = ir1*mpcoeff(0)
          + ir3*(mpcoeff(1)*y + mpcoeff(2)*z + mpcoeff(3)*x)
          + ir5*(mpcoeff(4)*xy + mpcoeff(5)*yz + (3.0*z2-r2)*mpcoeff(6)
               + mpcoeff(7)*zx + mpcoeff(8)*0.5*(x2-y2))
          + ir7*(y*tx2my2*mpcoeff(9) + x*x2mty2*mpcoeff(15)
               + xy*z*mpcoeff(10) + z*hx2my2*mpcoeff(14)
               + (5.0*z2-r2)*(y*mpcoeff(11) + x*mpcoeff(13))
               + z*(z2-3.0*r2)*mpcoeff(12))
          + ir9*(xy*hx2my2*mpcoeff(16) + 0.125*(x2*x2mty2-y2*tx2my2)*mpcoeff(24)
               + yz*tx2my2*mpcoeff(17) + zx*x2mty2*mpcoeff(23)
               + (7.0*z2-r2)*(xy*mpcoeff(18) + hx2my2*mpcoeff(22))
               + (7.0*z2-3.0*r2)*(yz*mpcoeff(19) + zx*mpcoeff(21))
               + (35.0*z2*z2-30.0*z2*r2+3.0*r2*r2)*mpcoeff(20));
        dst(0,k-1,j,i) = 2.0*phis - dst(0,k,j,i);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGMultipoleOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
//!            int is, int ie, int js, int je, int ks, int ke, int ngh,
//!            const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff, int mporder)
//!            const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
//!            const AthenaArray<Real> &mporigin, int mporder)
//! \brief Multipole boundary condition in the outer-X3 direction
//!  *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MGMultipoleOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
       int is, int ie, int js, int je, int ks, int ke, int ngh,
       const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
       const AthenaArray<Real> &mporigin, int mporder) {
  int k = ke;
  Real xorigin = mporigin(0);
  Real yorigin = mporigin(1);
  Real zorigin = mporigin(2);
  Real z = coord.x3f(k+1) - zorigin, z2 = z*z;
  if (mporder == 2) {
    for (int j = js; j <= je; ++j) {
      Real y = coord.x2v(j) - yorigin, y2 = y*y, yz = y*z;
#pragma ivdep
      for (int i = is; i <= ie; ++i) {
        Real x = coord.x1v(i) - xorigin, x2 = x*x, xy = x*y, zx = z*x;
        Real r2 = x2 + y2 + z2;
        Real ir2 = 1.0/r2, ir1 = std::sqrt(ir2);
        Real ir3 = ir2*ir1, ir5 = ir3*ir2;
        Real phis = ir1*mpcoeff(0)
          + ir3*(mpcoeff(1)*y + mpcoeff(2)*z + mpcoeff(3)*x)
          + ir5*(mpcoeff(4)*xy + mpcoeff(5)*yz + (3.0*z2-r2)*mpcoeff(6)
               + mpcoeff(7)*zx + mpcoeff(8)*0.5*(x2-y2));
        dst(0,k+1,j,i) = 2.0*phis - dst(0,k,j,i);
      }
    }
  } else if (mporder == 4) {
    for (int j = js; j <= je; ++j) {
      Real y = coord.x2v(j) - yorigin, y2 = y*y, yz = y*z;
#pragma ivdep
      for (int i = is; i <= ie; ++i) {
        Real x = coord.x1v(i) - xorigin, x2 = x*x, xy = x*y, zx = z*x;
        Real r2 = x2 + y2 + z2;
        Real ir2 = 1.0/r2, ir1 = std::sqrt(ir2);
        Real ir3 = ir2*ir1, ir5 = ir3*ir2, ir7 = ir5*ir2, ir9 = ir7*ir2;
        Real hx2my2 = 0.5*(x2-y2);
        Real x2mty2 = x2-3.0*y2;
        Real tx2my2 = 3.0*x2-y2;
        Real phis = ir1*mpcoeff(0)
          + ir3*(mpcoeff(1)*y + mpcoeff(2)*z + mpcoeff(3)*x)
          + ir5*(mpcoeff(4)*xy + mpcoeff(5)*yz + (3.0*z2-r2)*mpcoeff(6)
               + mpcoeff(7)*zx + mpcoeff(8)*0.5*(x2-y2))
          + ir7*(y*tx2my2*mpcoeff(9) + x*x2mty2*mpcoeff(15)
               + xy*z*mpcoeff(10) + z*hx2my2*mpcoeff(14)
               + (5.0*z2-r2)*(y*mpcoeff(11) + x*mpcoeff(13))
               + z*(z2-3.0*r2)*mpcoeff(12))
          + ir9*(xy*hx2my2*mpcoeff(16) + 0.125*(x2*x2mty2-y2*tx2my2)*mpcoeff(24)
               + yz*tx2my2*mpcoeff(17) + zx*x2mty2*mpcoeff(23)
               + (7.0*z2-r2)*(xy*mpcoeff(18) + hx2my2*mpcoeff(22))
               + (7.0*z2-3.0*r2)*(yz*mpcoeff(19) + zx*mpcoeff(21))
               + (35.0*z2*z2-30.0*z2*r2+3.0*r2*r2)*mpcoeff(20));
        dst(0,k+1,j,i) = 2.0*phis - dst(0,k,j,i);
      }
    }
  }
  return;
}

