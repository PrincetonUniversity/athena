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
#include "multigrid.hpp"


// constants for multipole expansion
static const Real c0  = 0.5/std::sqrt(PI);
static const Real c1  = std::sqrt(3.0/(4.0*PI))/3.0;
static const Real c2  = 0.25*std::sqrt(5.0/PI)/5.0;
static const Real c2a = 0.5*std::sqrt(15.0/PI)/5.0;
static const Real c30 = 0.25*std::sqrt(7.0/PI)/7.0;
static const Real c31 = 0.25*std::sqrt(21.0/TWO_PI)/7.0;
static const Real c32 = 0.5*std::sqrt(105.0/PI)/7.0;
static const Real c33 = 0.25*std::sqrt(35.0/TWO_PI)/7.0;
static const Real c40 = 0.1875/std::sqrt(PI)/9.0;
static const Real c41 = 0.75*std::sqrt(5.0/TWO_PI)/9.0;
static const Real c42 = 0.75*std::sqrt(5.0/PI)/9.0;
static const Real c43 = 0.75*std::sqrt(35.0/TWO_PI)/9.0;
static const Real c44 = 1.5*std::sqrt(35.0/PI)/9.0;


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::MGMultipole4InnerX1(AthenaArray<Real> &dst, Real time,
//                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                    const MGCoordinates &coord)
//  \brief Multipole4 boundary condition in the inner-X1 direction
//   *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MultigridDriver::MGMultipole4InnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord) {
  int i = is;
  Real x = coord.x1f(i+1), x2 = x*x;
  for (int k = ks; k <= ke; ++k) {
    Real z = coord.x3v(k);
    Real z2 = z*z, zx = z*x;
    for (int j = js; j <= je; ++j) {
      Real y = coord.x2v(j);
      Real y2 = y*y, yz = y*z, xy = x*y;
      Real r2 = x2 + y2 + z2;
      Real ir2 = 1.0/r2, ir = std::sqrt(ir2);
      Real ir3 = ir2*ir, ir5 = ir3*ir2;
      Real phis = ir1*mpcoeff_(0)
        + ir3*c1*(mpcoeff_(1)*y + mpcoeff_(2)*z + mpcoeff_(3)*x)
        + ir5*(c2a*(mpcoeff_(4)*xy + mpcoeff_(5)*yz + mpcoeff_(7)*zx
                   +mpcoeff_(8)*0.5*(x2-y2)) + c2*(3.0*z2-r2)*mpcoeff_(6));
      dst(0,k,j,i) = 2.0*phis - dst(0,k,j,i+1);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::MGMultipole4OuterX1(AthenaArray<Real> &dst, Real time,
//                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                    const MGCoordinates &coord)
//   *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MultigridDriver::MGMultipole4OuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord) {
  int i = ie;
  Real x = coord.x1f(i), x2 = x*x;
  for (int k = ks; k <= ke; ++k) {
    Real z = coord.x3v(k);
    Real z2 = z*z, zx = z*x;
    for (int j = js; j <= je; ++j) {
      Real y = coord.x2v(j);
      Real y2 = y*y, yz = y*z, xy = x*y;
      Real r2 = x2 + y2 + z2;
      Real ir2 = 1.0/r2, ir = std::sqrt(ir2);
      Real ir3 = ir2*ir, ir5 = ir3*ir2;
      Real phis = ir1*mpcoeff_(0)
        + ir3*c1*(mpcoeff_(1)*y + mpcoeff_(2)*z + mpcoeff_(3)*x)
        + ir5*(c2a*(mpcoeff_(4)*xy + mpcoeff_(5)*yz + mpcoeff_(7)*zx
                   +mpcoeff_(8)*0.5*(x2-y2)) + c2*(3.0*z2-r2)*mpcoeff_(6));
      dst(0,k,j,i) = 2.0*phis - dst(0,k,j,i-1);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::MGMultipole4InnerX2(AthenaArray<Real> &dst, Real time,
//                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                    const MGCoordinates &coord)
//  \brief Multipole4 boundary condition in the inner-X2 direction
//   *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MultigridDriver::MGMultipole4InnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord) {
  int j = js;
  Real y = coord.x2f(j+1), y2 = y*y;
  for (int k = ks; k <= ke; ++k) {
    Real z = coord.x3v(k);
    Real z2 = z*z, yz = y*z;
    for (int i = is; i <= ie; ++i) {
      Real x = coord.x1v(i), x2 = x*x, xy = x*y, zx = z*x;
      Real r2 = x2 + y2 + z2;
      Real ir2 = 1.0/r2, ir = std::sqrt(ir2);
      Real ir3 = ir2*ir, ir5 = ir3*ir2;
      Real phis = ir1*mpcoeff_(0)
        + ir3*c1*(mpcoeff_(1)*y + mpcoeff_(2)*z + mpcoeff_(3)*x)
        + ir5*(c2a*(mpcoeff_(4)*xy + mpcoeff_(5)*yz + mpcoeff_(7)*zx
                   +mpcoeff_(8)*0.5*(x2-y2)) + c2*(3.0*z2-r2)*mpcoeff_(6));
      dst(0,k,j,i) = 2.0*phis - dst(0,k,j+1,i);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::MGMultipole4OuterX2(AthenaArray<Real> &dst, Real time,
//                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                    const MGCoordinates &coord)
//  \brief Multipole4 boundary condition in the outer-X2 direction
//   *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MultigridDriver::MGMultipole4OuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord) {
  int j = je;
  Real y = coord.x2f(j), y2 = y*y;
  for (int k = ks; k <= ke; ++k) {
    Real z = coord.x3v(k);
    Real z2 = z*z, yz = y*z;
    for (int i = is; i <= ie; ++i) {
      Real x = coord.x1v(i), x2 = x*x, xy = x*y, zx = z*x;
      Real r2 = x2 + y2 + z2;
      Real ir2 = 1.0/r2, ir = std::sqrt(ir2);
      Real ir3 = ir2*ir, ir5 = ir3*ir2;
      Real phis = ir1*mpcoeff_(0)
        + ir3*c1*(mpcoeff_(1)*y + mpcoeff_(2)*z + mpcoeff_(3)*x)
        + ir5*(c2a*(mpcoeff_(4)*xy + mpcoeff_(5)*yz + mpcoeff_(7)*zx
                   +mpcoeff_(8)*0.5*(x2-y2)) + c2*(3.0*z2-r2)*mpcoeff_(6));
      dst(0,k,j,i) = 2.0*phis - dst(0,k,j-1,i);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::MGMultipole4InnerX3(AthenaArray<Real> &dst, Real time,
//                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                    const MGCoordinates &coord)
//  \brief Multipole4 boundary condition in the inner-X3 direction
//   *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MultigridDriver::MGMultipole4InnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord) {
  int k = ks;
  Real z = coord.x3f(k+1);
  Real z2 = z*z;
  for (int j = js; j <= je; ++j) {
    Real y = coord.x2v(j), y2 = y*y, yz = y*z;
    for (int i = is; i <= ie; ++i) {
      Real x = coord.x1v(i), x2 = x*x, xy = x*y, zx = z*x;
      Real r2 = x2 + y2 + z2;
      Real ir2 = 1.0/r2, ir = std::sqrt(ir2);
      Real ir3 = ir2*ir, ir5 = ir3*ir2;
      Real phis = ir1*mpcoeff_(0)
        + ir3*c1*(mpcoeff_(1)*y + mpcoeff_(2)*z + mpcoeff_(3)*x)
        + ir5*(c2a*(mpcoeff_(4)*xy + mpcoeff_(5)*yz + mpcoeff_(7)*zx
                   +mpcoeff_(8)*0.5*(x2-y2)) + c2*(3.0*z2-r2)*mpcoeff_(6));
      dst(0,k,j,i) = 2.0*phis - dst(0,k+1,j,i);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::MGMultipole4OuterX3(AthenaArray<Real> &dst, Real time,
//                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                    const MGCoordinates &coord)
//  \brief Multipole4 boundary condition in the outer-X3 direction
//   *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MultigridDriver::MGMultipole4OuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord) {
  int k = ke;
  Real z = coord.x3f(k);
  Real z2 = z*z;
  for (int j = js; j <= je; ++j) {
    Real y = coord.x2v(j), y2 = y*y, yz = y*z;
    for (int i = is; i <= ie; ++i) {
      Real x = coord.x1v(i), x2 = x*x, xy = x*y, zx = z*x;
      Real r2 = x2 + y2 + z2;
      Real ir2 = 1.0/r2, ir = std::sqrt(ir2);
      Real ir3 = ir2*ir, ir5 = ir3*ir2;
      Real phis = ir1*mpcoeff_(0)
        + ir3*c1*(mpcoeff_(1)*y + mpcoeff_(2)*z + mpcoeff_(3)*x)
        + ir5*(c2a*(mpcoeff_(4)*xy + mpcoeff_(5)*yz + mpcoeff_(7)*zx
                   +mpcoeff_(8)*0.5*(x2-y2)) + c2*(3.0*z2-r2)*mpcoeff_(6));
      dst(0,k,j,i) = 2.0*phis - dst(0,k-1,j,i);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::MGMultipole16InnerX1(AthenaArray<Real> &dst, Real time,
//                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                    const MGCoordinates &coord)
//  \brief Multipole16 boundary condition in the inner-X1 direction
//   *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MultigridDriver::MGMultipole16InnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       const MGCoordinates &coord) {
  int i = is;
  Real x = coord.x1f(i+1), x2 = x*x;
  for (int k = ks; k <= ke; ++k) {
    Real z = coord.x3v(k);
    Real z2 = z*z, zx = z*x;
    for (int j = js; j <= je; ++j) {
      Real y = coord.x2v(j);
      Real y2 = y*y, yz = y*z, xy = x*y;
      Real r2 = x2 + y2 + z2;
      Real ir2 = 1.0/r2, ir = std::sqrt(ir2);
      Real ir3 = ir2*ir, ir5 = ir3*ir2, ir7 = ir5*ir2, ir9 = ir7*ir2;
      Real hx2my2 = 0.5*(x2-y2);
      Real x2mty2 = x2-3.0*y2;
      Real tx2my2 = 3.0*x2-y2;
      Real phis = ir1*mpcoeff_(0)
        + ir3*c1*(mpcoeff_(1)*y + mpcoeff_(2)*z + mpcoeff_(3)*x)
        + ir5*(c2a*(mpcoeff_(4)*xy + mpcoeff_(5)*yz + mpcoeff_(7)*zx
                   +mpcoeff_(8)*hx2my2) + c2*(3.0*z2-r2)*mpcoeff_(6))
        + ir7*(c33*(y*tx2my2*mpcoeff_(9) + x*x2mty2*mpcoeff_(15))
              +c32*(xy*z*mpcoeff_(10) + z*hx2my2*mpcoeff_(14))
              +c31*(5.0*z2-r2)*(y*mpcoeff_(11) + x*mpcoeff_(13))
              +c30*z*(z2-3.0*r2)*mpcoeff_(12))
        + ir9*(c44*(xy*hx2my2*mpcoeff_(16) + 0.125*(x2*x2mty2-y2*tx2my2)*mpcoeff_(24))
              +c43*(yz*tx2my2*mpcoeff_(17) + zx*x2mty2*mpcoeff_(23))
              +c42*(7.0*z2-r2)*(xy*mpcoeff_(18) + hx2my2*mpcoeff_(22))
              +c41*(7.0*z2-3.0*r2)*(yz*mpcoeff_(19) + zx*mpcoeff_(21))
              *c40*(35.0*z2*z2-30.0*z2*r2+3.0*r2*r2)*mpcoeff_(20));
      dst(0,k,j,i) = 2.0*phis - dst(0,k,j,i+1);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::MGMultipole16OuterX1(AthenaArray<Real> &dst, Real time,
//                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                    const MGCoordinates &coord)
//  \brief Multipole16 boundary condition in the outer-X1 direction
//   *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MultigridDriver::MGMultipole16OuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                          int is, int ie, int js, int je, int ks, int ke, int ngh,
                          const MGCoordinates &coord) {
  int i = ie;
  Real x = coord.x1f(i), x2 = x*x;
  for (int k = ks; k <= ke; ++k) {
    Real z = coord.x3v(k);
    Real z2 = z*z, zx = z*x;
    for (int j = js; j <= je; ++j) {
      Real y = coord.x2v(j);
      Real y2 = y*y, yz = y*z, xy = x*y;
      Real r2 = x2 + y2 + z2;
      Real ir2 = 1.0/r2, ir = std::sqrt(ir2);
      Real ir3 = ir2*ir, ir5 = ir3*ir2, ir7 = ir5*ir2, ir9 = ir7*ir2;
      Real hx2my2 = 0.5*(x2-y2);
      Real x2mty2 = x2-3.0*y2;
      Real tx2my2 = 3.0*x2-y2;
      Real phis = ir1*mpcoeff_(0)
        + ir3*c1*(mpcoeff_(1)*y + mpcoeff_(2)*z + mpcoeff_(3)*x)
        + ir5*(c2a*(mpcoeff_(4)*xy + mpcoeff_(5)*yz + mpcoeff_(7)*zx
                   +mpcoeff_(8)*hx2my2) + c2*(3.0*z2-r2)*mpcoeff_(6))
        + ir7*(c33*(y*tx2my2*mpcoeff_(9) + x*x2mty2*mpcoeff_(15))
              +c32*(xy*z*mpcoeff_(10) + z*hx2my2*mpcoeff_(14))
              +c31*(5.0*z2-r2)*(y*mpcoeff_(11) + x*mpcoeff_(13))
              +c30*z*(z2-3.0*r2)*mpcoeff_(12))
        + ir9*(c44*(xy*hx2my2*mpcoeff_(16) + 0.125*(x2*x2mty2-y2*tx2my2)*mpcoeff_(24))
              +c43*(yz*tx2my2*mpcoeff_(17) + zx*x2mty2*mpcoeff_(23))
              +c42*(7.0*z2-r2)*(xy*mpcoeff_(18) + hx2my2*mpcoeff_(22))
              +c41*(7.0*z2-3.0*r2)*(yz*mpcoeff_(19) + zx*mpcoeff_(21))
              *c40*(35.0*z2*z2-30.0*z2*r2+3.0*r2*r2)*mpcoeff_(20));
      dst(0,k,j,i) = 2.0*phis - dst(0,k,j,i-1);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::MGMultipole16InnerX2(AthenaArray<Real> &dst, Real time,
//                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                    const MGCoordinates &coord)
//  \brief Multipole16 boundary condition in the inner-X2 direction
//   *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MultigridDriver::MGMultipole16InnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord) {
  int j = js;
  Real y = coord.x2f(j+1), y2 = y*y;
  for (int k = ks; k <= ke; ++k) {
    Real z = coord.x3v(k);
    Real z2 = z*z, yz = y*z;
    for (int i = is; i <= ie; ++i) {
      Real x = coord.x1v(i), x2 = x*x, xy = x*y, zx = z*x;
      Real r2 = x2 + y2 + z2;
      Real ir2 = 1.0/r2, ir = std::sqrt(ir2);
      Real ir3 = ir2*ir, ir5 = ir3*ir2, ir7 = ir5*ir2, ir9 = ir7*ir2;
      Real hx2my2 = 0.5*(x2-y2);
      Real x2mty2 = x2-3.0*y2;
      Real tx2my2 = 3.0*x2-y2;
      Real phis = ir1*mpcoeff_(0)
        + ir3*c1*(mpcoeff_(1)*y + mpcoeff_(2)*z + mpcoeff_(3)*x)
        + ir5*(c2a*(mpcoeff_(4)*xy + mpcoeff_(5)*yz + mpcoeff_(7)*zx
                   +mpcoeff_(8)*hx2my2) + c2*(3.0*z2-r2)*mpcoeff_(6))
        + ir7*(c33*(y*tx2my2*mpcoeff_(9) + x*x2mty2*mpcoeff_(15))
              +c32*(xy*z*mpcoeff_(10) + z*hx2my2*mpcoeff_(14))
              +c31*(5.0*z2-r2)*(y*mpcoeff_(11) + x*mpcoeff_(13))
              +c30*z*(z2-3.0*r2)*mpcoeff_(12))
        + ir9*(c44*(xy*hx2my2*mpcoeff_(16) + 0.125*(x2*x2mty2-y2*tx2my2)*mpcoeff_(24))
              +c43*(yz*tx2my2*mpcoeff_(17) + zx*x2mty2*mpcoeff_(23))
              +c42*(7.0*z2-r2)*(xy*mpcoeff_(18) + hx2my2*mpcoeff_(22))
              +c41*(7.0*z2-3.0*r2)*(yz*mpcoeff_(19) + zx*mpcoeff_(21))
              *c40*(35.0*z2*z2-30.0*z2*r2+3.0*r2*r2)*mpcoeff_(20));
      dst(0,k,j,i) = 2.0*phis - dst(0,k,j+1,i);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::MGMultipole16OuterX2(AthenaArray<Real> &dst, Real time,
//                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                    const MGCoordinates &coord)
//  \brief Multipole16 boundary condition in the outer-X2 direction
//   *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MultigridDriver::MGMultipole16OuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord) {
  int j = je;
  Real y = coord.x2f(j), y2 = y*y;
  for (int k = ks; k <= ke; ++k) {
    Real z = coord.x3v(k);
    Real z2 = z*z, yz = y*z;
    for (int i = is; i <= ie; ++i) {
      Real x = coord.x1v(i), x2 = x*x, xy = x*y, zx = z*x;
      Real r2 = x2 + y2 + z2;
      Real ir2 = 1.0/r2, ir = std::sqrt(ir2);
      Real ir3 = ir2*ir, ir5 = ir3*ir2, ir7 = ir5*ir2, ir9 = ir7*ir2;
      Real hx2my2 = 0.5*(x2-y2);
      Real x2mty2 = x2-3.0*y2;
      Real tx2my2 = 3.0*x2-y2;
      Real phis = ir1*mpcoeff_(0)
        + ir3*c1*(mpcoeff_(1)*y + mpcoeff_(2)*z + mpcoeff_(3)*x)
        + ir5*(c2a*(mpcoeff_(4)*xy + mpcoeff_(5)*yz + mpcoeff_(7)*zx
                   +mpcoeff_(8)*hx2my2) + c2*(3.0*z2-r2)*mpcoeff_(6))
        + ir7*(c33*(y*tx2my2*mpcoeff_(9) + x*x2mty2*mpcoeff_(15))
              +c32*(xy*z*mpcoeff_(10) + z*hx2my2*mpcoeff_(14))
              +c31*(5.0*z2-r2)*(y*mpcoeff_(11) + x*mpcoeff_(13))
              +c30*z*(z2-3.0*r2)*mpcoeff_(12))
        + ir9*(c44*(xy*hx2my2*mpcoeff_(16) + 0.125*(x2*x2mty2-y2*tx2my2)*mpcoeff_(24))
              +c43*(yz*tx2my2*mpcoeff_(17) + zx*x2mty2*mpcoeff_(23))
              +c42*(7.0*z2-r2)*(xy*mpcoeff_(18) + hx2my2*mpcoeff_(22))
              +c41*(7.0*z2-3.0*r2)*(yz*mpcoeff_(19) + zx*mpcoeff_(21))
              *c40*(35.0*z2*z2-30.0*z2*r2+3.0*r2*r2)*mpcoeff_(20));
      dst(0,k,j,i) = 2.0*phis - dst(0,k,j-1,i);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::MGMultipole16InnerX3(AthenaArray<Real> &dst, Real time,
//                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                    const MGCoordinates &coord)
//  \brief Multipole16 boundary condition in the inner-X3 direction
//   *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MultigridDriver::MGMultipole16InnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord) {
  int k = ks;
  Real z = coord.x3f(k+1);
  Real z2 = z*z;
  for (int j = js; j <= je; ++j) {
    Real y = coord.x2v(j), y2 = y*y, yz = y*z;
    for (int i = is; i <= ie; ++i) {
      Real x = coord.x1v(i), x2 = x*x, xy = x*y, zx = z*x;
      Real r2 = x2 + y2 + z2;
      Real ir2 = 1.0/r2, ir = std::sqrt(ir2);
      Real ir3 = ir2*ir, ir5 = ir3*ir2, ir7 = ir5*ir2, ir9 = ir7*ir2;
      Real hx2my2 = 0.5*(x2-y2);
      Real x2mty2 = x2-3.0*y2;
      Real tx2my2 = 3.0*x2-y2;
      Real phis = ir1*mpcoeff_(0)
        + ir3*c1*(mpcoeff_(1)*y + mpcoeff_(2)*z + mpcoeff_(3)*x)
        + ir5*(c2a*(mpcoeff_(4)*xy + mpcoeff_(5)*yz + mpcoeff_(7)*zx
                   +mpcoeff_(8)*hx2my2) + c2*(3.0*z2-r2)*mpcoeff_(6))
        + ir7*(c33*(y*tx2my2*mpcoeff_(9) + x*x2mty2*mpcoeff_(15))
              +c32*(xy*z*mpcoeff_(10) + z*hx2my2*mpcoeff_(14))
              +c31*(5.0*z2-r2)*(y*mpcoeff_(11) + x*mpcoeff_(13))
              +c30*z*(z2-3.0*r2)*mpcoeff_(12))
        + ir9*(c44*(xy*hx2my2*mpcoeff_(16) + 0.125*(x2*x2mty2-y2*tx2my2)*mpcoeff_(24))
              +c43*(yz*tx2my2*mpcoeff_(17) + zx*x2mty2*mpcoeff_(23))
              +c42*(7.0*z2-r2)*(xy*mpcoeff_(18) + hx2my2*mpcoeff_(22))
              +c41*(7.0*z2-3.0*r2)*(yz*mpcoeff_(19) + zx*mpcoeff_(21))
              *c40*(35.0*z2*z2-30.0*z2*r2+3.0*r2*r2)*mpcoeff_(20));
      dst(0,k,j,i) = 2.0*phis - dst(0,k+1,j,i);
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MultigridDriver::MGMultipole16OuterX3(AthenaArray<Real> &dst, Real time,
//                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
//                    const MGCoordinates &coord)
//  \brief Multipole16 boundary condition in the outer-X3 direction
//   *** Note ***: Currently this calculates the zeroth variable and nghost = 1 only.

void MultigridDriver::MGMultipole16OuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord) {
  int k = ke;
  Real z = coord.x3f(k);
  Real z2 = z*z;
  for (int j = js; j <= je; ++j) {
    Real y = coord.x2v(j), y2 = y*y, yz = y*z;
    for (int i = is; i <= ie; ++i) {
      Real x = coord.x1v(i), x2 = x*x, xy = x*y, zx = z*x;
      Real r2 = x2 + y2 + z2;
      Real ir2 = 1.0/r2, ir = std::sqrt(ir2);
      Real ir3 = ir2*ir, ir5 = ir3*ir2, ir7 = ir5*ir2, ir9 = ir7*ir2;
      Real hx2my2 = 0.5*(x2-y2);
      Real x2mty2 = x2-3.0*y2;
      Real tx2my2 = 3.0*x2-y2;
      Real phis = ir1*mpcoeff_(0)
        + ir3*c1*(mpcoeff_(1)*y + mpcoeff_(2)*z + mpcoeff_(3)*x)
        + ir5*(c2a*(mpcoeff_(4)*xy + mpcoeff_(5)*yz + mpcoeff_(7)*zx
                   +mpcoeff_(8)*hx2my2) + c2*(3.0*z2-r2)*mpcoeff_(6))
        + ir7*(c33*(y*tx2my2*mpcoeff_(9) + x*x2mty2*mpcoeff_(15))
              +c32*(xy*z*mpcoeff_(10) + z*hx2my2*mpcoeff_(14))
              +c31*(5.0*z2-r2)*(y*mpcoeff_(11) + x*mpcoeff_(13))
              +c30*z*(z2-3.0*r2)*mpcoeff_(12))
        + ir9*(c44*(xy*hx2my2*mpcoeff_(16) + 0.125*(x2*x2mty2-y2*tx2my2)*mpcoeff_(24))
              +c43*(yz*tx2my2*mpcoeff_(17) + zx*x2mty2*mpcoeff_(23))
              +c42*(7.0*z2-r2)*(xy*mpcoeff_(18) + hx2my2*mpcoeff_(22))
              +c41*(7.0*z2-3.0*r2)*(yz*mpcoeff_(19) + zx*mpcoeff_(21))
              *c40*(35.0*z2*z2-30.0*z2*r2+3.0*r2*r2)*mpcoeff_(20));
      dst(0,k,j,i) = 2.0*phis - dst(0,k-1,j,i);
    }
  }
  return;
}
