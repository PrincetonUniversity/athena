//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_orbital_advection.cpp
//! \brief functions of main calculations in orbital advection

// C/C++ headers
#include <algorithm>  // min()
#include <cfloat>     // FLT_MAX
#include <cmath>      // fabs
#include <cstring>    // memcpy
#include <iostream>   // cout, endl
#include <sstream>    //
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"

// this class header
#include "orbital_advection.hpp"


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::CalculateOrbitalAdvectionCC(const Real dt,
//!                            AthenaArray<Real> &u, AthenaArray<Real> &s)
//! \brief update hydro & passive scalars using orbital advection scheme

void OrbitalAdvection::CalculateOrbitalAdvectionCC(const Real dt,
                                AthenaArray<Real> &u, AthenaArray<Real> &s) {
  if (!orbital_advection_active) return;
  int is = pmb_->is, ie = pmb_->ie;
  int js = pmb_->js, je = pmb_->je;
  int ks = pmb_->ks, ke = pmb_->ke;
  if (orbital_direction ==1) { // cartesian or cylindrical
    if (orbital_uniform_mesh) { // uniform mesh
      for (int k=ks; k<=ke; ++k) {
        for (int i=is; i<=ie; ++i) {
          const Real epsilon = orc(k,i);
          int offset = ofc(k,i);
          const int osgn = (offset>0)?1:0;
          const int shift0 = osgn*onx - offset;
          const int shift  = shift0 + osgn;
          // calculate hydro flux
          for (int nph=0; nph<NHYDRO; ++nph) {
            hbuf.InitWithShallowSlice(orbital_cons, 4, nph, 1);
            if (xorder <= 2) {
              RemapFluxPlm(pflux, hbuf, epsilon, osgn, k, i, js, je+1, shift0);
            } else {
              RemapFluxPpm(pflux, hbuf, epsilon, osgn, k, i, js, je+1, shift0);
            }
            for (int j=js; j<=je; j++) {
              u(nph,k,j,i) = hbuf(k,i,j+shift) - (pflux(j+1) - pflux(j));
            }
          }
          // calculate passive scalar flux
          for (int nsc=0; nsc<NSCALARS; ++nsc) {
            hbuf.InitWithShallowSlice(orbital_scalar, 4, nsc, 1);
            if (xorder <= 2) {
              RemapFluxPlm(pflux, hbuf, epsilon, osgn, k, i, js, je+1, shift0);
            } else {
              RemapFluxPpm(pflux, hbuf, epsilon, osgn, k, i, js, je+1, shift0);
            }
            for (int j=js; j<=je; j++) {
              s(nsc,k,j,i) = hbuf(k,i,j+shift) - (pflux(j+1) - pflux(j));
            }
          }
        }
      }
    }
//    else { // non-uniform mesh
//    }
  } else if (orbital_direction == 2) { // spherical_polar
    if(orbital_uniform_mesh) { // uniform mesh
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          const Real epsilon = orc(j,i);
          int offset = ofc(j,i);
          const int osgn = (offset>0)?1:0;
          const int shift0 = osgn*onx - offset;
          const int shift  = shift0 + osgn;
          for (int nph=0; nph<NHYDRO; ++nph) {
            hbuf.InitWithShallowSlice(orbital_cons, 4, nph, 1);
            if (xorder <= 2) {
              RemapFluxPlm(pflux, hbuf, epsilon, osgn, j, i, ks, ke+1, shift0);
            } else {
              RemapFluxPpm(pflux, hbuf, epsilon, osgn, j, i, ks, ke+1, shift0);
            }
            for (int k=ks; k<=ke; k++) {
              u(nph,k,j,i) = hbuf(j,i,k+shift) - (pflux(k+1) - pflux(k));
            }
          }
          // calculate passive scalar flux
          for (int nsc=0; nsc<NSCALARS; ++nsc) {
            hbuf.InitWithShallowSlice(orbital_scalar, 4, nsc, 1);
            if (xorder <= 2) {
              RemapFluxPlm(pflux, hbuf, epsilon, osgn, j, i, ks, ke+1, shift0);
            } else {
              RemapFluxPpm(pflux, hbuf, epsilon, osgn, j, i, ks, ke+1, shift0);
            }
            for (int k=ks; k<=ke; k++) {
              s(nsc,k,j,i) = hbuf(j,i,k+shift) - (pflux(k+1) - pflux(k));
            }
          }
        }
      }
    }
//    else { // non-uniform mesh
//    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void OrbitalAdvection::CalculateOrbitalAdvectionFC(const Real dt, EdgeField &e)
//! \brief calculate field flux for orbital advection

void OrbitalAdvection::CalculateOrbitalAdvectionFC(const Real dt, EdgeField &e) {
  if (!orbital_advection_active) return;
  int is = pmb_->is, ie = pmb_->ie;
  int js = pmb_->js, je = pmb_->je;
  int ks = pmb_->ks, ke = pmb_->ke;
  if(orbital_direction ==1) { // cartesian or cylindrical
    if(orbital_uniform_mesh) { // uniform mesh
      // e1 = VK*(-B3)
      for (int k=ks; k<=ke+1; ++k) {
        for (int i=is; i<=ie  ; ++i) {
          const Real epsilon = orf[1](k,i);
          int offset   = off[1](k,i);
          const int osgn = (offset>0)?1:0;
          const int shift0 = osgn*onx - offset;
          AthenaArray<Real> &bbuf = orbital_b2;
          if (xorder<=2) {
            RemapFluxPlm(pflux, bbuf, epsilon, osgn, k, i, js, je+1, shift0);
          } else {
            RemapFluxPpm(pflux, bbuf, epsilon, osgn, k, i, js, je+1, shift0);
          }
          for (int jj=1; jj<offset; jj++) {
            const int shift = shift0+jj;
#pragma omp simd
            for (int j=js; j<=je+1; j++) {
              pflux(j) += bbuf(k,i,j+shift);
            }
          }
          for (int jj=0; jj>offset; jj--) {
            const int shift = shift0+jj-1;
#pragma omp simd
            for (int j=js; j<=je+1; j++) {
              pflux(j) -= bbuf(k,i,j+shift);
            }
          }
          Real len = pco_->h2v(i)*dx;
#pragma omp simd
          for (int j=js; j<=je+1; j++) {
            e.x1e(k,j,i) = pflux(j)*len;
          }
        }
      }

      // e2
      e.x2e.ZeroClear();

      // e3 = VK*B1
      for (int k=ks; k<=ke  ; ++k) {
        for (int i=is; i<=ie+1; ++i) {
          const Real epsilon = orf[0](k,i);
          int offset   = off[0](k,i);
          const int osgn = (offset>0)?1:0;
          const int shift0 = osgn*onx - offset;
          AthenaArray<Real> &bbuf = orbital_b1;
          if (xorder<=2) {
            RemapFluxPlm(pflux, bbuf, epsilon, osgn, k, i, js, je+1, shift0);
          } else {
            RemapFluxPpm(pflux, bbuf, epsilon, osgn, k, i, js, je+1, shift0);
          }
          for (int jj=1; jj<offset; jj++) {
            const int shift = shift0+jj;
#pragma omp simd
            for (int j=js; j<=je+1; j++) {
              pflux(j) += bbuf(k,i,j+shift);
            }
          }
          for (int jj=0; jj>offset; jj--) {
            const int shift = shift0+jj-1;
#pragma omp simd
            for (int j=js; j<=je+1; j++) {
              pflux(j) -= bbuf(k,i,j+shift);
            }
          }
          Real len = pco_->h2f(i)*dx;
#pragma omp simd
          for (int j=js; j<=je+1; j++) {
            e.x3e(k,j,i) = pflux(j)*len;
          }
        }
      }
    }
//    else { // non-uniform mesh
//    }
  } else if (orbital_direction ==2) { // spherical_polar
    if(orbital_uniform_mesh) { // uniform mesh
      // e1 = VK*B2
      for (int j=js; j<=je+1; ++j) {
        for (int i=is; i<=ie  ; ++i) {
          const Real epsilon = orf[1](j,i);
          int offset   = off[1](j,i);
          const int osgn = (offset>0)?1:0;
          const int shift0 = osgn*onx - offset;
          AthenaArray<Real> &bbuf = orbital_b2;
          if (xorder<=2) {
            RemapFluxPlm(pflux, bbuf, epsilon, osgn, j, i, ks, ke+1, shift0);
          } else {
            RemapFluxPpm(pflux, bbuf, epsilon, osgn, j, i, ks, ke+1, shift0);
          }
          for (int kk=1; kk<offset; kk++) {
            const int shift = shift0+kk;
#pragma omp simd
            for (int k=ks; k<=ke+1; k++) {
              pflux(k) += bbuf(j,i,k+shift);
            }
          }
          for (int kk=0; kk>offset; kk--) {
            const int shift = shift0+kk-1;
#pragma omp simd
            for (int k=ks; k<=ke+1; k++) {
              pflux(k) -= bbuf(j,i,k+shift);
            }
          }
          Real len = pco_->h2v(i)*pco_->h32f(j)*dx;
#pragma omp simd
          for (int k=ks; k<=ke+1; k++) {
            e.x1e(k,j,i) = pflux(k)*len;
          }
        }
      }

      // e2 = VK*(-B1)
      for (int j=js; j<=je  ; ++j) {
        for (int i=is; i<=ie+1; ++i) {
          const Real epsilon = orf[0](j,i);
          int offset   = off[0](j,i);
          const int osgn = (offset>0)?1:0;
          const int shift0 = osgn*onx - offset;
          AthenaArray<Real> &bbuf = orbital_b1;
          if (xorder<=2) {
            RemapFluxPlm(pflux, bbuf, epsilon, osgn, j, i, ks, ke+1, shift0);
          } else {
            RemapFluxPpm(pflux, bbuf, epsilon, osgn, j, i, ks, ke+1, shift0);
          }
          for (int kk=1; kk<offset; kk++) {
            const int shift = shift0+kk;
#pragma omp simd
            for (int k=ks; k<=ke+1; k++) {
              pflux(k) += bbuf(j,i,k+shift);
            }
          }
          for (int kk=0; kk>offset; kk--) {
            const int shift = shift0+kk-1;
#pragma omp simd
            for (int k=ks; k<=ke+1; k++) {
              pflux(k) -= bbuf(j,i,k+shift);
            }
          }
          Real len = pco_->h2f(i)*pco_->h32v(j)*dx;
#pragma omp simd
          for (int k=ks; k<=ke+1; k++) {
            e.x2e(k,j,i) = pflux(k)*len;
          }
        }
      }

      // e3
      e.x3e.ZeroClear();
    }
//    else { // non-uniform mesh
//    }
  }
  return;
}
