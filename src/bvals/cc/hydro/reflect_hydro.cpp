//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reflect_hydro.cpp
//! \brief implementation of reflecting BCs in each dimension

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "bvals_hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroBoundaryVariable::ReflectInnerX1(
//!         Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//! \brief REFLECTING boundary conditions, inner x1 boundary

void HydroBoundaryVariable::ReflectInnerX1(
    Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<=nu_; ++n) {
    if (n == (IVX)) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=1; i<=ngh; ++i) {
            (*var_cc)(IVX,k,j,il-i) = -(*var_cc)(IVX,k,j,(il+i-1));  // reflect 1-velocity
          }
        }
      }
    } else {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=1; i<=ngh; ++i) {
            (*var_cc)(n,k,j,il-i) = (*var_cc)(n,k,j,(il+i-1));
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroBoundaryVariable::ReflectOuterX1(
//!         Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//! \brief REFLECTING boundary conditions, outer x1 boundary

void HydroBoundaryVariable::ReflectOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v1
  for (int n=0; n<=nu_; ++n) {
    if (n == (IVX)) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=1; i<=ngh; ++i) {
            (*var_cc)(IVX,k,j,iu+i) = -(*var_cc)(IVX,k,j,(iu-i+1));  // reflect 1-velocity
          }
        }
      }
    } else {
      for (int k=kl; k<=ku; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=1; i<=ngh; ++i) {
            (*var_cc)(n,k,j,iu+i) = (*var_cc)(n,k,j,(iu-i+1));
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroBoundaryVariable::ReflectInnerX2(
//!         Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//! \brief REFLECTING boundary conditions, inner x2 boundary

void HydroBoundaryVariable::ReflectInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<=nu_; ++n) {
    if (n == (IVY)) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=1; j<=ngh; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            (*var_cc)(IVY,k,jl-j,i) = -(*var_cc)(IVY,k,jl+j-1,i);  // reflect 2-velocity
          }
        }
      }
    } else {
      for (int k=kl; k<=ku; ++k) {
        for (int j=1; j<=ngh; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            (*var_cc)(n,k,jl-j,i) = (*var_cc)(n,k,jl+j-1,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroBoundaryVariable::ReflectOuterX2(
//!         Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//! \brief REFLECTING boundary conditions, outer x2 boundary

void HydroBoundaryVariable::ReflectOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v2
  for (int n=0; n<=nu_; ++n) {
    if (n == (IVY)) {
      for (int k=kl; k<=ku; ++k) {
        for (int j=1; j<=ngh; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            (*var_cc)(IVY,k,ju+j,i) = -(*var_cc)(IVY,k,ju-j+1,i);  // reflect 2-velocity
          }
        }
      }
    } else {
      for (int k=kl; k<=ku; ++k) {
        for (int j=1; j<=ngh; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            (*var_cc)(n,k,ju+j,i) = (*var_cc)(n,k,ju-j+1,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroBoundaryVariable::ReflectInnerX3(
//!         Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh)
//! \brief REFLECTING boundary conditions, inner x3 boundary

void HydroBoundaryVariable::ReflectInnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {
  // copy hydro variables into ghost zones, reflecting v3
  for (int n=0; n<=nu_; ++n) {
    if (n == (IVZ)) {
      for (int k=1; k<=ngh; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            (*var_cc)(IVZ,kl-k,j,i) = -(*var_cc)(IVZ,kl+k-1,j,i);  // reflect 3-velocity
          }
        }
      }
    } else {
      for (int k=1; k<=ngh; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            (*var_cc)(n,kl-k,j,i) = (*var_cc)(n,kl+k-1,j,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void HydroBoundaryVariable::ReflectOuterX3(
//!         Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh)
//! \brief REFLECTING boundary conditions, outer x3 boundary

void HydroBoundaryVariable::ReflectOuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {
  // copy hydro variables into ghost zones, reflecting v3
  for (int n=0; n<=nu_; ++n) {
    if (n == (IVZ)) {
      for (int k=1; k<=ngh; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            (*var_cc)(IVZ,ku+k,j,i) = -(*var_cc)(IVZ,ku-k+1,j,i);  // reflect 3-velocity
          }
        }
      }
    } else {
      for (int k=1; k<=ngh; ++k) {
        for (int j=jl; j<=ju; ++j) {
#pragma omp simd
          for (int i=il; i<=iu; ++i) {
            (*var_cc)(n,ku+k,j,i) = (*var_cc)(n,ku-k+1,j,i);
          }
        }
      }
    }
  }
  return;
}
