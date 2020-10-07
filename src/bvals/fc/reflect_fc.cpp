//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reflect_fc.cpp
//! \brief implementation of reflecting BCs in each dimension

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "bvals_fc.hpp"

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ReflectInnerX1(
//!              Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//! \brief REFLECTING boundary conditions, inner x1 boundary

void FaceCenteredBoundaryVariable::ReflectInnerX1(
    Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh) {
  // copy face-centered magnetic fields into ghost zones, reflecting b1
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_fc).x1f(k,j,(il-i)) = -(*var_fc).x1f(k,j,(il+i  ));  // reflect 1-field
      }
    }
  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_fc).x2f(k,j,(il-i)) =  (*var_fc).x2f(k,j,(il+i-1));
      }
    }
  }

  for (int k=kl; k<=ku+1; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_fc).x3f(k,j,(il-i)) =  (*var_fc).x3f(k,j,(il+i-1));
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ReflectOuterX1(
//!              Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//! \brief REFLECTING boundary conditions, outer x1 boundary

void FaceCenteredBoundaryVariable::ReflectOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {
// copy face-centered magnetic fields into ghost zones, reflecting b1
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_fc).x1f(k,j,(iu+i+1)) = -(*var_fc).x1f(k,j,(iu-i+1));  // reflect 1-field
      }
    }
  }

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_fc).x2f(k,j,(iu+i  )) =  (*var_fc).x2f(k,j,(iu-i+1));
      }
    }
  }

  for (int k=kl; k<=ku+1; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_fc).x3f(k,j,(iu+i  )) =  (*var_fc).x3f(k,j,(iu-i+1));
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ReflectInnerX2(
//!              Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//! \brief REFLECTING boundary conditions, inner x2 boundary

void FaceCenteredBoundaryVariable::ReflectInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  // copy face-centered magnetic fields into ghost zones, reflecting b2
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        (*var_fc).x1f(k,(jl-j),i) =  (*var_fc).x1f(k,(jl+j-1),i);
      }
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x2f(k,(jl-j),i) = -(*var_fc).x2f(k,(jl+j  ),i);  // reflect 2-field
      }
    }
  }
  for (int k=kl; k<=ku+1; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x3f(k,(jl-j),i) =  (*var_fc).x3f(k,(jl+j-1),i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ReflectOuterX2(
//!              Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//! \brief REFLECTING boundary conditions, outer x2 boundary

void FaceCenteredBoundaryVariable::ReflectOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  // copy face-centered magnetic fields into ghost zones, reflecting b2
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        (*var_fc).x1f(k,(ju+j  ),i) =  (*var_fc).x1f(k,(ju-j+1),i);
      }
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x2f(k,(ju+j+1),i) = -(*var_fc).x2f(k,(ju-j+1),i);  // reflect 2-field
      }
    }
  }
  for (int k=kl; k<=ku+1; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x3f(k,(ju+j  ),i) =  (*var_fc).x3f(k,(ju-j+1),i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ReflectInnerX3(
//!              Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh)
//! \brief REFLECTING boundary conditions, inner x3 boundary

void FaceCenteredBoundaryVariable::ReflectInnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {
  // copy face-centered magnetic fields into ghost zones, reflecting b3
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        (*var_fc).x1f((kl-k),j,i) =  (*var_fc).x1f((kl+k-1),j,i);
      }
    }
  }
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x2f((kl-k),j,i) =  (*var_fc).x2f((kl+k-1),j,i);
      }
    }
  }
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x3f((kl-k),j,i) = -(*var_fc).x3f((kl+k  ),j,i);  // reflect 3-field
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::ReflectOuterX3(
//!              Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh)
//! \brief REFLECTING boundary conditions, outer x3 boundary

void FaceCenteredBoundaryVariable::ReflectOuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {
  // copy face-centered magnetic fields into ghost zones, reflecting b3
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        (*var_fc).x1f((ku+k  ),j,i) =  (*var_fc).x1f((ku-k+1),j,i);
      }
    }
  }
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x2f((ku+k  ),j,i) =  (*var_fc).x2f((ku-k+1),j,i);
      }
    }
  }
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x3f((ku+k+1),j,i) = -(*var_fc).x3f((ku-k+1),j,i);  // reflect 3-field
      }
    }
  }
  return;
}
