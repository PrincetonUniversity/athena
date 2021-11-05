//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file outflow_fc.cpp
//! \brief implementation of outflow BCs in each dimension

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "bvals_fc.hpp"

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::OutflowInnerX1(
//!              Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//! \brief OUTFLOW boundary conditions, inner x1 boundary

void FaceCenteredBoundaryVariable::OutflowInnerX1(
    Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh) {
  // copy face-centered magnetic fields into ghost zones
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_fc).x1f(k,j,(il-i)) = (*var_fc).x1f(k,j,il);
      }
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_fc).x2f(k,j,(il-i)) = (*var_fc).x2f(k,j,il);
      }
    }
  }
  for (int k=kl; k<=ku+1; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_fc).x3f(k,j,(il-i)) = (*var_fc).x3f(k,j,il);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::OutflowOuterX1(
//!              Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//! \brief OUTFLOW boundary conditions, outer x1 boundary

void FaceCenteredBoundaryVariable::OutflowOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy face-centered magnetic fields into ghost zones
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_fc).x1f(k,j,(iu+i+1)) = (*var_fc).x1f(k,j,(iu+1));
      }
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_fc).x2f(k,j,(iu+i)) = (*var_fc).x2f(k,j,iu);
      }
    }
  }
  for (int k=kl; k<=ku+1; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_fc).x3f(k,j,(iu+i)) = (*var_fc).x3f(k,j,iu);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::OutflowInnerX2(
//!              Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//! \brief OUTFLOW boundary conditions, inner x2 boundary

void FaceCenteredBoundaryVariable::OutflowInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  // copy face-centered magnetic fields into ghost zones
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        (*var_fc).x1f(k,(jl-j),i) = (*var_fc).x1f(k,jl,i);
      }
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x2f(k,(jl-j),i) = (*var_fc).x2f(k,jl,i);
      }
    }
  }
  for (int k=kl; k<=ku+1; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x3f(k,(jl-j),i) = (*var_fc).x3f(k,jl,i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::OutflowOuterX2(
//!              Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//! \brief OUTFLOW boundary conditions, outer x2 boundary

void FaceCenteredBoundaryVariable::OutflowOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  // copy face-centered magnetic fields into ghost zones
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        (*var_fc).x1f(k,(ju+j  ),i) = (*var_fc).x1f(k,(ju  ),i);
      }
    }
  }
  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x2f(k,(ju+j+1),i) = (*var_fc).x2f(k,(ju+1),i);
      }
    }
  }
  for (int k=kl; k<=ku+1; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x3f(k,(ju+j  ),i) = (*var_fc).x3f(k,(ju  ),i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::OutflowInnerX3(
//!              Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh)
//! \brief OUTFLOW boundary conditions, inner x3 boundary

void FaceCenteredBoundaryVariable::OutflowInnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {
  // copy face-centered magnetic fields into ghost zones
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        (*var_fc).x1f((kl-k),j,i) = (*var_fc).x1f(kl,j,i);
      }
    }
  }
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x2f((kl-k),j,i) = (*var_fc).x2f(kl,j,i);
      }
    }
  }
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x3f((kl-k),j,i) = (*var_fc).x3f(kl,j,i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FaceCenteredBoundaryVariable::OutflowOuterX3(
//!              Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh)
//! \brief OUTFLOW boundary conditions, outer x3 boundary

void FaceCenteredBoundaryVariable::OutflowOuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {
  // copy face-centered magnetic fields into ghost zones
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        (*var_fc).x1f((ku+k  ),j,i) = (*var_fc).x1f((ku  ),j,i);
      }
    }
  }
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x2f((ku+k  ),j,i) = (*var_fc).x2f((ku  ),j,i);
      }
    }
  }
  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        (*var_fc).x3f((ku+k+1),j,i) = (*var_fc).x3f((ku+1),j,i);
      }
    }
  }
  return;
}
