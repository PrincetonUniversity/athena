//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file outflow.cpp
//  \brief implementation of outflow BCs in each dimension

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"

//----------------------------------------------------------------------------------------
//! \fn void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions, inner x1 boundary

void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        prim(n,k,j,il-i) = prim(n,k,j,il);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x1f(k,j,(il-i)) = b.x1f(k,j,il);
      }
    }}

    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x2f(k,j,(il-i)) = b.x2f(k,j,il);
      }
    }}

    for (int k=kl; k<=ku+1; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x3f(k,j,(il-i)) = b.x3f(k,j,il);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                         FaceField &b, Real time, Real dt,
//                         int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions, outer x1 boundary

void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        prim(n,k,j,iu+i) = prim(n,k,j,iu);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x1f(k,j,(iu+i+1)) = b.x1f(k,j,(iu+1));
      }
    }}

    for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x2f(k,j,(iu+i)) = b.x2f(k,j,iu);
      }
    }}

    for (int k=kl; k<=ku+1; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        b.x3f(k,j,(iu+i)) = b.x3f(k,j,iu);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions, inner x2 boundary

void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        prim(n,k,jl-j,i) = prim(n,k,jl,i);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        b.x1f(k,(jl-j),i) = b.x1f(k,jl,i);
      }
    }}

    for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x2f(k,(jl-j),i) = b.x2f(k,jl,i);
      }
    }}

    for (int k=kl; k<=ku+1; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x3f(k,(jl-j),i) = b.x3f(k,jl,i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions, outer x2 boundary

void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        prim(n,k,ju+j,i) = prim(n,k,ju,i);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        b.x1f(k,(ju+j  ),i) = b.x1f(k,(ju  ),i);
      }
    }}

    for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x2f(k,(ju+j+1),i) = b.x2f(k,(ju+1),i);
      }
    }}

    for (int k=kl; k<=ku+1; ++k) {
    for (int j=1; j<=ngh; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x3f(k,(ju+j  ),i) = b.x3f(k,(ju  ),i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions, inner x3 boundary

void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        prim(n,kl-k,j,i) = prim(n,kl,j,i);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        b.x1f((kl-k),j,i) = b.x1f(kl,j,i);
      }
    }}

    for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju+1; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x2f((kl-k),j,i) = b.x2f(kl,j,i);
      }
    }}

    for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x3f((kl-k),j,i) = b.x3f(kl,j,i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
//                          FaceField &b, Real time, Real dt,
//                          int il, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions, outer x3 boundary

void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        prim(n,ku+k,j,i) = prim(n,ku,j,i);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu+1; ++i) {
        b.x1f((ku+k  ),j,i) = b.x1f((ku  ),j,i);
      }
    }}

    for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x2f((ku+k  ),j,i) = b.x2f((ku  ),j,i);
      }
    }}

    for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        b.x3f((ku+k+1),j,i) = b.x3f((ku+1),j,i);
      }
    }}
  }

  return;
}
