//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reflect.cpp
//  \brief implementation of reflecting BCs in each dimension

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../nr_radiation/radiation.hpp"
#include "./bvals_rad.hpp"


// The angular grid can change
// The angular octant is
//   1  |  0       5  |  4
//   -------      ---------
//   3  |  2       7  |  6
// in radiatin class, n_ang is angles per octant, noct is the number of octant
// radiation relfection means set specific intensitiy corresponding
// the value in the opposite octant

// Temporary function to copy intensity
void CopyIntensity(Real *iri, Real *iro, int li, int lo, int n_ang) {
  // here ir is only intensity for each cell and each frequency band
  for (int n=0; n<n_ang; ++n) {
    int angi = li * n_ang + n;
    int ango = lo * n_ang + n;
    iro[angi] = iri[ango];
    iro[ango] = iri[angi];
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectInnerX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x1 boundary

void RadBoundaryVariable::ReflectInnerX1(Real time, Real dt, int il, int jl, int ju,
                                         int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  const int& noct = pmy_block_->pnrrad->noct;
  int n_ang = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *iri = &var(k,j,(il+i-1),ifr*pmy_block_->pnrrad->nang);
          Real *iro = &var(k,j, il-i, ifr*pmy_block_->pnrrad->nang);
          CopyIntensity(iri, iro, 0, 1, n_ang);
          if (noct > 2) {
            CopyIntensity(iri, iro, 2, 3, n_ang);
          }
          if (noct > 3) {
            CopyIntensity(iri, iro, 4, 5, n_ang);
            CopyIntensity(iri, iro, 6, 7, n_ang);
          }
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectOuterX1(
//          Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x1 boundary

void RadBoundaryVariable::ReflectOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  const int& noct = pmy_block_->pnrrad->noct;
  int n_ang = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *iri = &var(k,j,(iu-i+1),ifr*pmy_block_->pnrrad->nang);
          Real *iro = &var(k,j, iu+i, ifr*pmy_block_->pnrrad->nang);
          CopyIntensity(iri, iro, 0, 1, n_ang);
          if (noct > 2) {
            CopyIntensity(iri, iro, 2, 3, n_ang);
          }
          if (noct > 3) {
            CopyIntensity(iri, iro, 4, 5, n_ang);
            CopyIntensity(iri, iro, 6, 7, n_ang);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectInnerX2(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x2 boundary

void RadBoundaryVariable::ReflectInnerX2(Real time, Real dt, int il, int iu, int jl,
                                         int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx

  const int& noct = pmy_block_->pnrrad->noct;
  int n_ang = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *iri = &var(k,jl+j-1,i,ifr*pmy_block_->pnrrad->nang);
          Real *iro = &var(k,jl-j,i, ifr*pmy_block_->pnrrad->nang);
          CopyIntensity(iri, iro, 0, 2, n_ang);
          CopyIntensity(iri, iro, 1, 3, n_ang);

          if (noct > 3) {
            CopyIntensity(iri, iro, 4, 6, n_ang);
            CopyIntensity(iri, iro, 5, 7, n_ang);
          }
        }
      }
    }
  }
  return;
}



//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectOuterX2(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x2 boundary

void RadBoundaryVariable::ReflectOuterX2(Real time, Real dt, int il, int iu, int ju,
                                         int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx

  const int& noct = pmy_block_->pnrrad->noct;
  int n_ang = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *iri = &var(k,ju-j+1,i,ifr*pmy_block_->pnrrad->nang);
          Real *iro = &var(k,ju+j,i, ifr*pmy_block_->pnrrad->nang);
          CopyIntensity(iri, iro, 0, 2, n_ang);
          CopyIntensity(iri, iro, 1, 3, n_ang);

          if (noct > 3) {
            CopyIntensity(iri, iro, 4, 6, n_ang);
            CopyIntensity(iri, iro, 5, 7, n_ang);
          }
        }
      }
    }
  }
  return;
}




//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectInnerX3(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x3 boundary

void RadBoundaryVariable::ReflectInnerX3(Real time, Real dt, int il, int iu, int jl,
                                         int ju, int kl, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx

  const int& noct = pmy_block_->pnrrad->noct;
  int n_ang = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *iri = &var(kl+k-1,j,i,ifr*pmy_block_->pnrrad->nang);
          Real *iro = &var(kl-k,j,i, ifr*pmy_block_->pnrrad->nang);
          CopyIntensity(iri, iro, 0, 4, n_ang);
          CopyIntensity(iri, iro, 1, 5, n_ang);
          CopyIntensity(iri, iro, 2, 6, n_ang);
          CopyIntensity(iri, iro, 3, 7, n_ang);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::ReflectOuterX3(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, outer x3 boundary

void RadBoundaryVariable::ReflectOuterX3(Real time, Real dt, int il, int iu, int jl,
                                         int ju, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  const int& noct = pmy_block_->pnrrad->noct;
  int n_ang = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *iri = &var(ku-k+1,j,i,ifr*pmy_block_->pnrrad->nang);
          Real *iro = &var(ku+k,j,i, ifr*pmy_block_->pnrrad->nang);
          CopyIntensity(iri, iro, 0, 4, n_ang);
          CopyIntensity(iri, iro, 1, 5, n_ang);
          CopyIntensity(iri, iro, 2, 6, n_ang);
          CopyIntensity(iri, iro, 3, 7, n_ang);
        }
      }
    }
  }
  return;
}
