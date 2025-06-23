//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file vacuum_rad.cpp
//  \brief implementation of vacuum BCs in each dimension

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

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::VacuumInnerX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief VACUUM boundary conditions, inner x1 boundary

void RadBoundaryVariable::VacuumInnerX1(Real time, Real dt, int il, int jl, int ju,
                                        int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  MeshBlock *pmb=pmy_block_;
  const int& nang = pmb->pnrrad->nang; // angles per octant
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for(int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
            const Real& miux=pmb->pnrrad->mu(0,k,j,il,n);
            if (miux < 0.0) {
              (*var_cc)(k,j,il-i,ang) = (*var_cc)(k,j,il,ang);
            } else {
              (*var_cc)(k,j,il-i,ang) = 0.0;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::VacuumOuterX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief VACUUM boundary conditions, outer x1 boundary

void RadBoundaryVariable::VacuumOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {

  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  MeshBlock *pmb=pmy_block_;
  const int& nang = pmb->pnrrad->nang; // angles per octant
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for(int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
            const Real& miux=pmb->pnrrad->mu(0,k,j,iu,n);
            if (miux > 0.0) {
              (*var_cc)(k,j,iu+i,ang) = (*var_cc)(k,j,iu,ang);
            } else {
              (*var_cc)(k,j,iu+i,ang) = 0.0;
            }
          }
        }
      }
    }
  }
  return;
}





//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::VacuumInnerX2(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief REFLECTING boundary conditions, inner x2 boundary

void RadBoundaryVariable::VacuumInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  MeshBlock *pmb=pmy_block_;
  const int& nang = pmb->pnrrad->nang; // angles per octant
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for(int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
            const Real& miuy=pmb->pnrrad->mu(1,k,jl,i,n);
            if (miuy < 0.0) {
              (*var_cc)(k,jl-j,i,ang) = (*var_cc)(k,jl,i,ang);
            } else {
              (*var_cc)(k,jl-j,i,ang) = 0.0;
            }
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::VacuumOuterX2(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief VACUUM boundary conditions, outer x2 boundary

void RadBoundaryVariable::VacuumOuterX2(
    Real time, Real dt, int il, int iu, int ju,  int kl, int ku, int ngh) {

  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  MeshBlock *pmb=pmy_block_;
  const int& nang = pmb->pnrrad->nang; // angles per octant
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for(int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
            const Real& miuy=pmb->pnrrad->mu(1,k,ju,i,n);
            if (miuy > 0.0) {
              (*var_cc)(k,ju+j,i,ang) = (*var_cc)(k,ju,i,ang);
            } else {
              (*var_cc)(k,ju+j,i,ang) = 0.0;
            }
          }
        }
      }
    }
  }
  return;
}






//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::VacuumInnerX2(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief VACUUM boundary conditions, inner x2 boundary

void RadBoundaryVariable::VacuumInnerX3(
    Real time, Real dt, int il, int iu, int jl,  int ju, int kl, int ngh) {

  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  MeshBlock *pmb=pmy_block_;
  const int& nang = pmb->pnrrad->nang; // angles per octant
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for(int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
            const Real& miuz=pmb->pnrrad->mu(2,kl,j,i,n);
            if (miuz < 0.0) {
              (*var_cc)(kl-k,j,i,ang) = (*var_cc)(kl,j,i,ang);
            } else {
              (*var_cc)(kl-k,j,i,ang) = 0.0;
            }
          }
        }
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void RadBoundaryVariable::VacuumOuterX3(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief VACUUM boundary conditions, outer X3 boundary

void RadBoundaryVariable::VacuumOuterX3(
    Real time, Real dt, int il, int iu, int jl,  int ju, int ku, int ngh) {

  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx
  MeshBlock *pmb=pmy_block_;
  const int& nang = pmb->pnrrad->nang; // angles per octant
  const int& nfreq = pmb->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for(int n=0; n<nang; ++n) {
            int ang=ifr*nang+n;
            const Real& miuz=pmb->pnrrad->mu(2,ku,j,i,n);
            if (miuz > 0.0) {
              (*var_cc)(ku+k,j,i,ang) = (*var_cc)(ku,j,i,ang);
            } else {
              (*var_cc)(ku+k,j,i,ang) = 0.0;
            }
          }
        }
      }
    }
  }
  return;
}
