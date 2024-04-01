//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file rotate_rad.cpp
//  \brief implementation of rotating radiation BCs in each dimension
//======================================================================================

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../nr_radiation/radiation.hpp"
#include "./bvals_rad.hpp"

// The angular octant ( in x-y plane) is
//   1  |  0       5  |  4
//   -------      ---------
//   3  |  2       7  |  6
// in X-Z plane, it is
//   5  |  4       7  |  6
//   -------      ---------
//   1  |  0       3  |  2


// in radiatin class, n_ang is angles per octant, noct is the number of octant


// Temporary function to copy intensity
// swap 0 and 3, 1 and 2
// change sign of both Fx and Fy
// only do this for 2D spherical polar
void CopyIntensity2(Real *ir, int n_ang, int direction) {
  // here ir is only intensity for each cell and each frequency band
  for (int n=0; n<n_ang; ++n) {
    // from 0 to 4, 4 to 5, 5 to 1, 1 to 0
    int ang1 = 0 * n_ang + n;
    int ang2 = 3 * n_ang + n;

    int ang3 = 1 * n_ang + n;
    int ang4 = 2 * n_ang + n;

    Real temp = ir[ang1];
    ir[ang1] = ir[ang2];
    ir[ang2] = temp;


    temp = ir[ang3];
    ir[ang3] = ir[ang4];
    ir[ang4] = temp;
  }
  return;
}

void CopyIntensity3(Real *ir, int n_ang, int direction) {
  // here ir is only intensity for each cell and each frequency band
  for (int n=0; n<n_ang; ++n) {
    // switch 1 and 2, switch 3 and 0
    int ang1 = 0 * n_ang + n;
    int ang2 = 3 * n_ang + n;
    int ang3 = 1 * n_ang + n;
    int ang4 = 2 * n_ang + n;
    Real temp = ir[ang1];
    ir[ang1] = ir[ang2];
    ir[ang2] = temp;

    temp = ir[ang4];
    ir[ang4] = ir[ang3];
    ir[ang3] = temp;
  }
  for (int n=0; n<n_ang; ++n) {
    // switch 5 and 6, switch 4 and 7
    int ang1 = 5 * n_ang + n;
    int ang2 = 6 * n_ang + n;
    int ang3 = 4 * n_ang + n;
    int ang4 = 7 * n_ang + n;
    Real temp = ir[ang1];
    ir[ang1] = ir[ang2];
    ir[ang2] = temp;
    temp = ir[ang4];
    ir[ang4] = ir[ang3];
    ir[ang3] = temp;
  }
  return;
}

// The angular octant ( in x-y plane) is
//   1  |  0       5  |  4
//   -------      ---------
//   3  |  2       7  |  6

void CopyIntensity4(Real *ir, int n_ang, int direction) {
  // here ir is only intensity for each cell and each frequency band
  for (int n=0; n<n_ang; ++n) {
    // from 0 to 1, 1 to 3, 3 to 2, 2 to 0
    int ang1 = 0 * n_ang + n;
    int ang2 = 1 * n_ang + n;
    int ang3 = 3 * n_ang + n;
    int ang4 = 2 * n_ang + n;
    Real temp = ir[ang1];
    ir[ang1] = ir[ang4];
    ir[ang4] = ir[ang3];
    ir[ang3] = ir[ang2];
    ir[ang2] = temp;
  }
  for (int n=0; n<n_ang; ++n) {
    // from 4 to 5, 5 to 7, 7 to 6, 6 to 1
    int ang1 = 4 * n_ang + n;
    int ang2 = 5 * n_ang + n;
    int ang3 = 7 * n_ang + n;
    int ang4 = 6 * n_ang + n;

    Real temp = ir[ang1];
    ir[ang1] = ir[ang4];
    ir[ang4] = ir[ang3];
    ir[ang3] = ir[ang2];
    ir[ang2] = temp;
  }
  return;
}

void CopyIntensity5(Real *ir, int n_ang, int direction) {
  // here ir is only intensity for each cell and each frequency band
  for (int n=0; n<n_ang; ++n) {
    // switch 0 and 4, switch 5 and 1
    int ang1 = 0 * n_ang + n;
    int ang2 = 4 * n_ang + n;
    int ang3 = 5 * n_ang + n;
    int ang4 = 1 * n_ang + n;
    Real temp = ir[ang1];
    ir[ang1] = ir[ang2];
    ir[ang2] = temp;

    temp = ir[ang4];
    ir[ang4] = ir[ang3];
    ir[ang3] = temp;
  }
  for (int n=0; n<n_ang; ++n) {
    // from 4 to 5, 5 to 7, 7 to 6, 6 to 4
    int ang1 = 2 * n_ang + n;
    int ang2 = 6 * n_ang + n;
    int ang3 = 7 * n_ang + n;
    int ang4 = 3 * n_ang + n;
    Real temp = ir[ang1];
    ir[ang1] = ir[ang2];
    ir[ang2] = temp;

    temp = ir[ang4];
    ir[ang4] = ir[ang3];
    ir[ang3] = temp;
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_InnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x2, inner x2 boundary by Pi/2

// anti-clockwise rotation
// This function is used after periodic copy

void RadBoundaryVariable::RotateHPi_InnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {
  // copy radiation variables into ghost zones
  const int& noct = pmy_block_->pnrrad->noct;
  int n_ang = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *ir = &var(k,jl-j,i,ifr*pmy_block_->pnrrad->nang);
          CopyIntensity5(ir, n_ang, 1);
        }
      }
    }
  }
  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_OuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  ROTATE boundary radiation conditions, outer x2 boundary

void RadBoundaryVariable::RotateHPi_OuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {
  const int& noct = pmy_block_->pnrrad->noct;
  int n_ang =pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=kl; k<=ku; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *ir = &var(k,ju+j,i,ifr*pmy_block_->pnrrad->nang);
          CopyIntensity5(ir, n_ang, -1);
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_InnerX3((MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                     int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x3, inner x3 boundary by Pi/2

// anti-clockwise rotation
// This function is used after periodic copy
void RadBoundaryVariable::RotateHPi_InnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {
  // copy radiation variables into ghost zones
  const int& noct = pmy_block_->pnrrad->noct;
  int n_ang = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq =pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *ir = &var(kl-k,j,i,ifr*pmy_block_->pnrrad->nang);
          CopyIntensity4(ir, n_ang, 1);
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_OuterX3((MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                     int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x3, outer x3 boundary by Pi/2

// anti-clockwise rotation
// This function is used after periodic copy
void RadBoundaryVariable::RotateHPi_OuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {
  // copy radiation variables into ghost zones
  const int& noct = pmy_block_->pnrrad->noct;
  int n_ang = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *ir = &var(ku+k,j,i,ifr*pmy_block_->pnrrad->nang);
          CopyIntensity4(ir, n_ang, 1);
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RotatePi_InnerX3((MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                     int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x3, inner x3 boundary by Pi

// anti-clockwise rotation
// This function is used after periodic copy

void RadBoundaryVariable::RotatePi_InnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {
  // copy radiation variables into ghost zones
  const int& noct = pmy_block_->pnrrad->noct;
  int n_ang = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *ir = &var(kl-k,j,i,ifr*pmy_block_->pnrrad->nang);
          CopyIntensity3(ir, n_ang, 1);
        }
      }
    }
  }
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RotatePi_OuterX3((MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                     int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x3, outer x3 boundary by Pi

// anti-clockwise rotation
// This function is used after periodic copy
void RadBoundaryVariable::RotatePi_OuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {
  // copy radiation variables into ghost zones
  const int& noct = pmy_block_->pnrrad->noct;
  int n_ang = pmy_block_->pnrrad->nang/noct; // angles per octant
  const int& nfreq = pmy_block_->pnrrad->nfreq; // number of frequency bands

  for (int k=1; k<=ngh; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        for (int ifr=0; ifr<nfreq; ++ifr) {
          AthenaArray<Real> &var = *var_cc;
          Real *ir = &var(ku+k,j,i,ifr*pmy_block_->pnrrad->nang);
          CopyIntensity3(ir, n_ang, 1);
        }
      }
    }
  }
  return;
}
