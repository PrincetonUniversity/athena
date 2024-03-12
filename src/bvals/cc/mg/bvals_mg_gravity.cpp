//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_mg_gravity.cpp
//! \brief

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy()
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../globals.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../multigrid/multigrid.hpp"
#include "../../../parameter_input.hpp"
#include "../../../utils/buffer_utils.hpp"
#include "bvals_mg.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class Multigrid;
class MultigridDriver;
class MGBoundaryValues;

//----------------------------------------------------------------------------------------
//! \fn int MGGravityBoundaryValues::LoadMultigridBoundaryBufferToCoarserFluxCons(
//!                                               Real *buf, const NeighborBlock& nb)
//! \brief Set Multigrid boundary buffers for sending to a block on the coarser level
//!        using the Mass Conservation formula for gravity

int MGGravityBoundaryValues::LoadMultigridBoundaryBufferToCoarserFluxCons(Real *buf,
                                                          const NeighborBlock& nb) {
  const AthenaArray<Real> &u = pmy_mg_->GetCurrentData();
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int cn = ngh, fs = ngh, cs = cn, fe = fs + nc - 1, ce = cs + nc/2 - 1;
  int p = 0;

  // take averages of four cells contacting the interface and pack
  if (nb.ni.ox1 != 0) { // x1 face
    int fi;
    if (nb.ni.ox1 < 0) fi = fs;
    else               fi = fe;
    for (int v=0; v<nvar; ++v) {
      for (int fk=fs; fk<=fe; fk+=2) {
#pragma ivdep
        for (int fj=fs; fj<=fe; fj+=2)
          buf[p++] = 0.25*((u(v, fk,   fj,   fi)+u(v, fk,   fj+1, fi))
                          +(u(v, fk+1, fj,   fi)+u(v, fk+1, fj+1, fi)));
      }
    }
  } else if (nb.ni.ox2 != 0) { // x2 face
    int fj;
    if (nb.ni.ox2 < 0) fj = fs;
    else               fj = fe;
    for (int v=0; v<nvar; ++v) {
      for (int fk=fs; fk<=fe; fk+=2) {
#pragma ivdep
        for (int fi=fs; fi<=fe; fi+=2)
          buf[p++] = 0.25*((u(v, fk,   fj, fi)+u(v, fk,   fj, fi+1))
                          +(u(v, fk+1, fj, fi)+u(v, fk+1, fj, fi+1)));
      }
    }
  } else { // x3 face
    int fk;
    if (nb.ni.ox3 < 0) fk = fs;
    else               fk = fe;
    for (int v=0; v<nvar; ++v) {
      for (int fj=fs; fj<=fe; fj+=2) {
#pragma ivdep
        for (int fi=fs; fi<=fe; fi+=2)
          buf[p++] = 0.25*((u(v, fk, fj,   fi)+u(v, fk, fj,   fi+1))
                          +(u(v, fk, fj+1, fi)+u(v, fk, fj+1, fi+1)));
      }
    }
  }

  return p;
}


//----------------------------------------------------------------------------------------
//! \fn int MGGravityBoundaryValues::LoadMultigridBoundaryBufferToFinerFluxCons(Real *buf,
//!                                                             const NeighborBlock& nb)
//! \brief Set Multigrid boundary buffers for sending to a block on the finer level
int MGGravityBoundaryValues::LoadMultigridBoundaryBufferToFinerFluxCons(Real *buf,
                                                               const NeighborBlock& nb) {
  const AthenaArray<Real> &u = pmy_mg_->GetCurrentData();
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int cn = ngh - 1, fs = ngh, fe = fs + nc - 1;
  int si = (nb.ni.ox1 > 0) ? (fe - cn) : fs;
  int ei = (nb.ni.ox1 < 0) ? (fs + cn) : fe;
  int sj = (nb.ni.ox2 > 0) ? (fe - cn) : fs;
  int ej = (nb.ni.ox2 < 0) ? (fs + cn) : fe;
  int sk = (nb.ni.ox3 > 0) ? (fe - cn) : fs;
  int ek = (nb.ni.ox3 < 0) ? (fs + cn) : fe;

  int nred = nc/2 - ngh;
  if (nb.ni.ox1 == 0) {
    if (nb.ni.fi1 == 1)   si += nred;
    else                  ei -= nred;
  }
  if (nb.ni.ox2 == 0) {
    if (nb.ni.ox1 != 0) {
      if (nb.ni.fi1 == 1) sj += nred;
      else                ej -= nred;
    } else {
      if (nb.ni.fi2 == 1) sj += nred;
      else                ej -= nred;
    }
  }
  if (nb.ni.ox3 == 0) {
    if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
      if (nb.ni.fi1 == 1) sk += nred;
      else                ek -= nred;
    } else {
      if (nb.ni.fi2 == 1) sk += nred;
      else                ek -= nred;
    }
  }

  int p = 0;
  BufferUtility::PackData(u, buf, 0, nvar-1, si, ei, sj, ej, sk, ek, p);

  return p;
}


//----------------------------------------------------------------------------------------
//! \fn void MGGravityBoundaryValues::SetMultigridBoundaryFromCoarserFluxCons(
//!                                      const Real *buf, const NeighborBlock& nb)
//! \brief Set hydro boundary received from a block on the same level

void MGGravityBoundaryValues::SetMultigridBoundaryFromCoarserFluxCons(const Real *buf,
                                                             const NeighborBlock& nb) {
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int si, sj, sk, ei, ej, ek;
  int cng = 1;
  int cs = cng, ce = cng + nc/2 -1;

  if (nb.ni.ox1 == 0) {
    si = cs, ei = ce;
    if ((loc.lx1 & 1LL) == 0LL) ei += cng;
    else                        si -= cng;
  } else if (nb.ni.ox1 > 0)  {
    si = ce + 1,   ei = ce + cng;
  } else {
    si = cs - cng, ei = cs - 1;
  }
  if (nb.ni.ox2 == 0) {
    sj = cs, ej = ce;
    if ((loc.lx2 & 1LL) == 0LL) ej += cng;
    else                        sj -= cng;
  } else if (nb.ni.ox2 > 0) {
    sj = ce + 1,   ej = ce + cng;
  } else {
    sj = cs - cng, ej = cs - 1;
  }
  if (nb.ni.ox3 == 0) {
    sk = cs, ek = ce;
    if ((loc.lx3 & 1LL) == 0LL) ek += cng;
    else                        sk -= cng;
  } else if (nb.ni.ox3 > 0)  {
    sk = ce + 1,   ek = ce + cng;
  } else {
    sk = cs - cng, ek = cs - 1;
  }

  int p = 0;
  BufferUtility::UnpackData(buf, cbuf_, 0, nvar-1, si, ei, sj, ej, sk, ek, p);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGGravityBoundaryValues::SetMultigridBoundaryFromFinerFluxCons(
//!                                      const Real *buf, const NeighborBlock& nb)

void MGGravityBoundaryValues::SetMultigridBoundaryFromFinerFluxCons(const Real *buf,
                                                             const NeighborBlock& nb) {
  AthenaArray<Real> &dst = pmy_mg_->GetCurrentData();
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  int fs = ngh, fe = fs + nc - 1;
  int si, sj, sk, ei, ej, ek;
  int oi = 0, oj = 0, ok = 0;

  // could reuse the index calculation for normal boundaries
  // but in general this depends on physics
  if (nb.ni.ox1 == 0) {
    si = fs, ei = fe;
    if (nb.ni.fi1 == 1)   si += nc/2;
    else                  ei -= nc/2;
  } else if (nb.ni.ox1 > 0) {
    si = fe + 1,   ei = fe + ngh, oi = -1;
  } else {
    si = fs - ngh, ei = fs - 1,   oi = 1;
  }
  if (nb.ni.ox2 == 0) {
    sj = fs, ej = fe;
    if (nb.ni.ox1 != 0) {
      if (nb.ni.fi1 == 1) sj += nc/2;
      else                ej -= nc/2;
    } else {
      if (nb.ni.fi2 == 1) sj += nc/2;
      else                ej -= nc/2;
    }
  } else if (nb.ni.ox2 > 0) {
    sj = fe + 1,   ej = fe + ngh, oj = -1;
  } else {
    sj = fs - ngh, ej = fs - 1,   oj = 1;
  }
  if (nb.ni.ox3 == 0) {
    sk = fs, ek = fe;
    if (nb.ni.ox1 != 0 && nb.ni.ox2 != 0) {
      if (nb.ni.fi1 == 1) sk += nc/2;
      else                ek -= nc/2;
    } else {
      if (nb.ni.fi2 == 1) sk += nc/2;
      else                ek -= nc/2;
    }
  } else if (nb.ni.ox3 > 0) {
    sk = fe + 1,   ek = fe + ngh, ok = -1;
  } else {
    sk = fs - ngh, ek = fs - 1,   ok = 1;
  }

  constexpr Real ot = 1.0/3.0;
  int p = 0;

  // correct the ghost values using the mass conservation formula
  for (int v=0; v<nvar; ++v) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
        for (int i=si; i<=ei; ++i) {
          dst(v,k,j,i) = ot * (4.0*buf[p++] - dst(v,k+ok,j+oj,i+oi));
        }
      }
    }
  }
  return;
}



//----------------------------------------------------------------------------------------
//! \fn void MGGravityBoundaryValues::ProlongateMultigridBoundariesFluxCons()
//! \brief prolongate boundaries for Multigrid

void MGGravityBoundaryValues::ProlongateMultigridBoundariesFluxCons() {
  const AthenaArray<Real> &u = pmy_mg_->GetCurrentData();
  AthenaArray<Real> &dst = pmy_mg_->GetCurrentData();
  int nc = pmy_mg_->GetCurrentNumberOfCells();
  int ngh = pmy_mg_->ngh_, nvar = pmy_mg_->nvar_;
  Real time = pmy_mesh_->time;
  int lev = pmy_mg_->GetCurrentLevel();
  int ll = pmy_mg_->nlevel_ - 1 - lev;
  MGCoordinates &coord = pmy_mg_->ccoord_[lev];
  constexpr Real ot = 1.0 / 3.0;
  int cn = 1, cs = cn, ce = cs + nc/2 -1;
  int fs = ngh, fe = fs + nc - 1;
  int flim = ngh + nc;

  // x1face
  for (int ox1=-1; ox1<=1; ox1+=2) {
    if (nblevel[1][1][ox1+1] == loc.level - 1) {
      int i, fi, fig;
      if (ox1 > 0) i = ce + 1, fi = fe, fig = fe + 1;
      else         i = cs - 1, fi = fs, fig = fs - 1;
      if (apply_bndry_fn_[BoundaryFace::inner_x2])
        DispatchBoundaryFunction(BoundaryFace::inner_x2, cbuf_, time, nvar,
                                 i, i, cs, ce, cs, ce, ngh, coord, false);
      if (apply_bndry_fn_[BoundaryFace::outer_x2])
        DispatchBoundaryFunction(BoundaryFace::outer_x2, cbuf_, time, nvar,
                                 i, i, cs, ce, cs, ce, ngh, coord, false);
      if (apply_bndry_fn_[BoundaryFace::inner_x3])
        DispatchBoundaryFunction(BoundaryFace::inner_x3, cbuf_, time, nvar,
                                 i, i, cs, ce, cs, ce, ngh, coord, false);
      if (apply_bndry_fn_[BoundaryFace::outer_x3])
        DispatchBoundaryFunction(BoundaryFace::outer_x3, cbuf_, time, nvar,
                                 i, i, cs, ce, cs, ce, ngh, coord, false);
      for (int v=0; v<nvar; ++v) {
        for(int k=cs; k<=ce; ++k) {
          int fk = (k - cs) * 2 + fs;
          for(int j=cs; j<=ce; ++j) {
            int fj = (j - cs) * 2 + fs;
            Real ccval = cbuf_(v, k, j, i);
            Real gx2c = 0.125*(cbuf_(v, k, j+1, i) -  cbuf_(v, k, j-1, i));
            Real gx3c = 0.125*(cbuf_(v, k+1, j, i) -  cbuf_(v, k-1, j, i));
            dst(v,fk  ,fj  ,fig) = ot*(2.0*(ccval - gx2c - gx3c) + u(v,fk  ,fj  ,fi));
            dst(v,fk  ,fj+1,fig) = ot*(2.0*(ccval + gx2c - gx3c) + u(v,fk  ,fj+1,fi));
            dst(v,fk+1,fj  ,fig) = ot*(2.0*(ccval - gx2c + gx3c) + u(v,fk+1,fj  ,fi));
            dst(v,fk+1,fj+1,fig) = ot*(2.0*(ccval + gx2c + gx3c) + u(v,fk+1,fj+1,fi));
          }
        }
      }
    }
  }

  // x2face
  for (int ox2=-1; ox2<=1; ox2+=2) {
    if (nblevel[1][ox2+1][1] == loc.level - 1) {
      int j, fj, fjg;
      if (ox2 > 0) j = ce + 1, fj = fe, fjg = fe + 1;
      else         j = cs - 1, fj = fs, fjg = fs - 1;
      if (apply_bndry_fn_[BoundaryFace::inner_x1])
        DispatchBoundaryFunction(BoundaryFace::inner_x1, cbuf_, time, nvar,
                                 cs, ce, j, j, cs, ce, ngh, coord, false);
      if (apply_bndry_fn_[BoundaryFace::outer_x1])
        DispatchBoundaryFunction(BoundaryFace::outer_x1, cbuf_, time, nvar,
                                 cs, ce, j, j, cs, ce, ngh, coord, false);
      if (apply_bndry_fn_[BoundaryFace::inner_x3])
        DispatchBoundaryFunction(BoundaryFace::inner_x3, cbuf_, time, nvar,
                                 cs, ce, j, j, cs, ce, ngh, coord, false);
      if (apply_bndry_fn_[BoundaryFace::outer_x3])
        DispatchBoundaryFunction(BoundaryFace::outer_x3, cbuf_, time, nvar,
                                 cs, ce, j, j, cs, ce, ngh, coord, false);
      for (int v=0; v<nvar; ++v) {
        for(int k=cs; k<=ce; ++k) {
          int fk = (k - cs) * 2 + fs;
          for(int i=cs; i<=ce; ++i) {
            int fi = (i - cs) * 2 + fs;
            Real ccval = cbuf_(v, k, j, i);
            Real gx1c = 0.125*(cbuf_(v, k, j, i+1) -  cbuf_(v, k, j, i-1));
            Real gx3c = 0.125*(cbuf_(v, k+1, j, i) -  cbuf_(v, k-1, j, i));
            dst(v,fk  ,fjg,fi  ) = ot*(2.0*(ccval - gx1c - gx3c) + u(v,fk,  fj,fi  ));
            dst(v,fk  ,fjg,fi+1) = ot*(2.0*(ccval + gx1c - gx3c) + u(v,fk,  fj,fi+1));
            dst(v,fk+1,fjg,fi  ) = ot*(2.0*(ccval - gx1c + gx3c) + u(v,fk+1,fj,fi  ));
            dst(v,fk+1,fjg,fi+1) = ot*(2.0*(ccval + gx1c + gx3c) + u(v,fk+1,fj,fi+1));
          }
        }
      }
    }
  }

  // x3face
  for (int ox3=-1; ox3<=1; ox3+=2) {
    if (nblevel[ox3+1][1][1] == loc.level - 1) {
      int k, fk, fkg;
      if (ox3 > 0) k = ce + 1, fk = fe, fkg = fe + 1;
      else         k = cs - 1, fk = fs, fkg = fs - 1;
      if (apply_bndry_fn_[BoundaryFace::inner_x1])
        DispatchBoundaryFunction(BoundaryFace::inner_x1, cbuf_, time, nvar,
                                 cs, ce, cs, ce, k, k, ngh, coord, false);
      if (apply_bndry_fn_[BoundaryFace::outer_x1])
        DispatchBoundaryFunction(BoundaryFace::outer_x1, cbuf_, time, nvar,
                                 cs, ce, cs, ce, k, k, ngh, coord, false);
      if (apply_bndry_fn_[BoundaryFace::inner_x2])
        DispatchBoundaryFunction(BoundaryFace::inner_x2, cbuf_, time, nvar,
                                 cs, ce, cs, ce, k, k, ngh, coord, false);
      if (apply_bndry_fn_[BoundaryFace::outer_x2])
        DispatchBoundaryFunction(BoundaryFace::outer_x2, cbuf_, time, nvar,
                                 cs, ce, cs, ce, k, k, ngh, coord, false);
      for (int v=0; v<nvar; ++v) {
        for(int j=cs; j<=ce; ++j) {
          int fj = (j - cs) * 2 + fs;
          for(int i=cs; i<=ce; ++i) {
            int fi = (i - cs) * 2 + fs;
            Real ccval = cbuf_(v, k, j, i);
            Real gx1c = 0.125*(cbuf_(v, k, j, i+1) -  cbuf_(v, k, j, i-1));
            Real gx2c = 0.125*(cbuf_(v, k, j+1, i) -  cbuf_(v, k, j-1, i));
            dst(v,fkg,fj  ,fi  ) = ot*(2.0*(ccval - gx1c - gx2c) + u(v,fk,fj  ,fi  ));
            dst(v,fkg,fj  ,fi+1) = ot*(2.0*(ccval + gx1c - gx2c) + u(v,fk,fj  ,fi+1));
            dst(v,fkg,fj+1,fi  ) = ot*(2.0*(ccval - gx1c + gx2c) + u(v,fk,fj+1,fi  ));
            dst(v,fkg,fj+1,fi+1) = ot*(2.0*(ccval + gx1c + gx2c) + u(v,fk,fj+1,fi+1));
          }
        }
      }
    }
  }

  return;
}
