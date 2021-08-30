//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file calculate_velocity_differences.cpp
//! \brief Calculate velocity differences for Carbuncle cure, called from RiemannSolver

// C headers

// C++ headers
#include <algorithm>   // min,max

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "hydro.hpp"

void Hydro::CalculateVelocityDifferences(const int k, const int j,
            const int il, const int iu, const int ivx,
            AthenaArray<Real> &dvn, AthenaArray<Real> &dvt) {
  if (pmy_block->pmy_mesh->f3) {
    if (ivx == IVX) {
      for (int i = il; i <= iu; ++i) {
        dvn(i) = w(IVX, k, j, i) - w(IVX, k, j, i-1);
        Real dvl = std::min(w(IVY, k, j+1, i-1) - w(IVY, k, j,   i-1),
                            w(IVY, k, j,   i-1) - w(IVY, k, j-1, i-1));
        Real dvr = std::min(w(IVY, k, j+1, i)   - w(IVY, k, j,   i),
                            w(IVY, k, j,   i)   - w(IVY, k, j-1, i));
        Real dwl = std::min(w(IVZ, k+1, j, i-1) - w(IVZ, k,   j, i-1),
                            w(IVZ, k,   j, i-1) - w(IVZ, k-1, j, i-1));
        Real dwr = std::min(w(IVZ, k+1, j, i)   - w(IVZ, k,   j, i),
                            w(IVZ, k,   j, i)   - w(IVZ, k-1, j, i));
        dvt(i) = std::min(std::min(dvl, dvr), std::min(dwl, dwr));
      }
    } else if (ivx == IVY) {
      // note: technically, we can reuse dvr/dwr as dvl/dwl in the next j-loop.
      for (int i = il; i <= iu; ++i) {
        dvn(i) = w(IVY, k, j, i) - w(IVY, k, j-1, i);
        Real dvl = std::min(w(IVZ, k+1, j-1, i) - w(IVZ, k,   j-1, i),
                            w(IVZ, k,   j-1, i) - w(IVZ, k-1, j-1, i));
        Real dvr = std::min(w(IVZ, k+1, j,   i) - w(IVZ, k,   j,   i),
                            w(IVZ, k,   j,   i) - w(IVZ, k-1, j,   i));
        Real dwl = std::min(w(IVX, k, j-1, i+1) - w(IVX, k, j-1, i),
                            w(IVX, k, j-1, i)   - w(IVX, k, j-1, i-1));
        Real dwr = std::min(w(IVX, k, j,   i+1) - w(IVX, k, j,   i),
                            w(IVX, k, j,   i)   - w(IVX, k, j,   i-1));
        dvt(i) = std::min(std::min(dvl, dvr), std::min(dwl, dwr));
      }
    } else { // (ivx == IVZ)
      for (int i = il; i <= iu; ++i) {
        dvn(i) = w(IVZ, k, j, i) - w(IVZ, k-1, j, i);
        Real dvl = std::min(w(IVX, k-1, j, i+1) - w(IVX, k-1, j, i),
                            w(IVX, k-1, j, i)   - w(IVX, k-1, j, i-1));
        Real dvr = std::min(w(IVX, k,   j, i+1) - w(IVX, k,   j, i),
                            w(IVX, k,   j, i)   - w(IVX, k,   j, i-1));
        Real dwl = std::min(w(IVY, k-1, j+1, i) - w(IVY, k-1, j,   i),
                            w(IVY, k-1, j,   i) - w(IVY, k-1, j-1, i));
        Real dwr = std::min(w(IVY, k,   j+1, i) - w(IVY, k,   j,   i),
                            w(IVY, k,   j,   i) - w(IVY, k,   j-1, i));
        dvt(i) = std::min(std::min(dvl, dvr), std::min(dwl, dwr));
      }
    }
  } else if (pmy_block->pmy_mesh->f2) {
    if (ivx == IVX) {
      for (int i = il; i <= iu; ++i) {
        dvn(i) = w(IVX, k, j, i) - w(IVX, k, j, i-1);
        Real dvl = std::min(w(IVY, k, j+1, i-1) - w(IVY, k, j,   i-1),
                            w(IVY, k, j,   i-1) - w(IVY, k, j-1, i-1));
        Real dvr = std::min(w(IVY, k, j+1, i)   - w(IVY, k, j,   i),
                            w(IVY, k, j,   i)   - w(IVY, k, j-1, i));
        dvt(i) = std::min(dvl, dvr);
      }
    } else { //ivx == IVY
      for (int i = il; i <= iu; ++i) {
        dvn(i) = w(IVY, k, j, i) - w(IVY, k, j-1, i);\
        Real dvl = std::min(w(IVX, k, j-1, i+1) - w(IVX, k, j-1, i),
                            w(IVX, k, j-1, i)   - w(IVX, k, j-1, i-1));
        Real dvr = std::min(w(IVX, k, j,   i+1) - w(IVX, k, j,   i),
                            w(IVX, k, j,   i)   - w(IVX, k, j,   i-1));
        dvt(i) = std::min(dvl, dvr);
      }
    }
  } else {
    for (int i = il; i <= iu; ++i)
      dvn(i) = w(IVX, k, j, i) - w(IVX, k, j, i-1);
  }
  return;
}
