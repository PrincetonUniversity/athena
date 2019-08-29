//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reflect_rad.cpp
//  \brief implementation of reflecting boundary conditions for radiation

// Athena++ headers
#include "bvals_rad.hpp"
#include "../../../athena.hpp"         // Real, indices
#include "../../../athena_arrays.hpp"  // AthenaArray

//----------------------------------------------------------------------------------------
// Radiation reflecting boundary at inner x^1-surface
// Inputs:
//   time: coordinate time of simulation
//   dt: coordinate timestep
//   il: x^1-index of surface
//   jl,ju: x^2-index bounds of surface
//   kl,ku: x^3-index bounds of surface
//   ngh: number of ghost cells to fill
// Outputs: (none)
// Notes:
//   Fills active and ghost angles in ghost cells.
//   Requires active and ghost angles set in active cells.
//   Reflection is achieved via simple bilinear interpolation from the angular grid to
//       itself, which may not exactly conserve any particular quantity.

void RadBoundaryVariable::ReflectInnerX1(Real time, Real dt, int il, int jl, int ju,
    int kl, int ku, int ngh) {
  for (int n = 0; n <= nu_; ++n) {
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int di = 0; di < ngh; ++di) {
          int i_g = il - ngh + di;
          int i_a = il + ngh - 1 - di;
          int ind_a = reflect_ind_ix1_(0,n,k,j,di);
          int ind_b = reflect_ind_ix1_(1,n,k,j,di);
          int ind_c = reflect_ind_ix1_(2,n,k,j,di);
          int ind_d = reflect_ind_ix1_(3,n,k,j,di);
          Real frac_a = reflect_frac_ix1_(0,n,k,j,di);
          Real frac_b = reflect_frac_ix1_(1,n,k,j,di);
          Real frac_c = reflect_frac_ix1_(2,n,k,j,di);
          Real frac_d = reflect_frac_ix1_(3,n,k,j,di);
          Real val_a = (*var_cc)(ind_a,k,j,i_a);
          Real val_b = (*var_cc)(ind_b,k,j,i_a);
          Real val_c = (*var_cc)(ind_c,k,j,i_a);
          Real val_d = (*var_cc)(ind_d,k,j,i_a);
          (*var_cc)(n,k,j,i_g) =
              frac_a * val_a + frac_b * val_b + frac_c * val_c + frac_d * val_d;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Radiation reflecting boundary at outer x^1-surface
// Inputs:
//   time: coordinate time of simulation
//   dt: coordinate timestep
//   iu: x^1-index of surface
//   jl,ju: x^2-index bounds of surface
//   kl,ku: x^3-index bounds of surface
//   ngh: number of ghost cells to fill
// Outputs: (none)
// Notes:
//   Fills active and ghost angles in ghost cells.
//   Requires active and ghost angles set in active cells.
//   Reflection is achieved via simple bilinear interpolation from the angular grid to
//       itself, which may not exactly conserve any particular quantity.

void RadBoundaryVariable::ReflectOuterX1(Real time, Real dt, int iu, int jl, int ju,
    int kl, int ku, int ngh) {
  for (int n = 0; n <= nu_; ++n) {
    for (int k = kl; k <= ku; ++k) {
      for (int j = jl; j <= ju; ++j) {
        for (int di = 0; di < ngh; ++di) {
          int i_g = iu + 1 + di;
          int i_a = iu - di;
          int ind_a = reflect_ind_ox1_(0,n,k,j,di);
          int ind_b = reflect_ind_ox1_(1,n,k,j,di);
          int ind_c = reflect_ind_ox1_(2,n,k,j,di);
          int ind_d = reflect_ind_ox1_(3,n,k,j,di);
          Real frac_a = reflect_frac_ox1_(0,n,k,j,di);
          Real frac_b = reflect_frac_ox1_(1,n,k,j,di);
          Real frac_c = reflect_frac_ox1_(2,n,k,j,di);
          Real frac_d = reflect_frac_ox1_(3,n,k,j,di);
          Real val_a = (*var_cc)(ind_a,k,j,i_a);
          Real val_b = (*var_cc)(ind_b,k,j,i_a);
          Real val_c = (*var_cc)(ind_c,k,j,i_a);
          Real val_d = (*var_cc)(ind_d,k,j,i_a);
          (*var_cc)(n,k,j,i_g) =
              frac_a * val_a + frac_b * val_b + frac_c * val_c + frac_d * val_d;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Radiation reflecting boundary at inner x^2-surface
// Inputs:
//   time: coordinate time of simulation
//   dt: coordinate timestep
//   il,iu: x^1-index bounds of surface
//   jl: x^2-index of surface
//   kl,ku: x^3-index bounds of surface
//   ngh: number of ghost cells to fill
// Outputs: (none)
// Notes:
//   Fills active and ghost angles in ghost cells.
//   Requires active and ghost angles set in active cells.
//   Reflection is achieved via simple bilinear interpolation from the angular grid to
//       itself, which may not exactly conserve any particular quantity.

void RadBoundaryVariable::ReflectInnerX2(Real time, Real dt, int il, int iu, int jl,
    int kl, int ku, int ngh) {
  for (int n = 0; n <= nu_; ++n) {
    for (int k = kl; k <= ku; ++k) {
      for (int dj = 0; dj < ngh; ++dj) {
        for (int i = il; i <= iu; ++i) {
          int j_g = jl - ngh + dj;
          int j_a = jl + ngh - 1 - dj;
          int ind_a = reflect_ind_ix2_(0,n,k,dj,i);
          int ind_b = reflect_ind_ix2_(1,n,k,dj,i);
          int ind_c = reflect_ind_ix2_(2,n,k,dj,i);
          int ind_d = reflect_ind_ix2_(3,n,k,dj,i);
          Real frac_a = reflect_frac_ix2_(0,n,k,dj,i);
          Real frac_b = reflect_frac_ix2_(1,n,k,dj,i);
          Real frac_c = reflect_frac_ix2_(2,n,k,dj,i);
          Real frac_d = reflect_frac_ix2_(3,n,k,dj,i);
          Real val_a = (*var_cc)(ind_a,k,j_a,i);
          Real val_b = (*var_cc)(ind_b,k,j_a,i);
          Real val_c = (*var_cc)(ind_c,k,j_a,i);
          Real val_d = (*var_cc)(ind_d,k,j_a,i);
          (*var_cc)(n,k,j_g,i) =
              frac_a * val_a + frac_b * val_b + frac_c * val_c + frac_d * val_d;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Radiation reflecting boundary at outer x^2-surface
// Inputs:
//   time: coordinate time of simulation
//   dt: coordinate timestep
//   il,iu: x^1-index bounds of surface
//   ju: x^2-index of surface
//   kl,ku: x^3-index bounds of surface
//   ngh: number of ghost cells to fill
// Outputs: (none)
// Notes:
//   Fills active and ghost angles in ghost cells.
//   Requires active and ghost angles set in active cells.
//   Reflection is achieved via simple bilinear interpolation from the angular grid to
//       itself, which may not exactly conserve any particular quantity.

void RadBoundaryVariable::ReflectOuterX2(Real time, Real dt, int il, int iu, int ju,
    int kl, int ku, int ngh) {
  for (int n = 0; n <= nu_; ++n) {
    for (int k = kl; k <= ku; ++k) {
      for (int dj = 0; dj < ngh; ++dj) {
        for (int i = il; i <= iu; ++i) {
          int j_g = ju + 1 + dj;
          int j_a = ju - dj;
          int ind_a = reflect_ind_ox2_(0,n,k,dj,i);
          int ind_b = reflect_ind_ox2_(1,n,k,dj,i);
          int ind_c = reflect_ind_ox2_(2,n,k,dj,i);
          int ind_d = reflect_ind_ox2_(3,n,k,dj,i);
          Real frac_a = reflect_frac_ox2_(0,n,k,dj,i);
          Real frac_b = reflect_frac_ox2_(1,n,k,dj,i);
          Real frac_c = reflect_frac_ox2_(2,n,k,dj,i);
          Real frac_d = reflect_frac_ox2_(3,n,k,dj,i);
          Real val_a = (*var_cc)(ind_a,k,j_a,i);
          Real val_b = (*var_cc)(ind_b,k,j_a,i);
          Real val_c = (*var_cc)(ind_c,k,j_a,i);
          Real val_d = (*var_cc)(ind_d,k,j_a,i);
          (*var_cc)(n,k,j_g,i) =
              frac_a * val_a + frac_b * val_b + frac_c * val_c + frac_d * val_d;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Radiation reflecting boundary at inner x^3-surface
// Inputs:
//   time: coordinate time of simulation
//   dt: coordinate timestep
//   il,iu: x^1-index bounds of surface
//   jl,ju: x^2-index bounds of surface
//   kl: x^3-index of surface
//   ngh: number of ghost cells to fill
// Outputs: (none)
// Notes:
//   Fills active and ghost angles in ghost cells.
//   Requires active and ghost angles set in active cells.
//   Reflection is achieved via simple bilinear interpolation from the angular grid to
//       itself, which may not exactly conserve any particular quantity.

void RadBoundaryVariable::ReflectInnerX3(Real time, Real dt, int il, int iu, int jl,
    int ju, int kl, int ngh) {
  for (int n = 0; n <= nu_; ++n) {
    for (int dk = 0; dk < ngh; ++dk) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = il; i <= iu; ++i) {
          int k_g = kl - ngh + dk;
          int k_a = kl + ngh - 1 - dk;
          int ind_a = reflect_ind_ix3_(0,n,dk,j,i);
          int ind_b = reflect_ind_ix3_(1,n,dk,j,i);
          int ind_c = reflect_ind_ix3_(2,n,dk,j,i);
          int ind_d = reflect_ind_ix3_(3,n,dk,j,i);
          Real frac_a = reflect_frac_ix3_(0,n,dk,j,i);
          Real frac_b = reflect_frac_ix3_(1,n,dk,j,i);
          Real frac_c = reflect_frac_ix3_(2,n,dk,j,i);
          Real frac_d = reflect_frac_ix3_(3,n,dk,j,i);
          Real val_a = (*var_cc)(ind_a,k_a,j,i);
          Real val_b = (*var_cc)(ind_b,k_a,j,i);
          Real val_c = (*var_cc)(ind_c,k_a,j,i);
          Real val_d = (*var_cc)(ind_d,k_a,j,i);
          (*var_cc)(n,k_g,j,i) =
              frac_a * val_a + frac_b * val_b + frac_c * val_c + frac_d * val_d;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// Radiation reflecting boundary at outer x^3-surface
// Inputs:
//   time: coordinate time of simulation
//   dt: coordinate timestep
//   il,iu: x^1-index bounds of surface
//   jl,ju: x^2-index bounds of surface
//   ku: x^3-index of surface
//   ngh: number of ghost cells to fill
// Outputs: (none)
// Notes:
//   Fills active and ghost angles in ghost cells.
//   Requires active and ghost angles set in active cells.
//   Reflection is achieved via simple bilinear interpolation from the angular grid to
//       itself, which may not exactly conserve any particular quantity.

void RadBoundaryVariable::ReflectOuterX3(Real time, Real dt, int il, int iu, int jl,
    int ju, int ku, int ngh) {
  for (int n = 0; n <= nu_; ++n) {
    for (int dk = 0; dk < ngh; ++dk) {
      for (int j = jl; j <= ju; ++j) {
        for (int i = il; i <= iu; ++i) {
          int k_g = ku + 1 + dk;
          int k_a = ku - dk;
          int ind_a = reflect_ind_ox3_(0,n,dk,j,i);
          int ind_b = reflect_ind_ox3_(1,n,dk,j,i);
          int ind_c = reflect_ind_ox3_(2,n,dk,j,i);
          int ind_d = reflect_ind_ox3_(3,n,dk,j,i);
          Real frac_a = reflect_frac_ox3_(0,n,dk,j,i);
          Real frac_b = reflect_frac_ox3_(1,n,dk,j,i);
          Real frac_c = reflect_frac_ox3_(2,n,dk,j,i);
          Real frac_d = reflect_frac_ox3_(3,n,dk,j,i);
          Real val_a = (*var_cc)(ind_a,k_a,j,i);
          Real val_b = (*var_cc)(ind_b,k_a,j,i);
          Real val_c = (*var_cc)(ind_c,k_a,j,i);
          Real val_d = (*var_cc)(ind_d,k_a,j,i);
          (*var_cc)(n,k_g,j,i) =
              frac_a * val_a + frac_b * val_b + frac_c * val_c + frac_d * val_d;
        }
      }
    }
  }
  return;
}
