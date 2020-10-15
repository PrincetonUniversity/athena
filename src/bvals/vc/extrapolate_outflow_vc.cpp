//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file extrapolate_outflow_vc.cpp
//  \brief implementation of extrapolated outflow BCs in each dimension for
//         vertex-centered AthenaArray
//
//  Notes:
//  Internal data is extrapolated out to ghost zones

// C, C++ headers
#include <iostream>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "bvals_vc.hpp"

#include "../../utils/interp_univariate.hpp"

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ExtrapolateOutflowInnerX1(
//          Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, inner x1 boundary

void VertexCenteredBoundaryVariable::ExtrapolateOutflowInnerX1(
    Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
// Leads to errors here
// #pragma omp simd
        for (int i = il-1; i >= il-ngh; --i) {
          (*var_vc)(n,k,j,i) = 0.;

          for (int ix_e=0; ix_e < NEXTRAPOLATE; ++ix_e) {
            Real const extr_coeff = \
              InterpolateLagrangeUniformBiasR<NEXTRAPOLATE>::coeff[ix_e];
            (*var_vc)(n,k,j,i) += extr_coeff * (*var_vc)(n,k,j,i+ix_e+1);
          }

        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ExtrapolateOutflowOuterX1(
//          Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, outer x1 boundary

void VertexCenteredBoundaryVariable::ExtrapolateOutflowOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
// #pragma omp simd
        for (int i=iu+1; i<=iu+ngh; ++i) {
          (*var_vc)(n,k,j,i) = 0.;

          for (int ix_e=0; ix_e < NEXTRAPOLATE; ++ix_e) {
            // Real const extr_coeff = \
            //   InterpolateLagrangeUniformBiasR<NEXTRAPOLATE>::coeff[ix_e];
            Real const extr_coeff = \
              InterpolateLagrangeUniformBiasL<NEXTRAPOLATE>::coeff[NEXTRAPOLATE-ix_e-1];
            (*var_vc)(n,k,j,i) += extr_coeff * (*var_vc)(n,k,j,i-ix_e-1);
          }

        }

      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ExtrapolateOutflowInnerX2(
//          Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, inner x2 boundary

void VertexCenteredBoundaryVariable::ExtrapolateOutflowInnerX2(
    Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) {

  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl-1; j>=jl-ngh; --j) {
// #pragma omp simd
        for (int i=il; i<=iu; ++i) {
          (*var_vc)(n,k,j,i) = 0.;

          for (int ix_e=0; ix_e < NEXTRAPOLATE; ++ix_e) {
            Real const extr_coeff = \
              InterpolateLagrangeUniformBiasR<NEXTRAPOLATE>::coeff[ix_e];
            (*var_vc)(n,k,j,i) += extr_coeff * (*var_vc)(n,k,j+ix_e+1,i);
          }

        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ExtrapolateOutflowOuterX2(
//          Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, outer x2 boundary

void VertexCenteredBoundaryVariable::ExtrapolateOutflowOuterX2(
    Real time, Real dt, int il, int iu, int ju, int kl, int ku, int ngh) {

  for (int n=0; n<=nu_; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=ju+1; j<=ju+ngh; ++j) {
// #pragma omp simd
        for (int i=il; i<=iu; ++i) {
          (*var_vc)(n,k,j,i) = 0.;

          for (int ix_e=0; ix_e < NEXTRAPOLATE; ++ix_e) {
            Real const extr_coeff = \
              InterpolateLagrangeUniformBiasL<NEXTRAPOLATE>::coeff[NEXTRAPOLATE-ix_e-1];
            (*var_vc)(n,k,j,i) += extr_coeff * (*var_vc)(n,k,j-ix_e-1,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ExtrapolateOutflowInnerX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, inner x3 boundary

void VertexCenteredBoundaryVariable::ExtrapolateOutflowInnerX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ngh) {

  for (int n=0; n<=nu_; ++n) {
    for (int k=kl-1; k>=kl-ngh; --k) {
      for (int j=jl; j<=ju; ++j) {
// #pragma omp simd
        for (int i=il; i<=iu; ++i) {
          (*var_vc)(n,k,j,i) = 0.;

          for (int ix_e=0; ix_e < NEXTRAPOLATE; ++ix_e) {
            Real const extr_coeff = \
              InterpolateLagrangeUniformBiasR<NEXTRAPOLATE>::coeff[ix_e];
            (*var_vc)(n,k,j,i) += extr_coeff * (*var_vc)(n,k+ix_e+1,j,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ExtrapolateOutflowOuterX3(
//          Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh)
//  \brief OUTFLOW boundary conditions with extrapolation, outer x3 boundary

void VertexCenteredBoundaryVariable::ExtrapolateOutflowOuterX3(
    Real time, Real dt, int il, int iu, int jl, int ju, int ku, int ngh) {

  for (int n=0; n<=nu_; ++n) {
    for (int k=ku+1; k<=ku+ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
// #pragma omp simd
        for (int i=il; i<=iu; ++i) {
          (*var_vc)(n,k,j,i) = 0.;

          for (int ix_e=0; ix_e < NEXTRAPOLATE; ++ix_e) {
            Real const extr_coeff = \
              InterpolateLagrangeUniformBiasL<NEXTRAPOLATE>::coeff[NEXTRAPOLATE-ix_e-1];
            (*var_vc)(n,k,j,i) += extr_coeff * (*var_vc)(n,k-ix_e-1,j,i);
          }
        }
      }
    }
  }

  return;
}
