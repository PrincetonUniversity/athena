//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file consistency_conditions_vc.cpp
//  \brief additive unpack of vertex-centered data needs averaging

// C headers

// C++ headers
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "bvals_vc.hpp"

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::AllocateNodeMult()
//  \brief Allocate node_mult array
void VertexCenteredBoundaryVariable::AllocateNodeMult() {
  if (!node_mult.IsAllocated())
    if (pmy_block_->pmy_mesh->ndim == 3)
      node_mult.NewAthenaArray(1, 7, 7, 7);
    else if (pmy_block_->pmy_mesh->ndim == 2) {
      c_kvs = c_kve = c_kpe = 0;
      node_mult.NewAthenaArray(1, 1, 7, 7);
    } else {
      c_kvs = c_kve = c_kpe = 0;
      c_jvs = c_jve = c_jpe = 0;
      node_mult.NewAthenaArray(1, 1, 1, 7);
    }
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::PrepareNodeMult()
//  \brief Prepare node_mult array allocating as required
void VertexCenteredBoundaryVariable::PrepareNodeMult() {
  AllocateNodeMult();
  node_mult.ZeroClear();

  for (int k = c_kvs; k <= c_kve; ++k)
    for (int j = c_jvs; j <= c_jve; ++j)
      for (int i = c_ivs; i <= c_ive; ++i)
        node_mult(0, k, j, i) = 1;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ZeroVertexGhosts()
//  \brief zero out ghost zones

void VertexCenteredBoundaryVariable::ZeroVertexGhosts() {
  if (DBGPR_CONSISTENCY_CONDITIONS_VC)
    coutYellow("VertexCenteredBoundaryVariable::ZeroVertexGhosts\n");

  // additive unpack used to populate ghosts entails that old ghost data needs
  // to be cleaned

  // BD: TODO coarse_buf needs to be zerod also...
  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_vc;
  if (pmb->block_size.nx3 > 1) {

    // fill bands
    for (int n=nl_; n<=nu_; ++n) {
      for (int k=0; k<pmb->ng; k++)
        for (int j=0; j<=pmb->jpe; j++)
#pragma omp simd
          for (int i=0; i<=pmb->ipe; i++) {
            var(n, k, j, i) = 0;
            var(n, pmb->kps + k, j, i) = 0;
          }

      for (int k=pmb->kvs; k<=pmb->kve; k++)
        for (int j=0; j<=pmb->jpe; j++)
#pragma omp simd
          for (int i=0; i<pmb->ng; i++) {
            var(n, k, j, i) = 0;
            var(n, k, j, pmb->ips + i) = 0;
          }

      for (int k=pmb->kvs; k<=pmb->kve; k++)
        for (int j=0; j<pmb->ng; j++)
#pragma omp simd
          for (int i=pmb->ivs; i<=pmb->ive; i++) {
            var(n, k, j, i) = 0;
            var(n, k, pmb->jps + j, i) = 0;
          }

    }

  } else if (pmb->block_size.nx2 > 1) {
    for (int n=nl_; n<=nu_; ++n) {
      // top and bottom bands
      for (int j=0; j<pmb->ng; j++) {
#pragma omp simd
        for (int i=0; i<=pmb->ipe; i++) {
          var(n, 0, j, i) = 0;
          var(n, 0, pmb->jps + j, i) = 0;
        }
      }
      // east and west truncated bands
      for (int j=pmb->jvs; j<=pmb->jve; j++) {
#pragma omp simd
        for (int i=0; i<pmb->ng; i++) {
          var(n, 0, j, i) = 0;
          var(n, 0, j, pmb->jps + i) = 0;
        }
      }
    }
  } else {
    for (int n=nl_; n<=nu_; ++n) {
#pragma omp simd
      for (int i=pmb->ims; i<=pmb->ime; ++i) {
        var(n, 0, 0, i) = 0;
      }
#pragma omp simd
      for (int i=pmb->ips; i<=pmb->ipe; ++i) {
        var(n, 0, 0, i) = 0;
      }
    }
  }

  return;
  }

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ApplyNodeMultiplicitesDim3(...)
//  \brief Division factors applied in three dimensions
inline void VertexCenteredBoundaryVariable::ApplyNodeMultiplicitesDim3(
  AthenaArray<Real> &var,
  int ims, int ivs, int ive, int ipe, int axis_half_size_x1,
  int jms, int jvs, int jve, int jpe, int axis_half_size_x2,
  int kms, int kvs, int kve, int kpe, int axis_half_size_x3) {

  // provide indices
  int const tk_c[7] = {kms, kvs, kvs + 1,
                       kvs  + axis_half_size_x3,
                       kve - 1, kve, kpe};
  int const tj_c[7] = {jms, jvs, jvs + 1,
                       jvs  + axis_half_size_x2,
                       jve - 1, jve, jpe};
  int const ti_c[7] = {ims, ivs, ivs + 1,
                       ivs + axis_half_size_x1,
                       ive - 1, ive, ipe};

  short int ufac;

  for (int K=0; K<7; ++K)
    for (int J=0; J<7; ++J)
      for (int I=0; I<7; ++I) {
        ufac = node_mult(0, K, J, I);

        if (ufac > 1) {

          int il, iu;
          if (I % 2 == 0) {
            if (I > 3)
              il = ti_c[I-1]+1, iu = ti_c[I];
            else
              il = ti_c[I], iu = ti_c[I+1] - 1;
          } else {
            il = iu = ti_c[I];
          }

          int jl, ju;
          if (J % 2 == 0) {
            if (J > 3)
              jl = tj_c[J-1]+1, ju = tj_c[J];
            else
              jl = tj_c[J], ju = tj_c[J+1] - 1;
          } else {
            jl = ju = tj_c[J];
          }

          int kl, ku;
          if (K % 2 == 0) {
            if (K > 3)
              kl = tk_c[K-1]+1, ku = tk_c[K];
            else
              kl = tk_c[K], ku = tk_c[K+1] - 1;
          } else {
            kl = ku = tk_c[K];
          }

          for (int n_=nl_; n_<=nu_; ++n_)
            for (int k=kl; k<=ku; ++k)
              for (int j=jl; j<=ju; ++j)
#pragma omp simd
                for (int i=il; i<=iu; ++i)
                  var(n_, k, j, i) /= ufac;


        }
      }


}


//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ApplyNodeMultiplicitesDim2(...)
//  \brief Division factors applied in two dimensions
inline void VertexCenteredBoundaryVariable::ApplyNodeMultiplicitesDim2(
  AthenaArray<Real> &var,
  int ims, int ivs, int ive, int ipe, int axis_half_size_x1,
  int jms, int jvs, int jve, int jpe, int axis_half_size_x2) {

  // provide indices
  int const tj_c[7] = {jms, jvs, jvs + 1,
                       jvs  + axis_half_size_x2,
                       jve - 1, jve, jpe};
  int const ti_c[7] = {ims, ivs, ivs + 1,
                       ivs + axis_half_size_x1,
                       ive - 1, ive, ipe};

  short int ufac;

  for (int J=0; J<7; ++J)
    for (int I=0; I<7; ++I) {
      ufac = node_mult(0, 0, J, I);

      if (ufac > 1) {

        int il, iu;
        if (I % 2 == 0) {
          if (I > 3)
            il = ti_c[I-1]+1, iu = ti_c[I];
          else
            il = ti_c[I], iu = ti_c[I+1] - 1;
        } else {
          il = iu = ti_c[I];
        }
        int jl, ju;
        if (J % 2 == 0) {
          if (J > 3)
            jl = tj_c[J-1]+1, ju = tj_c[J];
          else
            jl = tj_c[J], ju = tj_c[J+1] - 1;
        } else {
          jl = ju = tj_c[J];
        }

        for (int n_=nl_; n_<=nu_; ++n_)
          for (int j=jl; j<=ju; ++j)
#pragma omp simd
            for (int i=il; i<=iu; ++i)
              var(n_, 0, j, i) /= ufac;

      }
    }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::ApplyNodeMultiplicitesDim1(...)
//  \brief Division factors applied in one dimension
inline void VertexCenteredBoundaryVariable::ApplyNodeMultiplicitesDim1(
  AthenaArray<Real> &var,
  int ims, int ivs, int ive, int ipe, int axis_half_size_x1) {

  // provide indices
  int const ti_c[7] = {ims, ivs, ivs + 1,
                       ivs + axis_half_size_x1,
                       ive - 1, ive, ipe};

  short int ufac;

  for (int I=0; I<7; ++I) {
    ufac = node_mult(0, 0, 0, I);

    if (ufac > 1) {

      int il, iu;
      if (I % 2 == 0) {
        if (I > 3)
          il = ti_c[I-1]+1, iu = ti_c[I];
        else
          il = ti_c[I], iu = ti_c[I+1] - 1;
      } else {
        il = iu = ti_c[I];
      }

    for (int n_=nl_; n_<=nu_; ++n_)
#pragma omp simd
        for (int i=il; i<=iu; ++i)
          var(n_, 0, 0, i) /= ufac;

    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void VertexCenteredBoundaryVariable::FinalizeVertexConsistency()
//  \brief Apply division factors to shared vertices
void VertexCenteredBoundaryVariable::FinalizeVertexConsistency() {
  if (DBGPR_BVALS_VC)
    coutYellow("VertexCenteredBoundaryVariable::FinalizeVertexConsistency\n");

  MeshBlock *pmb = pmy_block_;
  AthenaArray<Real> &var = *var_vc;

  // different conditions in different dimensions
  if (pmb->block_size.nx3 > 1) {
    ApplyNodeMultiplicitesDim3(
      var, pmb->ims, pmb->ivs, pmb->ive, pmb->ipe,
      pmb->block_size.nx1 / 2,
      pmb->jms, pmb->jvs, pmb->jve, pmb->jpe,
      pmb->block_size.nx2 / 2,
      pmb->kms, pmb->kvs, pmb->kve, pmb->kpe,
      pmb->block_size.nx3 / 2);

    if (pmy_mesh_->multilevel) {
      AthenaArray<Real> &coarse_var = *coarse_buf;

      ApplyNodeMultiplicitesDim3(
        coarse_var, pmb->cims, pmb->civs, pmb->cive, pmb->cipe,
        pmb->block_size.nx1 / 4,
        pmb->cjms, pmb->cjvs, pmb->cjve, pmb->cjpe,
        pmb->block_size.nx2 / 4,
        pmb->ckms, pmb->ckvs, pmb->ckve, pmb->ckpe,
        pmb->block_size.nx3 / 4);
    }
  } else if (pmb->block_size.nx2 > 1) {
    ApplyNodeMultiplicitesDim2(
      var, pmb->ims, pmb->ivs, pmb->ive, pmb->ipe,
      pmb->block_size.nx1 / 2,
      pmb->jms, pmb->jvs, pmb->jve, pmb->jpe,
      pmb->block_size.nx2 / 2);

    if (pmy_mesh_->multilevel) {
      AthenaArray<Real> &coarse_var = *coarse_buf;

      ApplyNodeMultiplicitesDim2(
        coarse_var, pmb->cims, pmb->civs, pmb->cive, pmb->cipe,
        pmb->block_size.nx1 / 4,
        pmb->cjms, pmb->cjvs, pmb->cjve, pmb->cjpe,
        pmb->block_size.nx2 / 4);
    }
  } else {
    ApplyNodeMultiplicitesDim1(
      var, pmb->ims, pmb->ivs, pmb->ive, pmb->ipe,
      pmb->block_size.nx1 / 2);

    if (pmy_mesh_->multilevel) {
      AthenaArray<Real> &coarse_var = *coarse_buf;

      ApplyNodeMultiplicitesDim1(
        coarse_var, pmb->cims, pmb->civs, pmb->cive, pmb->cipe,
        pmb->block_size.nx1 / 4);
    }

  }

  return;
}
