#ifndef Z4c_MACRO_HPP
#define Z4c_MACRO_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file z4c_macro.hpp
//  \brief common macros for the Z4c class

#define IX_IL                                                                 \
  pmy_block->pz4c->mbi.il

#define IX_IU                                                                 \
  pmy_block->pz4c->mbi.iu

#define IX_JL                                                                 \
  pmy_block->pz4c->mbi.jl

#define IX_JU                                                                 \
  pmy_block->pz4c->mbi.ju

#define IX_KL                                                                 \
  pmy_block->pz4c->mbi.kl

#define IX_KU                                                                 \
  pmy_block->pz4c->mbi.ku

#define GSIZEI                                                                \
  (NGHOST)

#define GSIZEJ                                                                \
  ((pmy_block->block_size.nx2 > 1) ? (NGHOST) : (0))

#define GSIZEK                                                                \
  ((pmy_block->block_size.nx3 > 1) ? (NGHOST) : (0))

// 2D loop over k and j in the interior of the block
#define ILOOP2(k,j)                                                           \
  for(int k = IX_KL; k <= IX_KU; ++k)                         \
  for(int j = IX_JL; j <= IX_JU; ++j)

// 2D loop over k and j on the whole block
#define GLOOP2(k,j)                                                           \
  for(int k = IX_KL - GSIZEK; k <= IX_KU + GSIZEK; ++k)       \
  for(int j = IX_JL - GSIZEJ; j <= IX_JU + GSIZEJ; ++j)

// 1D loop over i in the interior of the block
#define ILOOP1(i)                                                             \
  _Pragma("omp simd")                                                         \
  for(int i = IX_IL; i <= IX_IU; ++i)

// 1D loop over i on the whole block
#define GLOOP1(i)                                                             \
  _Pragma("omp simd")                                                         \
  for(int i = IX_IL - GSIZEI; i <= IX_IU + GSIZEI; ++i)

// 3D loop over the interior of the block
#define ILOOP3(k,j,i)                                                         \
    ILOOP2(k,j)                                                               \
    ILOOP1(i)

// 3D loop over the whole block
#define GLOOP3(k,j,i)                                                         \
    GLOOP2(k,j)                                                               \
    GLOOP1(i)

#endif
