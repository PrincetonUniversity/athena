//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file buffer_utils.cpp
//  \brief namespace containing buffer utilities.

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "buffer_utils.hpp"

namespace BufferUtility {
//----------------------------------------------------------------------------------------
//! \fn template <typename T> void PackData(AthenaArray<T> &src, T *buf, int sn, int en,
//                     int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//  \brief pack a 4D AthenaArray into a one-dimensional buffer

template <typename T> void PackData(AthenaArray<T> &src, T *buf,
                                    int sn, int en,
                                    int si, int ei, int sj, int ej, int sk, int ek,
                                    int &offset) {
  for (int n=sn; n<=en; ++n) {
    for (int k=sk; k<=ek; k++) {
      for (int j=sj; j<=ej; j++) {
#pragma omp simd
        for (int i=si; i<=ei; i++)
          buf[offset++]=src(n,k,j,i);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn template <typename T> void PackData(AthenaArray<T> &src, T *buf,
//                      int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//  \brief pack a 3D AthenaArray into a one-dimensional buffer

template <typename T> void PackData(AthenaArray<T> &src, T *buf,
                                    int si, int ei, int sj, int ej, int sk, int ek,
                                    int &offset) {
  for (int k=sk; k<=ek; k++) {
    for (int j=sj; j<=ej; j++) {
#pragma omp simd
      for (int i=si; i<=ei; i++)
        buf[offset++]=src(k, j, i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn template <typename T> void UnpackData(T *buf, AthenaArray<T> &dst, int sn, int en,
//                        int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//  \brief unpack a one-dimensional buffer into a 4D AthenaArray

template <typename T> void UnpackData(T *buf, AthenaArray<T> &dst,
                                      int sn, int en,
                                      int si, int ei, int sj, int ej, int sk, int ek,
                                      int &offset) {
  for (int n=sn; n<=en; ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma omp simd
        for (int i=si; i<=ei; ++i)
          dst(n,k,j,i) = buf[offset++];
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn template <typename T> void UnpackData(T *buf, AthenaArray<T> &dst,
//                        int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//  \brief unpack a one-dimensional buffer into a 3D AthenaArray

template <typename T> void UnpackData(T *buf, AthenaArray<T> &dst,
                                      int si, int ei, int sj, int ej, int sk, int ek,
                                      int &offset) {
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma omp simd
      for (int i=si; i<=ei; ++i)
        dst(k,j,i) = buf[offset++];
    }
  }
  return;
}

} // end namespace BufferUtility
