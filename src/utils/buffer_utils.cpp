//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file buffer_utils.cpp
//! \brief namespace containing buffer utilities.

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "buffer_utils.hpp"

namespace BufferUtility {
//----------------------------------------------------------------------------------------
//! \fn template <typename T> void PackData(const AthenaArray<T> &src, T *buf,
//!     int sn, int en, int sm, int em,
//!     int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//! \brief pack a 5D AthenaArray into a one-dimensional buffer

template <typename T> void PackData(const AthenaArray<T> &src, T *buf,
         int sn, int en, int sm, int em,
         int si, int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n=sn; n<=en; ++n) {
    for (int k=sk; k<=ek; k++) {
      for (int j=sj; j<=ej; j++) {
        for (int i=si; i<=ei; i++) {
#pragma omp simd
          for (int m=sm; m<=em; m++) {
            buf[offset++] = src(n,k,j,i,m);
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn template <typename T> void PackData(const AthenaArray<T> &src, T *buf,
//!     int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//! \brief pack a 4D AthenaArray into a one-dimensional buffer

template <typename T> void PackData(const AthenaArray<T> &src, T *buf,
         int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n=sn; n<=en; ++n) {
    for (int k=sk; k<=ek; k++) {
      for (int j=sj; j<=ej; j++) {
        for (int i=si; i<=ei; i++)
          buf[offset++] = src(n,k,j,i);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn template <typename T> void PackData(const AthenaArray<T> &src, T *buf,
//!                     int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//! \brief pack a 3D AthenaArray into a one-dimensional buffer

template <typename T> void PackData(const AthenaArray<T> &src, T *buf,
                                    int si, int ei, int sj, int ej, int sk, int ek,
                                    int &offset) {
  for (int k=sk; k<=ek; k++) {
    for (int j=sj; j<=ej; j++) {
      for (int i=si; i<=ei; i++)
        buf[offset++] = src(k, j, i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn template <typename T> void UnpackData(const T *buf, AthenaArray<T> &dst,
//!     int sn, int en, int sm, int em,
//!     int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//! \brief unpack a one-dimensional buffer into a 5D AthenaArray

template <typename T> void UnpackData(const T *buf, AthenaArray<T> &dst,
         int sn, int en, int sm, int em,
         int si, int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n=sn; n<=en; ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
        for (int i=si; i<=ei; ++i) {
#pragma omp simd
          for (int m=sm; m<=em; ++m) {
            dst(n,k,j,i,m) = buf[offset++];
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn template <typename T> void UnpackData(const T *buf, AthenaArray<T> &dst,
//!     int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//! \brief unpack a one-dimensional buffer into a 4D AthenaArray

template <typename T> void UnpackData(const T *buf, AthenaArray<T> &dst,
         int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int n=sn; n<=en; ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
        for (int i=si; i<=ei; ++i)
          dst(n,k,j,i) = buf[offset++];
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn template <typename T> void UnpackData(const T *buf, AthenaArray<T> &dst,
//!                       int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//! \brief unpack a one-dimensional buffer into a 3D AthenaArray

template <typename T> void UnpackData(const T *buf, AthenaArray<T> &dst,
                           int si, int ei, int sj, int ej, int sk, int ek, int &offset) {
  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
      for (int i=si; i<=ei; ++i)
        dst(k,j,i) = buf[offset++];
    }
  }
  return;
}

// provide explicit instantiation definitions (C++03) to allow the template definitions to
// exist outside of header file (non-inline), but still provide the requisite instances
// for other TUs during linking time (~13x files include "buffer_utils.hpp")

// 13x files include buffer_utils.hpp
template void UnpackData<Real>(const Real *, AthenaArray<Real> &,
                               int, int, int, int, int, int, int, int, int, int, int &);
template void UnpackData<Real>(const Real *, AthenaArray<Real> &,
                               int, int, int, int, int, int, int, int, int &);
template void UnpackData<Real>(const Real *, AthenaArray<Real> &,
                               int, int, int, int, int, int, int &);

template void PackData<Real>(const AthenaArray<Real> &, Real *,
                             int, int, int, int, int, int, int, int, int, int, int &);
template void PackData<Real>(const AthenaArray<Real> &, Real *,
                             int, int, int, int, int, int, int, int, int &);
template void PackData<Real>(const AthenaArray<Real> &, Real *,
                             int, int, int, int, int, int, int &);

} // end namespace BufferUtility
