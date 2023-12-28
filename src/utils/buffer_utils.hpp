#ifndef UTILS_BUFFER_UTILS_HPP_
#define UTILS_BUFFER_UTILS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file buffer_utils.hpp
//! \brief prototypes of utility functions to pack/unpack buffers

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

namespace BufferUtility {
// 3x templated and overloaded functions
// 5D
template <typename T> void PackData(const AthenaArray<T> &src, T *buf,
         int sn, int en, int sm, int em,
         int si, int ei, int sj, int ej, int sk, int ek, int &offset);
// 4D
template <typename T> void PackData(const AthenaArray<T> &src, T *buf,
         int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset);
// 3D
template <typename T> void PackData(const AthenaArray<T> &src, T *buf,
                           int si, int ei, int sj, int ej, int sk, int ek, int &offset);
// 5D
template <typename T> void UnpackData(const T *buf, AthenaArray<T> &dst,
         int sn, int en, int sm, int em,
         int si, int ei, int sj, int ej, int sk, int ek, int &offset);
// 4D
template <typename T> void UnpackData(const T *buf, AthenaArray<T> &dst,
         int sn, int en, int si, int ei, int sj, int ej, int sk, int ek, int &offset);
// 3D
template <typename T> void UnpackData(const T *buf, AthenaArray<T> &dst,
                      int si, int ei, int sj, int ej, int sk, int ek, int &offset);
} // namespace BufferUtility
#endif // UTILS_BUFFER_UTILS_HPP_
