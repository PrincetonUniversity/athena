#ifndef UTILS_BUFFER_UTILS_HPP_
#define UTILS_BUFFER_UTILS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file buffer_utils.hpp
//  \brief prototypes of utility functions to pack/unpack buffers

#include "../athena.hpp"
#include "../athena_arrays.hpp"

namespace BufferUtility {
void Pack4DData(AthenaArray<Real> &src, Real *buf, int sn, int en,
               int si, int ei, int sj, int ej, int sk, int ek, int &offset);
void Unpack4DData(Real *buf, AthenaArray<Real> &dst, int sn, int en,
                  int si, int ei, int sj, int ej, int sk, int ek, int &offset);
void Pack3DData(AthenaArray<Real> &src, Real *buf,
               int si, int ei, int sj, int ej, int sk, int ek, int &offset);
void Unpack3DData(Real *buf, AthenaArray<Real> &dst,
                  int si, int ei, int sj, int ej, int sk, int ek, int &offset);
}
#endif // UTILS_BUFFER_UTILS_HPP_
