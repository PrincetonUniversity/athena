//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file buffer_utils.cpp
//  \brief namespace containing buffer utilities.

// C headers

// C++ headers
#include <iostream>  // for debug

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "buffer_utils.hpp"

// BD: debug
bool DBGPR_BUFFER_UTILS = false;
//----------

namespace BufferUtility {
//----------------------------------------------------------------------------------------
//! \fn template <typename T> void PackData(AthenaArray<T> &src, T *buf, int sn, int en,
//                     int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//  \brief pack a 4D AthenaArray into a one-dimensional buffer

template <typename T> void PackData(AthenaArray<T> &src, T *buf,
                                    int sn, int en,
                                    int si, int ei, int sj, int ej, int sk, int ek,
                                    int &offset) {
  if (DBGPR_BUFFER_UTILS) {
    coutBoldBlue("\npd4:");
    src.print_all();
  }

  for (int n=sn; n<=en; ++n) {
    for (int k=sk; k<=ek; k++) {
      for (int j=sj; j<=ej; j++) {
#pragma omp simd
        for (int i=si; i<=ei; i++)
          buf[offset++] = src(n,k,j,i);
      }
    }
  }

  if (DBGPR_BUFFER_UTILS) {
    coutBoldBlue("\nbuf: ");
    for (int kk=0; kk<offset; kk++) {
      if (kk < offset - 1)
        printf("%1.3f, ", buf[kk]);
      else
        printf("%1.3f\n", buf[kk]);
    }

    coutBoldBlue("\nix: ");
    printf("(p, sn, en, si, ei, sj, ej, sk, ek)="
          "(%d, %d, %d, %d, %d, %d, %d, %d, %d)\n",
          offset, sn, en, si, ei, sj, ej, sk, ek);
  }
}

//----------------------------------------------------------------------------------------
//! \fn template <typename T> void PackData(AthenaArray<T> &src, T *buf,
//                      int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//  \brief pack a 3D AthenaArray into a one-dimensional buffer

template <typename T> void PackData(AthenaArray<T> &src, T *buf,
                                    int si, int ei, int sj, int ej, int sk, int ek,
                                    int &offset) {
  if (DBGPR_BUFFER_UTILS) {
    coutBoldBlue("\npd3 ");
    src.print_all();
  }

  for (int k=sk; k<=ek; k++) {
    for (int j=sj; j<=ej; j++) {
#pragma omp simd
      for (int i=si; i<=ei; i++)
        buf[offset++] = src(k, j, i);
    }
  }

  if (DBGPR_BUFFER_UTILS) {
    coutBoldBlue("\nix: ");
    printf("(p, si, ei, sj, ej, sk, ek)="
          "(%d, %d, %d, %d, %d, %d, %d)\n",
          offset, si, ei, sj, ej, sk, ek);
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
  if (DBGPR_BUFFER_UTILS) {
    coutBoldBlue("\nud4 [pre]");
    dst.print_all();
  }

  for (int n=sn; n<=en; ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma omp simd
        for (int i=si; i<=ei; ++i)
          dst(n,k,j,i) = buf[offset++];
      }
    }
  }

  if (DBGPR_BUFFER_UTILS) {
    coutBoldBlue("\nud4 [post]");
    dst.print_all();

    coutBoldBlue("\nix: ");
    printf("(p, sn, en, si, ei, sj, ej, sk, ek)="
          "(%d, %d, %d, %d, %d, %d, %d, %d, %d)\n",
          offset, sn, en, si, ei, sj, ej, sk, ek);
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
  if (DBGPR_BUFFER_UTILS) {
    coutBoldBlue("\nud3 [pre]");
    dst.print_all();
  }

  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma omp simd
      for (int i=si; i<=ei; ++i)
        dst(k,j,i) = buf[offset++];
    }
  }

  if (DBGPR_BUFFER_UTILS) {
    coutBoldBlue("\nud3 [post]");
    dst.print_all();

    coutBoldBlue("\nix: ");
    printf("(p, si, ei, sj, ej, sk, ek)="
          "(%d, %d, %d, %d, %d, %d, %d)\n",
          offset, si, ei, sj, ej, sk, ek);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn template <typename T> void UnpackDataAdd(T *buf, AthenaArray<T> &dst, int sn, int en,
//                        int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//  \brief unpack a one-dimensional buffer into a 4D AthenaArray additively

template <typename T> void UnpackDataAdd(T *buf, AthenaArray<T> &dst,
                                         int sn, int en,
                                         int si, int ei, int sj, int ej, int sk, int ek,
                                         int &offset) {
  if (DBGPR_BUFFER_UTILS) {
    coutBoldBlue("\nud4avg [pre]");
    dst.print_all();
  }

  for (int n=sn; n<=en; ++n) {
    for (int k=sk; k<=ek; ++k) {
      for (int j=sj; j<=ej; ++j) {
#pragma omp simd
        for (int i=si; i<=ei; ++i)
          dst(n,k,j,i) += buf[offset++];
      }
    }
  }

  if (DBGPR_BUFFER_UTILS) {
    coutBoldBlue("\nud4avg [post]");
    dst.print_all();

    coutBoldBlue("\nix: ");
    printf("(p, sn, en, si, ei, sj, ej, sk, ek)="
          "(%d, %d, %d, %d, %d, %d, %d, %d, %d)\n",
          offset, sn, en, si, ei, sj, ej, sk, ek);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn template <typename T> void UnpackDataAdd(T *buf, AthenaArray<T> &dst,
//                        int si, int ei, int sj, int ej, int sk, int ek, int &offset)
//  \brief unpack a one-dimensional buffer into a 3D AthenaArray additively

template <typename T> void UnpackDataAdd(T *buf, AthenaArray<T> &dst,
                                         int si, int ei, int sj, int ej, int sk, int ek,
                                         int &offset) {
  if (DBGPR_BUFFER_UTILS) {
    coutBoldBlue("\nud3avg [pre]");
    dst.print_all();
  }

  for (int k=sk; k<=ek; ++k) {
    for (int j=sj; j<=ej; ++j) {
#pragma omp simd
      for (int i=si; i<=ei; ++i)
        dst(k,j,i) += buf[offset++];
    }
  }

  if (DBGPR_BUFFER_UTILS) {
    coutBoldBlue("\nud3avg [post]");
    dst.print_all();

    coutBoldBlue("\nix: ");
    printf("(p, si, ei, sj, ej, sk, ek)="
          "(%d, %d, %d, %d, %d, %d, %d)\n",
          offset, si, ei, sj, ej, sk, ek);
  }

  return;
}

// provide explicit instantiation definitions (C++03) to allow the template definitions to
// exist outside of header file (non-inline), but still provide the requisite instances
// for other TUs during linking time (~13x files include "buffer_utils.hpp")

// 13x files include buffer_utils.hpp
template void UnpackData<Real>(Real *, AthenaArray<Real> &, int, int, int, int, int, int,
                               int, int,
                               int &);
template void UnpackData<Real>(Real *, AthenaArray<Real> &, int, int, int, int, int, int,
                               int &);

template void UnpackDataAdd<Real>(Real *, AthenaArray<Real> &, int, int, int, int, int, int,
                                  int, int,
                                  int &);
template void UnpackDataAdd<Real>(Real *, AthenaArray<Real> &, int, int, int, int, int, int,
                                  int &);

template void PackData<Real>(AthenaArray<Real> &, Real *, int, int, int, int, int, int,
                             int, int,
                             int &);
template void PackData<Real>(AthenaArray<Real> &, Real *, int, int, int, int, int, int,
                             int &);

} // end namespace BufferUtility
