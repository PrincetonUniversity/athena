//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena headers
#include "../athena.hpp"
#include "./utils.hpp"


//======================================================================================
//! \file rotate_vectors.cpp
//======================================================================================


//--------------------------------------------------------------------------------------

void RotateVec(const Real sint, const Real cost,
               const Real sinp, const Real cosp,
               Real &v1, Real &v2, Real &v3) {
  // vel1, vel2, vel3 are input
  // v1, v2, v3 are output
  // The two rotation matrix
  //R_1=
  //[cos_p  sin_p 0]
  //[-sin_p cos_p 0]
  //[0       0    1]

  //R_2=
  //[sin_t  0 cos_t]
  //[0      1    0]
  //[-cos_t 0 sin_t]

  // First apply R1, then apply R2
  Real newv1 =  cosp * v1 + sinp * v2;
  v2 = -sinp * v1 + cosp * v2;

  // now apply R2
  v1 =  sint * newv1 + cost * v3;
  Real newv3 = -cost * newv1 + sint * v3;
  v3 = newv3;
  return;
}


void InvRotateVec(const Real sint, const Real cost,
                 const Real sinp, const Real cosp,
                 Real &v1, Real &v2, Real &v3) {
  // vel1, vel2, vel3 are input
  // v1, v2, v3 are output
  // The two rotation matrix
  //R_1^-1=
  //[cos_p  -sin_p 0]
  //[sin_p cos_p 0]
  //[0       0    1]

  //R_2^-1=
  //[sin_t  0 -cos_t]
  //[0      1    0]
  //[cos_t 0 sin_t]

  // First apply R2^-1, then apply R1^-1
  Real newv1 = sint * v1 - cost * v3;
  v3 = cost * v1 + sint * v3;

  // now apply R1^-1
  v1 = cosp * newv1 - sinp * v2;
  Real newv2 = sinp * newv1 + cosp * v2;
  v2 = newv2;
}
