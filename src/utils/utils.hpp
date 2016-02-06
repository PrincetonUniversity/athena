#ifndef UTILS_HPP
#define UTILS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file utils.hpp
//  \brief prototypes of utility functions in utils/*.cpp
//======================================================================================

void ChangeRunDir(const char *pdir);
double ran2(long int *idum);
void ShowConfig();


namespace WallTimeLimit {
  void InitWTLimit(void);
  void SendWTLimit(int nwtlimit);
  bool TestWTLimit(int &nwtlimit);
  void FinalizeWTLimit(int wtflag);
}

#endif // UTILS_HPP
