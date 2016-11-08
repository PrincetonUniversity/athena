#ifndef UTILS_HPP
#define UTILS_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file utils.hpp
//  \brief prototypes of functions and class definitions for utils/*.cpp files

#include <csignal>

void ChangeRunDir(const char *pdir);
double ran2(long int *idum);
void ShowConfig();

namespace WallTimeLimit {
  void InitWTLimit(void);
  void SendWTLimit(int nwtlimit);
  bool TestWTLimit(int &nwtlimit);
  void FinalizeWTLimit(int wtflag);
}

//----------------------------------------------------------------------------------------
//! SignalHandler
//  \brief static data and functions that implement a simple signal handling system

namespace SignalHandler {
  const int nsignal = 2;
  static volatile sig_atomic_t signalflag[nsignal];
  const int ITERM = 0, IINT = 1;
  void SignalHandlerInit(void);
  void SynchronizeSignalFlag(void);
  int GetSignalFlag(int s);
  void SetSignalFlag(int s);
}

#endif // UTILS_HPP
