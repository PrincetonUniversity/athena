#ifndef UTILS_UTILS_HPP_
#define UTILS_UTILS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file utils.hpp
//! \brief prototypes of functions and class definitions for utils/*.cpp files

// C headers

// C++ headers
#include <csignal>   // sigset_t POSIX C extension
#include <cstdint>   // std::int64_t

// Athena++ headers

void ChangeRunDir(const char *pdir);
double ran2(std::int64_t *idum);
void ShowConfig();

//----------------------------------------------------------------------------------------
//! \namespace SignalHandler
//! \brief static data and functions that implement a simple signal handling system

namespace SignalHandler {
const int nsignal = 3;
static volatile int signalflag[nsignal];
const int ITERM = 0, IINT = 1, IALRM = 2;
static sigset_t mask;
void SignalHandlerInit();
int CheckSignalFlags();
int GetSignalFlag(int s);
void SetSignalFlag(int s);
void SetWallTimeAlarm(int t);
void CancelWallTimeAlarm();
} // namespace SignalHandler

#endif // UTILS_UTILS_HPP_
