//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file signal_handler.cpp
//! \brief contains functions that implement a simple SignalHandler
//! These functions are based on TAG's signal handler written for Athena 8/19/2004

// C headers
#include <unistd.h>   // alarm() Unix OS utility; not in C standard --> no <cunistd>

// C++ headers
// first 2x macros and signal() are the only ISO C features; rest are POSIX C extensions
#include <csignal>    // SIGTERM, SIGINT, SIGALARM, signal(), sigemptyset(), ...
#include <iostream>

// Athena++ headers
#include "../defs.hpp"
#include "../globals.hpp"
#include "utils.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


namespace SignalHandler {

//----------------------------------------------------------------------------------------
//! \fn void SignalHandlerInit()
//! \brief install handlers for selected signals

void SignalHandlerInit() {
  for (int n=0; n<nsignal; n++) {
    signalflag[n]=0;
  }
  // C++11 standard guarantees that <csignal> places C-standard signal.h contents in std::
  // namespace. POSIX C extensions are likely only placed in global namespace (not std::)
  std::signal(SIGTERM, SetSignalFlag);
  std::signal(SIGINT,  SetSignalFlag);
  std::signal(SIGALRM, SetSignalFlag);

  // populate set of signals to block while the handler is running; prevent premption
  sigemptyset(&mask);
  sigaddset(&mask, SIGTERM);
  sigaddset(&mask, SIGINT);
  sigaddset(&mask, SIGALRM);
}

//----------------------------------------------------------------------------------------
//! \fn int CheckSignalFlags()
//! \brief Synchronize and check signal flags and return true if any of them is caught

int CheckSignalFlags() {
  // Currently, only checking for nonzero return code at the end of each timestep in
  // main.cpp; i.e. if an issue prevents a process from reaching the end of a cycle, the
  // signals will never be handled by that process / the solver may hang
  int ret = 0;
  sigprocmask(SIG_BLOCK, &mask, nullptr);
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE,
                const_cast<void *>(reinterpret_cast<volatile void *>(signalflag)),
                nsignal, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
  for (int n=0; n<nsignal; n++)
    ret += signalflag[n];
  sigprocmask(SIG_UNBLOCK, &mask, nullptr);
  return ret;
}

//----------------------------------------------------------------------------------------
//! \fn int GetSignalFlag(int s)
//! \brief Gets a signal flag assuming the signalflag array is already synchronized.
//!        Returns -1 if the specified signal is not handled.

int GetSignalFlag(int s) {
  int ret=-1;
  switch(s) {
    case SIGTERM:
      ret=signalflag[ITERM];
      break;
    case SIGINT:
      ret=signalflag[IINT];
      break;
    case SIGALRM:
      ret=signalflag[IALRM];
      break;
    default:
      // nothing
      break;
  }
  return ret;
}

//----------------------------------------------------------------------------------------
//! \fn void SetSignalFlag(int s)
//! \brief Sets signal flags and reinstalls the signal handler function.

void SetSignalFlag(int s) {
  // Signal handler functions must have C linkage; C++ linkage is implemantation-defined
  switch(s) {
    case SIGTERM:
      signalflag[ITERM]=1;
      signal(s, SetSignalFlag);
      break;
    case SIGINT:
      signalflag[IINT]=1;
      signal(s, SetSignalFlag);
      break;
    case SIGALRM:
      signalflag[IALRM]=1;
      signal(s, SetSignalFlag);
      break;
    default:
      // nothing
      break;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void SetWallTimeAlarm(int t)
//! \brief Set the wall time limit alarm

void SetWallTimeAlarm(int t) {
  alarm(t);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void CancelWallTimeAlarm()
//! \brief Cancel the wall time limit alarm

void CancelWallTimeAlarm() {
  alarm(0);
  return;
}

} // namespace SignalHandler
