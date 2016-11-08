//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file signal_handler.cpp
//  \brief contains functions that implement a simple SignalHandler
//  These functions are based on TAG's signal handler written for Athena 8/19/2004
 
// C++ headers
#include <csignal>

// Athena++ headers
#include "utils.hpp"

namespace SignalHandler {

//----------------------------------------------------------------------------------------
//! \fn void SignalHandlerInit(void)
//  \brief install handlers for selected signals

void SignalHandlerInit(void)
{
  for(int n=0; n<nsignal; n++) signalflag[n]=0;
  signal(SIGTERM, SetSignalFlag);
  signal(SIGINT,  SetSignalFlag);
}

//----------------------------------------------------------------------------------------
//! \fn void SynchronizeSignalFlag(void)
// \brief synchronize the signal flags across all MPI ranks

void SynchronizeSignalFlag(void)
{
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, signalflag, nsignal, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int GetSignalFlag(int s)
//  \brief Gets a signal flag assuming the signalflag array is already synchronized.
//         Returns -1 if the specified signal is not handled.

int GetSignalFlag(int s)
{
  int ret=-1;
  switch(s) {
  case SIGTERM:
    ret=signalflag[ITERM];
    break;
  case SIGINT:
    ret=signalflag[IINT];
    break;
  default:
    // nothing
    break;
  }
  return ret;
}

//----------------------------------------------------------------------------------------
//! \fn void SetSignalFlag(int s)
//  \brief Sets signal flags and reinstalls the signal handler function.

void SetSignalFlag(int s)
{
  switch(s) {
  case SIGTERM:
    signalflag[ITERM]=1;
    signal(s, SetSignalFlag);
    break;
  case SIGINT:
    signalflag[IINT]=1;
    signal(s, SetSignalFlag);
    break;
  default:
    // nothing
    break;
  }
  return;
}

} // namespace SignalHandler
