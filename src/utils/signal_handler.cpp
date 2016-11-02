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
//! \fn
//  \brief install handlers for selected signals

void SignalHandlerInit()
{
  sigterm_flag = 0;
  sigint_flag = 0;
  signal(SIGTERM, SetSigtermFlag);
  signal(SIGINT,  SetSigintFlag);
}

//----------------------------------------------------------------------------------------
//! \fn 
// \brief return specified signal flag, suitable synchronized across all MPI ranks

int GetSignalFlag(int sig)
{
  if (sig == SIGTERM) {
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &sigterm_flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    return sigterm_flag;

  } else if (sig == SIGINT) {
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &sigint_flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    return sigint_flag;

  } else {
    return 0;
  }
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief Sets flag for SIGTERM and reinstalls the signal handler function.

void SetSigtermFlag(int s){
  sigterm_flag = s;
  signal(s, SetSigtermFlag);   // Reinstall the signal handler function
  return;
}

//----------------------------------------------------------------------------------------
//! \fn
//  \brief Sets flag for SIGINT and reinstalls the signal handler function.

void SetSigintFlag(int s){
  sigint_flag = s;
  signal(s, SetSigintFlag);   // Reinstall the signal handler function
  return;
}

} // end of namespace
