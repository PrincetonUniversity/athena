//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file wtlimit.cpp 
//  \brief sharing the wall time limit between processes

// Athena headers
#include "../athena.hpp"
#include "../globals.hpp"
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

namespace WallTimeLimit {

#ifdef MPI_PARALLEL
MPI_Request *wtreq;
int nwtlbuf;
#endif

void InitWTLimit(void)
{
#ifdef MPI_PARALLEL
  if(Globals::my_rank==0) {
    wtreq=new MPI_Request[Globals::nranks-1];
    for(int i=0; i<Globals::nranks-1; i++) {
      wtreq[i]=MPI_REQUEST_NULL;
      MPI_Send_init(&nwtlbuf,1,MPI_INT,i+1,(int)TAG_WTLIM,MPI_COMM_WORLD,&wtreq[i]);
    }
  }
#endif
  return;
}


void SendWTLimit(int nwtlimit)
{
#ifdef MPI_PARALLEL
  if(Globals::my_rank==0) {
    nwtlbuf=nwtlimit;
    MPI_Startall(Globals::nranks-1,wtreq);
  }
#endif
  return;
}


bool TestWTLimit(int &nwtlimit)
{
#ifdef MPI_PARALLEL
  if(Globals::my_rank!=0) {
    int test;
    MPI_Iprobe(0,(int)TAG_WTLIM,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
    if(test==false) return false;
    MPI_Recv(&nwtlbuf,1,MPI_INT,0,(int)TAG_WTLIM,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    nwtlimit=nwtlbuf;
  }
#endif
  return true;
}

void FinalizeWTLimit(int wtflag)
{
#ifdef MPI_PARALLEL
  if(Globals::my_rank==0) {
    if(wtflag==2) // sent
      MPI_Waitall(Globals::nranks-1,wtreq,MPI_STATUSES_IGNORE);
    delete [] wtreq;
  }
#endif
  return;
}


}

