//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

// Athena headers
#include "../athena.hpp"
#include "../globals.hpp"
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//======================================================================================
//! \file wtlimit.cpp 
//  \brief sharing the wall time limit between processes
//======================================================================================

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
      MPI_Send_init(&nwtlbuf,1,MPI_INT,i+1,(int)tag_wtlimit,MPI_COMM_WORLD,&wtreq[i]);
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
    MPI_Iprobe(0,(int)tag_wtlimit,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
    if(test==false) return false;
    MPI_Recv(&nwtlbuf,1,MPI_INT,0,(int)tag_wtlimit,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
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

