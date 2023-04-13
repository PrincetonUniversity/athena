/* ----------------------------------------------------------------------
   fftMPI - library for computing 3d/2d FFTs in parallel
   http://fftmpi.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright 2018 National Technology & Engineering Solutions of
   Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.
   This software is distributed under the modified Berkeley Software
   Distribution (BSD) License.

   See the README file in the top-level fftMPI directory.
------------------------------------------------------------------------- */

// Error class

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "error.h"

using namespace FFTMPI_NS;

/* ---------------------------------------------------------------------- */

Error::Error(MPI_Comm world_caller)
{
  world = world_caller;
}

/* ----------------------------------------------------------------------
   must be called by all procs in world
   shuts down MPI and exits
------------------------------------------------------------------------- */

void Error::all(const char *str)
{
  MPI_Barrier(world);

  int me;
  MPI_Comm_rank(world,&me);
  if (me == 0) printf("ERROR: %s\n",str);
  MPI_Finalize();
  exit(1);
}

/* ----------------------------------------------------------------------
   called by one proc in world
   forces abort of entire world if any proc in world calls
------------------------------------------------------------------------- */

void Error::one(const char *str)
{
  int me;
  MPI_Comm_rank(world,&me);
  printf("ERROR on proc %d: %s\n",me,str);
  MPI_Abort(world,1);
}

/* ----------------------------------------------------------------------
   called by one proc in world
   only write to screen if non-NULL on this proc since could be file
------------------------------------------------------------------------- */

void Error::warning(const char *str)
{
  printf("WARNING: %s\n",str);
}
