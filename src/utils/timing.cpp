//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file timing.cpp
//! \brief utility function for timing

// C headers

// C++ headers
#include <ctime>      // clock(), CLOCKS_PER_SEC, clock_t
#include <fstream>
#include <iostream>
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

static bool newfile_ = true;
//----------------------------------------------------------------------------------------
//! \fn Real MarkTime()
//! \brief mark time using proper scheme

double MarkTime() {
#ifdef OPENMP_PARALLEL
  return omp_get_wtime();
#else
#ifdef MPI_PARALLEL
  return MPI_Wtime();
#else
  return static_cast<double> (clock())/static_cast<double> (CLOCKS_PER_SEC);
#endif
#endif
}

//----------------------------------------------------------------------------------------
//! \fn void OutputLoopTime(Real dt_array[])
//! \brief output loop time breakdown
void OutputLoopTime(const int ncycle, double dt_array[], std::string basename) {
#ifdef MPI_PARALLEL
  // pack array, MPI allreduce over array, then unpack into Mesh variables
  MPI_Allreduce(MPI_IN_PLACE, dt_array, 5, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  double time_per_step = dt_array[0] + dt_array[1] + dt_array[2]
                       + dt_array[3] + dt_array[4];
  if (Globals::my_rank == 0) {
    std::ofstream os;
    std::string fname;
    fname.assign(basename);
    fname.append(".loop_time.txt");
    // open 'loop_time.txt' file
    if (newfile_) {
      os.open(fname.c_str(), std::ofstream::out);
      newfile_ = false;
    } else {
      os.open(fname.c_str(), std::ofstream::app);
    }

    if (!os.is_open()) {
      std::cout << "### ERROR in function OutputLoopTime" << std::endl
                << "Cannot open " << fname << std::endl;
      return;
    }

    os << "ncycle=" << ncycle << ", time=" << time_per_step;
    os << ",Before=" << dt_array[0];
    os << ",TurbulenceDriver=" << dt_array[1];
    os << ",TimeIntegratorTaskList=" << dt_array[2];
    os << ",SelfGravity=" << dt_array[3];
    os << ",After=" << dt_array[4] << std::endl;

    os.close();
  }
}
