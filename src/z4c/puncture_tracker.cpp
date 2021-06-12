//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file puncture_tracker.cpp
//  \brief implementation of functions in the PunctureTracker classes

#include <cmath>
#include <sstream>
#include <unistd.h>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#include "puncture_tracker.hpp"

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/lagrange_interp.hpp"

#include "z4c.hpp"
#include "z4c_macro.hpp"


//----------------------------------------------------------------------------------------
PunctureTracker::PunctureTracker(Mesh * pmesh, ParameterInput * pin, int n):
    owns_puncture{false}, pos{NAN, NAN, NAN}, betap{NAN, NAN, NAN}, pmesh{pmesh} {
  ofname = pin->GetOrAddString("z4c", "filename", "puncture_");
  ofname += std::to_string(n) + ".txt";
  pos[0] = pin->GetOrAddReal("z4c", "bh_" + std::to_string(n) + "_x", 0.0);
  pos[1] = pin->GetOrAddReal("z4c", "bh_" + std::to_string(n) + "_y", 0.0);
  pos[2] = pin->GetOrAddReal("z4c", "bh_" + std::to_string(n) + "_z", 0.0);
  if (0 == Globals::my_rank) {
    // check if output file already exists
    if (access(ofname.c_str(), F_OK) == 0) {
      pofile = fopen(ofname.c_str(), "a");
    }
    else {
      pofile = fopen(ofname.c_str(), "w");
      if (NULL == pofile) {
        std::stringstream msg;
        msg << "### FATAL ERROR in PunctureTracker constructor" << std::endl;
        msg << "Could not open file '" << ofname << "' for writing!";
        throw std::runtime_error(msg.str().c_str());
      }
      fprintf(pofile, "# 1:iter 2:time 3:x 4:y 5:z 6:betax 7:betay 8:betaz\n");
    }
  }
}

//----------------------------------------------------------------------------------------
PunctureTracker::~PunctureTracker() {
  if (0 == Globals::my_rank) {
    fclose(pofile);
  }
}

//----------------------------------------------------------------------------------------
void PunctureTracker::InterpolateShift(MeshBlock * pmb, AthenaArray<Real> & u) {
  Z4c::Z4c_vars z4c;
  pmb->pz4c->SetZ4cAliases(u, z4c);

  if (pos[0] >= pmb->block_size.x1min && pos[0] < pmb->block_size.x1max &&
      pos[1] >= pmb->block_size.x2min && pos[1] < pmb->block_size.x2max &&
      pos[2] >= pmb->block_size.x3min && pos[2] < pmb->block_size.x3max) {
    owns_puncture = true;

    Real const origin[3] = {
      pmb->pcoord->x1f(0),
      pmb->pcoord->x1f(1),
      pmb->pcoord->x1f(2),
    };
    Real const delta[3] = {
      pmb->pcoord->dx1f(0),
      pmb->pcoord->dx2f(0),
      pmb->pcoord->dx3f(0),
    };
    int const size[3] = {
      pmb->nverts1,
      pmb->nverts2,
      pmb->nverts3,
    };

    LagrangeInterpND<2*NGHOST-1, 3> linterp(origin, delta, size, pos);
    for (int a = 0; a < NDIM; ++a) {
      Real & beta = z4c.beta_u(a, 0, 0, 0);
      betap[a] = linterp.eval(&beta);
    }
  }
}

//----------------------------------------------------------------------------------------
void PunctureTracker::EvolveTracker() {
#ifdef MPI_PARALLEL
#ifndef NDEBUG
  int count = owns_puncture;
  if (0 == Globals::my_rank) {
    MPI_Reduce(MPI_IN_PLACE, &count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  else {
    MPI_Reduce(&count, &count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  assert(count == 1);   // the puncture should be in exactly one location
#else
  int count = 1;
#endif // NDEBUG
#else
#ifndef NDEBUG
  int count = 1;
#else
  int count = owns_puncture;
  assert(count == 1);
#endif // NDEBUG
#endif // MPI_PARALLEL

  if (owns_puncture) {
    for (int a = 0; a < NDIM; ++a) {
      pos[a] -= pmesh->dt * betap[a];
    }
  }

#ifdef MPI_PARALLEL
  Real buf[2*NDIM] = {0., 0., 0., 0., 0., 0.};
  if (owns_puncture) {
    buf[0] = pos[0];
    buf[1] = pos[1];
    buf[2] = pos[2];
    buf[3] = betap[0];
    buf[4] = betap[1];
    buf[5] = betap[2];
  }
  MPI_Allreduce(MPI_IN_PLACE, buf, 2*NDIM, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  pos[0] = buf[0]/count;
  pos[1] = buf[1]/count;
  pos[2] = buf[2]/count;
  betap[0] = buf[3]/count;
  betap[1] = buf[3]/count;
  betap[2] = buf[5]/count;
#endif

  // After the puncture has moved it might have changed ownership
  owns_puncture = false;
}

//----------------------------------------------------------------------------------------
void PunctureTracker::WriteTracker(int iter, Real time) const {
  if (0 == Globals::my_rank) {
    fprintf(pofile, "%d %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
        iter, time, pos[0], pos[1], pos[2], betap[0], betap[1], betap[2]);
  }
}
