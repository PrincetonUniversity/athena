//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file binary_gravity.cpp
//  \brief Problem generator to test Multigrid Poisson's solver with Multipole Expansion

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstring>    // memset
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../gravity/fft_gravity.hpp"
#include "../gravity/gravity.hpp"
#include "../gravity/mg_gravity.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"
#include "../parameter_input.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#if SELF_GRAVITY_ENABLED != 2
#error "This problem generator requires Multigrid gravity solver."
#endif

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields."
#endif

Real four_pi_G;

void Mesh::InitUserMeshData(ParameterInput *pin) {
  four_pi_G = pin->GetReal("problem","four_pi_G");
  SetFourPiG(four_pi_G);
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for gravity from a binary
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real x0 = 12.0/1024.0, y0 = 0.0, z0 = 0.0, r = 6.0/1024.0;
  Real m1 = 2.0, m2 = 1.0;
  Real G = four_pi_G / (4.0 * PI);
  Real den1 = m1/(4.0*PI/3.0*r*r*r);
  Real den2 = m2/(4.0*PI/3.0*r*r*r);

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);

        Real r1 = sqrt(SQR(x-x0)+SQR(y-y0)+SQR(z-z0));
        Real r2 = sqrt(SQR(x+x0)+SQR(y+y0)+SQR(z+z0));

        if (r1 < r)
          phydro->u(IDN,k,j,i) = den1;
        else if (r2 < r)
          phydro->u(IDN,k,j,i) = den2;
        else
          phydro->u(IDN,k,j,i) = 1e-300;
        phydro->u(IEN,k,j,i) = phydro->u(IDN,k,j,i);
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = phydro->u(IEN,k,j,i);
      }
    }
  }

  // rescale the mass
  if (lid == pmy_mesh->nblocal - 1) {
    Real mass = 0.0;
    for (int b = 0; b < pmy_mesh->nblocal; ++b) {
      MeshBlock *pmb = pmy_mesh->my_blocks(b);
      Coordinates *pcoord = pmb->pcoord;
      Real vol = pcoord->dx1f(is)*pcoord->dx2f(js)*pcoord->dx3f(ks);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i)
            mass += pmb->phydro->u(IDN,k,j,i) * vol;
        }
      }
    }

#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &mass, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif

    if ((mass < (m1+m2)*0.7 || mass > (m1+m2)*1.3) && Globals::my_rank == 0)
      std::cout << "Too much or too little mass. Resolution is too low." << std::endl;

    Real fac = (m1+m2)/mass;
    for (int b = 0; b < pmy_mesh->nblocal; ++b) {
      MeshBlock *pmb = pmy_mesh->my_blocks(b);
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i)
            pmb->phydro->u(IDN,k,j,i) *= fac;
        }
      }
    }
  }
}


//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief post-processing error analysis
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  Real x0 = 12.0/1024.0, y0 = 0.0, z0 = 0.0, r = 6.0/1024.0;
  Real m1 = 2.0, m2 = 1.0;
  Real G = four_pi_G / (4.0 * PI);
  Real den1 = m1/(4.0*PI/3.0*r*r*r);
  Real den2 = m2/(4.0*PI/3.0*r*r*r);

  MeshBlock *pmb = my_blocks(0);
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  Real err = 0.0;

  for (int b = 0; b < nblocal; ++b) {
    pmb = my_blocks(b);
    Hydro *phydro = pmb->phydro;
    Gravity *pgrav = pmb->pgrav;
    Coordinates *pcoord = pmb->pcoord;
    Real vol = pcoord->dx1f(is)*pcoord->dx2f(js)*pcoord->dx3f(ks);
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
        for (int i = is; i <= ie; ++i) {
          Real x = pcoord->x1v(i);
          Real y = pcoord->x2v(j);
          Real z = pcoord->x3v(k);

          Real r1 = sqrt(SQR(x-x0)+SQR(y-y0)+SQR(z-z0));
          Real r2 = sqrt(SQR(x+x0)+SQR(y+y0)+SQR(z+z0));
          Real p1 = 0.0, p2 = 0.0;
          if (r1 > r)
            p1 = -G*m1/r1;
          else
            p1 = -G*PI*2.0/3.0*den1*(3.0*r*r-r1*r1);
          if (r2 > r)
            p2 = -G*m2/r2;
          else
            p2 = -G*PI*2.0/3.0*den2*(3.0*r*r-r2*r2);
          Real pot = p1 + p2;
          phydro->u(IM1,k,j,i) = pot;
          phydro->u(IM2,k,j,i) = std::abs((pot - pgrav->phi(k,j,i))/pot);
          err += phydro->u(IM2, k, j, i) * vol;
        }
      }
    }
  }
#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif

  Real x1size = mesh_size.x1max - mesh_size.x1min;
  Real x2size = mesh_size.x2max - mesh_size.x2min;
  Real x3size = mesh_size.x3max - mesh_size.x3min;
  Real tvol = x1size * x2size * x3size;
  err /= tvol;
  if (Globals::my_rank == 0) {
    std::cout << std::scientific
              << std::setprecision(std::numeric_limits<Real>::max_digits10 - 1);
    std::cout << "=====================================================" << std::endl;
    std::cout << "L1 : " << err << std::endl;
    std::cout << "=====================================================" << std::endl;
  }

  return;
}
