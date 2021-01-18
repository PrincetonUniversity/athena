//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file poisson.cpp
//! \brief Problem generator to test Poisson's solver

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
#include "../fft/athena_fft.hpp"
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

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

void Mesh::InitUserMeshData(ParameterInput *pin) {
  Real four_pi_G = pin->GetReal("problem","four_pi_G");
  Real eps = pin->GetOrAddReal("problem","grav_eps", 0.0);
  SetFourPiG(four_pi_G);
  SetGravityThreshold(eps);
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real x0=0.0, y0=0.0, z0=0.0;
  Real x, y, z;
  Real r2; // r
  Real den(0.0), phia(0.0);
  RegionSize &mesh_size = pmy_mesh->mesh_size;

  Real x1size = mesh_size.x1max - mesh_size.x1min;
  Real x2size = mesh_size.x2max - mesh_size.x2min;
  Real x3size = mesh_size.x3max - mesh_size.x3min;

  Real four_pi_G = pin->GetReal("problem","four_pi_G");
  Real gconst = four_pi_G / (4.0*PI);

  int iprob = pin->GetOrAddInteger("problem","iprob",1);
  int nlim = pin->GetInteger("time","nlim");

  int dim = 1;
  if (mesh_size.nx2 > 1) dim=2;
  if (mesh_size.nx3 > 1) dim=3;

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        x = pcoord->x1v(i);
        y = pcoord->x2v(j);
        z = pcoord->x3v(k);
        r2 = SQR(x - x0);
        if (dim > 1) r2 += SQR(y - y0);
        if (dim > 2) r2 += SQR(z - z0);
        if (iprob == 1) {
          den = std::sin(TWO_PI*x/x1size);
          if (dim > 1) den *= std::sin(TWO_PI*y/x2size);
          if (dim > 2) den *= std::sin(TWO_PI*z/x3size); // dkx=2*PI/Nx
          phia = SQR(TWO_PI/x1size);
          if (dim > 1) phia += SQR(TWO_PI/x2size);
          if (dim > 2) phia += SQR(TWO_PI/x3size);
          phia = -den*four_pi_G/phia;
          den+=2.0;
        } else if (iprob == 2) {
          Real M = pin->GetOrAddReal("problem","M",1.0);
          Real a0 = pin->GetOrAddReal("problem","a0",1.0);
          Real da = 3.0*M/(4*PI*a0*a0*a0);
          den = da*std::pow((1.0+r2/SQR(a0)),-2.5);
          phia = - gconst*M/std::sqrt(r2+SQR(a0));
        } else if (iprob == 3) {
          Real a0 = pin->GetOrAddReal("problem","a0",1.0);
          den = (4.0*SQR(a0)*r2-6.0*a0)*std::exp(-a0*r2);
          phia = four_pi_G*std::exp(-a0*r2);
        }

        if (nlim > 0) {
          phydro->u(IDN,k,j,i) = den;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          phydro->u(IM3,k,j,i) = 0.0;
        } else {
          phydro->u(IDN,k,j,i) = den;
          phydro->u(IM1,k,j,i) = den;
          phydro->u(IM2,k,j,i) = phia;
          phydro->u(IM3,k,j,i) = 0.0;
        }
      }
    }
  }
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop()
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void MeshBlock::UserWorkInLoop() {
  // do nothing
  return;
}
//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  MeshBlock *pmb = my_blocks(0);
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  int cnt = (ke-ks+1)*(je-js+1)*(ie-is+1);

  int nlim = pin->GetInteger("time","nlim");

  if (SELF_GRAVITY_ENABLED) {
    if (nlim == 0) {
      // timing measure after loop
      int ncycle = pin->GetInteger("problem","ncycle");
      if (Globals::my_rank == 0) {
        std::cout << "=====================================================" << std::endl;
        std::cout << "Call Poisson Solver  " << ncycle << " times          " << std::endl;
        std::cout << "=====================================================" << std::endl;
      }
      clock_t tstart = clock();

#ifdef OPENMP_PARALLEL
      double omp_start_time = omp_get_wtime();
#endif

      for (int n=0; n < ncycle; n++) {
        for (int b=0; b<nblocal; ++b) {
          pmb = my_blocks(b);
          std::memset(pmb->pgrav->phi.data(), 0, pmb->pgrav->phi.GetSizeInBytes());
        }
        if (SELF_GRAVITY_ENABLED == 1) pfgrd->Solve(1,1);
        else if (SELF_GRAVITY_ENABLED == 2) pmgrd->Solve(1);
      }

#ifdef OPENMP_PARALLEL
      double omp_time = omp_get_wtime() - omp_start_time;;
#endif
      clock_t tstop = clock();
      float cpu_time = (tstop>tstart ? static_cast<Real>(tstop-tstart) : 1.0) /
                       static_cast<Real>(CLOCKS_PER_SEC);
      std::int64_t zones = GetTotalCells();
      std::int64_t mb_zones = GetTotalCells()/nbtotal*nblist[Globals::my_rank];
      double zc_cpus = static_cast<double>(mb_zones*ncycle)/cpu_time;
      double zc_cpus2 =
          static_cast<double>(mb_zones
                              *static_cast<double>(std::log2(mb_zones))
                              *ncycle)/cpu_time;

      if (Globals::my_rank == 0) {
        std::cout << "Timing Possison Solver                               " << std::endl;
        std::cout << "number of zones in Mesh = " << zones << std::endl;
        std::cout << "Mesh configuration = " << mesh_size.nx1 << "x"
                  << mesh_size.nx2 << "x" << mesh_size.nx3 << std::endl;
        std::cout << "number of zones in MeshBlock = " << mb_zones << std::endl;
        std::cout << "MeshBlock configuration = " << pmb->block_size.nx1 << "x"
                  << pmb->block_size.nx2 << "x" << pmb->block_size.nx3 << std::endl;
        std::cout << "number of processors = " << Globals::nranks << std::endl;
        std::cout << "processor configuration = "
                  << nrbx1 << "x" << nrbx2 << "x" << nrbx3 << std::endl;
        std::cout << "cpu time used  = " << cpu_time << std::endl;
        std::cout << "cpu time used/cycle  = " << cpu_time/static_cast<Real>(ncycle)
                  << std::endl;
        std::cout << "zone-cycles/cpu_second = " << zc_cpus << std::endl;
        std::cout << "zone-cycles(NlogN)/cpu_second = " << zc_cpus2 << std::endl;
#ifdef OPENMP_PARALLEL
        std::cout << "=====================================================" << std::endl;
        float zc_omps = static_cast<Real>(zones*ncycle)/omp_time;
        std::cout << "omp number of threads = " << GetNumMeshThreads() << std::endl;
        std::cout << "omp wtime used = " << omp_time << std::endl;
        std::cout << "zone-cycles/omp_wsecond = " << zc_omps << std::endl;
#endif
        std::cout << "=====================================================" << std::endl;
      }

      Real err1 = 0.0, maxphi = 0.0; // err2 = 0.0
      for (int b=0; b<nblocal; ++b) {
        pmb = my_blocks(b);
        Hydro *phydro = pmb->phydro;
        Gravity *pgrav = pmb->pgrav;
        for (int k=ks; k<=ke; ++k) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              err1 += std::abs(pgrav->phi(k,j,i) - phydro->u(IM2,k,j,i));
              //          err2 += pgrav->phi(k,j,i)/phydro->u(IM2,k,j,i);
              maxphi = std::max(pgrav->phi(k,j,i),maxphi);
            }
          }
        }
      }
#ifdef MPI_PARALLEL
      MPI_Allreduce(MPI_IN_PLACE, &err1, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &maxphi, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif

      err1 = err1/(static_cast<Real>(cnt*nbtotal));
      // err2 = err2/cnt;

      Real x1size = mesh_size.x1max - mesh_size.x1min;
      Real x2size = mesh_size.x2max - mesh_size.x2min;
      Real x3size = mesh_size.x3max - mesh_size.x3min;
      Real four_pi_G = pin->GetReal("problem", "four_pi_G");
      Real phiamp = SQR(TWO_PI/x1size);
      phiamp += SQR(TWO_PI/x2size);
      phiamp += SQR(TWO_PI/x3size);
      phiamp = 1.0*four_pi_G/phiamp;

      if (Globals::my_rank == 0) {
        std::cout << std::scientific
                  << std::setprecision(std::numeric_limits<Real>::max_digits10 - 1);
        std::cout << "=====================================================" << std::endl;
        std::cout << "L1 : " << err1 <<" MaxPhi: " << maxphi
                  << " Amp: " << phiamp << std::endl;
        std::cout << "=====================================================" << std::endl;
      }
    }
  }

  return;
}
