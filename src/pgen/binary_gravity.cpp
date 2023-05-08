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
  Real x1 = 6.0/1024.0, x2 = -12.0/1024.0, y1 = 0.0, y2 = 0.0, z1 = 0.0, z2 = 0.0;
  Real r = 6.0/1024.0;
  Real m1 = 2.0, m2 = 1.0;
  Real G = four_pi_G / (4.0 * PI);
  Real den1 = m1/(4.0*PI/3.0*r*r*r);
  Real den2 = m2/(4.0*PI/3.0*r*r*r);
  Real dx = pcoord->dx1f(is);
  Real dd = 0.1*dx;
  Real dv = 1e-3;
  Real dr = 0.6*std::sqrt(3.0)*dx;

  for (int k = ks; k <= ke; ++k) {
    Real z = pcoord->x3v(k);
    for (int j = js; j <= je; ++j) {
     Real y = pcoord->x2v(j);
      for (int i = is; i <= ie; ++i) {
        Real x = pcoord->x1v(i);
        Real r1 = std::sqrt(SQR(x-x1)+SQR(y-y1)+SQR(z-z1));
        Real r2 = std::sqrt(SQR(x-x2)+SQR(y-y2)+SQR(z-z2));
        phydro->u(IDN,k,j,i) = 1e-300;

        if (r1 < r + dr) {
          if (r1 < r - dr) {
            phydro->u(IDN,k,j,i) = den1;
          } else {
            Real xf = pcoord->x1f(i);
            Real yf = pcoord->x2f(j);
            Real zf = pcoord->x3f(k);
            for (int kk = 0; kk < 10; ++kk) {
              Real zz = zf + (kk+0.5)*dd;
              for (int jj = 0; jj < 10; ++jj) {
                Real yy = yf + (jj+0.5)*dd;
                for (int ii = 0; ii < 10; ++ii) {
                  Real xx = xf + (ii+0.5)*dd;
                  Real rr = std::sqrt(SQR(xx-x1)+SQR(yy-y1)+SQR(zz-z1));
                  if (rr < r)
                    phydro->u(IDN,k,j,i) += dv*den1;
                }
              }
            }
          }
        }
        if (r2 < r + dr) {
          if (r2 < r - dr) {
            phydro->u(IDN,k,j,i) = den2;
          } else {
            Real xf = pcoord->x1f(i);
            Real yf = pcoord->x2f(j);
            Real zf = pcoord->x3f(k);
            for (int kk = 0; kk < 10; ++kk) {
              Real zz = zf + (kk+0.5)*dd;
              for (int jj = 0; jj < 10; ++jj) {
                Real yy = yf + (jj+0.5)*dd;
                for (int ii = 0; ii < 10; ++ii) {
                  Real xx = xf + (ii+0.5)*dd;
                  Real rr = std::sqrt(SQR(xx-x2)+SQR(yy-y2)+SQR(zz-z2));
                  if (rr < r)
                    phydro->u(IDN,k,j,i) += dv*den2;
                }
              }
            }
          }
        }
        phydro->u(IEN,k,j,i) = phydro->u(IDN,k,j,i);
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
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
  Real x1 = 6.0/1024.0, x2 = -12.0/1024.0, y1 = 0.0, y2 = 0.0, z1 = 0.0, z2 = 0.0;
  Real r = 6.0/1024.0;
  Real m1 = 2.0, m2 = 1.0;
  Real G = four_pi_G / (4.0 * PI);
  Real den1 = m1/(4.0*PI/3.0*r*r*r);
  Real den2 = m2/(4.0*PI/3.0*r*r*r);

  MeshBlock *pmb = my_blocks(0);
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  Real err1 = 0.0, err2 = 0.0;

  if (nlim == 0) {
    // timing measure after loop
    int ncycle = pin->GetOrAddInteger("problem", "ncycle", 1);
    if (Globals::my_rank == 0) {
      std::cout << "=====================================================" << std::endl;
      std::cout << "Call Poisson Solver  " << ncycle << " times          " << std::endl;
      std::cout << "=====================================================" << std::endl;
    }

#ifdef MPI_PARALLEL
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#ifdef OPENMP_PARALLEL
    double omp_start_time = omp_get_wtime();
#endif
    clock_t tstart = clock();

    AthenaArray<Real> nsol;
    nsol.NewAthenaArray(nblocal,pmb->ncells3,pmb->ncells2,pmb->ncells1);
    for (int b=0; b<nblocal; ++b) {
      pmb = my_blocks(b);
      for (int k=0; k< pmb->ncells3; ++k) {
        for (int j=0; j< pmb->ncells2; ++j) {
          for (int i=0; i< pmb->ncells1; ++i)
            nsol(b,k,j,i)=pmb->pgrav->phi(k,j,i);
        }
      }
    }
    for (int n=0; n < ncycle; n++) {
      for (int b=0; b<nblocal; ++b) {
        pmb = my_blocks(b);
        pmb->pgrav->phi.ZeroClear();
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
    double zc_cpus = static_cast<double>(zones*ncycle)/cpu_time;

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
#ifdef OPENMP_PARALLEL
      std::cout << "=====================================================" << std::endl;
      double zc_omps = static_cast<Real>(zones*ncycle)/omp_time;
      std::cout << "omp number of threads = " << GetNumMeshThreads() << std::endl;
      std::cout << "omp wtime used = " << omp_time << std::endl;
      std::cout << "zone-cycles/omp_wsecond = " << zc_omps << std::endl;
#endif
      std::cout << "=====================================================" << std::endl;
    }

    Real temax = 0.0, nemax = 0.0, defmax = 0.0;
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
            Real dx = pcoord->dx1f(i);
            Real dy = pcoord->dx2f(j);
            Real dz = pcoord->dx3f(k);

            Real r1 = std::sqrt(SQR(x-x1)+SQR(y-y1)+SQR(z-z1));
            Real r2 = std::sqrt(SQR(x-x2)+SQR(y-y2)+SQR(z-z2));

            Real ax = -(pgrav->phi(k,j,i+1)-pgrav->phi(k,j,i-1))/(2.0*dx);
            Real ay = -(pgrav->phi(k,j+1,i)-pgrav->phi(k,j-1,i))/(2.0*dy);
            Real az = -(pgrav->phi(k+1,j,i)-pgrav->phi(k-1,j,i))/(2.0*dz);

            Real p1, p2, pot0, ax1, ay1, az1, ax2, ay2, az2, ax0, ay0, az0;
            if (r1 > r) {
              p1 = -G*m1/r1;
              ax1 = -G*m1/(r1*r1*r1)*(x-x1);
              ay1 = -G*m1/(r1*r1*r1)*(y-y1);
              az1 = -G*m1/(r1*r1*r1)*(z-z1);
            } else {
              p1 = -G*PI*2.0/3.0*den1*(3.0*r*r-r1*r1);
              ax1 = -G*PI*4.0/3.0*den1*(x-x1);
              ay1 = -G*PI*4.0/3.0*den1*(y-y1);
              az1 = -G*PI*4.0/3.0*den1*(z-z1);
            }
            if (r2 > r) {
              p2 = -G*m2/r2;
              ax2 = -G*m2/(r2*r2*r2)*(x-x2);
              ay2 = -G*m2/(r2*r2*r2)*(y-y2);
              az2 = -G*m2/(r2*r2*r2)*(z-z2);
            } else {
              p2 = -G*PI*2.0/3.0*den2*(3.0*r*r-r2*r2);
              ax2 = -G*PI*4.0/3.0*den2*(x-x2);
              ay2 = -G*PI*4.0/3.0*den2*(y-y2);
              az2 = -G*PI*4.0/3.0*den2*(z-z2);
            }
            pot0 = p1 + p2;
            ax0 = ax1 + ax2;
            ay0 = ay1 + ay2;
            az0 = az1 + az2;

            Real perr = (pot0 - pgrav->phi(k,j,i))/pot0;
            Real perr2 = (nsol(b,k,j,i) - pgrav->phi(k,j,i))/pot0;
            Real aerr = std::sqrt((SQR(ax-ax0)+SQR(ay-ay0)+SQR(az-az0))
                      / (SQR(ax0)+SQR(ay0)+SQR(az0)));
            phydro->u(IM1,k,j,i) = perr;
            phydro->u(IM2,k,j,i) = perr2;
            phydro->u(IM3,k,j,i) = aerr;
            phydro->u(IEN,k,j,i) = pot0;
            err1 += std::abs(perr) * vol;
            err2 += aerr * vol;
            if (std::abs(z) < dz) {
              if (std::abs(perr) > temax)
                temax = std::abs(perr);
              if (std::abs(perr2) > nemax)
                nemax = std::abs(perr2);
              if (std::abs(pgrav->def(k,j,i)) > defmax)
                defmax = std::abs(pgrav->def(k,j,i));
            }
          }
        }
      }
    }
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &err1, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &err2, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &temax, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &nemax, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &defmax, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif

    Real x1size = mesh_size.x1max - mesh_size.x1min;
    Real x2size = mesh_size.x2max - mesh_size.x2min;
    Real x3size = mesh_size.x3max - mesh_size.x3min;
    Real tvol = x1size * x2size * x3size;
    err1 /= tvol;
    err2 /= tvol;
    err1 = std::sqrt(err1);
    err2 = std::sqrt(err2);
    if (Globals::my_rank == 0) {
      std::cout << std::scientific
                << std::setprecision(std::numeric_limits<Real>::max_digits10 - 1);
      std::cout << "=====================================================" << std::endl;
      std::cout << "Potential    L2       : " << err1 << std::endl;
      std::cout << "Acceleration L2       : " << err2 << std::endl;
      std::cout << "Max True Error        : " << temax << std::endl;
      std::cout << "Max Discretized Error : " << nemax << std::endl;
      std::cout << "Max Defect            : " << defmax << std::endl;
      std::cout << "=====================================================" << std::endl;
    }
  }

  return;
}
