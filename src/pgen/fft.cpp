//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file fft.cpp
//  \brief Problem generator for complex-to-complex FFT test.
//

// C++ headers
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <ctime>
#include <iomanip>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../fft/athena_fft.hpp"
#include "../mesh/mesh.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  Coordinates *pcoord = pblock->pcoord;
  Real x0=0.0, y0=0.0, z0=0.0;
  int is=pblock->is, ie=pblock->ie;
  int js=pblock->js, je=pblock->je;
  int ks=pblock->ks, ke=pblock->ke;

  AthenaArray<Real> src, dst;
  LogicalLocation &loc = pblock->loc;
  RegionSize &block_size = pblock->block_size;
  int nx1=block_size.nx1+2*NGHOST;
  int nx2=block_size.nx2+2*NGHOST;
  int nx3=block_size.nx3+2*NGHOST;

  src.NewAthenaArray(nx3,nx2,nx1);
  dst.NewAthenaArray(2,nx3,nx2,nx1);

  if (FFT_ENABLED) {
    FFTDriver *pfftd;
    pfftd = new FFTDriver(this, pin);
    pfftd->InitializeFFTBlock(true);
    pfftd->QuickCreatePlan();

    FFTBlock *pfft = pfftd->pmy_fb;
  // Repeating FFTs for timing
    if (Globals::my_rank == 0) {
      std::cout << "=====================================================" << std::endl;
      std::cout << "Initialize...                                        " << std::endl;
      std::cout << "=====================================================" << std::endl;
    }

    pfftd->QuickCreatePlan();
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          Real r2;
          if (COORDINATE_SYSTEM == "cartesian") {
            Real x = pcoord->x1v(i);
            Real y = pcoord->x2v(j);
            Real z = pcoord->x3v(k);
            r2 = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
          }
          src(k,j,i)= std::exp(-r2);
        }
      }
    }

    pfft->LoadSource(src,1,NGHOST,loc,block_size);

    if (Globals::my_rank == 0) {
      std::cout << "=====================================================" << std::endl;
      std::cout << "End Initialization...                                " << std::endl;
      std::cout << "=====================================================" << std::endl;
    }

    int ncycle = pin->GetOrAddInteger("problem","ncycle",100);
    if (Globals::my_rank == 0) {
      std::cout << "=====================================================" << std::endl;
      std::cout << "Execute FFT " << ncycle << "                         " << std::endl;
      std::cout << "=====================================================" << std::endl;
    }

    clock_t tstart = clock();
#ifdef OPENMP_PARALLEL
    double omp_start_time = omp_get_wtime();
#endif
    for (int n=0; n <= ncycle; n++) {
      pfft->ExecuteForward();
      pfft->ExecuteBackward();
    }
#ifdef OPENMP_PARALLEL
    double omp_time = omp_get_wtime() - omp_start_time;;
#endif
    clock_t tstop = clock();
    float cpu_time = (tstop>tstart ? static_cast<Real>(tstop-tstart) : 1.0) /
        static_cast<Real>(CLOCKS_PER_SEC);
    int64_t zones = GetTotalCells();
    float zc_cpus = static_cast<Real>(zones*ncycle)/cpu_time;

    if (Globals::my_rank == 0) {
      std::cout << std::endl << "cpu time used  = " << cpu_time << std::endl;
      std::cout << "zone-cycles/cpu_second = " << zc_cpus << std::endl;
#ifdef OPENMP_PARALLEL
      float zc_omps = static_cast<Real>(zones*ncycle)/omp_time;
      std::cout << std::endl << "omp wtime used = " << omp_time << std::endl;
      std::cout << "zone-cycles/omp_wsecond = " << zc_omps << std::endl;
#endif
    }
// Reset everything and do FFT once for error estimation
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      Real r2;
      if (COORDINATE_SYSTEM == "cartesian") {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        r2 = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
      }
      src(k,j,i) = std::exp(-r2);
    }}}

    pfft->LoadSource(src,1,NGHOST,loc,block_size);
    pfft->ExecuteForward();
    pfft->ApplyKernel(0);
    pfft->ExecuteBackward();
    pfft->RetrieveResult(dst,2,NGHOST,loc,block_size);

    Real err1=0.0,err2=0.0;
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      err1 += std::abs(dst(0,k,j,i) - src(k,j,i));
      err2 += std::abs(dst(1,k,j,i));
    }}}
    if (Globals::my_rank == 0) {
      std::cout << std::setprecision(15) << std::scientific;
      std::cout << "=====================================================" << std::endl;
      std::cout << "Error for Real: " << err1 <<" Imaginary: " << err2 << std::endl;
      std::cout << "=====================================================" << std::endl;
    }
  }

  return;
}
