//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file fft.cpp
//  \brief Problem generator for complex-to-complex FFT test.
//

// C headers

// C++ headers
#include <cmath>
#include <ctime>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../fft/athena_fft.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

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
  AthenaArray<Real> src, dst;
  MeshBlock *pmb = my_blocks(0);
  // TODO(changgoo): this does NOT assume 3D anymore, but need to check that 2D works
  src.NewAthenaArray(pmb->ncells3, pmb->ncells2, pmb->ncells1);
  dst.NewAthenaArray(2, pmb->ncells3, pmb->ncells2, pmb->ncells1);
#ifdef FFT
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
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          Real x = pcoord->x1v(i);
          Real y = pcoord->x2v(j);
          Real z = pcoord->x3v(k);
          r2 = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        }
        src(k,j,i) = std::exp(-r2);
      }
    }
  }

  pfft->LoadSource(src, 0, NGHOST, loc, block_size);

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
  std::int64_t zones = GetTotalCells();
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
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          Real x = pcoord->x1v(i);
          Real y = pcoord->x2v(j);
          Real z = pcoord->x3v(k);
          r2 = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        }
        src(k,j,i) = std::exp(-r2);
      }
    }
  }

  pfft->LoadSource(src, 0, NGHOST, loc, block_size);
  pfft->ExecuteForward();
  pfft->ApplyKernel(0);
  pfft->ExecuteBackward();
  pfft->RetrieveResult(dst, 1, NGHOST, loc, block_size);

  Real err1=0.0, err2=0.0;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        err1 += std::abs(dst(0,k,j,i) - src(k,j,i));
        err2 += std::abs(dst(1,k,j,i));
      }
    }
  }
  if (Globals::my_rank == 0) {
    std::cout << std::scientific
              << std::setprecision(std::numeric_limits<Real>::max_digits10 - 1);
    std::cout << "=====================================================" << std::endl;
    std::cout << "Error for Real: " << err1 <<" Imaginary: " << err2 << std::endl;
    std::cout << "=====================================================" << std::endl;
  }
#endif // FFT

  return;
}
