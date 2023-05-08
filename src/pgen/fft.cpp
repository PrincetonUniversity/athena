//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file fft.cpp
//! \brief Problem generator for complex-to-complex FFT test.
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
//! \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  Real x0=0.0, y0=0.0, z0=0.0;
  AthenaArray<Real> src,dst;

#ifdef FFT
  FFTDriver *pfftd;
  pfftd = new FFTDriver(this, pin);
  pfftd->InitializeFFTBlock(true);

  FFTBlock *pfft = pfftd->pmy_fb;
  // Repeating FFTs for timing
  if (Globals::my_rank == 0) {
    std::cout << "=====================================================" << std::endl;
    std::cout << "Initialize...                                        " << std::endl;
    std::cout << "=====================================================" << std::endl;
  }

  pfftd->QuickCreatePlan();

  MeshBlock *pmb = my_blocks(0);
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  src.NewAthenaArray(nblocal, pmb->ncells3, pmb->ncells2, pmb->ncells1);
  dst.NewAthenaArray(nblocal, 2, pmb->ncells3, pmb->ncells2, pmb->ncells1);

  for (int b=0; b<nblocal; ++b) {
    MeshBlock *pmb = my_blocks(b);
    Coordinates *pcoord = pmb->pcoord;
    LogicalLocation &loc = pmb->loc;
    RegionSize &block_size = pmb->block_size;

    AthenaArray<Real> src_;
    src_.InitWithShallowSlice(src, 4, b, 1);
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
          src_(k,j,i) = std::exp(-r2);
        }
      }
    }
    pfft->LoadSource(src_, 0, NGHOST, loc, block_size);
  }

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

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &zc_cpus, 1, MPI_ATHENA_REAL, MPI_MIN, MPI_COMM_WORLD);
#endif

  if (Globals::my_rank == 0) {
    std::cout << std::endl << "cpu time used  = " << cpu_time << std::endl;
    std::cout << "zone-cycles/cpu_second = " << zc_cpus << std::endl;
#ifdef OPENMP_PARALLEL
    float zc_omps = static_cast<Real>(zones*ncycle)/omp_time;
    std::cout << std::endl << "omp wtime used = " << omp_time << std::endl;
    std::cout << "zone-cycles/omp_wsecond = " << zc_omps << std::endl;
#endif
  }

  for (int b=0; b<nblocal; ++b) {
    MeshBlock *pmb = my_blocks(b);
    Coordinates *pcoord = pmb->pcoord;
    LogicalLocation &loc = pmb->loc;
    RegionSize &block_size = pmb->block_size;

    AthenaArray<Real> src_;
    src_.InitWithShallowSlice(src, 4, b, 1);
    pfft->LoadSource(src_, 0, NGHOST, loc, block_size);
  }

  // Do FFT once for error estimation
  pfft->ExecuteForward();
  pfft->ApplyKernel(0);
  pfft->ExecuteBackward();

  Real err1=0.0, err2=0.0;
  for (int b=0; b<nblocal; ++b) {
    MeshBlock *pmb = my_blocks(b);
    Coordinates *pcoord = pmb->pcoord;
    LogicalLocation &loc = pmb->loc;
    RegionSize &block_size = pmb->block_size;

    AthenaArray<Real> dst_,src_;
    src_.InitWithShallowSlice(src, 4, b, 1);
    dst_.InitWithShallowSlice(dst, 5, b, 1);
    pfft->RetrieveResult(dst_, 1, NGHOST, loc, block_size);

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          err1 += std::abs(dst_(0,k,j,i) - src_(k,j,i));
          err2 += std::abs(dst_(1,k,j,i));
        }
      }
    }
  }
  err1 = err1*pfft->norm_factor_;
  err2 = err2*pfft->norm_factor_;

#ifdef MPI_PARALLEL
  MPI_Allreduce(MPI_IN_PLACE, &err1, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &err2, 1, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif

  if (Globals::my_rank == 0) {
    std::cout << std::scientific
              << std::setprecision(std::numeric_limits<Real>::max_digits10 - 1);
    std::cout << "=====================================================" << std::endl;
    std::cout << "Error for Real: " << err1 <<" Imaginary: " << err2 << std::endl;
    std::cout << "=====================================================" << std::endl;

    // open output file and write out errors
    std::string fname;
    fname.assign("fft-errors.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = std::fopen(fname.c_str(),"r")) != nullptr) {
      if ((pfile = std::freopen(fname.c_str(),"a",pfile)) == nullptr) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Error output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }

      // The file does not exist -- open the file in write mode and add headers
    } else {
      if ((pfile = std::fopen(fname.c_str(),"w")) == nullptr) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Error output file could not be opened" <<std::endl;
        ATHENA_ERROR(msg);
      }
      std::fprintf(pfile,"# Nx1  Nx2  Nx3  Ncycle ZCS-min Real-Error Imag-Error\n");
    }

    // write errors
    std::fprintf(pfile,"%d  %d",mesh_size.nx1,mesh_size.nx2);
    std::fprintf(pfile,"  %d  %d",mesh_size.nx3,ncycle);
    std::fprintf(pfile,"  %e  %e  %e",zc_cpus,err1,err2);
    std::fprintf(pfile,"\n");
    std::fclose(pfile);
  }

#endif // FFT

  return;
}
