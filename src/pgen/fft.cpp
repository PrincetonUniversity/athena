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
#include "omp.h"
#endif

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  AthenaFFT *pfft = pblock->pfft;
  Coordinates *pcoord = pblock->pcoord;
  Real x0=0.0, y0=0.0, z0=0.0;
  int is=pblock->is, ie=pblock->ie;
  int js=pblock->js, je=pblock->je;
  int ks=pblock->ks, ke=pblock->ke;

  if(FFT_ENABLED){
    int ncycle =  pin->GetOrAddInteger("problem","ncycle",100);
  // Repeating FFTs for timing
    std::cout << "=====================================================" << std::endl;
    pfft->Initialize();
    pfft->fplan = pfft->QuickCreatePlan(AthenaFFTForward);
    pfft->bplan = pfft->QuickCreatePlan(AthenaFFTBackward);
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      Real r2;
      if (COORDINATE_SYSTEM == "cartesian") {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        r2 = sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
      }
      int idx=pfft->GetIndex(k-ks,j-js,i-is);
      pfft->work[idx][0] = std::exp(-r2);
      pfft->work[idx][1] = 0.0;
    }}}

    clock_t tstart = clock();
#ifdef OPENMP_PARALLEL
    double omp_start_time = omp_get_wtime();
#endif
    for (int n=0; n <= ncycle; n++) {
      pfft->Execute(pfft->fplan);
      pfft->Execute(pfft->bplan);
    }
#ifdef OPENMP_PARALLEL
    double omp_time = omp_get_wtime() - omp_start_time;;
#endif
    clock_t tstop = clock();
    float cpu_time = (tstop>tstart ? (float)(tstop-tstart) : 1.0)/(float)CLOCKS_PER_SEC;
    int64_t zones = GetTotalCells();
    float zc_cpus = (float)(zones*ncycle)/cpu_time;

    std::cout << std::endl << "cpu time used  = " << cpu_time << std::endl;
    std::cout << "zone-cycles/cpu_second = " << zc_cpus << std::endl;
#ifdef OPENMP_PARALLEL
    float zc_omps = (float)(zones*ncycle)/omp_time;
    std::cout << std::endl << "omp wtime used = " << omp_time << std::endl;
    std::cout << "zone-cycles/omp_wsecond = " << zc_omps << std::endl;
#endif
    std::cout << "=====================================================" << std::endl;
// Reset everything and do FFT once for error estimation
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      Real r2;
      if (COORDINATE_SYSTEM == "cartesian") {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        r2 = sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
      }
      int idx=pfft->GetIndex(k-ks,j-js,i-is);
      pfft->work[idx][0] = std::exp(-r2);
      pfft->work[idx][1] = 0.0;
    }}}


    pfft->Execute(pfft->fplan);
    pfft->Execute(pfft->bplan);

    Real err1=0.0,err2=0.0;
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      Real r2;
      if (COORDINATE_SYSTEM == "cartesian") {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        r2 = sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
      }
      int idx=pfft->GetIndex(k-ks,j-js,i-is);
      err1 += std::abs(pfft->work[idx][0]/pfft->gcnt - std::exp(-r2));
      err2 += std::abs(pfft->work[idx][1]/pfft->gcnt);
    }}}
    std::cout << "Error for Real: " << err1 <<" Imaginary: " << err2 << std::endl;
    std::cout << "=====================================================" << std::endl;
  }

  return;
}
