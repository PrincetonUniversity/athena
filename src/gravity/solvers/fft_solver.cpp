
//========================================================================================
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// C/C++ headers
#include <cmath>

// Athena++ headers
#include "../gravity.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../../bvals/bvals.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../fft/athena_fft.hpp"

//----------------------------------------------------------------------------------------
//! \fn  void Gravity::Initialize
//  \brief Initialize FFT gravity solver

void Gravity::Initialize(ParameterInput *pin)
{
  MeshBlock *pmb=pmy_block;
  AthenaFFT *pfft=pmb->pfft;

  pfft->Initialize();
  pfft->fplan = pfft->QuickCreatePlan(AthenaFFTForward);
  pfft->bplan = pfft->QuickCreatePlan(AthenaFFTBackward);
}

//----------------------------------------------------------------------------------------
//! \fn  void Gravity::Solver
//  \brief Calculate potential from density using FFT

void Gravity::Solver(const AthenaArray<Real> &den)
{
  MeshBlock *pmb=pmy_block;
  AthenaFFT *pfft=pmb->pfft;
  Coordinates *pcoord = pmb->pcoord;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real dx1=pcoord->dx1v(is);
  Real dx2=pcoord->dx2v(js);
  Real dx3=pcoord->dx3v(ks);
  Real dx1sq = SQR(dx1), dx2sq = SQR(dx2), dx3sq = SQR(dx3);

  Real pcoeff;
// Copy phi to phi_old

// Fill the work array

  for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=is; i<=ie; ++i){
        long int idx=pfft->GetIndex(k-ks,j-js,i-is);
        pfft->work[idx][0] = (den(k,j,i) - grav_mean_rho);
        pfft->work[idx][1] = 0.0;
      }
    }
  }

  pfft->Execute(pfft->fplan);

// Multiply kernel coefficient

  for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=is; i<=ie; ++i){
        long int gidx = pfft->GetGlobalIndex(k-ks,j-js,i-is); 
        if(gidx == 0){
          pcoeff = 0.0;
        } else { 
          pcoeff = ((2.0*std::cos(((i-is)+pfft->idisp)*pfft->dkx)-2.0)/dx1sq);
          if(pfft->dim > 1)
            pcoeff += ((2.0*std::cos(((j-js)+pfft->jdisp)*pfft->dky)-2.0)/dx2sq);
          else if (pfft->dim > 2)
            pcoeff += ((2.0*std::cos(((k-ks)+pfft->kdisp)*pfft->dkz)-2.0)/dx3sq);
        }
        long int idx=pfft->GetIndex(k-ks,j-js,i-is);
        pfft->work[idx][0] *= pcoeff;
        pfft->work[idx][1] *= pcoeff;
      }
    }
  }
  
  pfft->Execute(pfft->bplan);

// Return phi

  for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=is; i<=ie; ++i){
        long int idx=pfft->GetIndex(k-ks,j-js,i-is);
        phi(k,j,i) = four_pi_gconst * pfft->work[idx][0]/pfft->gcnt;
      }
    }
  }
}
