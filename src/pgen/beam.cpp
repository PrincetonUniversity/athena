//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <fstream>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../parameter_input.hpp"

// File scope variables
static int ang;
static int octnum;

int RefinementCondition(MeshBlock *pmb);

//======================================================================================
/*! \file beam.cpp
 *  \brief Beam test for the radiative transfer module
 *
 *====================================================================================*/

void TwoBeams(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
              const AthenaArray<Real> &w, FaceField &b,
              AthenaArray<Real> &ir, Real time, Real dt,
              int is, int ie, int js, int je, int ks, int ke, int ngh);

void TwoBeamHydro(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive)
    EnrollUserRefinementCondition(RefinementCondition);

  EnrollUserBoundaryFunction(BoundaryFace::inner_x2, TwoBeamHydro);
  if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED)
    EnrollUserRadBoundaryFunction(BoundaryFace::inner_x2, TwoBeams);
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  ang = pin->GetOrAddInteger("problem","ang",0);
  octnum = pin->GetOrAddInteger("problem","octnum",0);
  return;
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real gamma = peos->GetGamma();

  // Initialize hydro variable
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = 1.0/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  // Now initialize opacity and specific intensity
  if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
    int nfreq = pnrrad->nfreq;
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ifr=0; ifr < nfreq; ++ifr) {
            pnrrad->sigma_s(k,j,i,ifr) = 0.0;
            pnrrad->sigma_a(k,j,i,ifr) = 0.0;
            pnrrad->sigma_pe(k,j,i,ifr) = 0.0;
            pnrrad->sigma_p(k,j,i,ifr) = 0.0;
          }
          for (int n=0; n<pnrrad->n_fre_ang; ++n) {
              pnrrad->ir(k,j,i,n) = 0.0;
          }
        }
      }
    }
  }
  return;
}

void TwoBeams(MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
              const AthenaArray<Real> &w, FaceField &b,
              AthenaArray<Real> &ir,
              Real time, Real dt,
              int is, int ie, int js, int je, int ks, int ke, int ngh) {
  int nang=prad->nang;
  int noct=prad->noct;
  int nfreq=prad->nfreq;
  int ang_oct=nang/noct;

  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=ngh; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real const &x1 = pco->x1v(i);
        Real const &x2 = pco->x2v(js-j);
        for (int ifr=0; ifr<nfreq; ++ifr) {
          for (int l=0; l<noct; ++l) {
            for (int n=0; n<ang_oct; ++n) {
              int n_ang=l*ang_oct + n;

              Real slope1=-prad->mu(1,k,js-j,i,0)/prad->mu(0,k,js-j,i,0);
              Real slope2=prad->mu(1,k,js-j,i,0)/prad->mu(0,k,js-j,i,0);
              Real dis1=std::abs(slope1*(x1-0.1)+(x2+2.0));
              Real dis2=std::abs(slope2*(x1+0.1)+(x2+2.0));
              if (ifr == 0) {
                if (((l==0)&&(n==0)&&(dis1<pco->dx1v(i))) ||
                    ((l==1)&&(n==0)&&(dis2<pco->dx1v(i)))) {
                  ir(k,js-j,i,n_ang+ifr*nang) = 10.0;
                } else {
                  ir(k,js-j,i,n_ang+ifr*nang) = 0.0;
                }
              } else {
                if (((l==0)&&(n==1)&&(dis1<pco->dx1v(i))) ||
                    ((l==1)&&(n==1)&&(dis2<pco->dx1v(i)))) {
                  ir(k,js-j,i,n_ang+ifr*nang) = 10.0;
                } else {
                  ir(k,js-j,i,n_ang+ifr*nang) = 0.0;
                }
              }
            }
          }
        }
      }
    }
  }
  return;
}


void TwoBeamHydro(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh) {
  for (int n=0; n<NHYDRO; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=ngh; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,k,js-j,i) = prim(n,k,js,i);
        }
      }
    }
  }
}

// refinement condition: density curvature
int RefinementCondition(MeshBlock *pmb) {
  Coordinates *pco = pmb->pcoord;
  Real xmin = pco->x1f(pmb->is);
  Real xmax = pco->x1f(pmb->ie+1);
  Real ymin = pco->x2f(pmb->js);
  Real ymax = pco->x2f(pmb->je+1);
  Real zmin = pco->x3f(pmb->ks);
  Real zmax = pco->x3f(pmb->ke+1);
  // refine : delta rho > 0.9*amp
  if (ymin > -1 && ymax < 1) return 1;
  return 0;
}
