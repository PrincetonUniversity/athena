//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_crdiffusion.cpp
//! \brief create multigrid solver for gravity

// C headers

// C++ headers
#include <algorithm>
#include <iostream>
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"
#include "../parameter_input.hpp"
#include "../task_list/crdiffusion_task_list.hpp"
#include "mg_crdiffusion.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;

//----------------------------------------------------------------------------------------
//! \fn MGCRDiffusionDriver::MGCRDiffusionDriver(Mesh *pm, ParameterInput *pin)
//! \brief MGCRDiffusionDriver constructor

MGCRDiffusionDriver::MGCRDiffusionDriver(Mesh *pm, ParameterInput *pin)
    : MultigridDriver(pm, pm->MGCRDiffusionBoundaryFunction_,
                      pm->MGCRDiffusionSourceMaskFunction_, 1) {
  eps_ = pin->GetOrAddReal("crdiffusion", "threshold", -1.0);
  niter_ = pin->GetOrAddInteger("crdiffusion", "niteration", -1);
  ffas_ = pin->GetOrAddBoolean("crdiffusion", "fas", ffas_);
  std::string m = pin->GetOrAddString("crdiffusion", "mgmode", "none");
  std::transform(m.begin(), m.end(), m.begin(), ::tolower);
  if (m == "fmg") {
    mode_ = 0;
  } else if (m == "mgi") {
    mode_ = 1; // Iterative
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in MGCRDiffusionDriver::MGCRDiffusionDriver" << std::endl
        << "The \"mgmode\" parameter in the <gravity> block is invalid." << std::endl
        << "FMG: Full Multigrid + Multigrid iteration (default)" << std::endl
        << "MGI: Multigrid Iteration" << std::endl;
    ATHENA_ERROR(msg);
  }
  if (eps_ < 0.0 && niter_ < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MGCRDiffusionDriver::MGCRDiffusionDriver" << std::endl
        << "Either \"threshold\" or \"niteration\" parameter must be set "
        << "in the <crdiffusion> block." << std::endl
        << "When both parameters are specified, \"niteration\" is ignored." << std::endl
        << "Set \"threshold = 0.0\" for automatic convergence control." << std::endl;
    ATHENA_ERROR(msg);
  }
  mg_mesh_bcs_[inner_x1] =
              GetMGBoundaryFlag(pin->GetOrAddString("crdiffusion", "ix1_bc", "none"));
  mg_mesh_bcs_[outer_x1] =
              GetMGBoundaryFlag(pin->GetOrAddString("crdiffusion", "ox1_bc", "none"));
  mg_mesh_bcs_[inner_x2] =
              GetMGBoundaryFlag(pin->GetOrAddString("crdiffusion", "ix2_bc", "none"));
  mg_mesh_bcs_[outer_x2] =
              GetMGBoundaryFlag(pin->GetOrAddString("crdiffusion", "ox2_bc", "none"));
  mg_mesh_bcs_[inner_x3] =
              GetMGBoundaryFlag(pin->GetOrAddString("crdiffusion", "ix3_bc", "none"));
  mg_mesh_bcs_[outer_x3] =
              GetMGBoundaryFlag(pin->GetOrAddString("crdiffusion", "ox3_bc", "none"));
  CheckBoundaryFunctions();

  mgtlist_ = new MultigridTaskList(this);

  // Allocate the root multigrid
  mgroot_ = new MGCRDiffusion(this, nullptr);

  crtlist_ = new CRDiffusionBoundaryTaskList(pin, pm);
}


//----------------------------------------------------------------------------------------
//! \fn MGCRDiffusionDriver::~MGCRDiffusionDriver()
//! \brief MGCRDiffusionDriver destructor

MGCRDiffusionDriver::~MGCRDiffusionDriver() {
  delete crtlist_;
  delete mgroot_;
  delete mgtlist_;
}


//----------------------------------------------------------------------------------------
//! \fn MGCRDiffusion::MGCRDiffusion(MultigridDriver *pmd, MeshBlock *pmb)
//! \brief MGCRDiffusion constructor

MGCRDiffusion::MGCRDiffusion(MultigridDriver *pmd, MeshBlock *pmb)
  : Multigrid(pmd, pmb, 1, 1) {
  for (int l = 0; l < nlevel_; l++) {
    int ll=nlevel_-1-l;
    int ncx=(size_.nx1>>ll)+2*ngh_;
    int ncy=(size_.nx2>>ll)+2*ngh_;
    int ncz=(size_.nx3>>ll)+2*ngh_;
    D_[l].NewAthenaArray(NCOEFF,ncz,ncy,ncx);
    nlambda_[l].NewAthenaArray(ncz,ncy,ncx);
  }

  btype = btypef = BoundaryQuantity::mg;
  pmgbval = new MGBoundaryValues(this, mg_block_bcs_);
}


//----------------------------------------------------------------------------------------
//! \fn MGCRDiffusion::~MGCRDiffusion()
//! \brief MGCRDiffusion deconstructor

MGCRDiffusion::~MGCRDiffusion() {
  delete [] D_;
  delete [] nlambda_;
  delete pmgbval;
}


//----------------------------------------------------------------------------------------
//! \fn void MGCRDiffusionDriver::RestrictCoefficientsOctets()
//! \brief restrict coefficients in Octets for the CR diffusion equation

void MGCRDiffusionDriver::RestrictCoefficients() {
  if (nreflevel_ > 0) {
    const int &ngh = mgroot_->ngh_;
    for (int l = nreflevel_ - 1; l >= 1; --l) {  // fine octets to coarse octets
#pragma omp parallel for num_threads(nthreads_)
      for (int o = 0; o < noctets_[l]; ++o) {
        const LogicalLocation &loc = octets_[l][o].loc;
        LogicalLocation cloc;
        cloc.lx1 = (loc.lx1 >> 1);
        cloc.lx2 = (loc.lx2 >> 1);
        cloc.lx3 = (loc.lx3 >> 1);
        cloc.level = loc.level - 1;
        int oid = octetmap_[l-1][cloc];
        int oi = (static_cast<int>(loc.lx1) & 1) + ngh;
        int oj = (static_cast<int>(loc.lx2) & 1) + ngh;
        int ok = (static_cast<int>(loc.lx3) & 1) + ngh;
        
/*        for (int v = 0; v < NCOEFF; ++v)
          octets_[l-1][oid].D_(v, ok, oj, oi) = RestrictOne(octets_[l][o].D_,
                                                            v, ngh, ngh, ngh);
        octets_[l-1][oid].nlambda_(ok, oj, oi) = RestrictOne(octets_[l][o].nlambda_,
                                                             0, ngh, ngh, ngh);*/

      }
    }
#pragma omp parallel for num_threads(nthreads_)
    for (int o = 0; o < noctets_[0]; ++o) { // octets to the root grid
      const LogicalLocation &loc = octets_[0][o].loc;
      
//      for (int v = 0; v < nvar_; ++v)
//        Real r = RestrictOne(octets_[0][o].D_, v, ngh, ngh, ngh);//******
    }
  }

  return;

}


//----------------------------------------------------------------------------------------
//! \fn void MGCRDiffusionDriver::Solve(int stage)
//! \brief load the data and solve

void MGCRDiffusionDriver::Solve(int stage) {
  // Construct the Multigrid array
  vmg_.clear();
  for (int i = 0; i < pmy_mesh_->nblocal; ++i)
    vmg_.push_back(pmy_mesh_->my_blocks(i)->pcrdiff->pmg);

  // load the source
#pragma omp parallel for num_threads(nthreads_)
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    MGCRDiffusion *pmg = static_cast<MGCRDiffusion*>(*itr);
    // assume all the data are located on the same node
    pmg->LoadSource(pmg->pmy_block_->pcrdiff->ecr, 0, NGHOST, 1.0);
    if (mode_ == 1) // use the previous timestep data as the initial guess
      pmg->LoadFinestData(pmg->pmy_block_->pcrdiff->ecr, 0, NGHOST);
    pmg->LoadCoefficients(pmg->pmy_block_->pcrdiff->D,
                          pmg->pmy_block_->pcrdiff->nlambda, NGHOST);
  }

  SetupMultigrid();
  RestrictCoefficients();

  if (mode_ == 0) {
    SolveFMGCycle();
  } else {
    if (eps_ >= 0.0)
      SolveIterative();
    else
      SolveIterativeFixedTimes();
  }

  // Return the result
#pragma omp parallel for num_threads(nthreads_)
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    CRDiffusion *pcrdiff = pmg->pmy_block_->pcrdiff;
    pmg->RetrieveResult(pcrdiff->ecr, 0, NGHOST);
    if(pcrdiff->output_defect)
      pmg->RetrieveDefect(pcrdiff->def, 0, NGHOST);
  }

  crtlist_->DoTaskListOneStage(pmy_mesh_, stage);

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGCRDiffusion::RestrictCoefficients()
//! \brief restrict coefficients within a MGCRDiffusion object

void MGCRDiffusion::RestrictCoefficients() {
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  for (current_level_=nlevel_-1; current_level_>0; current_level_--) {
    int ll=nlevel_-current_level_;
    ie=is+(size_.nx1>>ll)-1, je=js+(size_.nx2>>ll)-1, ke=ks+(size_.nx3>>ll)-1;
    Restrict(D_[current_level_-1], D_[current_level_],
             NCOEFF, is, ie, js, je, ks, ke, false);
    Restrict(nlambda_[current_level_-1], nlambda_[current_level_],
             1, is, ie, js, je, ks, ke, false);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGCRDiffusion::Smooth(AthenaArray<Real> &u, const AthenaArray<Real> &src,
//!          int rlev, int il, int iu, int jl, int ju, int kl, int ku, int color, bool th)
//! \brief Implementation of the Red-Black Gauss-Seidel Smoother
//!        rlev = relative level from the finest level of this Multigrid block

void MGCRDiffusion::Smooth(AthenaArray<Real> &u, const AthenaArray<Real> &src, int rlev,
                int il, int iu, int jl, int ju, int kl, int ku, int color, bool th) {
  Real dx;
  if (rlev <= 0) dx = rdx_*static_cast<Real>(1<<(-rlev));
  else           dx = rdx_/static_cast<Real>(1<<rlev);
  Real dx2 = SQR(dx);
  Real isix = omega_/6.0;
  color ^= pmy_driver_->coffset_;
  if (th == true && (ku-kl) >=  minth_) {
#pragma omp parallel for num_threads(pmy_driver_->nthreads_)
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        int c = (color + k + j) & 1;
#pragma ivdep
        for (int i=il+c; i<=iu; i+=2);
      }
    }
  } else {
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        int c = (color + k + j) & 1;
#pragma ivdep
        for (int i=il+c; i<=iu; i+=2);
      }
    }
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGCRDiffusion::CalculateDefect(AthenaArray<Real> &def,
//!          const AthenaArray<Real> &u, const AthenaArray<Real> &src, int rlev,
//!          int il, int iu, int jl, int ju, int kl, int ku, bool th)
//! \brief Implementation of the Defect calculation
//!        rlev = relative level from the finest level of this Multigrid block

void MGCRDiffusion::CalculateDefect(AthenaArray<Real> &def, const AthenaArray<Real> &u,
                                const AthenaArray<Real> &src, int rlev,
                                int il, int iu, int jl, int ju, int kl, int ku, bool th) {
  Real dx;
  if (rlev <= 0) dx = rdx_*static_cast<Real>(1<<(-rlev));
  else           dx = rdx_/static_cast<Real>(1<<rlev);
  Real idx2 = 1.0/SQR(dx);
  if (th == true && (ku-kl) >=  minth_) {
#pragma omp parallel for num_threads(pmy_driver_->nthreads_)
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
#pragma omp simd
        for (int i=il; i<=iu; i++);
      }
    }
  } else {
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
#pragma omp simd
        for (int i=il; i<=iu; i++);
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGCRDiffusion::CalculateFASRHS(AthenaArray<Real> &src,
//!                         const AthenaArray<Real> &u, int rlev,
//!                         int il, int iu, int jl, int ju, int kl, int ku, bool th)
//! \brief Implementation of the RHS calculation for FAS
//!        rlev = relative level from the finest level of this Multigrid block

void MGCRDiffusion::CalculateFASRHS(AthenaArray<Real> &src, const AthenaArray<Real> &u,
                int rlev, int il, int iu, int jl, int ju, int kl, int ku, bool th) {
  Real dx;
  if (rlev <= 0) dx = rdx_*static_cast<Real>(1<<(-rlev));
  else           dx = rdx_/static_cast<Real>(1<<rlev);
  Real idx2 = 1.0/SQR(dx);
  if (th == true && (ku-kl) >=  minth_) {
#pragma omp parallel for num_threads(pmy_driver_->nthreads_)
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
#pragma omp simd
        for (int i=il; i<=iu; i++);
      }
    }
  } else {
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
#pragma omp simd
        for (int i=il; i<=iu; i++);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MGCRDiffusion::LoadCoefficients(const AthenaArray<Real> &D,
//!                                          const AthenaArray<Real> &nlambda, int ngh)
//! \brief Load coefficients of the diffusion and source terms

void MGCRDiffusion::LoadCoefficients(const AthenaArray<Real> &D,
                                     const AthenaArray<Real> &nlambda, int ngh) {
  AthenaArray<Real> &Dm=D_[nlevel_-1];
  AthenaArray<Real> &nl=nlambda_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+size_.nx1-1, je=js+size_.nx2-1, ke=ks+size_.nx3-1;
  for (int mk=ks; mk<=ke; ++mk) {
    int k = mk - ks + ngh;
    for (int mj=js; mj<=je; ++mj) {
      int j = mj - js + ngh;
#pragma omp simd
      for (int mi=is; mi<=ie; ++mi) {
        int i = mi - is + ngh;
        Dm(XX,mk,mj,mi) = D(XX,k,j,i);
        Dm(XY,mk,mj,mi) = D(XY,k,j,i);
        Dm(XZ,mk,mj,mi) = D(XZ,k,j,i);
        Dm(YY,mk,mj,mi) = D(YY,k,j,i);
        Dm(YZ,mk,mj,mi) = D(YZ,k,j,i);
        Dm(ZZ,mk,mj,mi) = D(ZZ,k,j,i);
        nl(mk,mj,mi) = nlambda(k,j,i);
      }
    }
  }
  return;
}

