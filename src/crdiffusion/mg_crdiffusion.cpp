//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_crdiffusion.cpp
//! \brief create multigrid solver for comic-ray diffusion

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
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../multigrid/multigrid.hpp"
#include "../parameter_input.hpp"
#include "../task_list/crdiffusion_task_list.hpp"
#include "crdiffusion.hpp"
#include "mg_crdiffusion.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;

namespace {
  AthenaArray<Real> *temp; // temporary data for the Jacobi iteration
}

//----------------------------------------------------------------------------------------
//! \fn MGCRDiffusionDriver::MGCRDiffusionDriver(Mesh *pm, ParameterInput *pin)
//! \brief MGCRDiffusionDriver constructor

MGCRDiffusionDriver::MGCRDiffusionDriver(Mesh *pm, ParameterInput *pin)
    : MultigridDriver(pm, pm->MGCRDiffusionBoundaryFunction_,
                          pm->MGCRDiffusionCoeffBoundaryFunction_,
                          pm->MGCRDiffusionSourceMaskFunction_,
                          pm->MGCRDiffusionCoeffMaskFunction_,
                          1, NCOEFF, NMATRIX) {
  eps_ = pin->GetOrAddReal("crdiffusion", "threshold", -1.0);
  niter_ = pin->GetOrAddInteger("crdiffusion", "niteration", -1);
  ffas_ = pin->GetOrAddBoolean("crdiffusion", "fas", ffas_);
  fsubtract_average_ = false;
  omega_ = pin->GetOrAddReal("crdiffusion", "omega", 1.0);
  fsteady_ = pin->GetOrAddBoolean("crdiffusion", "steady", false);
  npresmooth_ = pin->GetOrAddReal("crdiffusion", "npresmooth", 2);
  npostsmooth_ = pin->GetOrAddReal("crdiffusion", "npostsmooth", 2);
  fshowdef_ = pin->GetOrAddBoolean("crdiffusion", "show_defect", fshowdef_);
  std::string smoother = pin->GetOrAddString("crdiffusion", "smoother", "jacobi-rb");
  if (smoother == "jacobi-rb") {
    fsmoother_ = 1;
    redblack_ = true;
  } else if (smoother == "jacobi-double") {
    fsmoother_ = 0;
    redblack_ = true;
  } else { // jacobi
    fsmoother_ = 0;
    redblack_ = false;
  }
  std::string prol = pin->GetOrAddString("crdiffusion", "prolongation", "trilinear");
  if (prol == "tricubic")
    fprolongation_ = 1;

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

  int nth = 1;
#ifdef OPENMP_PARALLEL
  nth = omp_get_max_threads();
#endif
  temp = new AthenaArray<Real>[nth];
  int nx = std::max(pmy_mesh_->block_size.nx1, pmy_mesh_->nrbx1) + 2*mgroot_->ngh_;
  int ny = std::max(pmy_mesh_->block_size.nx2, pmy_mesh_->nrbx2) + 2*mgroot_->ngh_;
  int nz = std::max(pmy_mesh_->block_size.nx3, pmy_mesh_->nrbx3) + 2*mgroot_->ngh_;
  for (int n = 0; n < nth; ++n)
    temp[n].NewAthenaArray(nz, ny, nx);
}


//----------------------------------------------------------------------------------------
//! \fn MGCRDiffusionDriver::~MGCRDiffusionDriver()
//! \brief MGCRDiffusionDriver destructor

MGCRDiffusionDriver::~MGCRDiffusionDriver() {
  delete crtlist_;
  delete mgroot_;
  delete mgtlist_;
  delete [] temp;
}


//----------------------------------------------------------------------------------------
//! \fn MGCRDiffusion::MGCRDiffusion(MGCRDiffusionDriver *pmd, MeshBlock *pmb)
//! \brief MGCRDiffusion constructor

MGCRDiffusion::MGCRDiffusion(MGCRDiffusionDriver *pmd, MeshBlock *pmb)
  : Multigrid(pmd, pmb, 1), omega_(pmd->omega_), fsmoother_(pmd->fsmoother_) {
  btype = btypef = BoundaryQuantity::mg;
  pmgbval = new MGBoundaryValues(this, mg_block_bcs_);
}


//----------------------------------------------------------------------------------------
//! \fn MGCRDiffusion::~MGCRDiffusion()
//! \brief MGCRDiffusion deconstructor

MGCRDiffusion::~MGCRDiffusion() {
  delete pmgbval;
}


//----------------------------------------------------------------------------------------
//! \fn void MGCRDiffusionDriver::AddCRSource(const AthenaArray<Real> &src,
//!                                           int ngh, Real dt)
//! \brief Add the cosmic-ray source term

void MGCRDiffusion::AddCRSource(const AthenaArray<Real> &src, int ngh, Real dt) {
  AthenaArray<Real> &dst=src_[nlevel_-1];
  int is, ie, js, je, ks, ke;
  is=js=ks=ngh_;
  ie=is+size_.nx1-1, je=js+size_.nx2-1, ke=ks+size_.nx3-1;
  if (!(static_cast<MGCRDiffusionDriver*>(pmy_driver_)->fsteady_)) {
    for (int mk=ks; mk<=ke; ++mk) {
      int k = mk - ks + ngh;
      for (int mj=js; mj<=je; ++mj) {
        int j = mj - js + ngh;
#pragma omp simd
        for (int mi=is; mi<=ie; ++mi) {
          int i = mi - is + ngh;
          dst(mk,mj,mi) += dt * src(k,j,i);
        }
      }
    }
  } else {
    for (int mk=ks; mk<=ke; ++mk) {
      int k = mk - ks + ngh;
      for (int mj=js; mj<=je; ++mj) {
        int j = mj - js + ngh;
#pragma omp simd
        for (int mi=is; mi<=ie; ++mi) {
          int i = mi - is + ngh;
          dst(mk,mj,mi) = src(k,j,i);
        }
      }
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MGCRDiffusionDriver::Solve(int stage, Real dt)
//! \brief load the data and solve

void MGCRDiffusionDriver::Solve(int stage, Real dt) {
  // Construct the Multigrid array
  vmg_.clear();
  for (int i = 0; i < pmy_mesh_->nblocal; ++i)
    vmg_.push_back(pmy_mesh_->my_blocks(i)->pcrdiff->pmg);

  // load the source
#pragma omp parallel for num_threads(nthreads_)
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    MGCRDiffusion *pmg = static_cast<MGCRDiffusion*>(*itr);
    // assume all the data are located on the same node
    CRDiffusion *pcrdiff = pmg->pmy_block_->pcrdiff;
    Hydro *phydro = pmg->pmy_block_->phydro;
    Field *pfield = pmg->pmy_block_->pfield;
    pcrdiff->CalculateCoefficients(phydro->w, pfield->bcc);
    if (!fsteady_)
      pmg->LoadSource(pcrdiff->ecr, 0, NGHOST, 1.0);
    if (mode_ == 1) // load the current data as the initial guess
      pmg->LoadFinestData(pcrdiff->ecr, 0, NGHOST);
    pmg->LoadCoefficients(pcrdiff->coeff, NGHOST);
    pmg->AddCRSource(pcrdiff->source, NGHOST, dt);
  }

  if (dt > 0.0 || fsteady_) {
    SetupMultigrid(dt, false);
    if (mode_ == 0) {
      SolveFMGCycle();
    } else {
      if (eps_ >= 0.0)
        SolveIterative();
      else
        SolveIterativeFixedTimes();
    }
  } else { // just copy trivial solution and set boundaries
    SetupMultigrid(dt, true);
    if (mode_ != 1) {
#pragma omp parallel for num_threads(nthreads_)
      for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
        MGCRDiffusion *pmg = static_cast<MGCRDiffusion*>(*itr);
        AthenaArray<Real> &ecr = pmg->GetCurrentData();
        AthenaArray<Real> &ecr0 = pmg->GetCurrentSource();
        ecr = ecr0;
      }
    }
    mgtlist_->SetMGTaskListBoundaryCommunication();
    mgtlist_->DoTaskListOneStage(this);
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
//! \fn void MGCRDiffusion::Smooth(AthenaArray<Real> &u, const AthenaArray<Real> &src,
//!            const AthenaArray<Real> &coeff, const AthenaArray<Real> &matrix, int rlev,
//!            int il, int iu, int jl, int ju, int kl, int ku, int color, bool th)
//! \brief Implementation of the Red-Black Gauss-Seidel Smoother
//!        rlev = relative level from the finest level of this Multigrid block

void MGCRDiffusion::Smooth(AthenaArray<Real> &u, const AthenaArray<Real> &src,
         const AthenaArray<Real> &coeff, const AthenaArray<Real> &matrix, int rlev,
         int il, int iu, int jl, int ju, int kl, int ku, int color, bool th) {
  Real dx;
  if (rlev <= 0) dx = rdx_*static_cast<Real>(1<<(-rlev));
  else           dx = rdx_/static_cast<Real>(1<<rlev);
  Real dx2 = SQR(dx);
  Real isix = omega_/6.0;
  color ^= pmy_driver_->coffset_;
  if (fsmoother_ == 1) { // jacobi-rb
    if (th == true && (ku-kl) >=  minth_) {
      AthenaArray<Real> &work = temp[0];
#pragma omp parallel num_threads(pmy_driver_->nthreads_)
      {
#pragma omp for
        for (int k=kl; k<=ku; k++) {
          for (int j=jl; j<=ju; j++) {
            int c = (color + k + j) & 1;
#pragma ivdep
            for (int i=il+c; i<=iu; i+=2) {
              Real M = matrix(CCM,k,j,i)*u(k,j,i-1)   + matrix(CCP,k,j,i)*u(k,j,i+1)
                     + matrix(CMC,k,j,i)*u(k,j-1,i)   + matrix(CPC,k,j,i)*u(k,j+1,i)
                     + matrix(MCC,k,j,i)*u(k-1,j,i)   + matrix(PCC,k,j,i)*u(k+1,j,i)
                     + matrix(CMM,k,j,i)*u(k,j-1,i-1) + matrix(CMP,k,j,i)*u(k,j-1,i+1)
                     + matrix(CPM,k,j,i)*u(k,j+1,i-1) + matrix(CPP,k,j,i)*u(k,j+1,i+1)
                     + matrix(MCM,k,j,i)*u(k-1,j,i-1) + matrix(MCP,k,j,i)*u(k-1,j,i+1)
                     + matrix(PCM,k,j,i)*u(k+1,j,i-1) + matrix(PCP,k,j,i)*u(k+1,j,i+1)
                     + matrix(MMC,k,j,i)*u(k-1,j-1,i) + matrix(MPC,k,j,i)*u(k-1,j+1,i)
                     + matrix(PMC,k,j,i)*u(k+1,j-1,i) + matrix(PPC,k,j,i)*u(k+1,j+1,i);
                work(k,j,i) = (src(k,j,i) - M) / matrix(CCC,k,j,i);
            }
          }
        }
#pragma omp for
        for (int k=kl; k<=ku; k++) {
          for (int j=jl; j<=ju; j++) {
            int c = (color + k + j) & 1;
#pragma ivdep
            for (int i=il+c; i<=iu; i+=2)
              u(k,j,i) += omega_ * (work(k,j,i) - u(k,j,i));
          }
        }
      }
    } else {
      int t = 0;
#ifdef OPENMP_PARALLEL
      t = omp_get_thread_num();
#endif
      AthenaArray<Real> &work = temp[t];
      for (int k=kl; k<=ku; k++) {
        for (int j=jl; j<=ju; j++) {
          int c = (color + k + j) & 1;
#pragma ivdep
          for (int i=il+c; i<=iu; i+=2) {
            Real M = matrix(CCM,k,j,i)*u(k,j,i-1)   + matrix(CCP,k,j,i)*u(k,j,i+1)
                   + matrix(CMC,k,j,i)*u(k,j-1,i)   + matrix(CPC,k,j,i)*u(k,j+1,i)
                   + matrix(MCC,k,j,i)*u(k-1,j,i)   + matrix(PCC,k,j,i)*u(k+1,j,i)
                   + matrix(CMM,k,j,i)*u(k,j-1,i-1) + matrix(CMP,k,j,i)*u(k,j-1,i+1)
                   + matrix(CPM,k,j,i)*u(k,j+1,i-1) + matrix(CPP,k,j,i)*u(k,j+1,i+1)
                   + matrix(MCM,k,j,i)*u(k-1,j,i-1) + matrix(MCP,k,j,i)*u(k-1,j,i+1)
                   + matrix(PCM,k,j,i)*u(k+1,j,i-1) + matrix(PCP,k,j,i)*u(k+1,j,i+1)
                   + matrix(MMC,k,j,i)*u(k-1,j-1,i) + matrix(MPC,k,j,i)*u(k-1,j+1,i)
                   + matrix(PMC,k,j,i)*u(k+1,j-1,i) + matrix(PPC,k,j,i)*u(k+1,j+1,i);
            work(k,j,i) = (src(k,j,i) - M) / matrix(CCC,k,j,i);
          }
        }
      }
      for (int k=kl; k<=ku; k++) {
        for (int j=jl; j<=ju; j++) {
          int c = (color + k + j) & 1;
#pragma ivdep
          for (int i=il+c; i<=iu; i+=2)
            u(k,j,i) += omega_ * (work(k,j,i) - u(k,j,i));
        }
      }
    }
  } else { // jacobi
    if (th == true && (ku-kl) >=  minth_) {
      AthenaArray<Real> &work = temp[0];
#pragma omp parallel num_threads(pmy_driver_->nthreads_)
      {
#pragma omp for
        for (int k=kl; k<=ku; k++) {
          for (int j=jl; j<=ju; j++) {
#pragma ivdep
            for (int i=il; i<=iu; i++) {
              Real M = matrix(CCM,k,j,i)*u(k,j,i-1)   + matrix(CCP,k,j,i)*u(k,j,i+1)
                     + matrix(CMC,k,j,i)*u(k,j-1,i)   + matrix(CPC,k,j,i)*u(k,j+1,i)
                     + matrix(MCC,k,j,i)*u(k-1,j,i)   + matrix(PCC,k,j,i)*u(k+1,j,i)
                     + matrix(CMM,k,j,i)*u(k,j-1,i-1) + matrix(CMP,k,j,i)*u(k,j-1,i+1)
                     + matrix(CPM,k,j,i)*u(k,j+1,i-1) + matrix(CPP,k,j,i)*u(k,j+1,i+1)
                     + matrix(MCM,k,j,i)*u(k-1,j,i-1) + matrix(MCP,k,j,i)*u(k-1,j,i+1)
                     + matrix(PCM,k,j,i)*u(k+1,j,i-1) + matrix(PCP,k,j,i)*u(k+1,j,i+1)
                     + matrix(MMC,k,j,i)*u(k-1,j-1,i) + matrix(MPC,k,j,i)*u(k-1,j+1,i)
                     + matrix(PMC,k,j,i)*u(k+1,j-1,i) + matrix(PPC,k,j,i)*u(k+1,j+1,i);
                work(k,j,i) = (src(k,j,i) - M) / matrix(CCC,k,j,i);
            }
          }
        }
#pragma omp for
        for (int k=kl; k<=ku; k++) {
          for (int j=jl; j<=ju; j++) {
#pragma ivdep
            for (int i=il; i<=iu; i++)
              u(k,j,i) += omega_ * (work(k,j,i) - u(k,j,i));
          }
        }
      }
    } else {
      int t = 0;
#ifdef OPENMP_PARALLEL
      t = omp_get_thread_num();
#endif
      AthenaArray<Real> &work = temp[t];
      for (int k=kl; k<=ku; k++) {
        for (int j=jl; j<=ju; j++) {
#pragma ivdep
          for (int i=il; i<=iu; i++) {
            Real M = matrix(CCM,k,j,i)*u(k,j,i-1)   + matrix(CCP,k,j,i)*u(k,j,i+1)
                   + matrix(CMC,k,j,i)*u(k,j-1,i)   + matrix(CPC,k,j,i)*u(k,j+1,i)
                   + matrix(MCC,k,j,i)*u(k-1,j,i)   + matrix(PCC,k,j,i)*u(k+1,j,i)
                   + matrix(CMM,k,j,i)*u(k,j-1,i-1) + matrix(CMP,k,j,i)*u(k,j-1,i+1)
                   + matrix(CPM,k,j,i)*u(k,j+1,i-1) + matrix(CPP,k,j,i)*u(k,j+1,i+1)
                   + matrix(MCM,k,j,i)*u(k-1,j,i-1) + matrix(MCP,k,j,i)*u(k-1,j,i+1)
                   + matrix(PCM,k,j,i)*u(k+1,j,i-1) + matrix(PCP,k,j,i)*u(k+1,j,i+1)
                   + matrix(MMC,k,j,i)*u(k-1,j-1,i) + matrix(MPC,k,j,i)*u(k-1,j+1,i)
                   + matrix(PMC,k,j,i)*u(k+1,j-1,i) + matrix(PPC,k,j,i)*u(k+1,j+1,i);
            work(k,j,i) = (src(k,j,i) - M) / matrix(CCC,k,j,i);
          }
        }
      }
      for (int k=kl; k<=ku; k++) {
        for (int j=jl; j<=ju; j++) {
#pragma ivdep
          for (int i=il; i<=iu; i++)
            u(k,j,i) += omega_ * (work(k,j,i) - u(k,j,i));
        }
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGCRDiffusion::CalculateDefect(AthenaArray<Real> &def,
//!            const AthenaArray<Real> &u, const AthenaArray<Real> &src,
//!            const AthenaArray<Real> &coeff, const AthenaArray<Real> &matrix,
//!            int rlev, int il, int iu, int jl, int ju, int kl, int ku, bool th)
//! \brief Implementation of the Defect calculation
//!        rlev = relative level from the finest level of this Multigrid block

void MGCRDiffusion::CalculateDefect(AthenaArray<Real> &def, const AthenaArray<Real> &u,
                    const AthenaArray<Real> &src, const AthenaArray<Real> &coeff,
                    const AthenaArray<Real> &matrix, int rlev, int il, int iu,
                    int jl, int ju, int kl, int ku, bool th) {
  Real dx;
  if (rlev <= 0) dx = rdx_*static_cast<Real>(1<<(-rlev));
  else           dx = rdx_/static_cast<Real>(1<<rlev);
  Real idx2 = 1.0/SQR(dx);

#pragma omp parallel for num_threads(pmy_driver_->nthreads_) if (th && (ku-kl) >= minth_)
  for (int k=kl; k<=ku; k++) {
    for (int j=jl; j<=ju; j++) {
#pragma omp simd
      for (int i=il; i<=iu; i++) {
        Real M = matrix(CCC,k,j,i)*u(k,j,i)
               + matrix(CCM,k,j,i)*u(k,j,i-1)   + matrix(CCP,k,j,i)*u(k,j,i+1)
               + matrix(CMC,k,j,i)*u(k,j-1,i)   + matrix(CPC,k,j,i)*u(k,j+1,i)
               + matrix(MCC,k,j,i)*u(k-1,j,i)   + matrix(PCC,k,j,i)*u(k+1,j,i)
               + matrix(CMM,k,j,i)*u(k,j-1,i-1) + matrix(CMP,k,j,i)*u(k,j-1,i+1)
               + matrix(CPM,k,j,i)*u(k,j+1,i-1) + matrix(CPP,k,j,i)*u(k,j+1,i+1)
               + matrix(MCM,k,j,i)*u(k-1,j,i-1) + matrix(MCP,k,j,i)*u(k-1,j,i+1)
               + matrix(PCM,k,j,i)*u(k+1,j,i-1) + matrix(PCP,k,j,i)*u(k+1,j,i+1)
               + matrix(MMC,k,j,i)*u(k-1,j-1,i) + matrix(MPC,k,j,i)*u(k-1,j+1,i)
               + matrix(PMC,k,j,i)*u(k+1,j-1,i) + matrix(PPC,k,j,i)*u(k+1,j+1,i);
        def(k,j,i) = src(k,j,i) - M;
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGCRDiffusion::CalculateFASRHS(AthenaArray<Real> &src,
//!            const AthenaArray<Real> &u, const AthenaArray<Real> &coeff,
//!            const AthenaArray<Real> &matrix, int rlev, int il, int iu, int jl, int ju,
//!            int kl, int ku, bool th)
//! \brief Implementation of the RHS calculation for FAS
//!        rlev = relative level from the finest level of this Multigrid block

void MGCRDiffusion::CalculateFASRHS(AthenaArray<Real> &src, const AthenaArray<Real> &u,
                    const AthenaArray<Real> &coeff, const AthenaArray<Real> &matrix,
                    int rlev, int il, int iu, int jl, int ju, int kl, int ku, bool th) {
  Real dx;
  if (rlev <= 0) dx = rdx_*static_cast<Real>(1<<(-rlev));
  else           dx = rdx_/static_cast<Real>(1<<rlev);
  Real idx2 = 1.0/SQR(dx);
#pragma omp parallel for num_threads(pmy_driver_->nthreads_) if (th && (ku-kl) >= minth_)
  for (int k=kl; k<=ku; k++) {
    for (int j=jl; j<=ju; j++) {
#pragma omp simd
      for (int i=il; i<=iu; i++) {
        Real M = matrix(CCC,k,j,i)*u(k,j,i)
               + matrix(CCM,k,j,i)*u(k,j,i-1)   + matrix(CCP,k,j,i)*u(k,j,i+1)
               + matrix(CMC,k,j,i)*u(k,j-1,i)   + matrix(CPC,k,j,i)*u(k,j+1,i)
               + matrix(MCC,k,j,i)*u(k-1,j,i)   + matrix(PCC,k,j,i)*u(k+1,j,i)
               + matrix(CMM,k,j,i)*u(k,j-1,i-1) + matrix(CMP,k,j,i)*u(k,j-1,i+1)
               + matrix(CPM,k,j,i)*u(k,j+1,i-1) + matrix(CPP,k,j,i)*u(k,j+1,i+1)
               + matrix(MCM,k,j,i)*u(k-1,j,i-1) + matrix(MCP,k,j,i)*u(k-1,j,i+1)
               + matrix(PCM,k,j,i)*u(k+1,j,i-1) + matrix(PCP,k,j,i)*u(k+1,j,i+1)
               + matrix(MMC,k,j,i)*u(k-1,j-1,i) + matrix(MPC,k,j,i)*u(k-1,j+1,i)
               + matrix(PMC,k,j,i)*u(k+1,j-1,i) + matrix(PPC,k,j,i)*u(k+1,j+1,i);
        src(k,j,i) += M;
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGCRDiffusion::CalculateMatrix(AthenaArray<Real> &matrix,
//!                         const AthenaArray<Real> &coeff, int rlev, Real dt,
//!                         int il, int iu, int jl, int ju, int kl, int ku, bool th)
//! \brief calculate Matrix element for cosmic ray transport
//!        rlev = relative level from the finest level of this Multigrid block

void MGCRDiffusion::CalculateMatrix(AthenaArray<Real> &matrix,
                             const AthenaArray<Real> &coeff, Real dt, int rlev,
                             int il, int iu, int jl, int ju, int kl, int ku, bool th) {
  Real dx;
  if (rlev <= 0) dx = rdx_*static_cast<Real>(1<<(-rlev));
  else           dx = rdx_/static_cast<Real>(1<<rlev);
  if (!(static_cast<MGCRDiffusionDriver*>(pmy_driver_)->fsteady_)) { // time-dependent
    Real fac = dt/SQR(dx), efac = 0.125*fac;
#pragma omp parallel for num_threads(pmy_driver_->nthreads_) if (th && (ku-kl) >= minth_)
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
#pragma omp simd
        for (int i=il; i<=iu; i++) {
          // center
          matrix(CCC,k,j,i) = 1.0 + dt*coeff(NLAMBDA,k,j,i) + 0.5*fac*(
                              2.0*coeff(DXX,k,j,i)+coeff(DXX,k,j,i+1)+coeff(DXX,k,j,i-1)
                            + 2.0*coeff(DYY,k,j,i)+coeff(DYY,k,j+1,i)+coeff(DYY,k,j-1,i)
                            + 2.0*coeff(DZZ,k,j,i)+coeff(DZZ,k+1,j,i)+coeff(DZZ,k-1,j,i));
          // face
          matrix(CCM,k,j,i) = fac*(-0.5*(coeff(DXX,k,j,i) + coeff(DXX,k,j,i-1))
                                +0.125*((coeff(DYX,k,j+1,i)-coeff(DYX,k,j-1,i))
                                       +(coeff(DZX,k+1,j,i)-coeff(DZX,k-1,j,i))));
          matrix(CCP,k,j,i) = fac*(-0.5*(coeff(DXX,k,j,i+1)+coeff(DXX,k,j,i))
                                -0.125*((coeff(DYX,k,j+1,i)-coeff(DYX,k,j-1,i))
                                       +(coeff(DZX,k+1,j,i)-coeff(DZX,k-1,j,i))));
          matrix(CMC,k,j,i) = fac*(-0.5*(coeff(DYY,k,j,i) + coeff(DYY,k,j-1,i))
                                +0.125*((coeff(DXY,k,j,i+1)-coeff(DXY,k,j,i-1))
                                       +(coeff(DZY,k+1,j,i)-coeff(DZY,k-1,j,i))));
          matrix(CPC,k,j,i) = fac*(-0.5*(coeff(DYY,k,j+1,i)+coeff(DYY,k,j,i))
                                -0.125*((coeff(DXY,k,j,i+1)-coeff(DXY,k,j,i-1))
                                       +(coeff(DZY,k+1,j,i)-coeff(DZY,k-1,j,i))));
          matrix(MCC,k,j,i) = fac*(-0.5*(coeff(DZZ,k,j,i) + coeff(DZZ,k-1,j,i))
                                +0.125*((coeff(DYZ,k,j+1,i)-coeff(DYZ,k,j-1,i))
                                       +(coeff(DXZ,k,j,i+1)-coeff(DXZ,k,j,i-1))));
          matrix(PCC,k,j,i) = fac*(-0.5*(coeff(DZZ,k+1,j,i)+coeff(DZZ,k,j,i))
                                -0.125*((coeff(DYZ,k,j+1,i)-coeff(DYZ,k,j-1,i))
                                       +(coeff(DXZ,k,j,i+1)-coeff(DXZ,k,j,i-1))));
          // edge
          matrix(CMM,k,j,i) = -efac*(coeff(DXY,k,j,i) + coeff(DXY,k,j,i-1)
                                    +coeff(DYX,k,j,i) + coeff(DYX,k,j-1,i));
          matrix(CMP,k,j,i) =  efac*(coeff(DXY,k,j,i+1)+coeff(DXY,k,j,i)
                                    +coeff(DYX,k,j,i) + coeff(DYX,k,j-1,i));
          matrix(CPM,k,j,i) =  efac*(coeff(DXY,k,j,i) + coeff(DXY,k,j,i-1)
                                    +coeff(DYX,k,j+1,i)+coeff(DYX,k,j,i));
          matrix(CPP,k,j,i) = -efac*(coeff(DXY,k,j,i+1)+coeff(DXY,k,j,i)
                                    +coeff(DYX,k,j+1,i)+coeff(DYX,k,j,i));
          matrix(MCM,k,j,i) = -efac*(coeff(DZX,k,j,i) + coeff(DZX,k-1,j,i)
                                    +coeff(DXZ,k,j,i) + coeff(DXZ,k,j,i-1));
          matrix(MCP,k,j,i) =  efac*(coeff(DZX,k,j,i) + coeff(DZX,k-1,j,i)
                                    +coeff(DXZ,k,j,i+1)+coeff(DXZ,k,j,i));
          matrix(PCM,k,j,i) =  efac*(coeff(DZX,k+1,j,i)+coeff(DZX,k,j,i)
                                    +coeff(DXZ,k,j,i) + coeff(DXZ,k,j,i-1));
          matrix(PCP,k,j,i) = -efac*(coeff(DZX,k+1,j,i)+coeff(DZX,k,j,i)
                                    +coeff(DXZ,k,j,i+1)+coeff(DXZ,k,j,i));
          matrix(MMC,k,j,i) = -efac*(coeff(DYZ,k,j,i) + coeff(DYZ,k,j-1,i)
                                    +coeff(DZY,k,j,i) + coeff(DZY,k-1,j,i));
          matrix(MPC,k,j,i) =  efac*(coeff(DYZ,k,j+1,i)+coeff(DYZ,k,j,i)
                                    +coeff(DZY,k,j,i) + coeff(DZY,k-1,j,i));
          matrix(PMC,k,j,i) =  efac*(coeff(DYZ,k,j,i) + coeff(DYZ,k,j-1,i)
                                    +coeff(DZY,k+1,j,i)+coeff(DZY,k,j,i));
          matrix(PPC,k,j,i) = -efac*(coeff(DYZ,k,j+1,i)+coeff(DYZ,k,j,i)
                                    +coeff(DZY,k+1,j,i)+coeff(DZY,k,j,i));
        }
      }
    }
  } else { // steady state
    Real fac = 1.0/SQR(dx), efac = 0.125*fac;
#pragma omp parallel for num_threads(pmy_driver_->nthreads_) if (th && (ku-kl) >= minth_)
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
#pragma omp simd
        for (int i=il; i<=iu; i++) {
          // center
          matrix(CCC,k,j,i) = coeff(NLAMBDA,k,j,i) + 0.5*fac*(
                              2.0*coeff(DXX,k,j,i)+coeff(DXX,k,j,i+1)+coeff(DXX,k,j,i-1)
                            + 2.0*coeff(DYY,k,j,i)+coeff(DYY,k,j+1,i)+coeff(DYY,k,j-1,i)
                            + 2.0*coeff(DZZ,k,j,i)+coeff(DZZ,k+1,j,i)+coeff(DZZ,k-1,j,i));
          // face
          matrix(CCM,k,j,i) = fac*(-0.5*(coeff(DXX,k,j,i) + coeff(DXX,k,j,i-1))
                                +0.125*((coeff(DYX,k,j+1,i)-coeff(DYX,k,j-1,i))
                                       +(coeff(DZX,k+1,j,i)-coeff(DZX,k-1,j,i))));
          matrix(CCP,k,j,i) = fac*(-0.5*(coeff(DXX,k,j,i+1)+coeff(DXX,k,j,i))
                                -0.125*((coeff(DYX,k,j+1,i)-coeff(DYX,k,j-1,i))
                                       +(coeff(DZX,k+1,j,i)-coeff(DZX,k-1,j,i))));
          matrix(CMC,k,j,i) = fac*(-0.5*(coeff(DYY,k,j,i) + coeff(DYY,k,j-1,i))
                                +0.125*((coeff(DXY,k,j,i+1)-coeff(DXY,k,j,i-1))
                                       +(coeff(DZY,k+1,j,i)-coeff(DZY,k-1,j,i))));
          matrix(CPC,k,j,i) = fac*(-0.5*(coeff(DYY,k,j+1,i)+coeff(DYY,k,j,i))
                                -0.125*((coeff(DXY,k,j,i+1)-coeff(DXY,k,j,i-1))
                                       +(coeff(DZY,k+1,j,i)-coeff(DZY,k-1,j,i))));
          matrix(MCC,k,j,i) = fac*(-0.5*(coeff(DZZ,k,j,i) + coeff(DZZ,k-1,j,i))
                                +0.125*((coeff(DYZ,k,j+1,i)-coeff(DYZ,k,j-1,i))
                                       +(coeff(DXZ,k,j,i+1)-coeff(DXZ,k,j,i-1))));
          matrix(PCC,k,j,i) = fac*(-0.5*(coeff(DZZ,k+1,j,i)+coeff(DZZ,k,j,i))
                                -0.125*((coeff(DYZ,k,j+1,i)-coeff(DYZ,k,j-1,i))
                                       +(coeff(DXZ,k,j,i+1)-coeff(DXZ,k,j,i-1))));
          // edge
          matrix(CMM,k,j,i) = -efac*(coeff(DXY,k,j,i) + coeff(DXY,k,j,i-1)
                                    +coeff(DYX,k,j,i) + coeff(DYX,k,j-1,i));
          matrix(CMP,k,j,i) =  efac*(coeff(DXY,k,j,i+1)+coeff(DXY,k,j,i)
                                    +coeff(DYX,k,j,i) + coeff(DYX,k,j-1,i));
          matrix(CPM,k,j,i) =  efac*(coeff(DXY,k,j,i) + coeff(DXY,k,j,i-1)
                                    +coeff(DYX,k,j+1,i)+coeff(DYX,k,j,i));
          matrix(CPP,k,j,i) = -efac*(coeff(DXY,k,j,i+1)+coeff(DXY,k,j,i)
                                    +coeff(DYX,k,j+1,i)+coeff(DYX,k,j,i));
          matrix(MCM,k,j,i) = -efac*(coeff(DZX,k,j,i) + coeff(DZX,k-1,j,i)
                                    +coeff(DXZ,k,j,i) + coeff(DXZ,k,j,i-1));
          matrix(MCP,k,j,i) =  efac*(coeff(DZX,k,j,i) + coeff(DZX,k-1,j,i)
                                    +coeff(DXZ,k,j,i+1)+coeff(DXZ,k,j,i));
          matrix(PCM,k,j,i) =  efac*(coeff(DZX,k+1,j,i)+coeff(DZX,k,j,i)
                                    +coeff(DXZ,k,j,i) + coeff(DXZ,k,j,i-1));
          matrix(PCP,k,j,i) = -efac*(coeff(DZX,k+1,j,i)+coeff(DZX,k,j,i)
                                    +coeff(DXZ,k,j,i+1)+coeff(DXZ,k,j,i));
          matrix(MMC,k,j,i) = -efac*(coeff(DYZ,k,j,i) + coeff(DYZ,k,j-1,i)
                                    +coeff(DZY,k,j,i) + coeff(DZY,k-1,j,i));
          matrix(MPC,k,j,i) =  efac*(coeff(DYZ,k,j+1,i)+coeff(DYZ,k,j,i)
                                    +coeff(DZY,k,j,i) + coeff(DZY,k-1,j,i));
          matrix(PMC,k,j,i) =  efac*(coeff(DYZ,k,j,i) + coeff(DYZ,k,j-1,i)
                                    +coeff(DZY,k+1,j,i)+coeff(DZY,k,j,i));
          matrix(PPC,k,j,i) = -efac*(coeff(DYZ,k,j+1,i)+coeff(DYZ,k,j,i)
                                    +coeff(DZY,k+1,j,i)+coeff(DZY,k,j,i));
        }
      }
    }
  }
  return;
}
