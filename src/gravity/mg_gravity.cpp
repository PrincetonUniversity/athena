//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_gravity.cpp
//  \brief create multigrid solver for gravity

// C headers

// C++ headers
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
#include "gravity.hpp"
#include "mg_gravity.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

class MeshBlock;

//----------------------------------------------------------------------------------------
//! \fn MGGravityDriver::MGGravityDriver(Mesh *pm, ParameterInput *pin)
//  \brief MGGravityDriver constructor

MGGravityDriver::MGGravityDriver(Mesh *pm, ParameterInput *pin)
    : MultigridDriver(pm, pm->MGGravityBoundaryFunction_, 1) {
  four_pi_G_ = pmy_mesh_->four_pi_G_;
  eps_ = pmy_mesh_->grav_eps_;
  if (four_pi_G_==0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MGGravityDriver::MGGravityDriver" << std::endl
        << "Gravitational constant must be set in the Mesh::InitUserMeshData "
        << "using the SetGravitationalConstant or SetFourPiG function." << std::endl;
    ATHENA_ERROR(msg);
  }
  if (mode_>=2 && eps_<0.0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in MGGravityDriver::MGGravityDriver" << std::endl
        << "Convergence threshold must be set in the Mesh::InitUserMeshData "
        << "using the SetGravitatyThreshold for the iterative mode." << std::endl
        << "Set the threshold = 0.0 for automatic convergence control." << std::endl;
    ATHENA_ERROR(msg);
  }

  // Allocate the root multigrid
  mgroot_ = new MGGravity(this, nullptr);
}


//----------------------------------------------------------------------------------------
//! \fn MGGravityDriver::~MGGravityDriver()
//  \brief MGGravityDriver destructor
MGGravityDriver::~MGGravityDriver() {
  delete mgroot_;
}

//----------------------------------------------------------------------------------------
//! \fn void MGGravityDriver::Solve(int stage)
//  \brief load the data and solve

void MGGravityDriver::Solve(int stage) {
  // Construct the Multigrid array
  vmg_.clear();
  MeshBlock *pmb = pmy_mesh_->pblock;
  while (pmb != nullptr) {
    vmg_.push_back(pmb->pmg);
    pmb = pmb->next;
  }

  // load the source
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    // assume all the data are located on the same node
    pmg->LoadSource(pmg->pmy_block_->phydro->u, IDN, NGHOST, four_pi_G_);
    if (mode_ >= 2) // iterative mode - load initial guess
      pmg->LoadFinestData(pmg->pmy_block_->pgrav->phi, 0, NGHOST);
  }

  SetupMultigrid();
  Real mean_rho = 0.0;
  if (fsubtract_average_)
    mean_rho = last_ave_/four_pi_G_;

  if (mode_ <= 1) {
    SolveFMGCycle();
  } else {
    SolveIterative();
  }

  // Return the result
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    Gravity *pgrav = pmg->pmy_block_->pgrav;
    pmg->RetrieveResult(pgrav->phi,0,NGHOST);
    pgrav->grav_mean_rho = mean_rho;
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void MGGravity::Smooth(int color)
//  \brief Red-Black Gauss-Seidel Smoother
void MGGravity::Smooth(int color) {
  int c = color;
  AthenaArray<Real> &u = u_[current_level_];
  AthenaArray<Real> &src = src_[current_level_];
  int ll = nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is = js = ks = ngh_;
  ie = is+(size_.nx1>>ll)-1, je = js+(size_.nx2>>ll)-1, ke = ks+(size_.nx3>>ll)-1;
  Real dx = rdx_*static_cast<Real>(1<<ll);
  Real dx2 = SQR(dx);
  Real isix = omega_/6.0;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is+c; i<=ie; i+=2)
        u(0,k,j,i) -= ((6.0*u(0,k,j,i) - u(0,k+1,j,i) - u(0,k,j+1,i) - u(0,k,j,i+1)
                      - u(0,k-1,j,i) - u(0,k,j-1,i) - u(0,k,j,i-1))
                       + src(0,k,j,i)*dx2)*isix;
      c ^= 1;  // bitwise XOR assignment
    }
    c ^= 1;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void MGGravity::CalculateDefect()
//  \brief calculate the residual

void MGGravity::CalculateDefect() {
  AthenaArray<Real> &u = u_[current_level_];
  AthenaArray<Real> &src = src_[current_level_];
  AthenaArray<Real> &def = def_[current_level_];
  int ll = nlevel_-1-current_level_;
  int is, ie, js, je, ks, ke;
  is = js = ks = ngh_;
  ie = is+(size_.nx1>>ll)-1, je = js+(size_.nx2>>ll)-1, ke = ks+(size_.nx3>>ll)-1;
  Real dx = rdx_*static_cast<Real>(1<<ll);
  Real idx2 = 1.0/SQR(dx);
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        def(0,k,j,i) = (6.0*u(0,k,j,i) - u(0,k+1,j,i) - u(0,k,j+1,i) - u(0,k,j,i+1)
                        - u(0,k-1,j,i) - u(0,k,j-1,i) - u(0,k,j,i-1))*idx2
                       + src(0,k,j,i);
      }
    }
  }
  return;
}
