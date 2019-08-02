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

  if (mode_ <= 1)
    SolveFMGCycle();
  else
    SolveIterative();

  // Return the result
  for (auto itr = vmg_.begin(); itr < vmg_.end(); itr++) {
    Multigrid *pmg = *itr;
    Gravity *pgrav = pmg->pmy_block_->pgrav;
    pmg->RetrieveResult(pgrav->phi, 0, NGHOST);
    pgrav->grav_mean_rho = mean_rho;
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void MGGravity::Smooth(AthenaArray<Real> &u, const AthenaArray<Real> &src,
//            int rlev, int il, int iu, int jl, int ju, int kl, int ku, int color)
//  \brief Implementation of the Red-Black Gauss-Seidel Smoother
//         rlev = relative level from the finest level of this Multigrid block

void MGGravity::Smooth(AthenaArray<Real> &u, const AthenaArray<Real> &src, int rlev,
                       int il, int iu, int jl, int ju, int kl, int ku, int color) {
  int c = color;
  Real dx;
  if (rlev <= 0) dx = rdx_*static_cast<Real>(1<<(-rlev));
  else           dx = rdx_/static_cast<Real>(1<<rlev);
  Real dx2 = SQR(dx);
  Real isix = omega_/6.0;
  for (int k=kl; k<=ku; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=il+c; i<=iu; i+=2)
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
//! \fn  void MGGravity::CalculateDefect(AthenaArray<Real> &def, 
//                       const AthenaArray<Real> &u, const AthenaArray<Real> &src,
//                       int rlev, int il, int iu, int jl, int ju, int kl, int ku)
//  \brief Implementation of the Defect calculation
//         rlev = relative level from the finest level of this Multigrid block

void MGGravity::CalculateDefect(AthenaArray<Real> &def, const AthenaArray<Real> &u,
                                const AthenaArray<Real> &src, int rlev,
                                int il, int iu, int jl, int ju, int kl, int ku) {
  Real dx;
  if (rlev <= 0) dx = rdx_*static_cast<Real>(1<<(-rlev));
  else           dx = rdx_/static_cast<Real>(1<<rlev);
  Real idx2 = 1.0/SQR(dx);
  for (int k=kl; k<=ku; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=il; i<=iu; i++)
        def(0,k,j,i) = (6.0*u(0,k,j,i) - u(0,k+1,j,i) - u(0,k,j+1,i) - u(0,k,j,i+1)
                       - u(0,k-1,j,i) - u(0,k,j-1,i) - u(0,k,j,i-1))*idx2
                       + src(0,k,j,i);
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn  void MGGravity::CalculateFASRHS(AthenaArray<Real> &src,
//   const AthenaArray<Real> &u, int rlev, int il, int iu, int jl, int ju, int kl, int ku)
//  \brief Implementation of the RHS calculation for FAS
//         rlev = relative level from the finest level of this Multigrid block

void MGGravity::CalculateFASRHS(AthenaArray<Real> &src, const AthenaArray<Real> &u,
                         int rlev, int il, int iu, int jl, int ju, int kl, int ku) {
  Real dx;
  if (rlev <= 0) dx = rdx_*static_cast<Real>(1<<(-rlev));
  else           dx = rdx_/static_cast<Real>(1<<rlev);
  Real idx2 = 1.0/SQR(dx);
  for (int k=kl; k<=ku; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=il; i<=iu; i++)
        src(0,k,j,i) -= (6.0*u(0,k,j,i) - u(0,k+1,j,i) - u(0,k,j+1,i) - u(0,k,j,i+1)
                        - u(0,k-1,j,i) - u(0,k,j-1,i) - u(0,k,j,i-1))*idx2;
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void MGGravityDriver::SetOctetBoundaryFromCoarser(AthenaArray<Real> &dst,
//           const AthenaArray<Real> &un, LogicalLocation loc, int ox1, int ox2, int ox3)
//  \brief set an Octet boundary from a neighbor Octet on the coarser level

void MGGravityDriver::SetOctetBoundaryFromCoarser(AthenaArray<Real> &dst,
     const AthenaArray<Real> &un, LogicalLocation loc, int ox1, int ox2, int ox3) {
  constexpr Real itw = 1.0/12.0;
  const int ngh = mgroot_->ngh_;
  const int i = ngh, j = ngh, k = ngh;
  const AthenaArray<Real> &u = dst;
  int ci, cj, ck;
  if (loc.level == nrootlevel_ - 1) { // from root
    ci = loc.lx1;
    cj = loc.lx2;
    ck = loc.lx3;
  } else { // from a neighbor octet
    if (ox1 == 0)     ci = (static_cast<int>(loc.lx1) & 1) + ngh;
    else if (ox1 < 0) ci = ngh+1;
    else              ci = ngh;
    if (ox2 == 0)     cj = (static_cast<int>(loc.lx2) & 1) + ngh;
    else if (ox2 < 0) cj = ngh+1;
    else              cj = ngh; 
    if (ox3 == 0)     ck = (static_cast<int>(loc.lx3) & 1) + ngh;
    else if (ox3 < 0) ck = ngh+1;
    else              ck = ngh; 
  }

  if (ox1 != 0) {
    int ig, fi;
    if (ox1 < 0) ig = 0,       fi = ngh;
    else         ig = ngh + 2, fi = ngh + 1;
    Real cg = itw * (8.0*un(0, ck, cj, ci) + ((u(0, k,   j, fi) + u(0, k,   j+1, fi))
                                           +  (u(0, k+1, j, fi) + u(0, k+1, j+1, fi))));
    Real qdy = 0.25 * ((u(0, k+1, j+1, fi) - u(0, k+1, j,   fi))
                     + (u(0, k,   j+1, fi) - u(0, k,   j,   fi)));
    Real qdz = 0.25 * ((u(0, k+1, j+1, fi) - u(0, k,   j+1, fi))
                     + (u(0, k+1, j,   fi) - u(0, k,   j,   fi)));
    dst(0, k,   j,   ig) = cg - qdy - qdz;
    dst(0, k,   j+1, ig) = cg + qdy - qdz;
    dst(0, k+1, j,   ig) = cg - qdy + qdz;
    dst(0, k+1, j+1, ig) = cg + qdy + qdz;
  } else if (ox2 != 0) {
    int jg, fj;
    if (ox2 < 0) jg = 0,       fj = ngh;
    else         jg = ngh + 2, fj = ngh + 1;
    Real cg = itw * (8.0 * un(0, ck, cj, ci) + ((u(0, k,   fj, i) + u(0, k,   fj, i+1))
                                             +  (u(0, k+1, fj, i) + u(0, k+1, fj, i+1))));
    Real qdx = 0.25 * ((u(0, k+1, fj, i+1) - u(0, k+1, fj, i))
                     + (u(0, k,   fj, i+1) - u(0, k,   fj, i)));
    Real qdz = 0.25 * ((u(0, k+1, fj, i+1) - u(0, k,   fj, i+1))
                     + (u(0, k+1, fj, i  ) - u(0, k,   fj, i)));
    dst(0, k,   jg, i  ) = cg - qdx - qdz;
    dst(0, k,   jg, i+1) = cg + qdx - qdz;
    dst(0, k+1, jg, i  ) = cg - qdx + qdz;
    dst(0, k+1, jg, i+1) = cg + qdx + qdz;
  } else if (ox3 != 0) {
    int kg, fk;
    if (ox3 < 0) kg = 0,       fk = ngh;
    else         kg = ngh + 2, fk = ngh + 1;
    Real cg = itw * (8.0 * un(0, ck, cj, ci) + ((u(0, fk, j,   i) + u(0, fk, j,   i+1))
                                             +  (u(0, fk, j+1, i) + u(0, fk, j+1, i+1))));
    Real qdx = 0.25 * ((u(0, fk, j+1, i+1) - u(0, fk, j+1, i))
                     + (u(0, fk, j,   i+1) - u(0, fk, j,   i)));
    Real qdy = 0.25 * ((u(0, fk, j+1, i+1) - u(0, fk, j,   i+1))
                     + (u(0, fk, j+1, i  ) - u(0, fk, j,   i)));
    dst(0, kg, j,   i  ) = cg - qdx - qdy;
    dst(0, kg, j,   i+1) = cg + qdx - qdy;
    dst(0, kg, j+1, i  ) = cg - qdx + qdy;
    dst(0, kg, j+1, i+1) = cg + qdx + qdy;
  }

  return;
}
