#ifndef GRAVITY_MG_GRAVITY_HPP_
#define GRAVITY_MG_GRAVITY_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_gravity.hpp
//! \brief defines MGGravity class

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../multigrid/multigrid.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;
class Multigrid;
class GravityBoundaryTaskList;

//! \class MGGravity
//! \brief Multigrid gravity solver for each block

class MGGravity : public Multigrid {
 public:
  MGGravity(MultigridDriver *pmd, MeshBlock *pmb);
  ~MGGravity();

  void Smooth(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
              const AthenaArray<Real> &coeff, const AthenaArray<Real> &matrix, int rlev,
              int il, int iu, int jl, int ju, int kl, int ku, int color, bool th) final;
  void CalculateDefect(AthenaArray<Real> &def, const AthenaArray<Real> &u,
                const AthenaArray<Real> &src, const AthenaArray<Real> &coeff,
                const AthenaArray<Real> &matrix, int rlev, int il, int iu, int jl, int ju,
                int kl, int ku, bool th) final;
  void CalculateFASRHS(AthenaArray<Real> &def, const AthenaArray<Real> &src,
                const AthenaArray<Real> &coeff, const AthenaArray<Real> &matrix,
                int rlev, int il, int iu, int jl, int ju, int kl, int ku, bool th) final;
};


//! \class MGGravityDriver
//! \brief Multigrid gravity solver

class MGGravityDriver : public MultigridDriver {
 public:
  MGGravityDriver(Mesh *pm, ParameterInput *pin);
  ~MGGravityDriver();
  void Solve(int stage, Real dt = 0.0) final;
  void ProlongateOctetBoundariesFluxCons(AthenaArray<Real> &dst,
                 AthenaArray<Real> &cbuf, const AthenaArray<bool> &ncoarse) final;
  friend class MGGravity;
 private:
  Real four_pi_G_, omega_;
  GravityBoundaryTaskList *gtlist_;
};

#endif // GRAVITY_MG_GRAVITY_HPP_
