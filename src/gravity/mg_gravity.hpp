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

//! \class MGGravity
//! \brief Multigrid gravity solver for each block

class MGGravity : public Multigrid {
 public:
  MGGravity(MultigridDriver *pmd, MeshBlock *pmb);
  ~MGGravity();

  void Smooth(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
              int rlev, int il, int iu, int jl, int ju, int kl, int ku, int color) final;
  void CalculateDefect(AthenaArray<Real> &def, const AthenaArray<Real> &u,
                       const AthenaArray<Real> &src, int rlev,
                       int il, int iu, int jl, int ju, int kl, int ku) final;
  void CalculateFASRHS(AthenaArray<Real> &def, const AthenaArray<Real> &src,
                       int rlev, int il, int iu, int jl, int ju, int kl, int ku) final;

 private:
  static constexpr Real omega_ = 1.15;
};


//! \class MGGravityDriver
//! \brief Multigrid gravity solver

class MGGravityDriver : public MultigridDriver {
 public:
  MGGravityDriver(Mesh *pm, ParameterInput *pin);
  ~MGGravityDriver();
  void Solve(int stage) final;
  // void SolveCoarsestGrid() final;
  void ProlongateOctetBoundariesFluxCons(AthenaArray<Real> &dst) final;
 private:
  Real four_pi_G_;
};

#endif // GRAVITY_MG_GRAVITY_HPP_
