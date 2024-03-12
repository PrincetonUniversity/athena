#ifndef CRDIFFUSION_MG_CRDIFFUSION_HPP_
#define CRDIFFUSION_MG_CRDIFFUSION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file mg_crdiffusion.hpp
//! \brief defines MGCRDiffusion class

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../multigrid/multigrid.hpp"
#include "crdiffusion.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;
class Multigrid;
class CRDiffusionBoundaryTaskList;


//! \class MGCRDiffusion
//! \brief Multigrid CR Diffusion solver for each block

class MGCRDiffusion : public Multigrid {
 public:
  MGCRDiffusion(MGCRDiffusionDriver *pmd, MeshBlock *pmb);
  ~MGCRDiffusion();

  void AddCRSource(const AthenaArray<Real> &src, int ngh, Real dt);

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
  void CalculateMatrix(AthenaArray<Real> &matrix, const AthenaArray<Real> &coeff,
       Real dt, int rlev, int il, int iu, int jl, int ju, int kl, int ku, bool th) final;

  friend class MGCRDiffusionDriver;

 private:
  Real omega_;
  int fsmoother_;
};


//! \class MGCRDiffusionDriver
//! \brief Multigrid CR Diffusion solver

class MGCRDiffusionDriver : public MultigridDriver {
 public:
  MGCRDiffusionDriver(Mesh *pm, ParameterInput *pin);
  ~MGCRDiffusionDriver();
  void Solve(int stage, Real dt) final;

  friend class MGCRDiffusion;

 private:
  CRDiffusionBoundaryTaskList *crtlist_;
  Real omega_;
  int fsmoother_;
  bool fsteady_;
};

#endif // CRDIFFUSION_MG_CRDIFFUSION_HPP_
