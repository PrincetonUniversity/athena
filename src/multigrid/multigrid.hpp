#ifndef MULTIGRID_HPP
#define MULTIGRID_HPP

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity.hpp
//  \brief defines Gravity class which implements data and functions for gravitational potential

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;


//! \class MultigridDriver
//  \brief Driver class

class MultigridDriver {
public:
  MultigridDriver(Mesh *pm, ParameterInput *pin, int nvar);
  ~MultigridDriver();

  void SetupMultigrid(ParameterInput *pin);
  void Solve(const AthenaArray<Real> &u);

private:
  MeshBlock *pblock_;

  friend class Multigrid;
};


//! \class Multigrid
//  \brief gravitational block

class Multigrid {
public:
  Multigrid(MeshBlock *pmb, int invar);
  virtual ~Multigrid();

  MeshBlock *pmy_block;
  MultigridDriver *pmgd;

  void LoadFinestData(AthenaArray<Real> &src, int ns);
  void LoadFinestSource(AthenaArray<Real> &src, int ns, Real fac);
  void RetrieveResult(AthenaArray<real> &dst, int ns);
  void ZeroClearData(int lev);
  void Restrict(int clev);
  void ProlongateAndCorrect(int clev);
  void FMGProlongate(int clev);
  Real CalculateDefectNorm(int nrm);

  // pure virtual functions
  virtual void Smooth(int lev, int color) = 0;
  virtual void CalculateDefect(int lev) = 0;

protected:
  // the multigrid solver
  const int ngh_=1;
  Real dx_;
  AthenaArray<Real> *u_, *def_, *src_;
};
#endif // MULTIGRID_HPP
