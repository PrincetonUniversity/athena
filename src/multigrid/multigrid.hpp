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

//! \class Multigrid
//  \brief gravitational block

class Multigrid {
public:
  Multigrid(int invar, int nx, int ny, int nz, int ngh, Real dx);
  virtual ~Multigrid();

  MeshBlock *pmy_block;

  void LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh);
  void LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac);
  void RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh);
  void ZeroClearData(int lev);
  void ApplyPhysicalBoundary(int lev);
  void Restrict(int clev);
  void ProlongateAndCorrect(int clev);
  void FMGProlongate(int clev);
  Real CalculateDefectNorm(int n, int nrm);
  Real CalculateTotalSource(int n);
  void SubtractMeanSource(int n, Real ave);

  // pure virtual functions
  virtual void Smooth(int lev, int color) = 0;
  virtual void CalculateDefect(int lev) = 0;

protected:
  int nlev_, nx_, ny_, nz_, int ngh_;
  Real dx_;
  AthenaArray<Real> *u_, *def_, *src_;
};
#endif // MULTIGRID_HPP
