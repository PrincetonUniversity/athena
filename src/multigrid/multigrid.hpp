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
#include "../mesh/mesh.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;

typedef void (*MGBoundaryFunc_t)(AthenaArray<Real> &dst,Real time, Real dt,
             int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
             Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);

void MGPeriodicInnerX1(AthenaArray<Real> &dst,Real time, Real dt,
                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                    Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicOuterX1(AthenaArray<Real> &dst,Real time, Real dt,
                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                    Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicInnerX2(AthenaArray<Real> &dst,Real time, Real dt,
                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                    Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicOuterX2(AthenaArray<Real> &dst,Real time, Real dt,
                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                    Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicInnerX3(AthenaArray<Real> &dst,Real time, Real dt,
                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                    Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicOuterX3(AthenaArray<Real> &dst,Real time, Real dt,
                    int nvar, int is, int ie, int js, int je, int ks, int ke, int ngh,
                    Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);

//! \class Multigrid
//  \brief gravitational block

class Multigrid {
public:
  Multigrid(MeshBlock *pmb, int invar, int nx, int ny, int nz,
            RegionSize isize, MGBoundaryFunc_t *MGBoundary);
  virtual ~Multigrid();

  void LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh);
  void LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac);
  void RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh);
  void ZeroClearData(int lev);
  void ApplyPhysicalBoundaries(int lev);
  void Restrict(int clev);
  void ProlongateAndCorrect(int clev);
  void FMGProlongate(int clev);
  Real CalculateDefectNorm(int n, int nrm);
  Real CalculateTotalSource(int n);
  void SubtractMeanSource(int n, Real ave);

  // pure virtual functions
  virtual void Smooth(int lev, int color) = 0;
  virtual void CalculateDefect(int lev) = 0;

  friend class MultigridDriver;

protected:
  MeshBlock *pmy_block_;
  RegionSize size_;
  int nlev_, nx_, ny_, nz_, ngh_;
  int current_level_;
  Real rdx_, rdy_, rdz_;
  AthenaArray<Real> *u_, *def_, *src_;
  MGBoundaryFunc_t MGBoundaryFunction_[6];
};


//! \class MultigridDriver
//  \brief Multigrid driver

class MultigridDriver
{
public:
  MultigridDriver(Mesh *pm, MEshBlock *pmb, MGBoundaryFunc_t *MGBoundary,
                  ParameterInput *pin);
  virtual ~MultigridDriver();
  virtual Multigrid* GetMultigridBlock(MeshBlock *) = 0;
  int GetNumMeshBlocks(void) { return nblocks;};
  friend class Multigrid;

protected:
  Mesh *pmy_mesh_;
  MeshBlock *pblock_;
  int nblocks_;
  TaskState ts_;
  MGBoundaryFunc_t MGBoundaryFunction_[6];
};


#endif // MULTIGRID_HPP
