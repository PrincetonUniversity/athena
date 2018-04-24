#ifndef MULTIGRID_MULTIGRID_HPP_
#define MULTIGRID_MULTIGRID_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multigrid.hpp
//  \brief

// C/C++ headers
#include <iostream>
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../bvals/bvals_mg.hpp"
#include "../task_list/mg_task_list.hpp"

class Mesh;
class MeshBlock;
class ParameterInput;
class Coordinates;

void MGPeriodicInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGPeriodicOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);

//! \class Multigrid
//  \brief gravitational block

class Multigrid {
public:
  Multigrid(MultigridDriver *pmd, LogicalLocation iloc, int igid, int ilid,
     int invar, int nghost, RegionSize isize, MGBoundaryFunc_t *MGBoundary,
     enum BoundaryFlag *input_bcs, bool root);
  virtual ~Multigrid();

  MGBoundaryValues *pmgbval;
  enum BoundaryType btype, btypef;
  Multigrid *next, *prev;

  void LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh);
  void LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac);
  void RestrictFMGSource(void);
  void RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh);
  void ZeroClearData(void);
  void Restrict(void);
  void ProlongateAndCorrect(void);
  void FMGProlongate(void);
  void SetFromRootGrid(AthenaArray<Real> &src, int ci, int cj, int ck);
  Real CalculateDefectNorm(int n, int nrm);
  Real CalculateTotal(int type, int n);
  void SubtractAverage(int type, int n, Real ave);

  // small functions
  int GetCurrentNumberOfCells(void) { return 1<<current_level_; }
  int GetNumberOfLevels(void) { return nlevel_; }
  int GetCurrentLevel(void) { return current_level_; }
  AthenaArray<Real>& GetCurrentData(void) { return u_[current_level_]; }
  AthenaArray<Real>& GetCurrentSource(void) { return src_[current_level_]; }
  Real GetRootSource(int n) { return src_[0](n,ngh_,ngh_,ngh_); }

  // pure virtual functions
  virtual void Smooth(int color) = 0;
  virtual void CalculateDefect(void) = 0;

  friend class MultigridDriver;
  friend class MultigridTaskList;
  friend class MGBoundaryValues;
  friend class MGGravityDriver;

protected:
  int gid_, lid_;
  LogicalLocation loc_;
  MultigridDriver *pmy_driver_;
  RegionSize size_;
  int nlevel_, ngh_, nvar_, current_level_;
  Real rdx_, rdy_, rdz_;
  AthenaArray<Real> *u_, *def_, *src_;
private:
  bool root_flag_;
  TaskState ts_;
};


//! \class MultigridDriver
//  \brief Multigrid driver

class MultigridDriver {
public:
  MultigridDriver(Mesh *pm, MGBoundaryFunc_t *MGBoundary, int invar);
  virtual ~MultigridDriver();
  void AddMultigrid(Multigrid *nmg);
  void SubtractAverage(int type);
  void SetupMultigrid(void);
  void FillRootGridSource(void);
  void FMGProlongate(void);
  void TransferFromRootToBlocks(void);
  void OneStepToFiner(int nsmooth);
  void OneStepToCoarser(int nsmooth);
  void SolveVCycle(int npresmooth, int npostsmooth);
  void SolveFCycle(int npresmooth, int npostsmooth);
  void SolveFMGCycle(void);
  void SolveIterative(void);

  virtual void SolveCoarsestGrid(void);
  Real CalculateDefectNorm(int n, int nrm);
  Multigrid* FindMultigrid(int tgid);

  // small functions
  int GetNumMultigrids(void) { return nblist_[Globals::my_rank]; }

  // pure virtual functions
  virtual void Solve(int step) = 0;

  friend class Multigrid;
  friend class MultigridTaskList;
  friend class MGBoundaryValues;

protected:
  int nranks_, nvar_, nrootlevel_, nmblevel_, ntotallevel_, mode_;
  int current_level_;
  int *nslist_, *nblist_, *nvlist_, *nvslist_, *ranklist_;
  MGBoundaryFunc_t MGBoundaryFunction_[6];
  Mesh *pmy_mesh_;
  Multigrid *pmg_;
  Multigrid *mgroot_;
  bool fperiodic_;
  Real last_ave_;
  Real eps_;
private:
  MultigridTaskList *mgtlist_;
  Real *rootbuf_;
  AthenaArray<Real> rootsrc_;
#ifdef MPI_PARALLEL
  MPI_Comm MPI_COMM_MULTIGRID;
#endif
};


#endif // MULTIGRID_MULTIGRID_HPP_
