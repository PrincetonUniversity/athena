#ifndef MULTIGRID_MULTIGRID_HPP_
#define MULTIGRID_MULTIGRID_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multigrid.hpp
//  \brief

// C headers

// C++ headers
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/cc/mg/bvals_mg.hpp"
#include "../globals.hpp"
#include "../mesh/mesh.hpp"
#include "../task_list/mg_task_list.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

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
  Multigrid(MultigridDriver *pmd, MeshBlock *pmb, int invar, int nghost);
  virtual ~Multigrid();

  MGBoundaryValues *pmgbval;
  // KGF: Both btype=BoundaryQuantity::mggrav and btypef=BoundaryQuantity::mggrav_f (face
  // neighbors only) are passed to comm function calls in mg_task_list.cpp Only
  // BoundaryQuantity::mggrav is handled in a case in InitBoundaryData(). Passed directly
  // (not through btype) in MGBoundaryValues() ctor
  BoundaryQuantity btype, btypef;

  void LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh);
  void LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac);
  void RestrictFMGSource();
  void RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh);
  void ZeroClearData();
  void Restrict();
  void ProlongateAndCorrect();
  void FMGProlongate();
  void SetFromRootGrid(AthenaArray<Real> &src, int ci, int cj, int ck);
  Real CalculateDefectNorm(int n, int nrm);
  Real CalculateTotal(int type, int n);
  void SubtractAverage(int type, int n, Real ave);

  // small functions
  int GetCurrentNumberOfCells() { return 1<<current_level_; }
  int GetNumberOfLevels() { return nlevel_; }
  int GetCurrentLevel() { return current_level_; }
  AthenaArray<Real>& GetCurrentData() { return u_[current_level_]; }
  AthenaArray<Real>& GetCurrentSource() { return src_[current_level_]; }
  Real GetRootSource(int n) { return src_[0](n,ngh_,ngh_,ngh_); }

  // pure virtual functions
  virtual void Smooth(int color) = 0;
  virtual void CalculateDefect() = 0;

  friend class MultigridDriver;
  friend class MultigridTaskList;
  friend class MGBoundaryValues;
  friend class MGGravityDriver;

 protected:
  MultigridDriver *pmy_driver_;
  MeshBlock *pmy_block_;
  LogicalLocation loc_;
  RegionSize size_;
  int gid_, nlevel_, ngh_, nvar_, current_level_;
  Real rdx_, rdy_, rdz_;
  AthenaArray<Real> *u_, *def_, *src_;

 private:
  TaskStates ts_;
};


//! \class MultigridDriver
//  \brief Multigrid driver

class MultigridDriver {
 public:
  MultigridDriver(Mesh *pm, MGBoundaryFunc *MGBoundary, int invar);
  virtual ~MultigridDriver();
  void SubtractAverage(int type);
  void SetupMultigrid();
  void FillRootGridSource();
  void FMGProlongate();
  void TransferFromRootToBlocks();
  void OneStepToFiner(int nsmooth);
  void OneStepToCoarser(int nsmooth);
  void SolveVCycle(int npresmooth, int npostsmooth);
  void SolveFCycle(int npresmooth, int npostsmooth);
  void SolveFMGCycle();
  void SolveIterative();

  virtual void SolveCoarsestGrid();
  Real CalculateDefectNorm(int n, int nrm);
  Multigrid* FindMultigrid(int tgid);

  // small functions
  int GetNumMultigrids() { return nblist_[Globals::my_rank]; }

  // pure virtual functions
  virtual void Solve(int step) = 0;

  friend class Multigrid;
  friend class MultigridTaskList;
  friend class MGBoundaryValues;

 protected:
  int nranks_, nvar_, nrootlevel_, nmblevel_, ntotallevel_, mode_;
  int current_level_;
  int *nslist_, *nblist_, *nvlist_, *nvslist_, *ranklist_;
  int nrbx1_, nrbx2_, nrbx3_;
  MGBoundaryFunc MGBoundaryFunction_[6];
  Mesh *pmy_mesh_;
  std::vector<Multigrid*> vmg_;
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
  int mg_phys_id_;
#endif
};


#endif // MULTIGRID_MULTIGRID_HPP_
