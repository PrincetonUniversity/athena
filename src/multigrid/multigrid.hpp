#ifndef MULTIGRID_MULTIGRID_HPP_
#define MULTIGRID_MULTIGRID_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multigrid.hpp
//  \brief defines the Multigrid base class

// C headers

// C++ headers
#include <iostream>
#include <vector>

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


enum class MGVariable {src, u};
enum class MGNormType {max, l1, l2};

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
  void SetFromRootGrid(AthenaArray<Real> &src, int ck, int cj, int ci);
  Real CalculateDefectNorm(MGNormType nrm, int n);
  Real CalculateTotal(MGVariable type, int n);
  void SubtractAverage(MGVariable type, int n, Real ave);
  void StoreOldData();
  Real GetCoarsestData(MGVariable type, int n);
  void SetData(MGVariable type, int n, int k, int j, int i, Real v);

  // small functions
  int GetCurrentNumberOfCells() { return 1<<current_level_; }
  int GetNumberOfLevels() { return nlevel_; }
  int GetCurrentLevel() { return current_level_; }
  AthenaArray<Real>& GetCurrentData() { return u_[current_level_]; }
  AthenaArray<Real>& GetCurrentSource() { return src_[current_level_]; }

  // physics-dependent virtual functions
  virtual void Smooth(int color) = 0;
  virtual void CalculateDefect() = 0;
  virtual void CalculateFASRHS() {};

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
  AthenaArray<Real> *u_, *def_, *src_, *uold_;

  void RestrictOperation(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
                         int nx1, int nx2, int nx3);

 private:
  TaskStates ts_;
};


//! \class MultigridDriver
//  \brief Multigrid driver

class MultigridDriver {
 public:
  MultigridDriver(Mesh *pm, MGBoundaryFunc *MGBoundary, int invar);
  virtual ~MultigridDriver();
  void SubtractAverage(MGVariable type);
  void SetupMultigrid();
  void TransferFromBlocksToRoot(bool initflag = false);
  void FMGProlongate();
  void TransferFromRootToBlocks();
  void OneStepToFiner(int nsmooth);
  void OneStepToCoarser(int nsmooth);
  void SolveVCycle(int npresmooth, int npostsmooth);
  void SolveFCycle(int npresmooth, int npostsmooth);
  void SolveFMGCycle();
  void SolveIterative();

  virtual void SolveCoarsestGrid();
  Real CalculateDefectNorm(MGNormType nrm, int n);
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
  int current_level_, fmglevel_;
  int *nslist_, *nblist_, *nvlist_, *nvslist_, *nvlisti_, *nvslisti_, *ranklist_;
  int nrbx1_, nrbx2_, nrbx3_;
  MGBoundaryFunc MGBoundaryFunction_[6];
  Mesh *pmy_mesh_;
  std::vector<Multigrid*> vmg_;
  Multigrid *mgroot_;
  bool fsubtract_average_, ffas_;
  Real last_ave_;
  Real eps_;

 private:
  MultigridTaskList *mgtlist_;
  Real *rootbuf_;
#ifdef MPI_PARALLEL
  MPI_Comm MPI_COMM_MULTIGRID;
  int mg_phys_id_;
#endif
};

// Function Prototypes for Multigrid Boundary

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

void MGZeroGradientInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroGradientOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroGradientInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroGradientOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroGradientInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroGradientOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);

void MGZeroFixedInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroFixedOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroFixedInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroFixedOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroFixedInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
void MGZeroFixedOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);

#endif // MULTIGRID_MULTIGRID_HPP_
