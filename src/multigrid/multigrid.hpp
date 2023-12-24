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
#include <cstdint>  // std::int64_t
#include <cstdio> // std::size_t
#include <iostream>
#include <unordered_map>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals_interfaces.hpp"
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

enum class MGVariable {src, u, coeff};
enum class MGNormType {max, l1, l2};

constexpr int minth_ = 8;

//! \fn inline std::int64_t rotl(std::int64_t i, int s)
//  \brief left bit rotation function for 64bit integers (unsafe if s > 64)

inline std::int64_t rotl(std::int64_t i, int s) {
  return (i << s) | (i >> (64 - s));
}


//! \struct LogicalLocationHash
//  \brief Hash function object for LogicalLocation

struct LogicalLocationHash {
 public:
  std::size_t operator()(const LogicalLocation &l) const {
    return static_cast<std::size_t>(l.lx1^rotl(l.lx2,21)^rotl(l.lx3,42));
  }
};


//! \class MGCoordinates
//  \brief Minimum set of coordinate arrays for Multigrid

class MGCoordinates {
 public:
  void AllocateMGCoordinates(int nx, int ny, int nz);
  void CalculateMGCoordinates(const RegionSize &size, int ll, int ngh);

  AthenaArray<Real> x1f, x2f, x3f, x1v, x2v, x3v;
};


//! \class MGOctet
//  \brief structure containing eight (+ ghost) cells for Mesh Refinement

class MGOctet {
 public:
  void Allocate(int nvar, int ncoct, int nccoct, int ncoeff, int nmatrix);
  LogicalLocation loc;
  bool fleaf;
  AthenaArray<Real> u, def, src, uold, coeff, matrix;
  MGCoordinates coord, ccoord;
};


//! \class Multigrid
//  \brief Multigrid object containing each MeshBlock and/or the root block

class Multigrid {
 public:
  Multigrid(MultigridDriver *pmd, MeshBlock *pmb, int nghost);
  virtual ~Multigrid();

  MGBoundaryValues *pmgbval;
  // KGF: Both btype=BoundaryQuantity::mg and btypef=BoundaryQuantity::mg_faceonly (face
  // neighbors only) are passed to comm function calls in mg_task_list.cpp Only
  // BoundaryQuantity::mg is handled in a case in InitBoundaryData(). Passed directly
  // (not through btype) in MGBoundaryValues() ctor
  BoundaryQuantity btype, btypef;

  void LoadFinestData(const AthenaArray<Real> &src, int ns, int ngh);
  void LoadSource(const AthenaArray<Real> &src, int ns, int ngh, Real fac);
  void LoadCoefficients(const AthenaArray<Real> &coeff, int ngh);
  void ApplyMask();
  void RestrictFMGSource();
  void RetrieveResult(AthenaArray<Real> &dst, int ns, int ngh);
  void RetrieveDefect(AthenaArray<Real> &dst, int ns, int ngh);
  void ZeroClearData();
  void RestrictBlock();
  void RestrictCoefficients();
  void ProlongateAndCorrectBlock();
  void FMGProlongateBlock();
  void SmoothBlock(int color);
  void CalculateDefectBlock();
  void CalculateFASRHSBlock();
  void CalculateMatrixBlock(Real dt);
  void SetFromRootGrid(bool folddata);
  Real CalculateDefectNorm(MGNormType nrm, int n);
  Real CalculateTotal(MGVariable type, int n);
  void SubtractAverage(MGVariable type, int n, Real ave);
  void StoreOldData();
  Real GetCoarsestData(MGVariable type, int n);
  void SetData(MGVariable type, int n, int k, int j, int i, Real v);

  void CalculateMultipoleCoefficients(AthenaArray<Real> &mpcoeff);
  void CalculateCenterOfMass(AthenaArray<Real> &mpcoeff);

  // small functions
  int GetCurrentNumberOfCells() { return 1<<current_level_; }
  int GetNumberOfLevels() { return nlevel_; }
  int GetCurrentLevel() { return current_level_; }
  AthenaArray<Real>& GetCurrentData() { return u_[current_level_]; }
  AthenaArray<Real>& GetCurrentSource() { return src_[current_level_]; }
  AthenaArray<Real>& GetCurrentOldData() { return uold_[current_level_]; }
  AthenaArray<Real>& GetCurrentCoefficient() { return coeff_[current_level_]; }

  // actual implementations of Multigrid operations
  void Restrict(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
                int nvar, int il, int iu, int jl, int ju, int kl, int ku, bool th);
  void ProlongateAndCorrect(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
    int il, int iu, int jl, int ju, int kl, int ku, int fil, int fjl, int fkl, bool th);
  void FMGProlongate(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
    int il, int iu, int jl, int ju, int kl, int ku, int fil, int fjl, int fkl, bool th);

  // physics-dependent virtual functions
  virtual void Smooth(AthenaArray<Real> &dst, const AthenaArray<Real> &src,
                      const AthenaArray<Real> &coeff, const AthenaArray<Real> &matrx,
                      int rlev, int il, int iu, int jl, int ju, int kl, int ku,
                      int color, bool th) = 0;
  virtual void CalculateDefect(AthenaArray<Real> &def, const AthenaArray<Real> &u,
               const AthenaArray<Real> &src, const AthenaArray<Real> &coeff,
               const AthenaArray<Real> &matrix, int rlev, int il, int iu, int jl, int ju,
               int kl, int ku, bool th) = 0;
  virtual void CalculateFASRHS(AthenaArray<Real> &def, const AthenaArray<Real> &src,
                 const AthenaArray<Real> &coeff, const AthenaArray<Real> &matrix,
                 int rlev, int il, int iu, int jl, int ju, int kl, int ku, bool th) = 0;
  virtual void CalculateMatrix(AthenaArray<Real> &matrix, const AthenaArray<Real> &coeff,
          Real dt, int rlev, int il, int iu, int jl, int ju, int kl, int ku, bool th) {}

  friend class MultigridDriver;
  friend class MultigridTaskList;
  friend class MGBoundaryValues;
  friend class MGGravityBoundaryValues;
  friend class MGGravityDriver;
  friend class MGCRDiffusionDriver;

 protected:
  MultigridDriver *pmy_driver_;
  MeshBlock *pmy_block_;
  LogicalLocation loc_;
  RegionSize size_;
  BoundaryFlag mg_block_bcs_[6];
  int nlevel_, ngh_, nvar_, ncoeff_, nmatrix_, current_level_;
  Real rdx_, rdy_, rdz_;
  Real defscale_;
  AthenaArray<Real> *u_, *def_, *src_, *uold_, *coeff_, *matrix_;
  MGCoordinates *coord_, *ccoord_;

 private:
  TaskStates ts_;
};


//! \class MultigridDriver
//  \brief Multigrid driver

class MultigridDriver {
 public:
  MultigridDriver(Mesh *pm, MGBoundaryFunc *MGBoundary, MGBoundaryFunc *MGCoeffBoundary,
                  MGMaskFunc MGSourceMask, MGMaskFunc MGCoeffMask, int invar, int ncoeff,
                  int nmatrix);
  virtual ~MultigridDriver();

  // pure virtual function
  virtual void Solve(int step, Real dt = 0.0) = 0;

  friend class Multigrid;
  friend class MultigridTaskList;
  friend class MGGravity;
  friend class MGBoundaryValues;
  friend class MGGravityBoundaryValues;
  friend class MGCRDiffusion;

 protected:
  void CheckBoundaryFunctions();
  void SubtractAverage(MGVariable type);
  void SetupMultigrid(Real dt, bool ftrivial);
  void SetupCoefficients();
  void TransferFromBlocksToRoot(bool initflag);
  void TransferFromRootToBlocks(bool folddata);
  void TransferCoefficientFromBlocksToRoot();
  void OneStepToFiner(int nsmooth);
  void OneStepToCoarser(int nsmooth);
  void FMGProlongate();
  void SolveVCycle(int npresmooth, int npostsmooth);
  void SolveFMGCycle();
  void SolveIterative();
  void SolveIterativeFixedTimes();

  virtual void SolveCoarsestGrid();
  Real CalculateDefectNorm(MGNormType nrm, int n);
  void CalculateMatrix(Real dt);
  Multigrid* FindMultigrid(int tgid);

  // octet manipulation functions
  void CalculateOctetCoordinates();
  void RestrictFMGSourceOctets();
  void RestrictOctets();
  void ZeroClearOctets();
  void StoreOldDataOctets();
  void CalculateFASRHSOctets();
  void SmoothOctets(int color);
  void ProlongateAndCorrectOctets();
  void FMGProlongateOctets();
  void SetBoundariesOctets(bool fprolong, bool folddata, bool fcoeff);
  void SetOctetBoundarySameLevel(AthenaArray<Real> &dst, const AthenaArray<Real> &un,
                       AthenaArray<Real> &uold, const AthenaArray<Real> &unold,
                       AthenaArray<Real> &cbuf, AthenaArray<Real> &cbufold,
                       int nvar, int ox1, int ox2, int ox3, bool folddata);
  void SetOctetBoundaryFromCoarser(const AthenaArray<Real> &un,
                       const AthenaArray<Real> &unold, AthenaArray<Real> &cbuf,
                       AthenaArray<Real> &cbufold, int nvar, const LogicalLocation &loc,
                       int ox1, int ox2, int ox3, bool folddata);
  void ApplyPhysicalBoundariesOctet(AthenaArray<Real> &u, const LogicalLocation &loc,
                                    const MGCoordinates &coord, bool fcbuf, bool fcoeff);
  void ProlongateOctetBoundaries(AthenaArray<Real> &u, AthenaArray<Real> &uold,
                      AthenaArray<Real> &cbuf, AthenaArray<Real> &cbufold,
                      int nvar, const AthenaArray<bool> &ncoarse, bool folddata);
  void SetOctetBoundariesBeforeTransfer(bool folddata);
  void RestrictOctetsBeforeTransfer();
  virtual void ProlongateOctetBoundariesFluxCons(AthenaArray<Real> &dst,
               AthenaArray<Real> &cbuf, const AthenaArray<bool> &ncoarse);

  // for multipole expansion
  void AllocateMultipoleCoefficients();
  void CalculateMultipoleCoefficients();
  virtual void ScaleMultipoleCoefficients();
  void CalculateCenterOfMass();

  // small functions
  int GetNumMultigrids() { return nblist_[Globals::my_rank]; }

  int nranks_, nthreads_, nbtotal_, nvar_, ncoeff_, nmatrix_, mode_;
  int locrootlevel_, nrootlevel_, nmblevel_, ntotallevel_, nreflevel_, maxreflevel_;
  int current_level_, fmglevel_;
  int *nslist_, *nblist_, *nvlist_, *nvslist_, *nvlisti_, *nvslisti_,
                          *nclist_, *ncslist_, *ranklist_;
  int nrbx1_, nrbx2_, nrbx3_;
  BoundaryFlag mg_mesh_bcs_[6];
  MGBoundaryFunc MGBoundaryFunction_[6];
  MGBoundaryFunc MGCoeffBoundaryFunction_[6];
  MGMaskFunc srcmask_, coeffmask_;
  Mesh *pmy_mesh_;
  std::vector<Multigrid*> vmg_;
  Multigrid *mgroot_;
  bool fsubtract_average_, ffas_, redblack_, needinit_, fshowdef_;
  Real last_ave_;
  Real eps_;
  int niter_, npresmooth_, npostsmooth_;
  int os_, oe_;
  int coffset_;
  int fprolongation_;

  // for mesh refinement
  std::vector<MGOctet> *octets_;
  std::unordered_map<LogicalLocation, int, LogicalLocationHash> *octetmap_;
  std::vector<bool> *octetbflag_;
  int *noctets_, *pmaxnoct_;
  AthenaArray<Real> *cbuf_, *cbufold_;
  AthenaArray<bool> *ncoarse_;

  MultigridTaskList *mgtlist_;

  // for multipole expansion
  AthenaArray<Real> *mpcoeff_;
  int mporder_, nmpcoeff_;
  AthenaArray<Real> mpo_;
  bool autompo_, nodipole_;

 private:
  Real *rootbuf_;
  int nb_rank_;
#ifdef MPI_PARALLEL
  MPI_Comm MPI_COMM_MULTIGRID;
  int mg_phys_id_;
#endif
};


inline Real RestrictOne(const AthenaArray<Real> &src, int v, int fi, int fj, int fk) {
  return 0.125*(src(v, fk,   fj,   fi)+src(v, fk,   fj,   fi+1)
               +src(v, fk,   fj+1, fi)+src(v, fk,   fj+1, fi+1)
               +src(v, fk+1, fj,   fi)+src(v, fk+1, fj,   fi+1)
               +src(v, fk+1, fj+1, fi)+src(v, fk+1, fj+1, fi+1));
}


// Multigrid Boundary functions

void MGPeriodicInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       const MGCoordinates &coord);
void MGPeriodicOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       const MGCoordinates &coord);
void MGPeriodicInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       const MGCoordinates &coord);
void MGPeriodicOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       const MGCoordinates &coord);
void MGPeriodicInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       const MGCoordinates &coord);
void MGPeriodicOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                       int is, int ie, int js, int je, int ks, int ke, int ngh,
                       const MGCoordinates &coord);

void MGZeroGradientInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           const MGCoordinates &coord);
void MGZeroGradientOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           const MGCoordinates &coord);
void MGZeroGradientInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           const MGCoordinates &coord);
void MGZeroGradientOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           const MGCoordinates &coord);
void MGZeroGradientInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           const MGCoordinates &coord);
void MGZeroGradientOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                           int is, int ie, int js, int je, int ks, int ke, int ngh,
                           const MGCoordinates &coord);

void MGZeroFixedInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord);
void MGZeroFixedOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord);
void MGZeroFixedInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord);
void MGZeroFixedOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord);
void MGZeroFixedInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord);
void MGZeroFixedOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
                        int is, int ie, int js, int je, int ks, int ke, int ngh,
                        const MGCoordinates &coord);

void MGMultipoleInnerX1(AthenaArray<Real> &dst, Real time, int nvar,
       int is, int ie, int js, int je, int ks, int ke, int ngh,
       const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
       const AthenaArray<Real> &mporigin, int mporder);
void MGMultipoleOuterX1(AthenaArray<Real> &dst, Real time, int nvar,
       int is, int ie, int js, int je, int ks, int ke, int ngh,
       const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
       const AthenaArray<Real> &mporigin, int mporder);
void MGMultipoleInnerX2(AthenaArray<Real> &dst, Real time, int nvar,
       int is, int ie, int js, int je, int ks, int ke, int ngh,
       const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
       const AthenaArray<Real> &mporigin, int mporder);
void MGMultipoleOuterX2(AthenaArray<Real> &dst, Real time, int nvar,
       int is, int ie, int js, int je, int ks, int ke, int ngh,
       const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
       const AthenaArray<Real> &mporigin, int mporder);
void MGMultipoleInnerX3(AthenaArray<Real> &dst, Real time, int nvar,
       int is, int ie, int js, int je, int ks, int ke, int ngh,
       const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
       const AthenaArray<Real> &mporigin, int mporder);
void MGMultipoleOuterX3(AthenaArray<Real> &dst, Real time, int nvar,
       int is, int ie, int js, int je, int ks, int ke, int ngh,
       const MGCoordinates &coord, const AthenaArray<Real> &mpcoeff,
       const AthenaArray<Real> &mporigin, int mporder);

#endif // MULTIGRID_MULTIGRID_HPP_
