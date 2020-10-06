#ifndef BVALS_CC_MG_BVALS_MG_HPP_
#define BVALS_CC_MG_BVALS_MG_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_mg.hpp
//! \brief defines MGBoundaryValues class

// C headers

// C++ headers
#include <string>   // string

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../bvals.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// forward declarations
class Multigrid;
class MultigridDriver;
class MeshBlock;
class MeshBlockTree;
class ParameterInput;
class Coordinates;

//----------------------------------------------------------------------------------------
//! \class MGBoundaryValues
//! \brief BVals data and functions

class MGBoundaryValues : public BoundaryBase {
 public:
  MGBoundaryValues(Multigrid *pmg, BoundaryFlag *input_bcs);
  ~MGBoundaryValues();

  void InitBoundaryData(BoundaryQuantity type);
  void DestroyBoundaryData();

  void ApplyPhysicalBoundaries(int flag = 0);
  void StartReceivingMultigrid(BoundaryQuantity type, bool folddata);
  void ClearBoundaryMultigrid(BoundaryQuantity type);
  bool SendMultigridBoundaryBuffers(BoundaryQuantity type, bool folddata);
  bool ReceiveMultigridBoundaryBuffers(BoundaryQuantity type, bool folddata);
  void ProlongateMultigridBoundaries(bool folddata);
  void ProlongateMultigridBoundariesFluxCons();
  void CopyNeighborInfoFromMeshBlock();

 protected:
  Multigrid *pmy_mg_;
  MGBoundaryFunc MGBoundaryFunction_[6];
  BoundaryData<> bdata_;
  AthenaArray<Real> cbuf_, cbufold_;

#ifdef MPI_PARALLEL
  MPI_Comm mgcomm_;
#endif

  int LoadMultigridBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb,
                                           bool folddata);
  int LoadMultigridBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb,
                                           bool folddata);
  int LoadMultigridBoundaryBufferToFiner(Real *buf, const NeighborBlock& nb,
                                         bool folddata);
  void SetMultigridBoundarySameLevel(const Real *buf, const NeighborBlock& nb,
                                     bool folddata);
  void SetMultigridBoundaryFromCoarser(const Real *buf, const NeighborBlock& nb,
                                       bool folddata);
  void SetMultigridBoundaryFromFiner(const Real *buf, const NeighborBlock& nb,
                                     bool folddata);
  //!@{
  //! functions specific to physics
  virtual int LoadMultigridBoundaryBufferToCoarserFluxCons(Real *buf,
                                                           const NeighborBlock& nb) = 0;
  virtual int LoadMultigridBoundaryBufferToFinerFluxCons(Real *buf,
                                                         const NeighborBlock& nb) = 0;
  virtual void SetMultigridBoundaryFromCoarserFluxCons(const Real *buf,
                                                       const NeighborBlock& nb) = 0;
  virtual void SetMultigridBoundaryFromFinerFluxCons(const Real *buf,
                                                     const NeighborBlock& nb) = 0;
  //!@}
  friend class Multigrid;
  friend class MultigridDriver;
};


//----------------------------------------------------------------------------------------
//! \class MGGravityBoundaryValues
//! \brief BVals data and functions for Multigrid Gravity

class MGGravityBoundaryValues : public MGBoundaryValues {
 public:
  MGGravityBoundaryValues(Multigrid *pmg, BoundaryFlag *input_bcs)
    : MGBoundaryValues(pmg, input_bcs) {}
 private:
  int LoadMultigridBoundaryBufferToCoarserFluxCons(Real *buf,
                                                   const NeighborBlock& nb) final;
  int LoadMultigridBoundaryBufferToFinerFluxCons(Real *buf,
                                                 const NeighborBlock& nb) final;
  void SetMultigridBoundaryFromCoarserFluxCons(const Real *buf,
                                               const NeighborBlock& nb) final;
  void SetMultigridBoundaryFromFinerFluxCons(const Real *buf,
                                             const NeighborBlock& nb) final;

  friend class Multigrid;
  friend class MultigridDriver;
};

#endif // BVALS_CC_MG_BVALS_MG_HPP_
