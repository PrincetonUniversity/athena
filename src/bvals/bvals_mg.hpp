#ifndef BVALS_BVALS_MG_HPP_
#define BVALS_BVALS_MG_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_mg.hpp
//  \brief defines MGBoundaryValues class

// C++ headers
#include <string>   // string

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "bvals.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// forward declarations
class Multigrid;
class MeshBlock;
class MeshBlockTree;
class ParameterInput;
class Coordinates;

//! \struct MGBoundaryData
//  \brief structure storing multigrid boundary information
typedef struct MGBoundaryData {
  int nbmax;
  enum BoundaryStatus flag[56], sflag[56];
  Real *send[56], *recv[56];
#ifdef MPI_PARALLEL
  MPI_Request req_send[56], req_recv[56];
#endif
} MGBoundaryData;

//----------------------------------------------------------------------------------------
//! \class MGBoundaryValues
//  \brief BVals data and functions

class MGBoundaryValues : public BoundaryBase {
public:
  MGBoundaryValues(Multigrid *pmg, enum BoundaryFlag *input_bcs,
                   MGBoundaryFunc_t *MGBoundary);
  ~MGBoundaryValues();

  void InitBoundaryData(MGBoundaryData &bd, enum BoundaryType type);
  void DestroyBoundaryData(MGBoundaryData &bd);

  void ApplyPhysicalBoundaries(void);
  void StartReceivingMultigrid(int nc, enum BoundaryType type);
  void ClearBoundaryMultigrid(enum BoundaryType type);
  int LoadMultigridBoundaryBufferSameLevel(AthenaArray<Real> &src,
                   int nvar, int nc, int ngh, Real *buf, const NeighborBlock& nb);
  bool SendMultigridBoundaryBuffers(AthenaArray<Real> &src,
                                    int nc, enum BoundaryType type);
  void SetMultigridBoundarySameLevel(AthenaArray<Real> &dst,
                   int nvar, int nc, int ngh, Real *buf, const NeighborBlock& nb);
  bool ReceiveMultigridBoundaryBuffers(AthenaArray<Real> &dst,
                                       int nc, enum BoundaryType type);

private:
  Multigrid *pmy_mg_;
  MGBoundaryFunc_t MGBoundaryFunction_[6];
  MGBoundaryData bd_mggrav_;

#ifdef MPI_PARALLEL
  MPI_Comm mgcomm_;
#endif

  friend class Multigrid;
};

#endif // BVALS_BVALS_MG_HPP_
