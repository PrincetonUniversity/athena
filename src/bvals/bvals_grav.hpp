#ifndef BVALS_BVALS_GRAV_HPP_
#define BVALS_BVALS_GRAV_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_grav.hpp
//  \brief defines GravityBoundaryValues class

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
class FFTBlock;
class MeshBlock;
class MeshBlockTree;
class ParameterInput;
class Coordinates;

//! \struct GravityBoundaryData
//  \brief structure storing multigrid boundary information
typedef struct GravityBoundaryData {
  int nbmax;
  enum BoundaryStatus flag[56], sflag[56];
  Real *send[56], *recv[56];
#ifdef MPI_PARALLEL
  MPI_Request req_send[56], req_recv[56];
#endif
} GravityBoundaryData;

//----------------------------------------------------------------------------------------
//! \class GravityBoundaryValues
//  \brief BVals data and functions

class GravityBoundaryValues : public BoundaryBase {
public:
  GravityBoundaryValues(MeshBlock *pmb, enum BoundaryFlag *input_bcs);
  ~GravityBoundaryValues();

  void InitBoundaryData(GravityBoundaryData &bd);
  void DestroyBoundaryData(GravityBoundaryData &bd);

  void ApplyPhysicalBoundaries(void);
  void StartReceivingGravity(void);
  void ClearBoundaryGravity(void);
  int LoadGravityBoundaryBufferSameLevel(AthenaArray<Real> &src, Real *buf,
                                       const NeighborBlock& nb);
  bool SendGravityBoundaryBuffers(AthenaArray<Real> &src);
  void SetGravityBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                 const NeighborBlock& nb);
  bool ReceiveGravityBoundaryBuffers(AthenaArray<Real> &dst);

private:
  MeshBlock *pmy_block_;
  GravityBoundaryFunc_t GravityBoundaryFunction_[6];
  GravityBoundaryData bd_gravity_;
};

#endif // BVALS_BVALS_GRAV_HPP_
