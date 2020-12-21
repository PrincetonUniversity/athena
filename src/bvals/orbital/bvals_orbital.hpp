#ifndef BVALS_ORBITAL_BVALS_ORBITAL_HPP_
#define BVALS_ORBITAL_BVALS_ORBITAL_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_orbital.hpp
//! \brief handle orbital communication for orbital advection

// C headers

// C++ headers

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../bvals.hpp"
#include "../bvals_interfaces.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// forward declarations
class OrbitalAdvection;
using OrbitalBoundaryData = BoundaryData<4>;

//----------------------------------------------------------------------------------------
//! \class OrbitalBoundaryCommunication
//! \brief

class OrbitalBoundaryCommunication {
 public:
  explicit OrbitalBoundaryCommunication(OrbitalAdvection *porb);
  ~OrbitalBoundaryCommunication();
  // BoundaryCommunication:
  void SetupPersistentMPI();
  void ComputeOrbit(const Real dt);
  void StartReceiving(BoundaryCommSubset phase);
  void ClearBoundary(BoundaryCommSubset phase);

  // BoundaryBuffer public functions with shared implementations
  void SendBoundaryBuffersCC();
  void SendBoundaryBuffersFC();
  bool ReceiveBoundaryBuffersCC();
  bool ReceiveBoundaryBuffersFC();

  static constexpr int max_phys_id = 2;

 private:
  void InitBoundaryData(OrbitalBoundaryData &bd, BoundaryQuantity type);
  void DestroyBoundaryData(OrbitalBoundaryData &bd);

  void LoadHydroBufferSameLevel(Real *buf, int &p, const int nb);
  void LoadHydroBufferToCoarser(Real *buf, int &p, const int nb);
  void LoadHydroBufferToFiner(Real *buf, int &p, const int nb);
  void SetHydroBufferSameLevel(Real *buf, int &p, const int nb);
  void SetHydroBufferFromCoarser(Real *buf, int &p, const int nb);
  void SetHydroBufferFromFiner(Real *buf, int &p, const int nb);

  void LoadFieldBufferSameLevel(Real *buf, int &p, const int nb);
  void LoadFieldBufferToCoarser(Real *buf, int &p, const int nb);
  void LoadFieldBufferToFiner(Real *buf, int &p, const int nb);
  void SetFieldBufferSameLevel(Real *buf, int &p, const int nb);
  void SetFieldBufferFromCoarser(Real *buf, int &p, const int nb);
  void SetFieldBufferFromFiner(Real *buf, int &p, const int nb);

  void LoadScalarBufferSameLevel(Real *buf, int &p, const int nb);
  void LoadScalarBufferToCoarser(Real *buf, int &p, const int nb);
  void LoadScalarBufferToFiner(Real *buf, int &p, const int nb);
  void SetScalarBufferSameLevel(Real *buf, int &p, const int nb);
  void SetScalarBufferFromCoarser(Real *buf, int &p, const int nb);
  void SetScalarBufferFromFiner(Real *buf, int &p, const int nb);

  // BoundaryData objects for communications
  OrbitalBoundaryData orbital_bd_cc_[2], orbital_bd_fc_[2];
  SimpleNeighborBlock orbital_send_neighbor_[2][4], orbital_recv_neighbor_[2][4];
  int orbital_send_cc_count_[2][4], orbital_recv_cc_count_[2][4];
  int orbital_send_fc_count_[2][4], orbital_recv_fc_count_[2][4];
  int xgh;

  int *size_cc_send[2];  //same, coarser, fine*4
  int *size_cc_recv[2];  //same, coarser, fine*4
  int *size_fc_send[2];  //same, coarser, fine*4
  int *size_fc_recv[2];  //same, coarser, fine*4

  // ptr to MeshBlock, Mesh, BoundaryValues, OrbitalAdvection
  MeshBlock *pmy_block_;
  Mesh *pmy_mesh_;
  BoundaryValues *pbval_;
  OrbitalAdvection *pmy_orbital_;

#ifdef MPI_PARALLEL
  int orbital_advection_cc_phys_id_, orbital_advection_fc_phys_id_;
#endif
};

#endif // BVALS_ORBITAL_BVALS_ORBITAL_HPP_
