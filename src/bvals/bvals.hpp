#ifndef BOUNDARY_VALUES_HPP
#define BOUNDARY_VALUES_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file bvals.hpp
//  \brief defines BoundaryValues class used for setting BCs on all data types
//======================================================================================

// C++ headers
#include <string>   // string

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

// MPI headers
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// forward declarations
class MeshBlock;
class Hydro;
class Field;
class ParameterInput;
class Coordinates;
struct FaceField;
struct NeighborBlock;
struct PolarNeighborBlock;

// identifiers for all 6 faces of a MeshBlock on which boundary conditions are applied
enum BoundaryFace {FACE_UNDEF=-1, INNER_X1=0, OUTER_X1=1, INNER_X2=2, OUTER_X2=3, 
  INNER_X3=4, OUTER_X3=5};

// identifiers for boundary conditions
enum BoundaryFlag {BLOCK_BNDRY=-1, BNDRY_UNDEF=0, REFLECTING_BNDRY=1, OUTFLOW_BNDRY=2,
  USER_BNDRY=3, PERIODIC_BNDRY=4, POLAR_BNDRY=5};

//-------------------- prototypes for all BC functions ---------------------------------
void ReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void ReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void ReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);

void OutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void OutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
                    FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);

// function to return boundary flag given input string
enum BoundaryFlag GetBoundaryFlag(std::string input_string);


//! \class BoundaryValues
//  \brief BVals data and functions

class BoundaryValues {
public:
  BoundaryValues(MeshBlock *pmb, ParameterInput *pin);
  ~BoundaryValues();

  void Initialize(void);
  void StartReceivingForInit(void);
  void StartReceivingAll(void);

  void CheckBoundary(void);

  int LoadHydroBoundaryBufferSameLevel(AthenaArray<Real> &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadHydroBoundaryBufferToCoarser(AthenaArray<Real> &src, Real *buf,
                                       const NeighborBlock& nb, bool conserved_values);
  int LoadHydroBoundaryBufferToFiner(AthenaArray<Real> &src, Real *buf,
                                     const NeighborBlock& nb);
  void SendHydroBoundaryBuffers(AthenaArray<Real> &src, int step,
                                bool conserved_values = true);
  void SetHydroBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                 const NeighborBlock& nb);
  void SetHydroBoundaryFromCoarser(Real *buf, const NeighborBlock& nb,
                                   bool conserved_values);
  void SetHydroBoundaryFromFiner(AthenaArray<Real> &dst, Real *buf,
                                 const NeighborBlock& nb);
  bool ReceiveHydroBoundaryBuffers(AthenaArray<Real> &dst, int step);
  void ReceiveHydroBoundaryBuffersWithWait(AthenaArray<Real> &dst, int step);
  void RestrictHydro(AthenaArray<Real> &src,
                     int si, int ei, int sj, int ej, int sk, int ek);

  void SendFluxCorrection(int step);
  bool ReceiveFluxCorrection(int step);

  int LoadFieldBoundaryBufferSameLevel(FaceField &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadFieldBoundaryBufferToCoarser(FaceField &src, Real *buf,
                                       const NeighborBlock& nb);
  int LoadFieldBoundaryBufferToFiner(FaceField &src, Real *buf,
                                     const NeighborBlock& nb);
  void SendFieldBoundaryBuffers(FaceField &src, int step);
  void SetFieldBoundarySameLevel(FaceField &dst, Real *buf, const NeighborBlock& nb);
  void SetFieldBoundaryFromCoarser(Real *buf, const NeighborBlock& nb);
  void SetFieldBoundaryFromFiner(FaceField &dst, Real *buf, const NeighborBlock& nb);
  bool ReceiveFieldBoundaryBuffers(FaceField &dst, int step);
  void ReceiveFieldBoundaryBuffersWithWait(FaceField &dst, int step);

  int LoadEMFBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryBufferToCoarser(Real *buf, const NeighborBlock& nb);
  int LoadEMFBoundaryPolarBuffer(Real *buf, const PolarNeighborBlock &nb);
  void SendEMFCorrection(int step);
  void SetEMFBoundarySameLevel(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryFromFiner(Real *buf, const NeighborBlock& nb);
  void SetEMFBoundaryPolar(Real **buf_list, int num_bufs, bool north);
  void ClearCoarseEMFBoundary(void);
  void AverageEMFBoundary(void);
  void PolarSingleEMF(void);
  void PolarSingleHydro(AthenaArray<Real> &dst);
  void PolarSingleField(FaceField &dst);
  bool ReceiveEMFCorrection(int step);

  void ProlongateBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst, 
                            FaceField &bfdst, AthenaArray<Real> &bcdst);

  void ApplyPhysicalBoundaries(AthenaArray<Real> &pdst, AthenaArray<Real> &cdst,
                               FaceField &bfdst, AthenaArray<Real> &bcdst);

  void ClearBoundaryForInit(void);
  void ClearBoundaryAll(void);

private:
  MeshBlock *pmy_mblock_;  // ptr to MeshBlock containing this BVals

  BValFunc_t BoundaryFunction_[6];

  int nface_, nedge_;
  bool edge_flag_[12];
  int nedge_fine_[12];
  bool firsttime_[NSTEP];

  enum boundary_status hydro_flag_[NSTEP][56], field_flag_[NSTEP][56];
  enum boundary_status flcor_flag_[NSTEP][6][2][2];
  enum boundary_status emfcor_flag_[NSTEP][48];
  enum boundary_status *emf_north_flag_[NSTEP];
  enum boundary_status *emf_south_flag_[NSTEP];
  Real *hydro_send_[NSTEP][56],  *hydro_recv_[NSTEP][56];
  Real *field_send_[NSTEP][56],  *field_recv_[NSTEP][56];
  Real *flcor_send_[NSTEP][6],   *flcor_recv_[NSTEP][6][2][2];
  Real *emfcor_send_[NSTEP][48], *emfcor_recv_[NSTEP][48];
  Real **emf_north_send_[NSTEP], **emf_north_recv_[NSTEP];
  Real **emf_south_send_[NSTEP], **emf_south_recv_[NSTEP];
  AthenaArray<Real> sarea_[2];
  AthenaArray<Real> exc_;
  int num_north_polar_blocks_, num_south_polar_blocks_;

#ifdef MPI_PARALLEL
  MPI_Request req_hydro_send_[NSTEP][56],  req_hydro_recv_[NSTEP][56];
  MPI_Request req_field_send_[NSTEP][56],  req_field_recv_[NSTEP][56];
  MPI_Request req_flcor_send_[NSTEP][6],   req_flcor_recv_[NSTEP][6][2][2];
  MPI_Request req_emfcor_send_[NSTEP][48], req_emfcor_recv_[NSTEP][48];
  MPI_Request *req_emf_north_send_[NSTEP], *req_emf_north_recv_[NSTEP];
  MPI_Request *req_emf_south_send_[NSTEP], *req_emf_south_recv_[NSTEP];
#endif
  // temporary
  friend class Mesh;
};

unsigned int CreateBufferID(int ox1, int ox2, int ox3, int fi1, int fi2);
int BufferID(int dim, bool multilevel, bool face_only);
int FindBufferID(int ox1, int ox2, int ox3, int fi1, int fi2, int bmax);

typedef struct NeighborIndexes
{
  int ox1, ox2, ox3, fi1, fi2;
  enum neighbor_type type;
} NeighborIndexes;


#endif // BOUNDARY_VALUES_HPP
