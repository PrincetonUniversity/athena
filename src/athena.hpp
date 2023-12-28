#ifndef ATHENA_HPP_
#define ATHENA_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file athena.hpp
//! \brief contains Athena++ general purpose types, structures, enums, etc.

// C headers

// C++ headers
#include <cmath>
#include <cstdint>  // std::int64_t

// Athena++ headers
#include "athena_arrays.hpp"
#include "defs.hpp"

// primitive type alias that allows code to run with either floats or doubles
#if SINGLE_PRECISION_ENABLED
using Real = float;
#ifdef MPI_PARALLEL
#define MPI_ATHENA_REAL MPI_FLOAT
#endif
#else
using Real = double;
#ifdef MPI_PARALLEL
#define MPI_ATHENA_REAL MPI_DOUBLE
#endif
#endif

// for OpenMP 4.0 SIMD vectorization, control width of SIMD lanes
#if defined(__AVX512F__)
#define SIMD_WIDTH 8
#elif defined(__AVX__)
#define SIMD_WIDTH 4
#elif defined(__SSE2__)
#define SIMD_WIDTH 2
#else
#define SIMD_WIDTH 4
#endif

#define CACHELINE_BYTES 64

// forward declarations needed for function pointer type aliases
class MeshBlock;
class Coordinates;
class ParameterInput;
class HydroDiffusion;
class FieldDiffusion;
class MGCoordinates;
class OrbitalAdvection;
class NRRadiation;
class IMRadiation;
class CosmicRay;

//--------------------------------------------------------------------------------------
//! \struct LogicalLocation
//! \brief stores logical location and level of MeshBlock

struct LogicalLocation { // aggregate and POD type
  // These values can exceed the range of std::int32_t even if the root grid has only a
  // single MeshBlock if >30 levels of AMR are used, since the corresponding max index =
  // 1*2^31 > INT_MAX = 2^31 -1 for most 32-bit signed integer type impelementations
  std::int64_t lx1, lx2, lx3;
  int level;
  // comparison functions for sorting
  static bool Lesser(const LogicalLocation &left, const LogicalLocation &right) {
    return left.level < right.level;
  }
  static bool Greater(const LogicalLocation & left, const LogicalLocation &right) {
    return left.level > right.level;
  }
};

//! prototype for overloading the comparison operator (defined in meshblock_tree.cpp)
bool operator==(const LogicalLocation &l1, const LogicalLocation &l2);

//----------------------------------------------------------------------------------------
//! \struct RegionSize
//! \brief physical size and number of cells in a Mesh or a MeshBlock

struct RegionSize {  // aggregate and POD type; do NOT reorder member declarations:
  Real x1min, x2min, x3min;
  Real x1max, x2max, x3max;
  Real x1rat, x2rat, x3rat; // ratio of dxf(i)/dxf(i-1)
  // the size of the root grid or a MeshBlock should not exceed std::int32_t limits
  int nx1, nx2, nx3;        // number of active cells (not including ghost zones)
};

//---------------------------------------------------------------------------------------
//! \struct FaceField
//! \brief container for face-centered fields

struct FaceField {
  AthenaArray<Real> x1f, x2f, x3f;
  FaceField() = default;
  FaceField(int ncells3, int ncells2, int ncells1,
            AthenaArray<Real>::DataStatus init=AthenaArray<Real>::DataStatus::allocated) :
      x1f(ncells3, ncells2, ncells1+1, init), x2f(ncells3, ncells2+1, ncells1, init),
      x3f(ncells3+1, ncells2, ncells1, init) {}
};

//----------------------------------------------------------------------------------------
//! \struct EdgeField
//! \brief container for edge-centered fields

struct EdgeField {
  AthenaArray<Real> x1e, x2e, x3e;
  EdgeField() = default;
  EdgeField(int ncells3, int ncells2, int ncells1,
            AthenaArray<Real>::DataStatus init=AthenaArray<Real>::DataStatus::allocated) :
      x1e(ncells3+1, ncells2+1, ncells1, init), x2e(ncells3+1, ncells2, ncells1+1, init),
      x3e(ncells3, ncells2+1, ncells1+1, init) {}
};

//----------------------------------------------------------------------------------------
// enums used everywhere
// (not specifying underlying integral type (C++11) for portability & performance)

//! \todo (felker):
//! - C++ Core Guidelines Enum.5: Don’t use ALL_CAPS for enumerators
//!   (avoid clashes with preprocessor macros).
//! - Enumerated type definitions in this file and:
//!   athena_fft.hpp, io_wrapper.hpp, bvals.hpp, hydro_diffusion.hpp, field_diffusion.hpp,
//!   task_list.hpp, ???

//------------------
// named, weakly typed / unscoped enums:
//------------------

// enumerators only used for indexing AthenaArray and regular arrays; typename and
// explicitly specified enumerator values aare unnecessary, but provided for clarity:

//! array indices for conserved: density, momemtum, total energy
enum ConsIndex {IDN=0, IM1=1, IM2=2, IM3=3, IEN=4};
//! array indices for face-centered field
enum MagneticIndex {IB1=0, IB2=1, IB3=2};

//! array indices for 1D primitives: velocity, transverse components of field
enum PrimIndex {IVX=1, IVY=2, IVZ=3, IPR=4, IBY=(NHYDRO), IBZ=((NHYDRO)+1)};

//! array indices for face-centered electric fields returned by Riemann solver
enum ElectricIndex {X1E2=0, X1E3=1, X2E3=0, X2E1=1, X3E1=0, X3E2=1};

//! array indices for metric matrices in GR
enum MetricIndex {I00=0, I01=1, I02=2, I03=3, I11=4, I12=5, I13=6, I22=7, I23=8, I33=9,
                  NMETRIC=10};
//! array indices for triangular matrices in GR
enum TriangleIndex {T00=0, T10=1, T11=2, T20=3, T21=4, T22=5, T30=6, T31=7, T32=8, T33=9,
                    NTRIANGULAR=10};

// enumerator types that are used for variables and function parameters:

// needed for arrays dimensioned over grid directions
// enumerator type only used in Mesh::EnrollUserMeshGenerator()

//! array indices for grid directions
enum CoordinateDirection {X1DIR=0, X2DIR=1, X3DIR=2};

//------------------
// strongly typed / scoped enums (C++11):
//------------------
// KGF: Except for the 2x MG* enums, these may be unnessary w/ the new class inheritance
// Now, only passed to BoundaryVariable::InitBoundaryData(); could replace w/ bool switch
// TODO(tomo-ono): consider necessity of orbita_cc and orbital_fc
enum class BoundaryQuantity {cc, fc, cc_flcor, fc_flcor, mg, mg_faceonly, mg_coeff,
                             orbital_cc, orbital_fc};
enum class HydroBoundaryQuantity {cons, prim};
enum class BoundaryCommSubset {mesh_init, gr_amr, all, orbital, radiation, radhydro};
// TODO(felker): consider generalizing/renaming to QuantityFormulation
// TODO(Gong): currently disabled=background (with passive scalar advection),
// and fixed is without passive scalar advection.
enum class FluidFormulation {evolve, background, fixed, disabled};
enum class TaskType {op_split_before, main_int, op_split_after};
enum class UserHistoryOperation {sum, max, min};

//----------------------------------------------------------------------------------------
// function pointer prototypes for user-defined modules set at runtime

using BValFunc = void (*)(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt,
    int is, int ie, int js, int je, int ks, int ke, int ngh);
using AMRFlagFunc = int (*)(MeshBlock *pmb);
using MeshGenFunc = Real (*)(Real x, RegionSize rs);
using SrcTermFunc = void (*)(
    MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
    const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar);
using TimeStepFunc = Real (*)(MeshBlock *pmb);
using HistoryOutputFunc = Real (*)(MeshBlock *pmb, int iout);
using MetricFunc = void (*)(
    Real x1, Real x2, Real x3, ParameterInput *pin,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv,
    AthenaArray<Real> &dg_dx1, AthenaArray<Real> &dg_dx2, AthenaArray<Real> &dg_dx3);
using MGBoundaryFunc = void (*)(
    AthenaArray<Real> &dst, Real time, int nvar,
    int is, int ie, int js, int je, int ks, int ke, int ngh,
    const MGCoordinates &coord);
using ViscosityCoeffFunc = void (*)(
    HydroDiffusion *phdif, MeshBlock *pmb,
    const  AthenaArray<Real> &w, const AthenaArray<Real> &bc,
    int is, int ie, int js, int je, int ks, int ke);
using ConductionCoeffFunc = void (*)(
    HydroDiffusion *phdif, MeshBlock *pmb,
    const AthenaArray<Real> &w, const AthenaArray<Real> &bc,
    int is, int ie, int js, int je, int ks, int ke);
using FieldDiffusionCoeffFunc = void (*)(
    FieldDiffusion *pfdif, MeshBlock *pmb,
    const AthenaArray<Real> &w,
    const AthenaArray<Real> &bmag,
    int is, int ie, int js, int je, int ks, int ke);
using MGMaskFunc = void (*)(AthenaArray<Real> &dat,
    int is, int ie, int js, int je, int ks, int ke, const MGCoordinates &coord);
using OrbitalVelocityFunc = Real (*)(
    OrbitalAdvection *porb, Real x1, Real x2, Real x3);
using RadBoundaryFunc = void (*)(
     MeshBlock *pmb, Coordinates *pco, NRRadiation *prad,
     const AthenaArray<Real> &w, FaceField &b,
     AthenaArray<Real> &ir,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
using OpacityFunc = void (*)(MeshBlock *pmb, AthenaArray<Real> &prim);
using FrequencyFunc = void (*)(NRRadiation *prad);
using EmissionFunc = void(*)(NRRadiation *prad, Real tgas);
using CROpacityFunc = void (*)(MeshBlock *pmb, AthenaArray<Real> &u_cr,
                      AthenaArray<Real> &prim, AthenaArray<Real> &bcc);
using CRStreamingFunc = void (*)(MeshBlock *pmb, AthenaArray<Real> &u_cr,
                      AthenaArray<Real> &prim, AthenaArray<Real> &bcc,
                      AthenaArray<Real> &grad_pc, int k, int j, int is, int ie);
using SRJFunc = void (*)(IMRadiation *pimrad);

using CRBoundaryFunc = void (*)(
     MeshBlock *pmb, Coordinates *pco, CosmicRay *pcr,
     const AthenaArray<Real> &w, FaceField &b, AthenaArray<Real> &u_cr,
     Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
using CRSrcTermFunc = void (*)(
    MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim, FaceField &b, AthenaArray<Real> &u_cr);

#endif // ATHENA_HPP_
