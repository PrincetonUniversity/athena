#ifndef ATHENA_HPP_
#define ATHENA_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file athena.hpp
//  \brief contains Athena++ general purpose types, structures, enums, etc.

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

//--------------------------------------------------------------------------------------
//! \struct LogicalLocation
//  \brief stores logical location and level of MeshBlock

struct LogicalLocation {
  // These values can exceed the range of std::int32_t even if the root grid has only a
  // single MeshBlock if >30 levels of AMR are used, since the corresponding max index =
  // 1*2^31 > INT_MAX = 2^31 -1 for most 32-bit signed integer type impelementations
  std::int64_t lx1, lx2, lx3;
  int level;

  LogicalLocation() : lx1(-1), lx2(-1), lx3(-1), level(-1) {}

  // operators useful for sorting
  bool operator==(LogicalLocation &ll) {
    return ((ll.level==level) && (ll.lx1==lx1) && (ll.lx2==lx2) && (ll.lx3==lx3));
  }
  static bool Lesser(const LogicalLocation &left, const LogicalLocation &right) {
    return left.level < right.level;
  }
  static bool Greater(const LogicalLocation & left, const LogicalLocation &right) {
    return left.level > right.level;
  }
};

//----------------------------------------------------------------------------------------
//! \struct RegionSize
//  \brief physical size and number of cells in a Mesh or a MeshBlock

struct RegionSize {
  Real x1min, x2min, x3min;
  Real x1max, x2max, x3max;
  Real x1rat, x2rat, x3rat; // ratio of x(i)/x(i-1)
  // the size of the root grid or a MeshBlock should not exceed std::int32_t limits
  int nx1, nx2, nx3;        // number of active cells (not including ghost zones)
};

//---------------------------------------------------------------------------------------
//! \struct FaceField
//  \brief container for face-centered fields

struct FaceField {
  AthenaArray<Real> x1f,x2f,x3f;
};

//----------------------------------------------------------------------------------------
//! \struct EdgeField
//  \brief container for edge-centered fields

struct EdgeField {
  AthenaArray<Real> x1e,x2e,x3e;
};

//----------------------------------------------------------------------------------------
// enums used everywhere
// (not specifying underlying integral type (C++11) for portability & performance)

// TODO(felker): C++ Core Guidelines Enum.5: Don’t use ALL_CAPS for enumerators
// (avoid clashes with preprocessor macros). Enumerated type definitions in this file and:
// athena_fft.hpp, io_wrapper.hpp, bvals.hpp, hydro_diffusion.hpp, field_diffusion.hpp,
// task_list.hpp, ???

//------------------
// anonymous, weakly typed / unscoped enums:
//------------------

// TODO(felker): C++ Core Guidelines Enum.6: Avoid unnamed enumerations
// Either use "constexpr int idn=0;" (if the enumerators are unrelated)
// Or name it as an unscoped enum type (even if the type name is never used?)

// array indices for conserved: density, momemtum, total energy, face-centered field
enum {IDN, IM1, IM2, IM3, IEN};
enum {IB1, IB2, IB3};

// array indices for 1D primitives: velocity, transverse components of field
enum {IVX=1, IVY, IVZ, IPR, IBY=(NHYDRO), IBZ=((NHYDRO)+1)};

// array indices for face-centered electric fields returned by Riemann solver
enum {X1E2=0, X1E3=1, X2E3=0, X2E1=1, X3E1=0, X3E2=1};

// array indices for metric and triangular matrices in GR
enum {I00, I01, I02, I03, I11, I12, I13, I22, I23, I33, NMETRIC};
enum {T00, T10, T11, T20, T21, T22, T30, T31, T32, T33, NTRIANGULAR};

//------------------
// named, weakly typed / unscoped enums:
//------------------
// needed for arrays dimensioned over grid directions
enum CoordinateDirection {X1DIR, X2DIR, X3DIR};
// TODO(felker): probably remove the enum typename "CoordinateDirection" since it is
// only used in Mesh::EnrollUserMeshGenerator(enum CoordinateDirection,MeshGenFunc my_mg)

//------------------
// strongly typed / scoped enums (C++11):
//------------------
// needed wherever MPI communications are used.  Must be < 32 and unique
enum class AthenaTagMPI {hydro, field, rad, chem, hydflx, fldflx, radflx, chmflx, amr,
                         fldflx_pole, gravity, mggrav,
                         shbox_hydro, shbox_field, shbox_emf};

enum class BoundaryQuantity {hydro, field, gravity, mggrav, mggravf, flcor, emfcor};
enum class CCBoundaryQuantity {cons, prim};
enum class FluxCorrectionQuantity {hydro};

//----------------------------------------------------------------------------------------
// function pointer prototypes for user-defined modules set at runtime

using BValFunc = void (*)(
    MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
    Real time, Real dt,
    int is, int ie, int js, int je, int ks, int ke, int ngh);
using AMRFlagFunc = int (*)(MeshBlock *pmb);
using MeshGenFunc = Real (*)(Real x, RegionSize rs);
using SrcTermFunc = void (*)(
    MeshBlock *pmb, const Real time, const Real dt,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
using TimeStepFunc = Real (*)(MeshBlock *pmb);
using HistoryOutputFunc = Real (*)(MeshBlock *pmb, int iout);
using MetricFunc = void (*)(
    Real x1, Real x2, Real x3, ParameterInput *pin,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv,
    AthenaArray<Real> &dg_dx1, AthenaArray<Real> &dg_dx2, AthenaArray<Real> &dg_dx3);
using MGBoundaryFunc = void (*)(
    AthenaArray<Real> &dst,Real time, int nvar,
    int is, int ie, int js, int je, int ks, int ke, int ngh,
    Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);
using GravityBoundaryFunc = void (*)(
    MeshBlock *pmb, Coordinates *pco,
    AthenaArray<Real> &dst, Real time, Real dt,
    int is, int ie, int js, int je, int ks, int ke);
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

#endif // ATHENA_HPP_
