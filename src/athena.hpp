#ifndef ATHENA_HPP
#define ATHENA_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file athena.hpp
//  \brief contains Athena++ general purpose types, structures, enums, etc.
//======================================================================================

#include "defs.hpp"
#include "athena_arrays.hpp"
#include <math.h>

// typedefs that allow code to run with either floats or doubles
typedef double Real;
#ifdef MPI_PARALLEL
#define MPI_ATHENA_REAL MPI_DOUBLE
#endif

class MeshBlock;
class Coordinates;
struct RegionSize;

//--------------------------------------------------------------------------------------
//! \struct FaceField
//  \brief container for face-centered fields

typedef struct FaceField {
  AthenaArray<Real> x1f,x2f,x3f;
} FaceField;

//--------------------------------------------------------------------------------------
//! \struct EdgeField
//  \brief container for edge-centered fields

typedef struct EdgeField {
  AthenaArray<Real> x1e,x2e,x3e;
} EdgeField;

//--------------------------------------------------------------------------------------
// enums used everywhere

// array indices for conserved: density, momemtum, total energy, face-centered field
enum {IDN=0, IM1=1, IM2=2, IM3=3, IEN=4};
enum {IB1=0, IB2=1, IB3=2};

// array indices for 1D primitives: velocity, transverse components of field
enum {IVX=1, IVY=2, IVZ=3, IPR=4, IBY=(NHYDRO), IBZ=((NHYDRO)+1)};

// array indices for face-centered electric fields returned by Riemann solver
enum {X1E2=0, X1E3=1, X2E3=0, X2E1=1, X3E1=0, X3E2=1};

// array indices for metric in GR
enum {I00, I01, I02, I03, I11, I12, I13, I22, I23, I33, NMETRIC};

// needed for arrays dimensioned over grid directions
enum CoordinateDirection {X1DIR=0, X2DIR=1, X3DIR=2};

// needed wherever MPI communications are used.  Must be < 16 and unique
// [JMSHI
enum Athena_MPI_Tag {TAG_HYDRO=0, TAG_FIELD=1, TAG_RAD=2, TAG_CHEM=3, TAG_HYDFLX=4,
  TAG_FLDFLX=5, TAG_RADFLX=6, TAG_CHMFLX=7, TAG_AMR=8, TAG_FLDFLX_POLE=9, TAG_WTLIM=10,
  TAG_SHBOX_HYDRO=13,TAG_SHBOX_FIELD=14,TAG_SHBOX_EMF=15};
//enum Athena_MPI_Tag {TAG_HYDRO=0, TAG_FIELD=1, TAG_RAD=2, TAG_CHEM=3, TAG_HYDFLX=4,
//  TAG_FLDFLX=5, TAG_RADFLX=6, TAG_CHMFLX=7, TAG_AMR=8, TAG_FLDFLX_POLE=9, TAG_WTLIM=10};
//JMSHI]
//
//--------------------------------------------------------------------------------------
// function pointer prototypes for user-defined modules set at runtime

typedef void (*BValFunc_t)(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &buf,
  FaceField &buf2, int is, int ie, int js, int je, int ks, int ke);
typedef int (*AMRFlagFunc_t)(MeshBlock *pmb);
typedef Real (*MeshGenFunc_t)(Real x, RegionSize rs);
typedef void (*SrcTermFunc_t)(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

#endif // ATHENA_HPP
