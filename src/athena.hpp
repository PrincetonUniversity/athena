#ifndef ATHENA_HPP
#define ATHENA_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file athena.hpp
//  \brief contains Athena++ general purpose types, structures, enums, etc.

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
class ParameterInput;
struct RegionSize;


//---------------------------------------------------------------------------------------
//! \struct FaceField
//  \brief container for face-centered fields

typedef struct FaceField {
  AthenaArray<Real> x1f,x2f,x3f;
} FaceField;

//----------------------------------------------------------------------------------------
//! \struct EdgeField
//  \brief container for edge-centered fields

typedef struct EdgeField {
  AthenaArray<Real> x1e,x2e,x3e;
} EdgeField;

//----------------------------------------------------------------------------------------
// enums used everywhere

// array indices for conserved: density, momemtum, total energy, face-centered field 
enum {IDN=0, IM1=1, IM2=2, IM3=3, IEN=4};
enum {IB1=0, IB2=1, IB3=2};

// array indices for 1D primitives: velocity, transverse components of field
enum {IVX=1, IVY=2, IVZ=3, IPR=4, IBY=(NHYDRO), IBZ=((NHYDRO)+1)};

// array indices for face-centered electric fields returned by Riemann solver
enum {X1E2=0, X1E3=1, X2E3=0, X2E1=1, X3E1=0, X3E2=1};

// array indices for metric and triangular matrices in GR
enum {I00, I01, I02, I03, I11, I12, I13, I22, I23, I33, NMETRIC};
enum {T00, T10, T11, T20, T21, T22, T30, T31, T32, T33, NTRIANGULAR};

// needed for arrays dimensioned over grid directions
enum CoordinateDirection {X1DIR=0, X2DIR=1, X3DIR=2};

// needed wherever MPI communications are used.  Must be < 32 and unique
enum Athena_MPI_Tag {TAG_HYDRO=0, TAG_FIELD=1, TAG_RAD=2, TAG_CHEM=3, TAG_HYDFLX=4,
  TAG_FLDFLX=5, TAG_RADFLX=6, TAG_CHMFLX=7, TAG_AMR=8, TAG_FLDFLX_POLE=9, TAG_GRAVITY=11,
  TAG_MGGRAV=12};

enum BoundaryType {BNDRY_HYDRO=0, BNDRY_FIELD=1, BNDRY_GRAVITY=2, BNDRY_MGGRAV=3,
                   BNDRY_MGGRAVF=4, BNDRY_FLCOR=5, BNDRY_EMFCOR=6};
enum CCBoundaryType {HYDRO_CONS=0, HYDRO_PRIM=1};
enum FluxCorrectionType {FLUX_HYDRO=0};

//----------------------------------------------------------------------------------------
// function pointer prototypes for user-defined modules set at runtime

typedef void (*BValFunc_t)(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
  FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
typedef int (*AMRFlagFunc_t)(MeshBlock *pmb);
typedef Real (*MeshGenFunc_t)(Real x, RegionSize rs);
typedef void (*SrcTermFunc_t)(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
typedef Real (*TimeStepFunc_t)(MeshBlock *pmb);
typedef Real (*HistoryOutputFunc_t)(MeshBlock *pmb, int iout);
typedef void (*MetricFunc_t)(Real x1, Real x2, Real x3, ParameterInput *pin,
             AthenaArray<Real> &g, AthenaArray<Real> &g_inv, AthenaArray<Real> &dg_dx1,
             AthenaArray<Real> &dg_dx2, AthenaArray<Real> &dg_dx3);
typedef void (*MGBoundaryFunc_t)(AthenaArray<Real> &dst,Real time, int nvar,
             int is, int ie, int js, int je, int ks, int ke, int ngh,
             Real x0, Real y0, Real z0, Real dx, Real dy, Real dz);


#endif // ATHENA_HPP
