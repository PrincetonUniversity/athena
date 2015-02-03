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

typedef double Real;
#ifdef MPI_PARALLEL
#define MPI_ATHENA_REAL MPI_DOUBLE
#endif

typedef struct InterfaceField {
  AthenaArray<Real> x1f,x2f,x3f;
} InterfaceField;

typedef struct EdgeField {
  AthenaArray<Real> x1e,x2e,x3e;
} EdgeField;


//! \struct NeighborBlock
//  \brief neighbor rank, level, and ids

typedef struct NeighborBlock {
  int rank, level, gid, lid;
  NeighborBlock() : rank(-1), level(-1), gid(-1), lid(-1) {};
} NeighborBlock;

enum {IDN=0, IM1=1, IM2=2, IM3=3, IEN=4};
enum {IVX=1, IVY=2, IVZ=3, IBY=(NFLUID), IBZ=((NFLUID)+1)};
enum {IB1=0, IB2=1, IB3=2};
enum {X1E2=0, X1E3=1, X2E3=0, X2E1=1, X3E1=0, X3E2=1};
enum {I00, I01, I02, I03, I11, I12, I13, I22, I23, I33, NMETRIC};

enum direction {inner_x1=0, outer_x1=1, inner_x2=2, outer_x2=3, inner_x3=4, outer_x3=5};
enum face {x1face=0, x2face=1, x3face=2};
enum rwmode {readmode,writemode};
enum mpitag {tag_fluid=0, tag_field=1, tag_eflux=2}; // mpitag must be < 16 and unique

enum task {
  none=0, 

  primitives_0=1,
  fluid_integrate_sendx1_0=2,
  eflux_recv_0=3,
  field_integrate_sendx1_0=4,
  fluid_recvx1_0=5, field_recvx1_0=6, // for 1D
  fluid_recvx1_sendx2_0=5, field_recvx1_sendx2_0=6, 
  fluid_recvx2_0=7, field_recvx2_0=8, // for2D
  fluid_recvx2_sendx3_0=7, field_recvx2_sendx3_0=8,
  fluid_recvx3_0=9, field_recvx3_0=10,

  primitives_1=11,
  fluid_integrate_sendx1_1=12,
  eflux_recv_1=13,
  field_integrate_sendx1_1=14,
  fluid_recvx1_1=15, field_recvx1_1=16, 
  fluid_recvx1_sendx2_1=15, field_recvx1_sendx2_1=16,
  fluid_recvx2_1=17, field_recvx2_1=18,
  fluid_recvx2_sendx3_1=17, field_recvx2_sendx3_1=18,
  fluid_recvx3_1=19, field_recvx3_1=20,

  new_blocktimestep=21
};

enum tlstatus { running, stuck, complete, nothing };

extern int myrank, nproc;

#endif
