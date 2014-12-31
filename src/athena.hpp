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
enum mpitag {tag_fluid=0, tag_field=1}; // mpitag must be < 16 and unique

enum task {  // reserved for up to 4-step integrators
  none=0L, new_blocktimestep=1L<<0,
  primitives_0=1L<<1, primitives_1=1L<<2, primitives_2=1L<<3, primitives_3=1L<<4,
  fluid_integrate_sendx1_0=1L<<5, fluid_integrate_sendx1_1=1L<<6,
  fluid_integrate_sendx1_2=1L<<7, fluid_integrate_sendx1_3=1L<<8,
  field_integrate_sendx1_0=1L<<9, field_integrate_sendx1_1=1L<<10,
  field_integrate_sendx1_2=1L<<11, field_integrate_sendx1_3=1L<<12,
  fluid_startrecv_0=1L<<13, fluid_startrecv_1=1L<<14,
  fluid_startrecv_2=1L<<15, fluid_startrecv_3=1L<<16,
  field_startrecv_0=1L<<17, field_startrecv_1=1L<<18,
  field_startrecv_2=1L<<19, field_startrecv_3=1L<<20,
  fluid_recvx1_0=1L<<21, fluid_recvx1_1=1L<<22, // for 1D
  fluid_recvx1_2=1L<<23, fluid_recvx1_3=1L<<24,
  field_recvx1_0=1L<<25, field_recvx1_1=1L<<26,
  field_recvx1_2=1L<<27, field_recvx1_3=1L<<28,
  fluid_recvx1_sendx2_0=1L<<21, fluid_recvx1_sendx2_1=1L<<22,
  fluid_recvx1_sendx2_2=1L<<23, fluid_recvx1_sendx2_3=1L<<24,
  field_recvx1_sendx2_0=1L<<25, field_recvx1_sendx2_1=1L<<26,
  field_recvx1_sendx2_2=1L<<27, field_recvx1_sendx2_3=1L<<28,
  fluid_recvx2_0=1L<<29, fluid_recvx2_1=1L<<30, // for 2D
  fluid_recvx2_2=1L<<31, fluid_recvx2_3=1L<<32,
  field_recvx2_0=1L<<33, field_recvx2_1=1L<<34,
  field_recvx2_2=1L<<35, field_recvx2_3=1L<<36,
  fluid_recvx2_sendx3_0=1L<<29, fluid_recvx2_sendx3_1=1L<<30,
  fluid_recvx2_sendx3_2=1L<<31, fluid_recvx2_sendx3_3=1L<<32,
  field_recvx2_sendx3_0=1L<<33, field_recvx2_sendx3_1=1L<<34,
  field_recvx2_sendx3_2=1L<<35, field_recvx2_sendx3_3=1L<<36,
  fluid_recvx3_0=1L<<37, fluid_recvx3_1=1L<<38,
  fluid_recvx3_2=1L<<39, fluid_recvx3_3=1L<<40,
  field_recvx3_0=1L<<41, field_recvx3_1=1L<<42,
  field_recvx3_2=1L<<43, field_recvx3_3=1L<<44
};

enum tlstatus { running, stuck, complete, nothing };

extern int myrank, nproc, tag_shift;

#endif
