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
enum mpitag {tag_fluid=0, tag_field=1, tag_eflux=2}; // mpitag must be < 16

enum ActionOnBlock
  {pgen,          primitives_n, primitives_nhalf, new_blocktimestep,
  fluid_predict, fluid_correct,   field_predict, field_correct,
  fluid_start_recv_n, fluid_start_recv_nhalf, 
  fluid_loadsend_bcsx1_n, fluid_loadsend_bcsx2_n, fluid_loadsend_bcsx3_n,
  fluid_recvset_bcsx1_n, fluid_recvset_bcsx2_n, fluid_recvset_bcsx3_n,
  fluid_loadsend_bcsx1_nhalf, fluid_loadsend_bcsx2_nhalf, fluid_loadsend_bcsx3_nhalf,
  fluid_recvset_bcsx1_nhalf, fluid_recvset_bcsx2_nhalf, fluid_recvset_bcsx3_nhalf,
  fluid_waitsend_bcsx1, fluid_waitsend_bcsx2, fluid_waitsend_bcsx3,
  field_start_recv_n, field_start_recv_nhalf,
  field_loadsend_bcsx1_n, field_loadsend_bcsx2_n, field_loadsend_bcsx3_n,
  field_recvset_bcsx1_n, field_recvset_bcsx2_n, field_recvset_bcsx3_n,
  field_loadsend_bcsx1_nhalf, field_loadsend_bcsx2_nhalf, field_loadsend_bcsx3_nhalf,
  field_recvset_bcsx1_nhalf, field_recvset_bcsx2_nhalf, field_recvset_bcsx3_nhalf,
  field_waitsend_bcsx1, field_waitsend_bcsx2, field_waitsend_bcsx3,
  eflux_start_recv_n, eflux_start_recv_nhalf,
  eflux_loadsend_bcs_n, eflux_loadsend_bcs_nhalf, 
  eflux_recvset_bcs_n, eflux_recvset_bcs_nhalf, 
  eflux_waitsend
  };

extern int myrank, nproc;

#endif
