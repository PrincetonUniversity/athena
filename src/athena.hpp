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

enum {IDN=0, IM1=1, IM2=2, IM3=3, IEN=4};
enum {IVX=1, IVY=2, IVZ=3, IBY=(NFLUID), IBZ=((NFLUID)+1)};
enum {IB1=0, IB2=1, IB3=2};
enum {X1E2=0, X1E3=1, X2E3=0, X2E1=1, X3E1=0, X3E2=1};
enum {I00, I01, I02, I03, I11, I12, I13, I22, I23, I33, NMETRIC};

enum direction {inner_x1=0, outer_x1=1, inner_x2=2, outer_x2=3, inner_x3=4, outer_x3=5};
enum face {x1face=0, x2face=1, x3face=2};
enum rwmode {readmode,writemode};
enum mpitag {tag_fluid=0, tag_field=1}; // mpitag must be < 16 and unique
enum neighbor_type {neighbor_none, neighbor_face, neighbor_edge, neighbor_corner};

enum task {
  none=0, 

  fluid_integrate_1=1L<<0, field_integrate_1=1L<<1,
  fluid_send_1=1L<<2, fluid_recv_1=1L<<3,
  flux_correction_send_1=1L<<4, flux_correction_recv_1=1L<<5,
  fluid_prolongation_1=1L<<6, fluid_boundary_1=1L<<7,
  field_send_1=1L<<8, field_recv_1=1L<<9,
  emf_correction_send_1=1L<<10, emf_correction_recv_1=1L<<11,
  field_prolongation_1=1L<<12, field_boundary_1=1L<<13,
  primitives_1=1L<<14,

  fluid_integrate_0=1L<<15, field_integrate_0=1L<<16,
  fluid_send_0=1L<<17, fluid_recv_0=1L<<18,
  flux_correction_send_0=1L<<19, flux_correction_recv_0=1L<<20,
  fluid_prolongation_0=1L<<21, fluid_boundary_0=1L<<22,
  field_send_0=1L<<23, field_recv_0=1L<<24,
  emf_correction_send_0=1L<<25, emf_correction_recv_0=1L<<26,
  field_prolongation_0=1L<<27, field_boundary_0=1L<<28,
  primitives_0=1L<<29,

  new_blocktimestep=1L<<30
};

enum tasklist_status { tl_running, tl_stuck, tl_complete, tl_nothing };

enum task_status { task_failure, task_success, task_donext};

enum boundary_status {boundary_waiting, boundary_arrived, boundary_completed};

extern int myrank, nproc;

#endif
