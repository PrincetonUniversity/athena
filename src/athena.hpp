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

enum direction {dir_undefined=-1, inner_x1=0, outer_x1=1, inner_x2=2, outer_x2=3, inner_x3=4, outer_x3=5};
enum direction {edgeid_undefined = -1, em2m1=0, em2p1=1, ep2m1=2, ep2p2=3, 
                em3m1=4, em3p1=5, ep3m1=6, ep3p2=7, em3m2=8, em3p2=9, ep3m2=10, ep3p2=11};
enum face {x1face=0, x2face=1, x3face=2};
enum rwmode {readmode,writemode};
enum mpitag {tag_fluid=0, tag_field=1, tag_flcor=2, tag_emfcor_face=3, tag_emfcor_edge=4}; // mpitag must be < 16 and unique
enum neighbor_type {neighbor_none, neighbor_face, neighbor_edge, neighbor_corner};

enum task {
  none=0, 

  fluid_integrate_1=1L<<0, calculate_emf_1=1L<<1, field_integrate=1L<<2,
  fluid_send_1=1L<<3, fluid_recv_1=1L<<4,
  flux_correct_send_1=1L<<5, flux_correct_recv_1=1L<<6,
  fluid_prolong_1=1L<<7, fluid_boundary_1=1L<<8,
  field_send_1=1L<<9, field_recv_1=1L<<10,
  emf_correct_send_1=1L<<11, emf_correct_recv_1=1L<<12,
  field_prolong_1=1L<<13, field_boundary_1=1L<<14,
  primitives_1=1L<<15,

  fluid_integrate_0=1L<<16, calculate_emf_0=1L<<17, field_integrate=1L<<18,
  fluid_send_0=1L<<19, fluid_recv_0=1L<<20,
  flux_correct_send_0=1L<<21, flux_correct_recv_0=1L<<22,
  fluid_prolong_0=1L<<23, fluid_boundary_0=1L<<24,
  field_send_0=1L<<25, field_recv_0=1L<<26,
  emf_correct_send_0=1L<<27, emf_correct_recv_0=1L<<28,
  field_prolong_0=1L<<29, field_boundary_0=1L<<30,
  primitives_0=1L<<31,

  new_blocktimestep=1L<<32
};

enum tasklist_status { tl_running, tl_stuck, tl_complete, tl_nothing };

enum task_status { task_failure, task_success, task_donext};

enum boundary_status {boundary_waiting, boundary_arrived, boundary_completed};

extern int myrank, nproc;

#endif
