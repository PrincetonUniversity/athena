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

//! \struct LogicalLocation
//  \brief logical location and level of meshblocks
typedef struct LogicalLocation {
  long int lx1, lx2, lx3;
  int level;
  LogicalLocation() : lx1(-1), lx2(-1), lx3(-1), level(-1) {};
} LogicalLocation;



enum {IDN=0, IM1=1, IM2=2, IM3=3, IEN=4};
enum {IVX=1, IVY=2, IVZ=3, IBY=(NHYDRO), IBZ=((NHYDRO)+1)};
enum {IB1=0, IB2=1, IB3=2};
enum {X1E2=0, X1E3=1, X2E3=0, X2E1=1, X3E1=0, X3E2=1};
enum {I00, I01, I02, I03, I11, I12, I13, I22, I23, I33, NMETRIC};

enum direction {dir_undefined=-1, inner_x1=0, outer_x1=1, inner_x2=2, outer_x2=3, inner_x3=4, outer_x3=5};
enum edgeid {edgeid_undefined = -1, em2m1=0, em2p1=1, ep2m1=2, ep2p2=3, 
                em3m1=4, em3p1=5, ep3m1=6, ep3p1=7, em3m2=8, em3p2=9, ep3m2=10, ep3p2=11};
enum face {x1face=0, x2face=1, x3face=2};
enum mpitag {tag_hydro=0, tag_field=1, tag_flcor=2, tag_emfcor=3}; // mpitag must be < 16 and unique
enum neighbor_type {neighbor_none, neighbor_face, neighbor_edge, neighbor_corner};

enum boundary_status {boundary_waiting, boundary_arrived, boundary_completed};

#endif // define ATHENA_HPP
