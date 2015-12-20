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

//! \struct FaceField
//  \brief container for face-centered fields
typedef struct FaceField {
  AthenaArray<Real> x1f,x2f,x3f;
} FaceField;

//! \struct EdgeField
//  \brief container for edge-centered fields
typedef struct EdgeField {
  AthenaArray<Real> x1e,x2e,x3e;
} EdgeField;

//! \struct LogicalLocation
//  \brief logical location and level of meshblocks
typedef struct LogicalLocation {
  long int lx1, lx2, lx3;
  int level;
  LogicalLocation() : lx1(-1), lx2(-1), lx3(-1), level(-1) {};
  // for sort from the finest level
  bool operator==(LogicalLocation &rloc) { return ((rloc.level==level) &&
                 (rloc.lx1==lx1) && (rloc.lx2==lx2) && (rloc.lx3==lx3)); }
  static bool Less(const LogicalLocation & lloc, const LogicalLocation &rloc)
  { return lloc.level < rloc.level; };
  static bool Greater(const LogicalLocation & lloc, const LogicalLocation &rloc)
  { return lloc.level > rloc.level; };
} LogicalLocation;

// array indices for conserved: density, momemtum, total energy, face-centered field 
enum {IDN=0, IM1=1, IM2=2, IM3=3, IEN=4};
enum {IB1=0, IB2=1, IB3=2};

// array indices for 1D primitives: velocity, transverse components of field
enum {IVX=1, IVY=2, IVZ=3, IBY=(NHYDRO), IBZ=((NHYDRO)+1)};

// array indices for face-centered electric fields returned by Riemann solver
enum {X1E2=0, X1E3=1, X2E3=0, X2E1=1, X3E1=0, X3E2=1};

// array indices for metric in GR
enum {I00, I01, I02, I03, I11, I12, I13, I22, I23, I33, NMETRIC};

// flags to denote which components of hydro and B-field arrays to flip at polar bndrys
static bool flip_across_pole_hydro[] = {false, false, true, true, false};
static bool flip_across_pole_field[] = {false, true, true};

enum edgeid {edgeid_undefined = -1, em2m1=0, em2p1=1, ep2m1=2, ep2p2=3, 
                em3m1=4, em3p1=5, ep3m1=6, ep3p1=7, em3m2=8, em3p2=9, ep3m2=10, ep3p2=11};
enum face {x1face=0, x2face=1, x3face=2};
enum mpitag {tag_hydro=0, tag_field=1, tag_flcor=2, tag_emfcor=3, tag_amr=4}; // mpitag must be < 16 and unique
enum neighbor_type {neighbor_none, neighbor_face, neighbor_edge, neighbor_corner};

enum mbtflag {mbt_node=0, mbt_refined=1, mbt_deref=2, mbt_newr=3, mbt_newd=4, mbt_leaf=5};

enum boundary_status {boundary_waiting, boundary_arrived, boundary_completed};

#endif // define ATHENA_HPP
