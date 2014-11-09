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
#include <math.h>

typedef double Real;

typedef struct ThreeVector {
  Real x1, x2, x3;
} ThreeVector;

enum {IDN=0, IM1=1, IM2=2, IM3=3, IEN=4};
enum {IVX=1, IVY=2, IVZ=3, IBY=(NFLUID), IBZ=(NFLUID+1)};

enum {I00, I01, I02, I03, I11, I12, I13, I22, I23, I33, NMETRIC};

enum ActionOnDomain 
  {pgen,          primitives_n, primitives_nhalf, new_timestep,
   fluid_predict, fluid_correct,   bfield_predict, bfield_correct,
   fluid_bcs_n,   fluid_bcs_nhalf, bfield_bcs_n,   bfield_bcs_nhalf};

#endif
