#ifndef ATHENA_HPP
#define ATHENA_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file athena.hpp
 *  \brief contains Athena++ specific types, structures, macros, etc.
 *====================================================================================*/

#define NGHOST 2
#define NVAR 5
#define PI 3.14159265358979323846
#define TINY_NUMBER 1.0e-20
#define HUGE_NUMBER 1.0e+36

#define COORDINATE_SYSTEM cartesian_coordinates
#define RIEMANN_SOLVER hllc_hydro_solver
#define RECONSTRUCTION_ALGORITHM piecewise_linear_method

typedef double Real;
enum {IDN=0, IM1=1, IM2=2, IM3=3, IEN=4};
enum {IVX=1, IVY=2, IVZ=3};
enum UpdateAction {fluid_predict, fluid_correct,   bfield_predict, bfield_correct,
                   fluid_bcs_n,   fluid_bcs_nhalf, bfield_bcs_n,   bfield_bcs_nhalf,
                   convert_vars_n, convert_vars_nhalf, new_timestep, data_output};
enum QuantityToBeInit {bfield_bcs, fluid_bcs, geometry, initial_conditions};

#endif
