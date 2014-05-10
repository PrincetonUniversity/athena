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

typedef double Real;
enum {IDN=0, IM1=1, IM2=2, IM3=3, IEN=4};
enum {IVX=1, IVY=2, IVZ=3};
enum AlgorithmSteps {fluid_predict, fluid_correct, bfield_predict, bfield_correct,
                     fluid_bvals_n, fluid_bvals_nhalf, bfield_bvals_n, bfield_bvals_nhalf,
                     convert_vars_n, convert_vars_nhalf, data_output};
enum QuantityToBeInitialized {fluid};

#endif
