#ifndef CONVERT_VARIABLES_HPP
#define CONVERT_VARIABLES_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file convert_var.hpp
 *  \brief defines class ConvertVariables
 *  Handles conversion from conserved to primitive variables
 *====================================================================================*/

class Fluid;

//! \class ConvertVariables
//  \brief variable conversion data and functions

class ConvertVariables {
public:
  ConvertVariables(Fluid *pf);
  ~ConvertVariables();

// converts conserved to primitive variables
  void ComputePrimitives(AthenaArray<Real> &c, AthenaArray<Real> &p);

private:
  Fluid *pparent_fluid_;  // ptr to parent Fluid

};
#endif
