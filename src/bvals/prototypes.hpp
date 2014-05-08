#ifndef BVAL_PROTOTYPES_HPP
#define BVAL_PROTOTYPES_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file prototypes.hpp
 *  \brief prototypes of all BC functions
 *====================================================================================*/

// reflecting BC functions

  void ReflectInnerX1(Fluid *pf);
  void ReflectOuterX1(Fluid *pf);
  void ReflectInnerX2(Fluid *pf);
  void ReflectOuterX2(Fluid *pf);
  void ReflectInnerX3(Fluid *pf);
  void ReflectOuterX3(Fluid *pf);

#endif
