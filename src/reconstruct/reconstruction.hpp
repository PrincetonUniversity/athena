#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP
//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 * See LICENSE file for full public license information.
 *====================================================================================*/
/*! \file reconstruction.hpp
 *  \brief defines class Reconstruction
 *  implements algorithms for reconstruction in the primitive variables
 *====================================================================================*/

class Fluid;

//! \class Reconstruction
//  \brief reconstruction data and functions

class Reconstruction {
public:
  Reconstruction(Fluid *pf);
  ~Reconstruction();

  void PiecewiseLinear(const int k,const int j,const int il,const int iu,const int dir,
      AthenaArray<Real> &w, AthenaArray<Real> &wl, AthenaArray<Real> &wr);

private:
  Fluid *pmy_fluid_;  // pointer to parent Fluid object

};
#endif
