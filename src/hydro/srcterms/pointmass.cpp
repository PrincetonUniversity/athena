//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief Adds source terms due to point mass AT ORIGIN
// C/C++ headers
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "hydro_srcterms.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../hydro.hpp"

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::PointMass
//  \brief Adds source terms due to point mass AT ORIGIN

void HydroSourceTerms::PointMass(const Real dt, const AthenaArray<Real> *flux,
  const AthenaArray<Real> &prim, AthenaArray<Real> &cons)
{

  MeshBlock *pmb = pmy_hydro_->pmy_block;
//  for (int k=pmb->ks; k<=pmb->ke; ++k) {
//#pragma omp parallel for schedule(static)
//    for (int j=pmb->js; j<=pmb->je; ++j) {
//#pragma simd
//      for (int i=pmb->is; i<=pmb->ie; ++i) {
//        Real den = prim(IDN,k,j,i);
//        Real src = dt*den*pmb->pcoord->coord_src1_i_(i)*gm_/pmb->pcoord->x1v(i);
//        cons(IM1,k,j,i) -= src;
//        if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) -=
//          dt*0.5*(pmb->pcoord->phy_src1_i_(i)*flux[X1DIR](IDN,k,j,i)*gm_
//                 +pmb->pcoord->phy_src2_i_(i)*flux[X1DIR](IDN,k,j,i+1)*gm_);
//      }
//    }
//  }
//[JMSHI
  if (COORDINATE_SYSTEM == "cylindrical" ||
      COORDINATE_SYSTEM == "spherical_polar") {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real den = prim(IDN,k,j,i);
          Real src = dt*den*pmb->pcoord->coord_src1_i_(i)*gm_/pmb->pcoord->x1v(i);
          cons(IM1,k,j,i) -= src;
          if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) -=
            dt*0.5*(pmb->pcoord->phy_src1_i_(i)*flux[X1DIR](IDN,k,j,i)*gm_
                   +pmb->pcoord->phy_src2_i_(i)*flux[X1DIR](IDN,k,j,i+1)*gm_);
        }
    }}
  } else if (COORDINATE_SYSTEM == "cartesian") {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real den = prim(IDN,k,j,i);
          Real x1  = pmb->pcoord->x1v(i);
          Real x2  = pmb->pcoord->x2v(j);
          Real x3  = pmb->pcoord->x3v(k);
          Real rad = sqrt(SQR(x1)+SQR(x2)+SQR(x3)+SQR(rsoft_));
          Real dpt = gm_/SQR(rad)/rad;
          Real src = dt*den*dpt;
          cons(IM1,k,j,i) -= src*x1;
          cons(IM2,k,j,i) -= src*x2;
          cons(IM3,k,j,i) -= src*x3;
          if (NON_BAROTROPIC_EOS) {
            cons(IEN,k,j,i) -= 0.5*dt*dpt*x1*(flux[X1DIR](IDN,k,j,i)+
                                              flux[X1DIR](IDN,k,j,i+1));
            if (pmb->block_size.nx2 > 1)
              cons(IEN,k,j,i) -= 0.5*dt*dpt*x2*(flux[X2DIR](IDN,k,j,i)+
                                                flux[X2DIR](IDN,k,j+1,i));
            if (pmb->block_size.nx3 > 1)
              cons(IEN,k,j,i) -= 0.5*dt*dpt*x3*(flux[X3DIR](IDN,k,j,i)+
                                                flux[X3DIR](IDN,k+1,j,i));
          }
          //Real src = dt*den*gm_/SQR(rad)/rad;
          //cons(IM1,k,j,i) -= src*x1;
          //cons(IM2,k,j,i) -= src*x2;
          //cons(IM3,k,j,i) -= src*x3;
          //if (NON_BAROTROPIC_EOS) cons(IEN,k,j,i) -=
          //  src*(x1*prim(IVX,k,j,i)+x2*prim(IVY,k,j,i)+x3*prim(IVZ,k,j,i));
    }}}
  } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in HydroSourceTerms::PointMass" << std::endl
          << "The point source gravity works only in spherical polar coordinates "
          << ", cylindrical coordinates or cartesian coordinates." << std::endl;
      throw std::runtime_error(msg.str().c_str());

  }
  //JMSHI]
  return;
}
