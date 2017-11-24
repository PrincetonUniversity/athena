//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file reconstruction.cpp
//  \brief

// C/C++ headers
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"

// constructor

Reconstruction::Reconstruction(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block_ = pmb;

  // read and set type of spatial reconstruction
  characteristic_reconstruction = false;
  std::string input_recon = pin->GetOrAddString("time","xorder","2");
  if (input_recon == "1") {
    xorder = 1;
  } else if (input_recon == "2") {
    xorder = 2;
  } else if (input_recon == "2c") {
    xorder = 2;
    characteristic_reconstruction = true;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in Reconstruction constructor" << std::endl
        << "xorder=" << input_recon << " not valid choice for reconstruction"<< std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // set function pointers for reconstruction functions in each direction
  // First-order (donor cell) reconstruction
  if (xorder == 1) {
    ReconstructFuncX1 = DonorCellX1;
    ReconstructFuncX2 = DonorCellX2;
    ReconstructFuncX3 = DonorCellX3;

  // Second-order (piecewise linear) reconstruction
  } else if (xorder == 2) {
      ReconstructFuncX1 = PiecewiseLinearX1;
      ReconstructFuncX2 = PiecewiseLinearX2;
      ReconstructFuncX3 = PiecewiseLinearX3;

  // Third/Fourth-order (piecewise parabolic) reconstruction
  } else if (xorder == 3) {
    if (pmb->block_size.x1rat == 1.0) {
      ReconstructFuncX1 = PPMUniformX1;
    } else {
      ReconstructFuncX1 = PPMX1;
    }

    if (pmb->block_size.x2rat == 1.0) {
      ReconstructFuncX2 = PPMUniformX2;
    } else {
      ReconstructFuncX2 = PPMX2;
    }

    if (pmb->block_size.x3rat == 1.0) {
      ReconstructFuncX3 = PPMUniformX3;
    } else {
      ReconstructFuncX3 = PPMX3;
    }
  // Error; unknown order
  } else {
    std:: stringstream msg;
    msg << "### FATAL ERROR in function [Reconstruction constructor]" << std::endl
        << "spatial order xorder= " << xorder << " not supported" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Allocate memory for scratch arrays
//  int ncells1 = ((pmy_block->ie)-(pmy_block->is) + 1) + 2*(NGHOST);
//  dwl_.NewAthenaArray(NWAVE,ncells1);
//  dwr_.NewAthenaArray(NWAVE,ncells1);
//  dw2_.NewAthenaArray(NWAVE,ncells1);
//  dwm_.NewAthenaArray(NWAVE,ncells1);
//  wc_.NewAthenaArray(NWAVE,ncells1);
}

// destructor

Reconstruction::~Reconstruction()
{
//  dwl_.DeleteAthenaArray();
//  dwr_.DeleteAthenaArray();
//  dw2_.DeleteAthenaArray();
//  dwm_.DeleteAthenaArray();
//  wc_.DeleteAthenaArray();
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::VectorRightEigenmatrixProduct()
//  \brief Computes inner-product of vector and right-eigenmatrix of Roe's matrix A in the
//  primitive variables.  Only terms involving non-zero matrix elements are included to
//  improve performance.  This operation converts characteristic to primitive variables.
//
// REFERENCES:
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper. 

void Reconstruction::VectorRightEigenmatrixProduct(MeshBlock *pmb, const int ivx,
  const AthenaArray<Real> &w, const int il, const int iu, AthenaArray<Real> &vect)
{
  // permute components of output primitive vector depending on direction
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Adiabatic hydrodynamics (eq. A3)
  if (NON_BAROTROPIC_EOS) {
    Real gamma = pmb->peos->GetGamma();
#pragma simd
    for (int i=il; i<=iu; ++i) {
      Real asq = gamma*w(IEN,i)/w(IDN,i);
      Real a   = sqrt(asq);

      Real v_0 = vect(0,i) + vect(1,i) + vect(4,i);
      Real v_1 = a*(vect(4,i) - vect(0,i))/w(IDN,i);
      Real v_2 = vect(2,i);
      Real v_3 = vect(3,i);
      Real v_4 = asq*(vect(0,i) + vect(4,i));

      vect(IDN,i) = v_0; 
      vect(ivx,i) = v_1; 
      vect(ivy,i) = v_2; 
      vect(ivz,i) = v_3; 
      vect(IEN,i) = v_4; 
    }

  // Isothermal hydrodynamics (eq. A3)
  } else {
    Real iso_cs = pmb->peos->GetIsoSoundSpeed();
#pragma simd
    for (int i=il; i<=iu; ++i) {
      Real v_0 = vect(0,i) + vect(3,i);
      Real v_1 = iso_cs*(vect(3,i) - vect(0,i))/w(IDN,i);
      Real v_2 = vect(1,i);
      Real v_3 = vect(2,i);

      vect(IDN,i) = v_0; 
      vect(ivx,i) = v_1; 
      vect(ivy,i) = v_2; 
      vect(ivz,i) = v_3; 
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::LeftEigenmatrixVectorProduct()
//  \brief Computes inner-product of left-eigenmatrix of Roe's matrix A in the primitive
//  variables and vector.  Only terms involving non-zero matrix elements are included to
//  improve performance.  This operation converts primitive to characteristic variables.
//
// REFERENCES:
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper. 

void Reconstruction::LeftEigenmatrixVectorProduct(MeshBlock *pmb, const int ivx,
  const AthenaArray<Real> &w, const int il, const int iu, AthenaArray<Real> &vect)
{
  // permute components of input primitive vector depending on direction
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Adiabatic hydrodynamics (eq. A4)
  if (NON_BAROTROPIC_EOS) {
    Real gamma = pmb->peos->GetGamma();
#pragma simd
    for (int i=il; i<=iu; ++i) {
      Real asq = gamma*w(IEN,i)/w(IDN,i);
      Real a   = sqrt(asq);

      Real v_0 = 0.5*(vect(4,i)/asq - w(IDN,i)*vect(ivx,i)/a);
      Real v_1 = vect(0,i) - vect(4,i)/asq;
      Real v_2 = vect(ivy,i);
      Real v_3 = vect(ivz,i);
      Real v_4 = 0.5*(vect(4,i)/asq + w(IDN,i)*vect(ivx,i)/a);

      vect(0,i) = v_0; 
      vect(1,i) = v_1; 
      vect(2,i) = v_2; 
      vect(3,i) = v_3; 
      vect(4,i) = v_4; 
    }

  // Isothermal hydrodynamics (eq. A7)
  } else {
    Real iso_cs = pmb->peos->GetIsoSoundSpeed();
#pragma simd
    for (int i=il; i<=iu; ++i) {
      Real v_0 = 0.5*(vect(0,i) - w(IDN,i)*vect(ivx,i)/iso_cs);
      Real v_1 = vect(ivy,i);
      Real v_2 = vect(ivz,i);
      Real v_3 = 0.5*(vect(0,i) + w(IDN,i)*vect(ivx,i)/iso_cs);

      vect(0,i) = v_0; 
      vect(1,i) = v_1; 
      vect(2,i) = v_2; 
      vect(3,i) = v_3; 
    }
  }
  return;
}

