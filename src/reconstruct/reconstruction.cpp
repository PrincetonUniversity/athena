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
//! \fn Reconstruction::VectorDotRightEigenmatrix()
//  \brief Computes inner-product of vector and right-eigenmatrix of Roe's matrix A in the
//  primitive variables.  Only terms involving non-zero matrix elements are included to
//  improve performance.  This operation converts characteristic to primitive variables.
//
// REFERENCES:
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper. 

void Reconstruction::VectorDotRightEigenmatrix(MeshBlock *pmb, const int ivx,
  const int il, const int iu, const AthenaArray<Real> &b1, const AthenaArray<Real> &w,
  AthenaArray<Real> &vect)
{
  // permute components of output primitive vector depending on direction
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  if (MAGNETIC_FIELDS_ENABLED) {
    // Adiabatic MHD (eq. A3)
    if (NON_BAROTROPIC_EOS) {
      Real gamma = pmb->peos->GetGamma();
#pragma simd
      for (int i=il; i<=iu; ++i) {
        Real d = w(IDN,i);
        Real sqrtd = sqrt(d);

        Real btsq  = SQR(w(IBY,i)) + SQR(w(IBZ,i));
        Real ct2   = btsq/d;
        Real vaxsq = b1(i)*b1(i)/d;
        Real asq   = gamma*w(IEN,i)/d;
        Real a = sqrt(asq);

        // Compute fast- and slow-magnetosonic speeds (eq. A10)
        Real tdif = vaxsq + ct2 - asq;
        Real cf2_cs2 = sqrt(tdif*tdif + 4.0*asq*ct2);

        Real cfsq = 0.5*(vaxsq + ct2 + asq + cf2_cs2);
        Real cf = sqrt(cfsq);

        Real cssq = asq*vaxsq/cfsq;
        Real cs = sqrt(cssq);

        // Compute beta(s) (eq A17)
        Real bt  = sqrt(btsq);
        Real bet2 = 1.0;
        Real bet3 = 0.0;
        if (bt != 0.0) {
          bet2 = w(IBY,i)/bt;
          bet3 = w(IBZ,i)/bt;
        }

        // Compute alpha(s) (eq A16)
        Real alpha_f, alpha_s;
        if (cf2_cs2 == 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else if ( (asq - cssq) <= 0.0) {
          alpha_f = 0.0;
          alpha_s = 1.0;
        } else if ( (cfsq - asq) <= 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else {
          alpha_f = sqrt((asq - cssq)/cf2_cs2);
          alpha_s = sqrt((cfsq - asq)/cf2_cs2);
        }

        // Compute Q(s) and A(s) (eq. A14-15), etc.
        Real s = SIGN(b1(i));
        Real qf = cf*alpha_f*s;
        Real qs = cs*alpha_s*s;
        Real af = a*alpha_f*sqrtd;
        Real as = a*alpha_s*sqrtd;

        Real v_0 = d*(alpha_f*(vect(0,i) + vect(6,i)) + alpha_s*(vect(2,i) + vect(4,i)));
        Real v_1 = cf*alpha_f*(vect(6,i)-vect(0,i)) + cs*alpha_s*(vect(4,i)-vect(2,i));
        Real v_2 = bet2*(qs*(vect(0,i) - vect(6,i)) + qf*(vect(4,i) - vect(2,i)))
                 + bet3*(vect(5,i) - vect(1,i)); 
        Real v_3 = bet3*(qs*(vect(0,i) - vect(6,i)) + qf*(vect(4,i) - vect(2,i)))
                 + bet2*(vect(1,i) - vect(5,i)); 
        Real v_4 = d*asq*(alpha_f*(vect(0,i)+vect(6,i)) + alpha_s*(vect(2,i)+vect(4,i))); 
        Real v_5 = bet2*(as*(vect(0,i) + vect(6,i)) - af*(vect(2,i) + vect(4,i)))
                 - bet3*s*sqrtd*(vect(5,i) + vect(1,i)); 
        Real v_6 = bet3*(as*(vect(0,i) + vect(6,i)) - af*(vect(2,i) + vect(4,i)))
                 + bet2*s*sqrtd*(vect(5,i) + vect(1,i)); 

        vect(IDN,i) = v_0; 
        vect(ivx,i) = v_1; 
        vect(ivy,i) = v_2; 
        vect(ivz,i) = v_3; 
        vect(IEN,i) = v_4; 
        vect(IBY,i) = v_5; 
        vect(IBZ,i) = v_6; 
      }

    // Isothermal MHD (eq. A12)
    } else {
      Real iso_cs = pmb->peos->GetIsoSoundSpeed();
      Real iso_cs2 = SQR(iso_cs);
#pragma simd
      for (int i=il; i<=iu; ++i) {
        Real d = w(IDN,i);
        Real sqrtd = sqrt(d);

        Real btsq  = SQR(w(IBY,i)) + SQR(w(IBZ,i));
        Real ct2   = btsq/d;
        Real vaxsq = b1(i)*b1(i)/d;

        // Compute fast- and slow-magnetosonic speeds (eq. A10)
        Real tdif = vaxsq + ct2 - iso_cs2;
        Real cf2_cs2 = sqrt(tdif*tdif + 4.0*iso_cs2*ct2);

        Real cfsq = 0.5*(vaxsq + ct2 + iso_cs2 + cf2_cs2);
        Real cf = sqrt(cfsq);

        Real cssq = iso_cs2*vaxsq/cfsq;
        Real cs = sqrt(cssq);

        // Compute beta(s) (eq A17)
        Real bt  = sqrt(btsq);
        Real bet2 = 1.0;
        Real bet3 = 0.0;
        if (bt != 0.0) {
          bet2 = w(IBY,i)/bt;
          bet3 = w(IBZ,i)/bt;
        }

        // Compute alpha(s) (eq A16)
        Real alpha_f, alpha_s;
        if (cf2_cs2 == 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else if ( (iso_cs2 - cssq) <= 0.0) {
          alpha_f = 0.0;
          alpha_s = 1.0;
        } else if ( (cfsq - iso_cs2) <= 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else {
          alpha_f = sqrt((iso_cs2 - cssq)/cf2_cs2);
          alpha_s = sqrt((cfsq - iso_cs2)/cf2_cs2);
        }

        // Compute Q(s) and A(s) (eq. A14-15), etc.
        Real s = SIGN(b1(i));
        Real qf = cf*alpha_f*s;
        Real qs = cs*alpha_s*s;
        Real af = iso_cs*alpha_f*sqrtd;
        Real as = iso_cs*alpha_s*sqrtd;

        Real v_0 = d*(alpha_f*(vect(0,i) + vect(5,i)) + alpha_s*(vect(2,i) + vect(3,i)));
        Real v_1 = cf*alpha_f*(vect(5,i) - vect(0,i)) + cs*alpha_s*(vect(3,i)-vect(2,i));
        Real v_2 = bet2*(qs*(vect(0,i) - vect(5,i)) + qf*(vect(3,i) - vect(2,i)))
                 + bet3*(vect(4,i) - vect(1,i)); 
        Real v_3 = bet3*(qs*(vect(0,i) - vect(5,i)) + qf*(vect(3,i) - vect(2,i)))
                 + bet2*(vect(1,i) - vect(4,i)); 
        Real v_4 = bet2*(as*(vect(0,i) + vect(5,i)) - af*(vect(2,i) + vect(3,i)))
                 - bet3*s*sqrtd*(vect(4,i) + vect(1,i)); 
        Real v_5 = bet3*(as*(vect(0,i) + vect(5,i)) - af*(vect(2,i) + vect(3,i)))
                 + bet2*s*sqrtd*(vect(4,i) + vect(1,i)); 

        vect(IDN,i) = v_0; 
        vect(ivx,i) = v_1; 
        vect(ivy,i) = v_2; 
        vect(ivz,i) = v_3; 
        vect(IBY,i) = v_4; 
        vect(IBZ,i) = v_5; 
      }

    }

  } else {
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
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::LeftEigenmatrixDotVector()
//  \brief Computes inner-product of left-eigenmatrix of Roe's matrix A in the primitive
//  variables and vector.  Only terms involving non-zero matrix elements are included to
//  improve performance.  This operation converts primitive to characteristic variables.
//
// REFERENCES:
// - J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon "Athena: A new code for
//   astrophysical MHD", ApJS, (2008), Appendix A.  Equation numbers refer to this paper. 

void Reconstruction::LeftEigenmatrixDotVector(MeshBlock *pmb, const int ivx,
  const int il, const int iu, const AthenaArray<Real> &b1, const AthenaArray<Real> &w,
  AthenaArray<Real> &vect)
{
  // permute components of input primitive vector depending on direction
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  if (MAGNETIC_FIELDS_ENABLED) {
    // Adiabatic MHD (eq. A3)
    if (NON_BAROTROPIC_EOS) {
      Real gamma = pmb->peos->GetGamma();
#pragma simd
      for (int i=il; i<=iu; ++i) {
        Real d = w(IDN,i);
        Real id = 1.0/d;
        Real sqrtd = sqrt(d);
        Real isqrtd = 1.0/sqrtd;

        Real btsq  = SQR(w(IBY,i)) + SQR(w(IBZ,i));
        Real ct2   = btsq*id;
        Real vaxsq = b1(i)*b1(i)*id;
        Real asq   = gamma*w(IEN,i)*id;
        Real a = sqrt(asq);

        // Compute fast- and slow-magnetosonic speeds (eq. A10)
        Real tdif = vaxsq + ct2 - asq;
        Real cf2_cs2 = sqrt(tdif*tdif + 4.0*asq*ct2);

        Real cfsq = 0.5*(vaxsq + ct2 + asq + cf2_cs2);
        Real cf = sqrt(cfsq);

        Real cssq = asq*vaxsq/cfsq;
        Real cs = sqrt(cssq);

        // Compute beta(s) (eq A17)
        Real bt  = sqrt(btsq);
        Real bet2 = 1.0;
        Real bet3 = 0.0;
        if (bt != 0.0) {
          bet2 = w(IBY,i)/bt;
          bet3 = w(IBZ,i)/bt;
        }

        // Compute alpha(s) (eq A16)
        Real alpha_f, alpha_s;
        if (cf2_cs2 == 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else if ( (asq - cssq) <= 0.0) {
          alpha_f = 0.0;
          alpha_s = 1.0;
        } else if ( (cfsq - asq) <= 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else {
          alpha_f = sqrt((asq - cssq)/cf2_cs2);
          alpha_s = sqrt((cfsq - asq)/cf2_cs2);
        }

        // Compute Q(s) and A(s) (eq. A14-15), etc.
        Real s = SIGN(b1(i));
        Real na = 0.5/asq;
        Real qf = na*cf*alpha_f*s;
        Real qs = na*cs*alpha_s*s;
        Real af_prime = 0.5*alpha_f/(a*sqrtd);
        Real as_prime = 0.5*alpha_s/(a*sqrtd);

        Real v_0 = na*alpha_f*(vect(IEN,i)*id - cf*vect(ivx,i)) + 
                      qs*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) +
                as_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));
        Real v_1 = 0.5*(bet2*(vect(IBZ,i)*s*isqrtd + vect(ivz,i)) - 
                        bet3*(vect(IBY,i)*s*isqrtd + vect(ivy,i)));
        Real v_2 = na*alpha_s*(vect(IEN,i)*id - cs*vect(ivx,i)) - 
                      qf*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) -
                af_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));
        Real v_3 = vect(IDN,i) - vect(IEN,i)/asq; 
        Real v_4 = na*alpha_s*(vect(IEN,i)*id + cs*vect(ivx,i)) + 
                      qf*(bet2*vect(ivy,i) - bet3*vect(ivz,i)) +
                af_prime*(bet2*vect(IBY,i) - bet3*vect(IBZ,i));
        Real v_5 = 0.5*(bet2*(vect(IBZ,i)*s*isqrtd - vect(ivz,i)) - 
                        bet3*(vect(IBY,i)*s*isqrtd - vect(ivy,i)));
        Real v_6 = na*alpha_f*(vect(IEN,i)*id + cf*vect(ivx,i)) - 
                      qs*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) +
                as_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));

        vect(0,i) = v_0; 
        vect(1,i) = v_1; 
        vect(2,i) = v_2; 
        vect(3,i) = v_3; 
        vect(4,i) = v_4; 
        vect(5,i) = v_5; 
        vect(6,i) = v_6; 
      }

    // Isothermal MHD (eq. A3)
    } else {
      Real iso_cs = pmb->peos->GetIsoSoundSpeed();
      Real iso_cs2 = SQR(iso_cs);
#pragma simd
      for (int i=il; i<=iu; ++i) {
        Real d = w(IDN,i);
        Real id = 1.0/d;
        Real sqrtd = sqrt(d);
        Real isqrtd = 1.0/sqrtd;

        Real btsq  = SQR(w(IBY,i)) + SQR(w(IBZ,i));
        Real ct2   = btsq*id;
        Real vaxsq = b1(i)*b1(i)*id;

        // Compute fast- and slow-magnetosonic speeds (eq. A10)
        Real tdif = vaxsq + ct2 - iso_cs2;
        Real cf2_cs2 = sqrt(tdif*tdif + 4.0*iso_cs2*ct2);

        Real cfsq = 0.5*(vaxsq + ct2 + iso_cs2 + cf2_cs2);
        Real cf = sqrt(cfsq);

        Real cssq = iso_cs2*vaxsq/cfsq;
        Real cs = sqrt(cssq);

        // Compute beta(s) (eq A17)
        Real bt  = sqrt(btsq);
        Real bet2 = 1.0;
        Real bet3 = 0.0;
        if (bt != 0.0) {
          bet2 = w(IBY,i)/bt;
          bet3 = w(IBZ,i)/bt;
        }

        // Compute alpha(s) (eq A16)
        Real alpha_f, alpha_s;
        if (cf2_cs2 == 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else if ( (iso_cs2 - cssq) <= 0.0) {
          alpha_f = 0.0;
          alpha_s = 1.0;
        } else if ( (cfsq - iso_cs2) <= 0.0) {
          alpha_f = 1.0;
          alpha_s = 0.0;
        } else {
          alpha_f = sqrt((iso_cs2 - cssq)/cf2_cs2);
          alpha_s = sqrt((cfsq - iso_cs2)/cf2_cs2);
        }

        // Compute Q(s) and A(s) (eq. A14-15), etc.
        Real s = SIGN(b1(i));
        Real qf = 0.5*(cf*alpha_f*s)/iso_cs2;
        Real qs = 0.5*(cs*alpha_s*s)/iso_cs2;
        Real af_prime = 0.5*alpha_f/(iso_cs*sqrtd);
        Real as_prime = 0.5*alpha_s/(iso_cs*sqrtd);

        Real v_0 = 0.5*alpha_f*(vect(IDN,i)*id - cf*vect(ivx,i)/iso_cs2) + 
                      qs*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) +
                as_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));
        Real v_1 = 0.5*(bet2*(vect(IBZ,i)*s*isqrtd + vect(ivz,i)) - 
                        bet3*(vect(IBY,i)*s*isqrtd + vect(ivy,i)));
        Real v_2 = 0.5*alpha_s*(vect(IDN,i)*id - cs*vect(ivx,i)/iso_cs2) - 
                      qf*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) -
                af_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));
        Real v_3 = 0.5*alpha_s*(vect(IDN,i)*id - cs*vect(ivx,i)/iso_cs2) - 
                      qf*(bet2*vect(ivy,i) - bet3*vect(ivz,i)) +
                af_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));
        Real v_4 = 0.5*(bet2*(vect(IBZ,i)*s*isqrtd - vect(ivz,i)) - 
                        bet3*(vect(IBY,i)*s*isqrtd - vect(ivy,i)));
        Real v_5 = 0.5*alpha_f*(vect(IDN,i)*id - cf*vect(ivx,i)/iso_cs2) - 
                      qs*(bet2*vect(ivy,i) + bet3*vect(ivz,i)) +
                as_prime*(bet2*vect(IBY,i) + bet3*vect(IBZ,i));

        vect(0,i) = v_0; 
        vect(1,i) = v_1; 
        vect(2,i) = v_2; 
        vect(3,i) = v_3; 
        vect(4,i) = v_4; 
        vect(5,i) = v_5; 
      }
    }

  } else {
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
  }
  return;
}

