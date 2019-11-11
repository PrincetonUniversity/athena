//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adm_z4c.cpp
//  \brief implementation of functions in the Z4c class related to ADM decomposition

// C++ standard headers
#include <cmath> // pow
#include <iostream>
#include <fstream>

// Athena++ headers
#include "z4c.hpp"
#include "z4c_macro.hpp"
#include "../mesh/mesh.hpp"

//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMToZ4c(AthenaArray<Real> & u_adm, AthenaArray<Real> & u)
// \brief Compute Z4c variables from ADM variables
//
// p  = detgbar^(-1/3)
// p0 = psi^(-4)
//
// gtilde_ij = p gbar_ij
// Ktilde_ij = p p0 K_ij
//
// phi = - log(p) / 4
// K   = gtildeinv^ij Ktilde_ij
// Atilde_ij = Ktilde_ij - gtilde_ij K / 3
//
// G^i = - del_j gtildeinv^ji
//
// BAM: Z4c_init()
// https://git.tpi.uni-jena.de/bamdev/z4
// https://git.tpi.uni-jena.de/bamdev/z4/blob/master/z4_init.m
//
// The Z4c variables will be set on the whole MeshBlock with the exception of
// the Gamma's that can only be set in the interior of the MeshBlock.

void Z4c::ADMToZ4c(AthenaArray<Real> & u_adm, AthenaArray<Real> & u)
{
  ADM_vars adm;
  SetADMAliases(u_adm, adm);
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);

  //--------------------------------------------------------------------------------------
  // Conformal factor, conformal metric, and trace of extrinsic curvature
  //
  GLOOP2(k,j) {
    
    // Conformal factor
    GLOOP1(i) {
      detg(i)   = SpatialDet(adm.g_dd, k, j, i);
      oopsi4(i) = pow(detg(i), -1./3.);
      z4c.chi(k,j,i) = pow(detg(i), 1./12.*opt.chi_psi_power);
    }
    
    // Conformal metric and extrinsic curvature
    for(int a = 0; a < NDIM; ++a)
    for(int b = a; b < NDIM; ++b) {
      GLOOP1(i) {
        z4c.g_dd(a,b,k,j,i) = oopsi4(i) * adm.g_dd(a,b,k,j,i);
        Kt_dd(a,b,i)        = oopsi4(i) * adm.K_dd(a,b,k,j,i);
      }
    }

    // Determinant of the conformal metric and trace of conf. extr. curvature
    GLOOP1(i) {
      detg(i) = SpatialDet(z4c.g_dd, k, j, i);
      z4c.Khat(k,j,i) = Trace(1.0/detg(i),
          z4c.g_dd(0,0,k,j,i), z4c.g_dd(0,1,k,j,i), z4c.g_dd(0,2,k,j,i),
          z4c.g_dd(1,1,k,j,i), z4c.g_dd(1,2,k,j,i), z4c.g_dd(2,2,k,j,i),
          Kt_dd(0,0,i), Kt_dd(0,1,i), Kt_dd(0,2,i),
          Kt_dd(1,1,i), Kt_dd(1,2,i), Kt_dd(2,2,i));
    }
    
    // Conformal traceless extrinsic curvatore
    for(int a = 0; a < NDIM; ++a)
    for(int b = a; b < NDIM; ++b) {
      GLOOP1(i) {
          z4c.A_dd(a,b,k,j,i) = Kt_dd(a,b,i) - (1./3.) * z4c.Khat(k,j,i) * z4c.g_dd(a,b,k,j,i);
      }
    }
  }

  //--------------------------------------------------------------------------------------
  // Gamma's
  //
  // Allocate temporary memory for the inverse conformal metric
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if(pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if(pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> g_uu;
  g_uu.NewAthenaTensor(ncells3, ncells2, ncells1);

  // Inverse conformal metric
  GLOOP3(k,j,i) {
    detg(i) = SpatialDet(z4c.g_dd, k, j, i);
    SpatialInv(1.0/detg(i),
        z4c.g_dd(0,0,k,j,i), z4c.g_dd(0,1,k,j,i), z4c.g_dd(0,2,k,j,i),
        z4c.g_dd(1,1,k,j,i), z4c.g_dd(1,2,k,j,i), z4c.g_dd(2,2,k,j,i),
        &g_uu(0,0,k,j,i),    &g_uu(0,1,k,j,i),    &g_uu(0,2,k,j,i),
        &g_uu(1,1,k,j,i),    &g_uu(1,2,k,j,i),    &g_uu(2,2,k,j,i));
  }

  // Compute Gamma's
  z4c.Gam_u.ZeroClear();
  ILOOP2(k,j) {
    for(int a = 0; a < NDIM; ++a)
    for(int b = 0; b < NDIM; ++b) {
      ILOOP1(i) {
        z4c.Gam_u(a,k,j,i) -= FD.Dx(b, g_uu(b,a,k,j,i)); // Is it ba or ab like in the pseudocode? Is the contraction correct?
      }
    }
  }

  g_uu.DeleteAthenaTensor();

  //--------------------------------------------------------------------------------------
  // Theta
  //
  z4c.Theta.ZeroClear();

  //--------------------------------------------------------------------------------------
  // Algebraic constraints enforcement
  //
  AlgConstr(u);

}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::Z4cToADM(AthenaArray<Real> & u, AthenaArray<Real> & u_adm)
// \brief Compute ADM Psi4, g_ij, and K_ij from Z4c variables
//
// This sets the ADM variables everywhere in the MeshBlock

void Z4c::Z4cToADM(AthenaArray<Real> & u, AthenaArray<Real> & u_adm)
{
  ADM_vars adm;
  SetADMAliases(u_adm, adm);
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);

  GLOOP2(k,j) {
    // psi4
    GLOOP1(i) {
      adm.psi4(k,j,i) = std::pow(z4c.chi(k,j,i), 4./opt.chi_psi_power);
    }
    // g_ab
    for(int a = 0; a < NDIM; ++a)
    for(int b = a; b < NDIM; ++b) {
      GLOOP1(i) {
        adm.g_dd(a,b,k,j,i) = adm.psi4(k,j,i) * z4c.g_dd(a,b,k,j,i);
      }
    }
    // K_ab
    for(int a = 0; a < NDIM; ++a)
    for(int b = a; b < NDIM; ++b) {
      GLOOP1(i) {
        adm.K_dd(a,b,k,j,i) = adm.psi4(k,j,i) * z4c.A_dd(a,b,k,j,i) +
          (1./3.) * (z4c.Khat(k,j,i) + 2.*z4c.Theta(k,j,i)) * adm.g_dd(a,b,k,j,i);
      }
    }
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMConstraints(AthenaArray<Real> & u_adm, AthenaArray<Real> & u_mat)
// \brief compute constraints ADM vars
//
// Note: we are assuming that u_adm has been initialized with the correct
// metric and matter quantities
//
// BAM: adm_constraints_N()
// https://git.tpi.uni-jena.de/bamdev/adm
// https://git.tpi.uni-jena.de/bamdev/adm/blob/master/adm_constraints_N.m
//
// The constraints are set only in the MeshBlock interior, because derivatives
// of the ADM quantities are neded to compute them.

void Z4c::ADMConstraints(AthenaArray<Real> & u_con, AthenaArray<Real> & u_adm,
                         AthenaArray<Real> & u_mat, AthenaArray<Real> & u_z4c)
{
  u_con.ZeroClear();

  Constraint_vars con;
  SetConstraintAliases(u_con, con);

  ADM_vars adm;
  SetADMAliases(u_adm, adm);

  Matter_vars mat;
  SetMatterAliases(u_mat, mat);

  Z4c_vars z4c;
  SetZ4cAliases(u_z4c, z4c);

  ILOOP2(k,j) {
    // -----------------------------------------------------------------------------------
    // derivatives
    //
    // first derivatives of g and K
    for(int c = 0; c < NDIM; ++c)
    for(int a = 0; a < NDIM; ++a)
    for(int b = a; b < NDIM; ++b) {
      ILOOP1(i) {
        dg_ddd(c,a,b,i) = FD.Dx(c, adm.g_dd(a,b,k,j,i));
        dK_ddd(c,a,b,i) = FD.Dx(c, adm.K_dd(a,b,k,j,i));
      }
    }
    // second derivatives of g
    for(int a = 0; a < NDIM; ++a)
    for(int b = a; b < NDIM; ++b)
    for(int c = 0; c < NDIM; ++c)
    for(int d = c; d < NDIM; ++d) {
      if(a == b) {
        ILOOP1(i) {
          ddg_dddd(a,a,c,d,i) = FD.Dxx(a, adm.g_dd(c,d,k,j,i));
        }
      }
      else {
        ILOOP1(i) {
          ddg_dddd(a,b,c,d,i) = FD.Dxy(a, b, adm.g_dd(c,d,k,j,i));
        }
      }
    }

    // -----------------------------------------------------------------------------------
    // inverse metric
    //
    ILOOP1(i) {
      detg(i) = SpatialDet(adm.g_dd,k,j,i);
      SpatialInv(1./detg(i),
          adm.g_dd(0,0,k,j,i), adm.g_dd(0,1,k,j,i), adm.g_dd(0,2,k,j,i),
          adm.g_dd(1,1,k,j,i), adm.g_dd(1,2,k,j,i), adm.g_dd(2,2,k,j,i),
          &g_uu(0,0,i), &g_uu(0,1,i), &g_uu(0,2,i),
          &g_uu(1,1,i), &g_uu(1,2,i), &g_uu(2,2,i));
    }

    // -----------------------------------------------------------------------------------
    // Christoffel symbols
    //
    for(int c = 0; c < NDIM; ++c)
    for(int a = 0; a < NDIM; ++a)
    for(int b = a; b < NDIM; ++b) {
      ILOOP1(i) {
        Gamma_ddd(c,a,b,i) = 0.5*(dg_ddd(a,b,c,i) + dg_ddd(b,a,c,i) - dg_ddd(c,a,b,i));
      }
    }

    Gamma_udd.ZeroClear();
    for(int c = 0; c < NDIM; ++c)
    for(int a = 0; a < NDIM; ++a)
    for(int b = a; b < NDIM; ++b)
    for(int d = 0; d < NDIM; ++d) {
      ILOOP1(i) {
        Gamma_udd(c,a,b,i) += g_uu(c,d,i)*Gamma_ddd(d,a,b,i);
      }
    }

    Gamma_u.ZeroClear();
    for(int a = 0; a < NDIM; ++a)
    for(int b = 0; b < NDIM; ++b)
    for(int c = 0; c < NDIM; ++c) {
      ILOOP1(i) {
        Gamma_u(a,i) += g_uu(b,c,i)*Gamma_udd(a,b,c,i);
      }
    }

    // -----------------------------------------------------------------------------------
    // Ricci tensor and Ricci scalar
    //
    R.ZeroClear();
    R_dd.ZeroClear();
    for(int a = 0; a < NDIM; ++a)
    for(int b = a; b < NDIM; ++b) {
      for(int c = 0; c < NDIM; ++c)
      for(int d = 0; d < NDIM; ++d) {
        // Part with the Christoffel symbols
        for(int e = 0; e < NDIM; ++e) {
          ILOOP1(i) {
            R_dd(a,b,i) += g_uu(c,d,i) * Gamma_udd(e,a,c,i) * Gamma_ddd(e,b,d,i);
            R_dd(a,b,i) -= g_uu(c,d,i) * Gamma_udd(e,a,b,i) * Gamma_ddd(e,c,d,i);
          }
        }
        // Wave operator part of the Ricci
        ILOOP1(i) {
          R_dd(a,b,i) += 0.5*g_uu(c,d,i)*(
              - ddg_dddd(c,d,a,b,i) - ddg_dddd(a,b,c,d,i) +
                ddg_dddd(a,c,b,d,i) + ddg_dddd(b,c,a,d,i));
        }
      }
      ILOOP1(i) {
        R(i) += g_uu(a,b,i) * R_dd(a,b,i);
      }
    }

    // -----------------------------------------------------------------------------------
    // Extrinsic curvature: traces and derivatives
    //
    K.ZeroClear();
    K_ud.ZeroClear();
    for(int a = 0; a < NDIM; ++a) {
      for(int b = a; b < NDIM; ++b) {
        for(int c = 0; c < NDIM; ++c) {
          ILOOP1(i) {
            K_ud(a,b,i) += g_uu(a,c,i) * adm.K_dd(c,b,k,j,i);
          }
        }
      }
      ILOOP1(i) {
        K(i) += K_ud(a,a,i);
      }
    }
    // K^a_b K^b_a
    KK.ZeroClear();
    for(int a = 0; a < NDIM; ++a)
    for(int b = 0; b < NDIM; ++b) {
      ILOOP1(i) {
        KK(i) += K_ud(a,b,i) * K_ud(b,a,i);
      }
    }
    // Covariant derivative of K
    for(int a = 0; a < NDIM; ++a)
    for(int b = 0; b < NDIM; ++b)
    for(int c = b; c < NDIM; ++c) {
      ILOOP1(i) {
        DK_ddd(a,b,c,i) = dK_ddd(a,b,c,i);
      }
      for(int d = 0; d < NDIM; ++d) {
        ILOOP1(i) {
          DK_ddd(a,b,c,i) -= Gamma_udd(d,a,b,i) * adm.K_dd(d,c,k,j,i);
          DK_ddd(a,b,c,i) -= Gamma_udd(d,a,c,i) * adm.K_dd(b,d,k,j,i);
        }
      }
    }
    DK_udd.ZeroClear();
    for(int a = 0; a < NDIM; ++a)
    for(int b = 0; b < NDIM; ++b)
    for(int c = b; c < NDIM; ++c)
    for(int d = 0; d < NDIM; ++d) {
      ILOOP1(i) {
        DK_udd(a,b,c,i) += g_uu(a,d,i) * DK_ddd(d,b,c,i);
      }
    }

    // -----------------------------------------------------------------------------------
    // Actual constraints
    //
    // Hamiltonian constraint
    //
    ILOOP1(i) {
      con.H(k,j,i) = R(i) + SQR(K(i)) - KK(i) - 16*M_PI * mat.rho(k,j,i);
    }
    // Momentum constraint (contravariant)
    //
    M_u.ZeroClear();
    for(int a = 0; a < NDIM; ++a)
    for(int b = 0; b < NDIM; ++b) {
      ILOOP1(i) {
        M_u(a,i) -= 8*M_PI * g_uu(a,b,i) * mat.S_d(b,k,j,i);
      }
      for(int c = 0; c < NDIM; ++c) {
        ILOOP1(i) {
          M_u(a,i) += g_uu(a,b,i) * DK_udd(c,b,c,i);
          M_u(a,i) -= g_uu(b,c,i) * DK_udd(a,b,c,i);
        }
      }
    }
    // Momentum constraint (covariant)
    for(int a = 0; a < NDIM; ++a)
    for(int b = 0; b < NDIM; ++b) {
      ILOOP1(i) {
        con.M_d(a,k,j,i) += adm.g_dd(a,b,k,j,i) * M_u(b,i);
      }
    }
    // Momentum constraint (norm squared)
    for(int a = 0; a < NDIM; ++a)
    for(int b = 0; b < NDIM; ++b) {
      ILOOP1(i) {
        con.M(k,j,i) += adm.g_dd(a,b,k,j,i) * M_u(a,i) * M_u(b,i);
      }
    }
    // Constraint violation Z (norm squared)
    for(int a = 0; a < NDIM; ++a)
    for(int b = 0; b < NDIM; ++b) {
      ILOOP1(i) {
        con.Z(k,j,i) += 0.25*adm.g_dd(a,b,k,j,i)*(z4c.Gam_u(a,k,j,i) - Gamma_u(a,i))
                                                *(z4c.Gam_u(b,k,j,i) - Gamma_u(b,i));
      }
    }
    // Constraint violation monitor C^2
    ILOOP1(i) {
      con.C(k,j,i) = SQR(con.H(k,j,i)) + con.M(k,j,i) + SQR(z4c.Theta(k,j,i)) + 4.0*con.Z(k,j,i);
    }
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::ADMMinkowski(AthenaArray<Real> & u)
// \brief Initialize ADM vars to Minkowski

void Z4c::ADMMinkowski(AthenaArray<Real> & u_adm)
{
  ADM_vars adm;
  SetADMAliases(u_adm, adm);
  adm.psi4.Fill(1.);
  adm.K_dd.ZeroClear();

  GLOOP3(k,j,i) {
    for(int a = 0; a < NDIM; ++a)
    for(int b = a; b < NDIM; ++b) {
      adm.g_dd(a,b,k,j,i) = (a == b ? 1. : 0.);
    }
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::GaugeGeodesic(AthenaArray<Real> & u)
// \brief Initialize lapse to 1 and shift to 0

void Z4c::GaugeGeodesic(AthenaArray<Real> & u)
{
  Z4c_vars z4c;
  SetZ4cAliases(u, z4c);
  z4c.alpha.Fill(1.);
  z4c.beta_u.ZeroClear();
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::MatterVacuum(AthenaArray<Real> & u_mat)
// \brief Initialize ADM vars to vacuum

void Z4c::MatterVacuum(AthenaArray<Real> & u_mat)
{
  u_mat.ZeroClear();
}
