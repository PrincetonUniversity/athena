//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file shock_tube.cpp
//! \brief Problem generator for shock tube problems.
//!
//! Problem generator for shock tube (1-D Riemann) problems. Initializes plane-parallel
//! shock along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <cstdio>     // fopen(), freopen(), fprintf(), fclose()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

//========================================================================================
//! \fn Real press(Real rho, Real T)
//! \brief Calculate pressure as a function of density and temperature for H EOS.
//========================================================================================

Real press(Real rho, Real T) {
  // Ionization fraction
  Real x = 2. /(1 + std::sqrt(1 + 4. * rho * std::exp(1. / T) * std::pow(T, -1.5)));
  return rho * T * (1. + x);
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//! \brief Calculate L1 errors in Sod (hydro) and RJ2a (MHD) tests
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  MeshBlock *pmb = my_blocks(0);

  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;

  // Read shock direction and set array indices
  int shk_dir = pin->GetInteger("problem","shock_dir");
  int im1,im2,im3,ib1,ib2,ib3;
  if (shk_dir == 1) {
    im1 = IM1; im2 = IM2; im3 = IM3;
    ib1 = IB1; ib2 = IB2; ib3 = IB3;
  } else if (shk_dir == 2) {
    im1 = IM2; im2 = IM3; im3 = IM1;
    ib1 = IB2; ib2 = IB3; ib3 = IB1;
  } else {
    im1 = IM3; im2 = IM1; im3 = IM2;
    ib1 = IB3; ib2 = IB1; ib3 = IB2;
  }

  // Initialize errors to zero
  Real err[NHYDRO+NFIELD];
  for (int i=0; i<(NHYDRO+NFIELD); ++i) err[i]=0.0;

  // Errors in RJ2a test (Dai & Woodward 1994 Tables Ia and Ib)
  if (MAGNETIC_FIELDS_ENABLED) {
    Real xfp = 2.2638*tlim;
    Real xrp = (0.53432 + 1.0/std::sqrt(PI*1.309))*tlim;
    Real xsp = (0.53432 + 0.48144/1.309)*tlim;
    Real xc = 0.57538*tlim;
    Real xsm = (0.60588 - 0.51594/1.4903)*tlim;
    Real xrm = (0.60588 - 1.0/std::sqrt(PI*1.4903))*tlim;
    Real xfm = (1.2 - 2.3305/1.08)*tlim;
    Real gm1;
    if (!GENERAL_EOS) gm1 = pmb->peos->GetGamma() - 1.0;
    for (int k=pmb->ks; k<=pmb->ke; k++) {
      for (int j=pmb->js; j<=pmb->je; j++) {
        for (int i=pmb->is; i<=pmb->ie; i++) {
          Real r(0.0), d0, mx, my, mz, e0, bx, by, bz;
          if (shk_dir == 1) r = pmb->pcoord->x1v(i);
          if (shk_dir == 2) r = pmb->pcoord->x2v(j);
          if (shk_dir == 3) r = pmb->pcoord->x3v(k);

          bx = 2.0/std::sqrt(4.0*PI);
          if (r > xfp) {
            d0 = 1.0;
            mx = 0.0;
            my = 0.0;
            mz = 0.0;
            by = 4.0/std::sqrt(4.0*PI);
            bz = 2.0/std::sqrt(4.0*PI);
            if (GENERAL_EOS) {
              e0 = pmb->peos->EgasFromRhoP(d0, 1.0);
            } else {
              e0 = 1.0/gm1;
            }
            e0 += 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
          } else if (r > xrp) {
            d0 = 1.3090;
            mx = 0.53432*d0;
            my = -0.094572*d0;
            mz = -0.047286*d0;
            by = 5.3452/std::sqrt(4.0*PI);
            bz = 2.6726/std::sqrt(4.0*PI);
            if (GENERAL_EOS) {
              e0 = pmb->peos->EgasFromRhoP(d0, 1.5844);
            } else {
              e0 = 1.5844/gm1;
            }
            e0 += 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
          } else if (r > xsp) {
            d0 = 1.3090;
            mx = 0.53432*d0;
            my = -0.18411*d0;
            mz = 0.17554*d0;
            by = 5.7083/std::sqrt(4.0*PI);
            bz = 1.7689/std::sqrt(4.0*PI);
            if (GENERAL_EOS) {
              e0 = pmb->peos->EgasFromRhoP(d0, 1.5844);
            } else {
              e0 = 1.5844/gm1;
            }
            e0 += 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
          } else if (r > xc) {
            d0 = 1.4735;
            mx = 0.57538*d0;
            my = 0.047601*d0;
            mz = 0.24734*d0;
            by = 5.0074/std::sqrt(4.0*PI);
            bz = 1.5517/std::sqrt(4.0*PI);
            if (GENERAL_EOS) {
              e0 = pmb->peos->EgasFromRhoP(d0, 1.9317);
            } else {
              e0 = 1.9317/gm1;
            }
            e0 += 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
          } else if (r > xsm) {
            d0 = 1.6343;
            mx = 0.57538*d0;
            my = 0.047601*d0;
            mz = 0.24734*d0;
            by = 5.0074/std::sqrt(4.0*PI);
            bz = 1.5517/std::sqrt(4.0*PI);
            e0 = 1.9317/gm1 + 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
          } else if (r > xrm) {
            d0 = 1.4903;
            mx = 0.60588*d0;
            my = 0.22157*d0;
            mz = 0.30125*d0;
            by = 5.5713/std::sqrt(4.0*PI);
            bz = 1.7264/std::sqrt(4.0*PI);
            if (GENERAL_EOS) {
              e0 = pmb->peos->EgasFromRhoP(d0, 1.6558);
            } else {
              e0 = 1.6558/gm1;
            }
            e0 += 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
          } else if (r > xfm) {
            d0 = 1.4903;
            mx = 0.60588*d0;
            my = 0.11235*d0;
            mz = 0.55686*d0;
            by = 5.0987/std::sqrt(4.0*PI);
            bz = 2.8326/std::sqrt(4.0*PI);
            if (GENERAL_EOS) {
              e0 = pmb->peos->EgasFromRhoP(d0, 1.6558);
            } else {
              e0 = 1.6558/gm1;
            }
            e0 += 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
          } else {
            d0 = 1.08;
            mx = 1.2*d0;
            my = 0.01*d0;
            mz = 0.5*d0;
            by = 3.6/std::sqrt(4.0*PI);
            bz = 2.0/std::sqrt(4.0*PI);
            if (GENERAL_EOS) {
              e0 = pmb->peos->EgasFromRhoP(d0, 0.95);
            } else {
              e0 = 0.95/gm1;
            }
            e0 += 0.5*((mx*mx+my*my+mz*mz)/d0 + (bx*bx+by*by+bz*bz));
          }

          err[IDN] += std::abs(d0 - pmb->phydro->u(IDN,k,j,i));
          err[im1] += std::abs(mx - pmb->phydro->u(im1,k,j,i));
          err[im2] += std::abs(my - pmb->phydro->u(im2,k,j,i));
          err[im3] += std::abs(mz - pmb->phydro->u(im3,k,j,i));
          err[IEN] += std::abs(e0 - pmb->phydro->u(IEN,k,j,i));
          err[NHYDRO + ib1] += std::abs(bx - pmb->pfield->bcc(ib1,k,j,i));
          err[NHYDRO + ib2] += std::abs(by - pmb->pfield->bcc(ib2,k,j,i));
          err[NHYDRO + ib3] += std::abs(bz - pmb->pfield->bcc(ib3,k,j,i));
        }
      }
    }

    // Errors in Sod solution
  } else {
    // Positions of shock, contact, head and foot of rarefaction for Sod test
    Real xs = 1.7522*tlim;
    Real xc = 0.92745*tlim;
    Real xf = -0.07027*tlim;
    Real xh = -1.1832*tlim;

    for (int k=pmb->ks; k<=pmb->ke; k++) {
      for (int j=pmb->js; j<=pmb->je; j++) {
        for (int i=pmb->is; i<=pmb->ie; i++) {
          Real r(0.0), d0, m0, e0;
          if (shk_dir == 1) r = pmb->pcoord->x1v(i);
          if (shk_dir == 2) r = pmb->pcoord->x2v(j);
          if (shk_dir == 3) r = pmb->pcoord->x3v(k);

          if (r > xs) {
            d0 = 0.125;
            m0 = 0.0;
            e0 = 0.25;
          } else if (r > xc) {
            d0 = 0.26557;
            m0 = 0.92745*d0;
            e0 = 0.87204;
          } else if (r > xf) {
            d0 = 0.42632;
            m0 = 0.92745*d0;
            e0 = 0.94118;
          } else if (r > xh) {
            Real v0 = 0.92745*(r-xh)/(xf-xh);
            d0 = 0.42632*std::pow((1.0+0.20046*(0.92745-v0)),5);
            m0 = v0*d0;
            e0 = (0.30313*std::pow((1.0+0.20046*(0.92745-v0)),7))/0.4 + 0.5*d0*v0*v0;
          } else {
            d0 = 1.0;
            m0 = 0.0;
            e0 = 2.5;
          }
          err[IDN] += std::abs(d0  - pmb->phydro->u(IDN,k,j,i));
          err[im1] += std::abs(m0  - pmb->phydro->u(im1,k,j,i));
          err[im2] += std::abs(0.0 - pmb->phydro->u(im2,k,j,i));
          err[im3] += std::abs(0.0 - pmb->phydro->u(im3,k,j,i));
          err[IEN] += std::abs(e0  - pmb->phydro->u(IEN,k,j,i));
        }
      }
    }
  }

  // normalize errors by number of cells, compute RMS
  for (int i=0; i<(NHYDRO+NFIELD); ++i) {
    err[i] = err[i]/static_cast<Real>(GetTotalCells());
  }
  Real rms_err = 0.0;
  for (int i=0; i<(NHYDRO+NFIELD); ++i) rms_err += SQR(err[i]);
  rms_err = std::sqrt(rms_err);

  // open output file and write out errors
  std::string fname;
  fname.assign("shock-errors.dat");
  std::stringstream msg;
  FILE *pfile;

  // The file exists -- reopen the file in append mode
  if ((pfile = std::fopen(fname.c_str(),"r")) != nullptr) {
    if ((pfile = std::freopen(fname.c_str(),"a",pfile)) == nullptr) {
      msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
          << std::endl << "Error output file could not be opened" <<std::endl;
      ATHENA_ERROR(msg);
    }

    // The file does not exist -- open the file in write mode and add headers
  } else {
    if ((pfile = std::fopen(fname.c_str(),"w")) == nullptr) {
      msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
          << std::endl << "Error output file could not be opened" <<std::endl;
      ATHENA_ERROR(msg);
    }
    std::fprintf(pfile,"# Nx1  Nx2  Nx3  Ncycle  RMS-Error  d  M1  M2  M3  E");
    if (MAGNETIC_FIELDS_ENABLED) std::fprintf(pfile,"  B1c  B2c  B3c");
    std::fprintf(pfile,"\n");
  }

  // write errors
  std::fprintf(pfile,"%d  %d",pmb->block_size.nx1,pmb->block_size.nx2);
  std::fprintf(pfile,"  %d  %d  %e",pmb->block_size.nx3,ncycle,rms_err);
  std::fprintf(pfile,"  %e  %e  %e  %e  %e",err[IDN],err[IM1],err[IM2],err[IM3],err[IEN]);
  if (MAGNETIC_FIELDS_ENABLED) {
    std::fprintf(pfile,"  %e  %e  %e",err[NHYDRO+IB1],err[NHYDRO+IB2],err[NHYDRO+IB3]);
  }
  std::fprintf(pfile,"\n");
  std::fclose(pfile);

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Problem Generator for the shock tube tests
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  std::stringstream msg;

  // parse shock direction: {1,2,3} -> {x1,x2,x3}
  int shk_dir = pin->GetInteger("problem","shock_dir");

  // parse shock location (must be inside grid)
  Real xshock = pin->GetReal("problem","xshock");
  if (shk_dir == 1 && (xshock < pmy_mesh->mesh_size.x1min ||
                       xshock > pmy_mesh->mesh_size.x1max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x1 domain for shkdir=" << shk_dir << std::endl;
    ATHENA_ERROR(msg);
  }
  if (shk_dir == 2 && (xshock < pmy_mesh->mesh_size.x2min ||
                       xshock > pmy_mesh->mesh_size.x2max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x2 domain for shkdir=" << shk_dir << std::endl;
    ATHENA_ERROR(msg);
  }
  if (shk_dir == 3 && (xshock < pmy_mesh->mesh_size.x3min ||
                       xshock > pmy_mesh->mesh_size.x3max)) {
    msg << "### FATAL ERROR in Problem Generator" << std::endl << "xshock="
        << xshock << " lies outside x3 domain for shkdir=" << shk_dir << std::endl;
    ATHENA_ERROR(msg);
  }

  // Parse left state read from input file: dl,ul,vl,wl,[pl]
  Real wl[NHYDRO+NFIELD];
  wl[IDN] = pin->GetReal("problem","dl");
  wl[IVX] = pin->GetReal("problem","ul");
  wl[IVY] = pin->GetReal("problem","vl");
  wl[IVZ] = pin->GetReal("problem","wl");
  if (NON_BAROTROPIC_EOS) {
    if (pin->DoesParameterExist("problem","Tl"))
      wl[IPR] = press(wl[IDN], pin->GetReal("problem","Tl"));
    else
      wl[IPR] = pin->GetReal("problem","pl");
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    wl[NHYDRO  ] = pin->GetReal("problem","bxl");
    wl[NHYDRO+1] = pin->GetReal("problem","byl");
    wl[NHYDRO+2] = pin->GetReal("problem","bzl");
  }

  // Parse right state read from input file: dr,ur,vr,wr,[pr]
  Real wr[NHYDRO+NFIELD];
  wr[IDN] = pin->GetReal("problem","dr");
  wr[IVX] = pin->GetReal("problem","ur");
  wr[IVY] = pin->GetReal("problem","vr");
  wr[IVZ] = pin->GetReal("problem","wr");
  if (NON_BAROTROPIC_EOS) {
    if (pin->DoesParameterExist("problem","Tr"))
      wr[IPR] = press(wr[IDN], pin->GetReal("problem","Tr"));
    else
      wr[IPR] = pin->GetReal("problem","pr");
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    wr[NHYDRO  ] = pin->GetReal("problem","bxr");
    wr[NHYDRO+1] = pin->GetReal("problem","byr");
    wr[NHYDRO+2] = pin->GetReal("problem","bzr");
  }

  // Initialize the discontinuity in the Hydro variables ---------------------------------

  switch(shk_dir) {
    //--- shock in 1-direction
    case 1:
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            if (pcoord->x1v(i) < xshock) {
              phydro->u(IDN,k,j,i) = wl[IDN];
              phydro->u(IM1,k,j,i) = wl[IVX]*wl[IDN];
              phydro->u(IM2,k,j,i) = wl[IVY]*wl[IDN];
              phydro->u(IM3,k,j,i) = wl[IVZ]*wl[IDN];
              if (NON_BAROTROPIC_EOS) {
                if (GENERAL_EOS) {
                  phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(wl[IDN], wl[IPR]);
                } else {
                  phydro->u(IEN,k,j,i) = wl[IPR]/(peos->GetGamma() - 1.0);
                }
                phydro->u(IEN,k,j,i) += 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY]
                                                     + wl[IVZ]*wl[IVZ]);
              }
            } else {
              phydro->u(IDN,k,j,i) = wr[IDN];
              phydro->u(IM1,k,j,i) = wr[IVX]*wr[IDN];
              phydro->u(IM2,k,j,i) = wr[IVY]*wr[IDN];
              phydro->u(IM3,k,j,i) = wr[IVZ]*wr[IDN];
              if (NON_BAROTROPIC_EOS) {
                if (GENERAL_EOS) {
                  phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(wr[IDN], wr[IPR]);
                } else {
                  phydro->u(IEN,k,j,i) = wr[IPR]/(peos->GetGamma() - 1.0);
                }
                phydro->u(IEN,k,j,i) += 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY]
                                                     + wr[IVZ]*wr[IVZ]);
              }
            }
          }
        }
      }
      break;
      //--- shock in 2-direction
    case 2:
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          if (pcoord->x2v(j) < xshock) {
            for (int i=is; i<=ie; ++i) {
              phydro->u(IDN,k,j,i) = wl[IDN];
              phydro->u(IM2,k,j,i) = wl[IVX]*wl[IDN];
              phydro->u(IM3,k,j,i) = wl[IVY]*wl[IDN];
              phydro->u(IM1,k,j,i) = wl[IVZ]*wl[IDN];
              if (NON_BAROTROPIC_EOS) {
                if (GENERAL_EOS) {
                  phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(wl[IDN], wl[IPR]);
                } else {
                  phydro->u(IEN,k,j,i) = wl[IPR]/(peos->GetGamma() - 1.0);
                }
                phydro->u(IEN,k,j,i) += 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY]
                                                     + wl[IVZ]*wl[IVZ]);
              }
            }
          } else {
            for (int i=is; i<=ie; ++i) {
              phydro->u(IDN,k,j,i) = wr[IDN];
              phydro->u(IM2,k,j,i) = wr[IVX]*wr[IDN];
              phydro->u(IM3,k,j,i) = wr[IVY]*wr[IDN];
              phydro->u(IM1,k,j,i) = wr[IVZ]*wr[IDN];
              if (NON_BAROTROPIC_EOS) {
                if (GENERAL_EOS) {
                  phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(wr[IDN], wr[IPR]);
                } else {
                  phydro->u(IEN,k,j,i) = wr[IPR]/(peos->GetGamma() - 1.0);
                }
                phydro->u(IEN,k,j,i) += 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY]
                                                     + wr[IVZ]*wr[IVZ]);
              }
            }
          }
        }
      }
      break;

      //--- shock in 3-direction
    case 3:
      for (int k=ks; k<=ke; ++k) {
        if (pcoord->x3v(k) < xshock) {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->u(IDN,k,j,i) = wl[IDN];
              phydro->u(IM3,k,j,i) = wl[IVX]*wl[IDN];
              phydro->u(IM1,k,j,i) = wl[IVY]*wl[IDN];
              phydro->u(IM2,k,j,i) = wl[IVZ]*wl[IDN];
              if (NON_BAROTROPIC_EOS) {
                if (GENERAL_EOS) {
                  phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(wl[IDN], wl[IPR]);
                } else {
                  phydro->u(IEN,k,j,i) = wl[IPR]/(peos->GetGamma() - 1.0);
                }
                phydro->u(IEN,k,j,i) += 0.5*wl[IDN]*(wl[IVX]*wl[IVX] + wl[IVY]*wl[IVY]
                                                     + wl[IVZ]*wl[IVZ]);
              }
            }
          }
        } else {
          for (int j=js; j<=je; ++j) {
            for (int i=is; i<=ie; ++i) {
              phydro->u(IDN,k,j,i) = wr[IDN];
              phydro->u(IM3,k,j,i) = wr[IVX]*wr[IDN];
              phydro->u(IM1,k,j,i) = wr[IVY]*wr[IDN];
              phydro->u(IM2,k,j,i) = wr[IVZ]*wr[IDN];
              if (NON_BAROTROPIC_EOS) {
                if (GENERAL_EOS) {
                  phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(wr[IDN], wr[IPR]);
                } else {
                  phydro->u(IEN,k,j,i) = wr[IPR]/(peos->GetGamma() - 1.0);
                }
                phydro->u(IEN,k,j,i) += 0.5*wr[IDN]*(wr[IVX]*wr[IVX] + wr[IVY]*wr[IVY]
                                                     + wr[IVZ]*wr[IVZ]);
              }
            }
          }
        }
      }
      break;

    default:
      msg << "### FATAL ERROR in Problem Generator" << std::endl
          << "shock_dir=" << shk_dir << " must be either 1,2, or 3" << std::endl;
      ATHENA_ERROR(msg);
  }

  // now set face-centered (interface) magnetic fields -----------------------------------

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          if (shk_dir==1 && pcoord->x1v(i) < xshock) {
            pfield->b.x1f(k,j,i) = wl[NHYDRO  ];
            pfield->b.x2f(k,j,i) = wl[NHYDRO+1];
            pfield->b.x3f(k,j,i) = wl[NHYDRO+2];
          } else if (shk_dir==2 && pcoord->x2v(j) < xshock) {
            pfield->b.x1f(k,j,i) = wl[NHYDRO+2];
            pfield->b.x2f(k,j,i) = wl[NHYDRO  ];
            pfield->b.x3f(k,j,i) = wl[NHYDRO+1];
          } else if (shk_dir==3 && pcoord->x3v(k) < xshock) {
            pfield->b.x1f(k,j,i) = wl[NHYDRO+1];
            pfield->b.x2f(k,j,i) = wl[NHYDRO+2];
            pfield->b.x3f(k,j,i) = wl[NHYDRO];
          }

          if (shk_dir==1 && pcoord->x1v(i) >= xshock) {
            pfield->b.x1f(k,j,i) = wr[NHYDRO  ];
            pfield->b.x2f(k,j,i) = wr[NHYDRO+1];
            pfield->b.x3f(k,j,i) = wr[NHYDRO+2];
          } else if (shk_dir==2 && pcoord->x2v(j) >= xshock) {
            pfield->b.x1f(k,j,i) = wr[NHYDRO+2];
            pfield->b.x2f(k,j,i) = wr[NHYDRO  ];
            pfield->b.x3f(k,j,i) = wr[NHYDRO+1];
          } else if (shk_dir==3 && pcoord->x3v(k) >= xshock)  {
            pfield->b.x1f(k,j,i) = wr[NHYDRO+1];
            pfield->b.x2f(k,j,i) = wr[NHYDRO+2];
            pfield->b.x3f(k,j,i) = wr[NHYDRO];
          }
          if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) += 0.5*(SQR(pfield->b.x1f(k,j,i))
                                         + SQR(pfield->b.x2f(k,j,i))
                                         + SQR(pfield->b.x3f(k,j,i)));
          }
        }
      }
    }
    // end by adding bi.x1 at ie+1, bi.x2 at je+1, and bi.x3 at ke+1
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        pfield->b.x1f(k,j,ie+1) = pfield->b.x1f(k,j,ie);
      }
    }
    for (int k=ks; k<=ke; ++k) {
      for (int i=is; i<=ie; ++i) {
        pfield->b.x2f(k,je+1,i) = pfield->b.x2f(k,je,i);
      }
    }
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        pfield->b.x3f(ke+1,j,i) = pfield->b.x3f(ke,j,i);
      }
    }
  }

  // uniformly fill all scalars to have equal concentration
  // mass fraction? or concentration?
  constexpr int scalar_norm = NSCALARS > 0 ? NSCALARS : 1.0;
  if (NSCALARS > 0) {
    for (int n=0; n<NSCALARS; ++n) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            pscalars->s(n,k,j,i) = 1.0/scalar_norm*phydro->u(IDN,k,j,i);
          }
        }
      }
    }
  }
  return;
}
