//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file linear_wave.c
//  \brief Linear wave problem generator for 1D/2D/3D problems.
//
// In 1D, the problem is setup along one of the three coordinate axes (specified by
// setting [ang_2,ang_3] = 0.0 or PI/2 in the input file).  In 2D/3D this routine
// automatically sets the wavevector along the domain diagonal.
//========================================================================================

// C++ headers
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <algorithm>  // min, max
#include <cmath>

// Athena++ headers
#include "../globals.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../eos/eos.hpp"
#include "../coordinates/coordinates.hpp"
#include "../gravity/gravity.hpp"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator does not support magnetic fields"
#endif

// with functions A1,2,3 which compute vector potentials
static Real ang_2, ang_3; // Rotation angles about the y and z' axis
static Real sin_a2, cos_a2, sin_a3, cos_a3;
static Real amp, njeans, lambda, kwave; // amplitude, Wavelength, 2*PI/wavelength
static Real cs2,gam,gm1,omega,omega2, gconst;
static Real ev[NWAVE], rem[NWAVE][NWAVE], lem[NWAVE][NWAVE];
static Real d0,p0,v0,u0,w0,va,b0;


void Mesh::InitUserMeshData(ParameterInput *pin) {
  Real x1size = mesh_size.x1max - mesh_size.x1min;
  Real x2size = mesh_size.x2max - mesh_size.x2min;
  Real x3size = mesh_size.x3max - mesh_size.x3min;
  amp = pin->GetReal("problem","amp");
  njeans = pin->GetReal("problem","njeans");
  ang_2 = pin->GetOrAddReal("problem","ang_2",-999.9);
  ang_3 = pin->GetOrAddReal("problem","ang_3",-999.9);
  // User should never input -999.9 in angles
  if (ang_3 == -999.9) ang_3 = atan(x1size/x2size);
  sin_a3 = sin(ang_3);
  cos_a3 = cos(ang_3);
  if (ang_2 == -999.9) ang_2 = atan(0.5*(x1size*cos_a3 + x2size*sin_a3)/x3size);
  sin_a2 = sin(ang_2);
  cos_a2 = cos(ang_2);
  Real x1 = x1size*cos_a2*cos_a3;
  Real x2 = x2size*cos_a2*sin_a3;
  Real x3 = x3size*sin_a2;
  njeans = pin->GetReal("problem","njeans");
  // For lambda choose the smaller of the 3
  lambda = x1;
  if (mesh_size.nx2 > 1 && ang_3 != 0.0) lambda = std::min(lambda,x2);
  if (mesh_size.nx3 > 1 && ang_2 != 0.0) lambda = std::min(lambda,x3);

  d0 = 1.0, p0 = 1.0;
  u0 = 0.0, v0 = 0.0, w0 = 0.0;
  va = 0.0, b0 = 0.0;

  if (NON_BAROTROPIC_EOS) {
    gam = pin->GetReal("hydro","gamma");
    p0 = 1.0/gam;
    gm1 = gam-1.0;
    cs2 = gam*p0/d0;
  } else {
    Real iso_cs = pin->GetReal("hydro","iso_sound_speed");
    cs2 = SQR(iso_cs);
  }
  gconst = cs2*PI*njeans*njeans/(d0*lambda*lambda);

  kwave = 2.0*PI/lambda;
  omega2 = SQR(kwave)*cs2*(1.0 - SQR(njeans));
  omega = std::sqrt(fabs(omega2));

  if (SELF_GRAVITY_ENABLED) {
    SetGravitationalConstant(gconst);
    Real eps = pin->GetOrAddReal("problem","grav_eps", 0.0);
    SetGravityThreshold(eps);
    SetMeanDensity(d0);
  }

  if (Globals::my_rank==0) {
    //moved print statements here from Meshblock::ProblemGenerator
    std::cout << "four_pi_G " << gconst*4.0*PI << std::endl;
    std::cout << "lambda " << lambda << std::endl;
    std::cout << "period " << (2*PI/omega) << std::endl;
    std::cout << "angle2 " << ang_2*180./PI << " "
              << sin_a2 << " " << cos_a2 << std::endl;
    std::cout << "angle3 " << ang_3*180./PI << " "
              << sin_a3 << " " << cos_a3 << std::endl;
  }




  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real x, sinkx, coskx;

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
    Real x = cos_a2*(pcoord->x1v(i)*cos_a3 + pcoord->x2v(j)*sin_a3)
           + pcoord->x3v(k)*sin_a2;
    sinkx = sin(x*kwave);
    coskx = cos(x*kwave);

    phydro->u(IDN,k,j,i) = d0*(1.0+amp*sinkx+amp*amp*sin(pcoord->x1v(i)*kwave));

    //when unstable initial v omega/kwave*amp*coskx
    //when stable initial v 0
    Real m = (omega2 < 0) ? d0*(omega/kwave)*amp*coskx:0.0;

    phydro->u(IM1,k,j,i) = m*cos_a3*cos_a2;
    phydro->u(IM2,k,j,i) = m*sin_a3*cos_a2;
    phydro->u(IM3,k,j,i) = m*sin_a2;

    if (NON_BAROTROPIC_EOS) {
      phydro->u(IEN,k,j,i) = p0/gm1*(1.0 + gam*amp*sinkx);
      phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
      phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
      phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
    }
  }}}

//  pmy_mesh->tlim=pin->SetReal("time","tlim",2.0*PI/omega*2.0);

}

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;
  if (omega2 < 0) {
    if (Globals::my_rank==0)
      std::cout << "This problem is Jeans unstable, njeans = " << njeans << std::endl;
    //    return;
  }

  MeshBlock *pmb = pblock;
  // Initialize errors to zero
  Real l1_err[NHYDRO+NFIELD],max_err[NHYDRO+NFIELD];
  for (int i=0; i<(NHYDRO+NFIELD); ++i) {
    l1_err[i]=0.0;
    max_err[i]=0.0;
  }

  Gravity *pgrav = pblock->pgrav;
  Hydro *phydro = pblock->phydro;
  Coordinates *pcoord = pblock->pcoord;
  Real sinkx, coskx, sinot, cosot;
  int is=pblock->is, ie=pblock->ie;
  int js=pblock->js, je=pblock->je;
  int ks=pblock->ks, ke=pblock->ke;

  Real tlim = time;
  while (pmb != NULL) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real x = cos_a2*(pcoord->x1v(i)*cos_a3 + pcoord->x2v(j)*sin_a3)
               + pcoord->x3v(k)*sin_a2;
        sinkx = sin(x*kwave);
        coskx = cos(x*kwave);
        if (omega2 < 0) {
          sinot = -exp(omega*tlim);//time dependent factor of vel
          // unstable case v = amp*omega/k * coskx * e^omega*t
          // minus sign counters minus sign in m
          cosot = exp(omega*tlim);//time dependent factor of rho
        } else {
          sinot = sin(omega*tlim);//time dependent factor of vel
          cosot = cos(omega*tlim);//time dependent factor of rho
        }
        Real den=d0*(1.0+amp*sinkx*cosot);
        l1_err[IDN] += fabs(den - phydro->u(IDN,k,j,i));
        max_err[IDN] = std::max(static_cast<Real>(fabs(den - phydro->u(IDN,k,j,i))),
                                max_err[IDN]);

        Real m = -den*(omega/kwave)*amp*coskx*sinot;
        Real m1 = m*cos_a3*cos_a2;
        Real m2 = m*sin_a3*cos_a2;
        Real m3 = m*sin_a2;

        l1_err[IM1] += fabs(m1-phydro->u(IM1,k,j,i));
        l1_err[IM2] += fabs(m2-phydro->u(IM2,k,j,i));
        l1_err[IM3] += fabs(m3-phydro->u(IM3,k,j,i));
        max_err[IM1] = std::max(static_cast<Real>(fabs(m1-phydro->u(IM1,k,j,i))),
                                max_err[IM1]);
        max_err[IM2] = std::max(static_cast<Real>(fabs(m2-phydro->u(IM2,k,j,i))),
                                max_err[IM2]);
        max_err[IM3] = std::max(static_cast<Real>(fabs(m3-phydro->u(IM3,k,j,i))),
                                max_err[IM3]);
        if (NON_BAROTROPIC_EOS) {
          Real e0 = p0*(1 + gam*amp*sinkx*cosot);///gm1 + 0.5*m*m/den;
          l1_err[IEN] += fabs(e0 - phydro->w(IEN,k,j,i));
          max_err[IEN] = std::max(static_cast<Real>(fabs(e0 - phydro->w(IEN,k,j,i))),
                                  max_err[IEN]);
        }
      }
    }}
    pmb=pmb->next;
  }

  for (int i=0; i<(NHYDRO+NFIELD); ++i) {
    l1_err[i] = l1_err[i]/static_cast<Real>(GetTotalCells());
  }
  Real rms_err = 0.0, max_max_over_l1=0.0;
#ifdef MPI_PARALLEL
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE,&l1_err,(NHYDRO+NFIELD),MPI_ATHENA_REAL,MPI_SUM,0,
               MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE,&max_err,(NHYDRO+NFIELD),MPI_ATHENA_REAL,MPI_MAX,0,
               MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&l1_err,&l1_err,(NHYDRO+NFIELD),MPI_ATHENA_REAL,MPI_SUM,0,
               MPI_COMM_WORLD);
    MPI_Reduce(&max_err,&max_err,(NHYDRO+NFIELD),MPI_ATHENA_REAL,MPI_MAX,0,
               MPI_COMM_WORLD);
  }
#endif

  // only the root process outputs the data
  if (Globals::my_rank == 0) {
    // compute rms error
    for (int i=0; i<(NHYDRO+NFIELD); ++i) {
       rms_err += SQR(l1_err[i]);
       max_max_over_l1 = std::max(max_max_over_l1, (max_err[i]/l1_err[i]));
    }
    rms_err = std::sqrt(rms_err);

    // open output file and write out errors
    std::string fname;
    fname.assign("jeans-errors.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = fopen(fname.c_str(),"r")) != NULL) {
      if ((pfile = freopen(fname.c_str(),"a",pfile)) == NULL) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Error output file could not be opened" <<std::endl;
        throw std::runtime_error(msg.str().c_str());
      }

    // The file does not exist -- open the file in write mode and add headers
    } else {
      if ((pfile = fopen(fname.c_str(),"w")) == NULL) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Error output file could not be opened" <<std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      fprintf(pfile,"# Nx1  Nx2  Nx3  Ncycle  ");
      fprintf(pfile,"RMS-L1-Error  d_L1  M1_L1  M2_L1  M3_L1");
      if (NON_BAROTROPIC_EOS) fprintf(pfile,"  E_L1 ");
      if (MAGNETIC_FIELDS_ENABLED) fprintf(pfile,"  B1c_L1  B2c_L1  B3c_L1");
      fprintf(pfile,"  Largest-Max/L1  d_max  M1_max  M2_max  M3_max");
      if (NON_BAROTROPIC_EOS) fprintf(pfile,"  E_max ");
      if (MAGNETIC_FIELDS_ENABLED) fprintf(pfile,"  B1c_max  B2c_max  B3c_max");
      fprintf(pfile,"\n");
    }

    // write errors
    fprintf(pfile,"%d  %d",mesh_size.nx1,mesh_size.nx2);
    fprintf(pfile,"  %d  %d",mesh_size.nx3,ncycle);
    fprintf(pfile,"  %e  %e",rms_err,l1_err[IDN]);
    fprintf(pfile,"  %e  %e  %e",l1_err[IM1],l1_err[IM2],l1_err[IM3]);
    if (NON_BAROTROPIC_EOS) fprintf(pfile,"  %e",l1_err[IEN]);
    if (MAGNETIC_FIELDS_ENABLED) {
      fprintf(pfile,"  %e",l1_err[NHYDRO+IB1]);
      fprintf(pfile,"  %e",l1_err[NHYDRO+IB2]);
      fprintf(pfile,"  %e",l1_err[NHYDRO+IB3]);
    }
    fprintf(pfile,"  %e  %e  ",max_max_over_l1,max_err[IDN]);
    fprintf(pfile,"%e  %e  %e",max_err[IM1],max_err[IM2],max_err[IM3]);
    if (NON_BAROTROPIC_EOS) fprintf(pfile,"  %e",max_err[IEN]);
    if (MAGNETIC_FIELDS_ENABLED) {
      fprintf(pfile,"  %e",max_err[NHYDRO+IB1]);
      fprintf(pfile,"  %e",max_err[NHYDRO+IB2]);
      fprintf(pfile,"  %e",max_err[NHYDRO+IB3]);
    }
    fprintf(pfile,"\n");
    fclose(pfile);
  }

  return;
}
