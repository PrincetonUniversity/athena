//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cpaw.c
//  \brief Circularly polarized Alfven wave (CPAW) for 1D/2D/3D problems
//
// In 1D, the problem is setup along one of the three coordinate axes (specified by
// setting [ang_2,ang_3] = 0.0 or PI/2 in the input file).  In 2D/3D this routine
// automatically sets the wavevector along the domain diagonal.
//
// Can be used for [standing/traveling] waves [(problem/v_par=1.0)/(problem/v_par=0.0)]
//
// REFERENCE: G. Toth,  "The div(B)=0 constraint in shock capturing MHD codes", JCP,
//   161, 605 (2000)

// C/C++ headers
#include <algorithm>
#include <cmath>      // std::sqrt(), std::fabs()
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../reconstruct/reconstruction.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

// Parameters which define initial solution -- made global so that they can be shared
// with functions A1,2,3 which compute vector potentials
static Real den, pres, gm1, b_par, b_perp, v_perp, v_par;
static Real ang_2, ang_3; // Rotation angles about the y and z' axis
static bool ang_2_vert, ang_3_vert; // Switches to set ang_2 and/or ang_3 to pi/2
static Real fac, sin_a2, cos_a2, sin_a3, cos_a3;
static Real lambda, k_par; // Wavelength, 2*PI/wavelength

// functions to compute vector potential to initialize the solution
static Real A1(const Real x1, const Real x2, const Real x3);
static Real A2(const Real x1, const Real x2, const Real x3);
static Real A3(const Real x1, const Real x2, const Real x3);

// edge-averaged values
static Real A1_ave(const Real x1f, const Real x1f_ip1, const Real x2, const Real x3);
static Real A2_ave(const Real x1, const Real x2f, const Real x2f_jp1, const Real x3);
static Real A3_ave(const Real x1, const Real x2, const Real x3f, const Real x3f_kp1);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
// Initialize magnetic field parameters
// For wavevector along coordinate axes, set desired values of ang_2/ang_3.
//    For example, for 1D problem use ang_2 = ang_3 = 0.0
//    For wavevector along grid diagonal, do not input values for ang_2/ang_3.
// Code below will automatically calculate these imposing periodicity and exactly one
// wavelength along each grid direction
  b_par = pin->GetReal("problem","b_par");
  b_perp = pin->GetReal("problem","b_perp");
  v_par = pin->GetReal("problem","v_par");
  ang_2 = pin->GetOrAddReal("problem","ang_2",-999.9);
  ang_3 = pin->GetOrAddReal("problem","ang_3",-999.9);
  ang_2_vert = pin->GetOrAddBoolean("problem", "ang_2_vert", false);
  ang_3_vert = pin->GetOrAddBoolean("problem", "ang_3_vert", false);
  Real dir = pin->GetOrAddReal("problem","dir",1); // right(1)/left(2) polarization
  if (NON_BAROTROPIC_EOS) {
    Real gam   = pin->GetReal("hydro","gamma");
    gm1 = (gam - 1.0);
  }
  pres = pin->GetReal("problem","pres");
  den = 1.0;

  Real x1size = mesh_size.x1max - mesh_size.x1min;
  Real x2size = mesh_size.x2max - mesh_size.x2min;
  Real x3size = mesh_size.x3max - mesh_size.x3min;

  // User should never input -999.9 in angles
  if (ang_3 == -999.9) ang_3 = atan(x1size/x2size);
  sin_a3 = sin(ang_3);
  cos_a3 = cos(ang_3);

  // Override ang_3 input and hardcode vertical (along x2 axis) wavevector
  if (ang_3_vert == true) {
    sin_a3 = 1.0;
    cos_a3 = 0.0;
    ang_3 = 0.5*M_PI;
  }

  if (ang_2 == -999.9) ang_2 = atan(0.5*(x1size*cos_a3 + x2size*sin_a3)/x3size);
  sin_a2 = sin(ang_2);
  cos_a2 = cos(ang_2);

  // Override ang_2 input and hardcode vertical (along x3 axis) wavevector
  if (ang_2_vert == true) {
    sin_a2 = 1.0;
    cos_a2 = 0.0;
    ang_2 = 0.5*M_PI;
  }

  Real x1 = x1size*cos_a2*cos_a3;
  Real x2 = x2size*cos_a2*sin_a3;
  Real x3 = x3size*sin_a2;

  // For lambda choose the smaller of the 3
  lambda = x1;
  if (mesh_size.nx2 > 1 && ang_3 != 0.0) lambda = std::min(lambda,x2);
  if (mesh_size.nx3 > 1 && ang_2 != 0.0) lambda = std::min(lambda,x3);

  // If cos_a2 or cos_a3 = 0, need to override lambda
  if (ang_3_vert == true)
    lambda = x2;
  if (ang_2_vert == true)
    lambda = x3;

  // Initialize k_parallel
  k_par = 2.0*(PI)/lambda;
  v_perp = b_perp/std::sqrt(den);

  if (dir == 1) // right polarization
    fac = 1.0;
  else          // left polarization
    fac = -1.0;
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Compute L1 error in CPAW and output to file
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;

  // Initialize errors to zero
  Real err[NHYDRO+NFIELD];
  for (int i=0; i<(NHYDRO+NFIELD); ++i) err[i]=0.0;

  MeshBlock *pmb = pblock;
  BoundaryValues *pbval;
  while (pmb != NULL) {
    pbval=pmb->pbval;
    int il=pmb->is, iu=pmb->ie, jl=pmb->js, ju=pmb->je, kl=pmb->ks, ku=pmb->ke;
    // adjust loop limits for fourth order error calculation
    //------------------------------------------------
    if (pmb->precon->correct_err == true) {
      // Expand loop limits on all sides by one
      if (pbval->nblevel[1][1][0]!=-1) il-=1;
      if (pbval->nblevel[1][1][2]!=-1) iu+=1;
      if (pbval->nblevel[1][0][1]!=-1) jl-=1;
      if (pbval->nblevel[1][2][1]!=-1) ju+=1;
      if (pbval->nblevel[0][1][1]!=-1) kl-=1;
      if (pbval->nblevel[2][1][1]!=-1) ku+=1;
    }
    // Save analytic solution of conserved variables in 4D scratch array
    AthenaArray<Real> cons_;
    int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
    int ncells2 = 1, ncells3 = 1;
    if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
    if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);
    // Even for MHD, there are only cell-centered mesh variables
    int ncells4 = NHYDRO + NFIELD;
    int nl = 0;
    int nu = ncells4-1;
    cons_.NewAthenaArray(ncells4, ncells3, ncells2, ncells1);

    //  Compute errors at cell centers
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=il; i<=iu; i++) {
          Real x = cos_a2*(pmb->pcoord->x1v(i)*cos_a3 + pmb->pcoord->x2v(j)*sin_a3)
              + pmb->pcoord->x3v(k)*sin_a2;
          Real sn = sin(k_par*x);
          Real cs = fac*cos(k_par*x);

          Real mx = den*v_par;
          Real my = -fac*den*v_perp*sn;
          Real mz = -fac*den*v_perp*cs;
          Real m1 = mx*cos_a2*cos_a3 - my*sin_a3 - mz*sin_a2*cos_a3;
          Real m2 = mx*cos_a2*sin_a3 + my*cos_a3 - mz*sin_a2*sin_a3;
          Real m3 = mx*sin_a2                    + mz*cos_a2;

          // Store analytic solution at cell-centers
          cons_(IDN,k,j,i) = den;
          cons_(IM1,k,j,i) = m1;
          cons_(IM2,k,j,i) = m2;
          cons_(IM3,k,j,i) = m3;

          Real bx = b_par;
          Real by = b_perp*sn;
          Real bz = b_perp*cs;
          Real b1 = bx*cos_a2*cos_a3 - by*sin_a3 - bz*sin_a2*cos_a3;
          Real b2 = bx*cos_a2*sin_a3 + by*cos_a3 - bz*sin_a2*sin_a3;
          Real b3 = bx*sin_a2                    + bz*cos_a2;
          cons_(NHYDRO+IB1,k,j,i) = b1;
          cons_(NHYDRO+IB2,k,j,i) = b2;
          cons_(NHYDRO+IB3,k,j,i) = b3;

          if (NON_BAROTROPIC_EOS) {
            Real e0 = pres/gm1 + 0.5*(m1*m1 + m2*m2 + m3*m3)/den
                + 0.5*(b1*b1+b2*b2+b3*b3);
            cons_(IEN,k,j,i) = e0;
          }
        }
      }
    }
    // begin fourth-order error correction
    // -------------------------------
    if (pmb->precon->correct_err == true) {
      // Restore loop limits to real cells only
      il=pmb->is, iu=pmb->ie, jl=pmb->js, ju=pmb->je, kl=pmb->ks, ku=pmb->ke;

      // Compute and store Laplacian of cell-centered conserved variables, Hydro and Bcc
      AthenaArray<Real> delta_cons_;
      delta_cons_.NewAthenaArray(ncells4, ncells3, ncells2, ncells1);
      pmb->pcoord->Laplacian(cons_, delta_cons_, il, iu, jl, ju, kl, ku, nl, nu);

      // TODO(kfelker): assuming uniform mesh with dx1f=dx2f=dx3f, so this factors out
      // TODO(kfelker): also, this may need to be dx1v, since Laplacian is cell-centered
      Real h = pmb->pcoord->dx1f(il);  // pco->dx1f(i); inside loop
      Real C = (h*h)/24.0;

      // Compute fourth-order approximation to cell-averaged conserved variables
      for (int n=nl; n<=nu; ++n) {
        for (int k=kl; k<=ku; ++k) {
          for (int j=jl; j<=ju; ++j) {
            for (int i=il; i<=iu; ++i) {
              cons_(n,k,j,i) = cons_(n,k,j,i) + C*delta_cons_(n,k,j,i);
            }
          }
        }
      }
    } // end if (pmb->precon->correct_err == true)
    // ------- end fourth-order error calculation

    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          // Load cell-averaged <U>, either midpoint approx. or fourth-order approx
          Real den = cons_(IDN,k,j,i);
          Real m1 = cons_(IM1,k,j,i);
          Real m2 = cons_(IM2,k,j,i);
          Real m3 = cons_(IM3,k,j,i);
          // Weight l1 error by cell volume
          Real vol = pmb->pcoord->GetCellVolume(k, j, i);

          err[IDN] += std::fabs(den - pmb->phydro->u(IDN,k,j,i))*vol;
          err[IM1] += std::fabs(m1 - pmb->phydro->u(IM1,k,j,i))*vol;
          err[IM2] += std::fabs(m2 - pmb->phydro->u(IM2,k,j,i))*vol;
          err[IM3] += std::fabs(m3 - pmb->phydro->u(IM3,k,j,i))*vol;

          Real b1 = cons_(NHYDRO+IB1,k,j,i);
          Real b2 = cons_(NHYDRO+IB2,k,j,i);
          Real b3 = cons_(NHYDRO+IB3,k,j,i);
          err[NHYDRO + IB1] += std::fabs(b1 - pmb->pfield->bcc(IB1,k,j,i))*vol;
          err[NHYDRO + IB2] += std::fabs(b2 - pmb->pfield->bcc(IB2,k,j,i))*vol;
          err[NHYDRO + IB3] += std::fabs(b3 - pmb->pfield->bcc(IB3,k,j,i))*vol;

          if (NON_BAROTROPIC_EOS) {
            Real e0 = cons_(IEN,k,j,i);
            err[IEN] += std::fabs(e0 - pmb->phydro->u(IEN,k,j,i))*vol;
          }
        }
      }
    }
    // TODO(kfelker): don't de/allocate per MeshBlock
    cons_.DeleteAthenaArray();
    pmb=pmb->next;
  }

  // TODO(kfelker): extend to MPI and l1, lmax errors like linear_wave.cpp
  // normalize errors by volume
  Real vol= (mesh_size.x1max-mesh_size.x1min)*(mesh_size.x2max-mesh_size.x2min)
      *(mesh_size.x3max-mesh_size.x3min);
  for (int i=0; i<(NHYDRO+NFIELD); ++i) err[i] = err[i]/vol;

  // compute RMS error
  Real rms_err = 0.0;
  for (int i=0; i<(NHYDRO+NFIELD); ++i) rms_err += SQR(err[i]);
  rms_err = std::sqrt(rms_err);

  // open output file and write out errors
  std::string fname;
  fname.assign("cpaw-errors.dat");
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
    fprintf(pfile,"# Nx1  Nx2  Nx3  Ncycle  RMS-Error  d  M1  M2  M3");
    if (NON_BAROTROPIC_EOS) fprintf(pfile,"  E");
    fprintf(pfile,"  B1c  B2c  B3c");
    fprintf(pfile,"\n");
  }

  // write errors
  fprintf(pfile,"%d  %d",mesh_size.nx1,mesh_size.nx2);
  fprintf(pfile,"  %d  %d  %e",mesh_size.nx3,ncycle,rms_err);
  fprintf(pfile,"  %e  %e  %e  %e",err[IDN],err[IM1],err[IM2],err[IM3]);
  if (NON_BAROTROPIC_EOS) fprintf(pfile,"  %e",err[IEN]);
  fprintf(pfile,"  %e  %e  %e",err[NHYDRO+IB1],err[NHYDRO+IB2],err[NHYDRO+IB3]);
  fprintf(pfile,"\n");
  fclose(pfile);

  return;
}

//========================================================================================
//! \fn ProblemGenerator
//  \brief circularly polarized Alfven wave problem generator for 1D/2D/3D problems.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  AthenaArray<Real> a1,a2,a3;
  int nx1 = (ie-is)+1 + 2*(NGHOST);
  int nx2 = (je-js)+1 + 2*(NGHOST);
  int nx3 = (ke-ks)+1 + 2*(NGHOST);
  a1.NewAthenaArray(nx3,nx2,nx1);
  a2.NewAthenaArray(nx3,nx2,nx1);
  a3.NewAthenaArray(nx3,nx2,nx1);

  int level=loc.level;
  // Initialize components of the vector potential
  if (block_size.nx3 > 1) {
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          if ((pbval->nblevel[1][0][1]>level && j==js)
           || (pbval->nblevel[1][2][1]>level && j==je+1)
           || (pbval->nblevel[0][1][1]>level && k==ks)
           || (pbval->nblevel[2][1][1]>level && k==ke+1)
           || (pbval->nblevel[0][0][1]>level && j==js   && k==ks)
           || (pbval->nblevel[0][2][1]>level && j==je+1 && k==ks)
           || (pbval->nblevel[2][0][1]>level && j==js   && k==ke+1)
           || (pbval->nblevel[2][2][1]>level && j==je+1 && k==ke+1)) {
            Real x1l = pcoord->x1f(i)+0.25*pcoord->dx1f(i);
            Real x1r = pcoord->x1f(i)+0.75*pcoord->dx1f(i);
            a1(k,j,i) = 0.5*(A1(x1l, pcoord->x2f(j), pcoord->x3f(k)) +
                             A1(x1r, pcoord->x2f(j), pcoord->x3f(k)));
          } else {
            // hardcode ang_2=0 initialization of <A>
            if (i != ie+1)
              a1(k,j,i) = A1_ave(pcoord->x1f(i), pcoord-> x1f(i+1), pcoord->x2f(j),
                                 pcoord->x3f(k));
            //a1(k,j,i) = A1(pcoord->x1v(i), pcoord->x2f(j), pcoord->x3f(k));
          }

          if ((pbval->nblevel[1][1][0]>level && i==is)
           || (pbval->nblevel[1][1][2]>level && i==ie+1)
           || (pbval->nblevel[0][1][1]>level && k==ks)
           || (pbval->nblevel[2][1][1]>level && k==ke+1)
           || (pbval->nblevel[0][1][0]>level && i==is   && k==ks)
           || (pbval->nblevel[0][1][2]>level && i==ie+1 && k==ks)
           || (pbval->nblevel[2][1][0]>level && i==is   && k==ke+1)
           || (pbval->nblevel[2][1][2]>level && i==ie+1 && k==ke+1)) {
            Real x2l = pcoord->x2f(j)+0.25*pcoord->dx2f(j);
            Real x2r = pcoord->x2f(j)+0.75*pcoord->dx2f(j);
            a2(k,j,i) = 0.5*(A2(pcoord->x1f(i), x2l, pcoord->x3f(k)) +
                             A2(pcoord->x1f(i), x2r, pcoord->x3f(k)));
          } else {
            // hardcode ang_2=0 initialization of <A>
            if (j != je+1)
              a2(k,j,i) = A2_ave(pcoord->x1f(i), pcoord->x2f(j), pcoord->x2f(j+1),
                                 pcoord->x3f(k));
            //a2(k,j,i) = A2(pcoord->x1f(i), pcoord->x2v(j), pcoord->x3f(k));
          }

          if ((pbval->nblevel[1][1][0]>level && i==is)
           || (pbval->nblevel[1][1][2]>level && i==ie+1)
           || (pbval->nblevel[1][0][1]>level && j==js)
           || (pbval->nblevel[1][2][1]>level && j==je+1)
           || (pbval->nblevel[1][0][0]>level && i==is   && j==js)
           || (pbval->nblevel[1][0][2]>level && i==ie+1 && j==js)
           || (pbval->nblevel[1][2][0]>level && i==is   && j==je+1)
           || (pbval->nblevel[1][2][2]>level && i==ie+1 && j==je+1)) {
            Real x3l = pcoord->x3f(k)+0.25*pcoord->dx3f(k);
            Real x3r = pcoord->x3f(k)+0.75*pcoord->dx3f(k);
            a3(k,j,i) = 0.5*(A3(pcoord->x1f(i), pcoord->x2f(j), x3l) +
                             A3(pcoord->x1f(i), pcoord->x2f(j), x3r));
          } else {
            // hardcode ang_2=0 initialization of <A>
            if (k != ke+1)
              a3(k,j,i) = A3_ave(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3f(k),
                                 pcoord->x3f(k+1));
            //a3(k,j,i) = A3(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3v(k));
          }
        }
      }
    }
  } else { // 2D or 1D:
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          if (i != ie+1)
            a1(k,j,i) = A1_ave(pcoord->x1f(i), pcoord-> x1f(i+1), pcoord->x2f(j),
                               pcoord->x3f(k));
          if (j != je+1)
            a2(k,j,i) = A2_ave(pcoord->x1f(i), pcoord->x2f(j), pcoord->x2f(j+1),
                               pcoord->x3f(k));
          if (k != ke+1)
            a3(k,j,i) = A3_ave(pcoord->x1f(i), pcoord->x2f(j), pcoord->x3f(k),
                               pcoord->x3f(k+1));
        }
      }
    }
  }

  // Initialize interface fields
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      pfield->b.x1f(k,j,i) = (a3(k  ,j+1,i) - a3(k,j,i))/pcoord->dx2f(j) -
                             (a2(k+1,j  ,i) - a2(k,j,i))/pcoord->dx3f(k);
    }
  }}

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x2f(k,j,i) = (a1(k+1,j,i  ) - a1(k,j,i))/pcoord->dx3f(k) -
                             (a3(k  ,j,i+1) - a3(k,j,i))/pcoord->dx1f(i);
    }
  }}

  for (int k=ks; k<=ke+1; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
     pfield->b.x3f(k,j,i) = (a2(k,j  ,i+1) - a2(k,j,i))/pcoord->dx1f(i) -
                            (a1(k,j+1,i  ) - a1(k,j,i))/pcoord->dx2f(j);
    }
  }}
  a1.DeleteAthenaArray();
  a2.DeleteAthenaArray();
  a3.DeleteAthenaArray();

  // Now initialize rest of the cell centered quantities
  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      Real x = cos_a2*(pcoord->x1v(i)*cos_a3 + pcoord->x2v(j)*sin_a3) +
               pcoord->x3v(k)*sin_a2;
      Real sn = sin(k_par*x);
      Real cs = fac*cos(k_par*x);

      phydro->u(IDN,k,j,i) = den;

      Real mx = den*v_par;
      Real my = -fac*den*v_perp*sn;
      Real mz = -fac*den*v_perp*cs;

      phydro->u(IM1,k,j,i) = mx*cos_a2*cos_a3 - my*sin_a3 - mz*sin_a2*cos_a3;
      phydro->u(IM2,k,j,i) = mx*cos_a2*sin_a3 + my*cos_a3 - mz*sin_a2*sin_a3;
      phydro->u(IM3,k,j,i) = mx*sin_a2                    + mz*cos_a2;

      if (NON_BAROTROPIC_EOS) {
        phydro->u(IEN,k,j,i) = pres/gm1 + 0.5*(b_par*b_par + b_perp*b_perp) +
            // TODO(kfelker): evaluate the impact of this change for 2nd order
            // 0.5*(SQR(0.5*(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1))) +
            //      SQR(0.5*(pfield->b.x2f(k,j,i) + pfield->b.x2f(k,j+1,i))) +
            //      SQR(0.5*(pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i)))) +
            (0.5/den)*(SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i)) +
                       SQR(phydro->u(IM3,k,j,i)));
      }
    }
  }}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn static Real A1(const Real x1,const Real x2,const Real x3)
//  \brief A1: 1-component of vector potential, using a gauge such that Ax = 0, and Ay,
//  Az are functions of x and y alone.

static Real A1(const Real x1, const Real x2, const Real x3) {
  Real x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real Ay = fac*(b_perp/k_par)*sin(k_par*(x));
  Real Az = (b_perp/k_par)*cos(k_par*(x)) + b_par*y;

  return -Ay*sin_a3 - Az*sin_a2*cos_a3;
}

//----------------------------------------------------------------------------------------
//! \fn static Real A2(const Real x1,const Real x2,const Real x3)
//  \brief A2: 2-component of vector potential

static Real A2(const Real x1, const Real x2, const Real x3) {
  Real x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real Ay = fac*(b_perp/k_par)*sin(k_par*(x));
  Real Az = (b_perp/k_par)*cos(k_par*(x)) + b_par*y;

  return Ay*cos_a3 - Az*sin_a2*sin_a3;
}

//----------------------------------------------------------------------------------------
//! \fn static Real A3(const Real x1,const Real x2,const Real x3)
//  \brief A3: 3-component of vector potential

static Real A3(const Real x1, const Real x2, const Real x3) {
  Real x =  x1*cos_a2*cos_a3 + x2*cos_a2*sin_a3 + x3*sin_a2;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real Az = (b_perp/k_par)*cos(k_par*(x)) + b_par*y;

  return Az*cos_a2;
}

//----------------------------------------------------------------------------------------
//! \fn static Real A1_ave(const Real x1,const Real x2,const Real x3)
//  \brief A1_ave: 1-component of vector potential averaged along x1 edge
// TODO(kfelker): consider switching dx1f to pco->dx1f, etc

static Real A1_ave(const Real x1f, const Real x1f_ip1, const Real x2, const Real x3) {
  Real x =  cos_a2*cos_a3*x1f + cos_a2*sin_a3*x2 + sin_a2*x3;
  Real x_ip1 =  cos_a2*cos_a3*x1f_ip1 + cos_a2*sin_a3*x2 + sin_a2*x3;
  Real dx1f = x1f_ip1 - x1f;
  Real Ay, Az;

  if (cos_a3 != 0.0) {
	Ay =  (fac*b_perp/(SQR(k_par)*cos_a3))*(-cos(k_par*x_ip1) + cos(k_par*x));
  } else { // vertical (+x2) coordinate aligned wave--- Ay uniform on the x1 edge
	Ay = dx1f*((fac*b_perp/k_par)*sin(k_par*(x))); // x  = x_ip1
  }
  if (cos_a2 != 0.0) {
    // cancelling cos_a3 factor
    Az = b_perp/(SQR(k_par)*cos_a2)*(sin(k_par*x_ip1) - sin(k_par*x));
  } else { // polar (+x3) coordinate aligned wave: only linear term in Az changes along x1
    Az = dx1f*((cos_a3*b_perp/k_par)*cos(k_par*x)); // x = x3
  }
  Az += b_par*cos_a3*(-sin_a3*0.5*(SQR(x1f_ip1) - SQR(x1f)) + dx1f*cos_a3*x2);
  return (-Ay*sin_a3 - Az*sin_a2)/dx1f;
}

//----------------------------------------------------------------------------------------
//! \fn static Real A2(const Real x1,const Real x2,const Real x3)
//  \brief A2_ave: 2-component of vector potential averaged along x2 edge

static Real A2_ave(const Real x1, const Real x2f, const Real x2f_jp1, const Real x3) {
  Real x =  cos_a2*cos_a3*x1 + cos_a2*sin_a3*x2f + sin_a2*x3;
  Real x_jp1 =  cos_a2*cos_a3*x1 + cos_a2*sin_a3*x2f_jp1 + sin_a2*x3;
  Real dx2f = x2f_jp1 - x2f;
  Real Ay, Az;

  if (sin_a3 != 0.0) {
	Ay =  (fac*b_perp/(SQR(k_par)*sin_a3))*(-cos(k_par*x_jp1) + cos(k_par*x));
  } else { // horizontal (+x1) coordinate aligned wave--- Ay uniform on the x2 edge
	Ay = dx2f*((fac*b_perp/k_par)*sin(k_par*(x))); // x = x_jp1
  }
  if (cos_a2 != 0.0) {
    // cancelling sin_a3 factor
    Az = b_perp/(SQR(k_par)*cos_a2)*(sin(k_par*x_jp1) - sin(k_par*x));
  } else { // polar (+x3) coordinate aligned wave: only linear term in Az changes along x2
    Az = dx2f*((sin_a3*b_perp/k_par)*cos(k_par*x)); // x = x3
  }
  Az += b_par*cos_a3*(-dx2f*sin_a3*x1 + cos_a3*0.5*(SQR(x2f_jp1) - SQR(x2f)));
  return (Ay*cos_a3 - Az*sin_a2)/dx2f;
}

//----------------------------------------------------------------------------------------
//! \fn static Real A3_ave(const Real x1,const Real x2,const Real x3f,const Real x3f_kp1)
//  \brief A3_ave: 3-component of vector potential averaged along x3 edge

static Real A3_ave(const Real x1, const Real x2, const Real x3f, const Real x3f_kp1) {
  Real x =  cos_a2*cos_a3*x1 + cos_a2*sin_a3*x2 + sin_a2*x3f;
  Real x_kp1 =  cos_a2*cos_a3*x1 + cos_a2*sin_a3*x2 + sin_a2*x3f_kp1;
  Real y = -x1*sin_a3        + x2*cos_a3;
  Real dx3f = x3f_kp1 - x3f;
  Real Az;

  if (sin_a2 != 0.0) {
    Az = cos_a2*b_perp/(SQR(k_par)*sin_a2)*(sin(k_par*x_kp1) - sin(k_par*x));
  } else { // wave propagates in x1-x2 plane
    Az = dx3f*((b_perp/k_par)*cos(k_par*(x)) + b_par*y); // x, y do not depend on x3
  }
  return Az/dx3f;
}
