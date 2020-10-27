//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hgb.cpp
//
//  \brief Problem generator for 3D shearing sheet.
//
// PURPOSE:  Problem generator for 3D shearing sheet.  Based on the initial
//   conditions described in "Local Three-dimensional Magnetohydrodynamic
//   Simulations of Accretion Disks" by Hawley, Gammie & Balbus, or HGB.
//
// Several different field configurations and perturbations are possible:
//
//- ifield = 0 - uses field set by choice of ipert flag
//- ifield = 1 - Bz=B0std::sin(kx*x1) field with zero-net-flux [default] (kx input)
//- ifield = 2 - uniform Bz
//- ifield = 3 - B=(0,B0std::cos(kx*x1),B0std::sin(kx*x1))= zero-net flux w helicity
//- ifield = 4 - B=(0,B0/std::sqrt(2),B0/std::sqrt(2))= net toroidal+vertical field
//- ifield = 5 - uniform By
//
//- ipert = 1 - random perturbations to P and V [default, used by HGB]
//- ipert = 2 - uniform Vx=amp (epicyclic wave test)
//- ipert = 3 - nonlinear density wave test of Fromang & Papaloizou (2007)
//- ipert = 4 - nonlinear shearing wave test of Heinemann & Papaloizou (2009)
//- ipert = 5 - MHD shwave test of JGG (2008) -- their figure 11
//
// To run simulations of stratified disks (including vertical gravity), use the
// strat.c problem generator.
//
// Code must be configured using -shear
//
// REFERENCE: Hawley, J. F. & Balbus, S. A., ApJ 400, 595-609 (1992).
//            Johnson, Guan, & Gammie, ApJSupp, (2008)
//============================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
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
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp" // ran2()

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

namespace {
Real iso_cs, gm1, d0, p0;
Real nwx, nwy, nwz; // Wavenumbers
Real Lx, Ly, Lz; // root grid size, global to share with output functions
Real Omega_0, qshear;
Real amp, beta;
int ipert, shboxcoord;

Real HistoryBxBy(MeshBlock *pmb, int iout);
Real HistorydVxVy(MeshBlock *pmb, int iout);
Real HistorydBy(MeshBlock *pmb, int iout);
} // namespace

// ===================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  ipert = pin->GetOrAddInteger("problem","ipert", 1);
  if (ipert == 5) { // ipert = 5
    AllocateUserHistoryOutput(3);
    EnrollUserHistoryOutput(0, HistoryBxBy, "-BxBy");
    EnrollUserHistoryOutput(1, HistorydVxVy, "dVxVy");
    EnrollUserHistoryOutput(2, HistorydBy, "dBy");
  } else if (ipert > 0 || ipert < 4) { // ipert = 1-4
    AllocateUserHistoryOutput(2);
    EnrollUserHistoryOutput(0, HistoryBxBy, "-BxBy");
    EnrollUserHistoryOutput(1, HistorydVxVy, "dVxVy");
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
        << "ipert should be chosen from 1 to 5." << std::endl;
    ATHENA_ERROR(msg);
  }

  if (!shear_periodic) {
    std::stringstream msg;
    msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box" << std::endl;
    ATHENA_ERROR(msg);
  }

  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief mhd shearing waves and unstratified disk problem generator for
//  3D problems.
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real SumRvx=0.0, SumRvy=0.0, SumRvz=0.0;
  if (pmy_mesh->mesh_size.nx2 == 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
        << "This problem generator works only in 2D or 3D." << std::endl;
    ATHENA_ERROR(msg);
  }

  // gamma, press, sound speed
  Real gamma = 1.0;
  d0    = pin->GetOrAddReal("problem","d0", 1.0);
  if (NON_BAROTROPIC_EOS) {
    p0    = pin->GetReal("problem","p0");
    gamma  = peos->GetGamma();
    iso_cs = std::sqrt(gamma*p0/d0);
  } else {
    iso_cs = peos->GetIsoSoundSpeed();
    p0 = d0*SQR(iso_cs);
  }

  // shearing box parameter
  shboxcoord = pin->GetOrAddInteger("problem","shboxcoord",1);
  if (shboxcoord != 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box in x-y plane." << std::endl;
    ATHENA_ERROR(msg);
  }
  Omega_0 = porb->Omega0;
  qshear  = porb->qshear;

  // Read problem parameters for initial conditions
  Real amp = pin->GetReal("problem","amp");

  Real dir_sgn;
  int ifield, Bdir;
  beta    = pin->GetReal("problem","beta");
  ifield  = pin->GetOrAddInteger("problem","ifield", 1);
  if ((ipert < 5)
      && (ifield <1 || ifield >5)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
        << "ifield should be from 1 to 5 for ipert = 1-4." << std::endl;
    ATHENA_ERROR(msg);
  }
  if ((ipert == 2) && (ifield != 0)) {
    std::stringstream msg;
    msg << "### WARNING in hb3.cpp ProblemGenerator" << std::endl
        << "For ipert = 5, ifield should be set at 0." << std::endl;
    ifield = 0;
  }

  // For net-flux calc, provide the direction of the B field
  Bdir    = pin->GetOrAddInteger("problem","Bdir",1);
  if (Bdir > 0)
    dir_sgn = 1.0;
  else
    dir_sgn = -1.0;

  // Compute pressure based on the EOS.
  if (NON_BAROTROPIC_EOS) {
    gamma  = peos->GetGamma();
    p0     = pin->GetReal("problem","p0");
  } else {
    iso_cs = peos->GetIsoSoundSpeed();
    p0   = d0*SQR(iso_cs);
  }
  // Compute field strength based on beta.
  Real B0  = std::sqrt(static_cast<Real>(2.0*p0/beta));

  // Ensure a different initial random seed for each meshblock.
  std::int64_t iseed = -1 - gid;

  // Initialize boxsize
  Lx = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
  Ly = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
  Lz = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;

  // initialize wavenumbers
  nwx = pin->GetOrAddInteger("problem","nwx",1);
  nwy = pin->GetOrAddInteger("problem","nwy",1);
  nwz = pin->GetOrAddInteger("problem","nwz",1);
  Real kx = (TWO_PI/Lx)*(static_cast<Real>(nwx));// nxw=-ve for leading wave
  Real ky = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/Lz)*(static_cast<Real>(nwz));
  if (pmy_mesh->mesh_size.nx3 == 1) kz = 0.0;

  // Rescale amp to sound speed for ipert 2
  if (ipert == 2) {
    if (NON_BAROTROPIC_EOS)
      amp *= std::sqrt(gamma*p0/d0);
    else
      amp *= iso_cs;
  }

  // Initialize perturbations
  // ipert = 1 - random perturbations to P and V [default, used by HGB]
  // ipert = 2 - uniform Vx=amp (epicyclic wave test)
  // ipert = 3 - nonlinear density wave test of Fromang & Papaloizou (2007)
  // ipert = 4 - nonlinear shearing wave test of Heinemann & Papaloizou (2009)
  // ipert = 5 - MHD shwave test of JGG (2008) -- their figure 11
  if (ipert < 5) { // ipert = 1-4
    Real x1;
    // hydro
    if (ipert == 1) {
      Real rd(0.0), rp(0.0), rvx(0.0), rvy(0.0), rvz(0.0);
      Real rval;
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            x1  = pcoord->x1v(i);
            // In HGB, amp = 2.5e-2
            rval = amp*(ran2(&iseed) - 0.5);
            if (NON_BAROTROPIC_EOS) {
              rp = p0*(1.0 + 2.0*rval);
              rd = d0;
            } else {
              rd = d0; //den*(1.0 + 2.0*rval);
            }
            // Following HGB: the perturbations to V/Cs are
            // (1/5)amp/std::sqrt(gamma)
            rval = amp*(ran2(&iseed) - 0.5);
            rvx = (0.4/std::sqrt(3.0)) *rval*1e-3/std::sqrt(gamma);
            //rvx  = 0.4*rval*iso_cs/std::sqrt(gamma);
            SumRvx += rvx;

            rval = amp*(ran2(&iseed) - 0.5);
            rvy = (0.4/std::sqrt(3.0)) *rval*1e-3/std::sqrt(gamma);
            //rvy  = 0.4*rval*iso_cs/std::sqrt(gamma);
            SumRvy += rvy;

            rval = amp*(ran2(&iseed) - 0.5);
            rvz = (0.4/std::sqrt(3.0)) *rval*1e-3/std::sqrt(gamma);
            //rvz  = 0.4*rval*iso_cs/std::sqrt(gamma);
            SumRvz += rvz;

            // Initialize (d, M, P)
            phydro->u(IDN,k,j,i) = rd;
            phydro->u(IM1,k,j,i) = rd*rvx;
            phydro->u(IM2,k,j,i) = rd*rvy;
            if(!porb->orbital_advection_defined)
              phydro->u(IM2,k,j,i) -= rd*qshear*Omega_0*x1;
            phydro->u(IM3,k,j,i) = rd*rvz;
            if (NON_BAROTROPIC_EOS) {
              phydro->u(IEN,k,j,i) = rp/(gamma-1.0)
                                     + 0.5*(SQR(phydro->u(IM1,k,j,i))
                                            + SQR(phydro->u(IM2,k,j,i))
                                            + SQR(phydro->u(IM3,k,j,i)))/rd;
            }
          }
        }
      }
      // For random perturbations as in HGB, ensure net momentum is zero by
      // subtracting off mean of perturbations
      int cell_num = block_size.nx1*block_size.nx2*block_size.nx3;
      SumRvx /= cell_num;
      SumRvy /= cell_num;
      SumRvz /= cell_num;
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            phydro->u(IM1,k,j,i) -= rd*SumRvx;
            phydro->u(IM2,k,j,i) -= rd*SumRvy;
            phydro->u(IM3,k,j,i) -= rd*SumRvz;
          }
        }
      }
    } else if (ipert == 2) {
      Real x1;
      Real rd(0.0), rp(0.0), rvx(0.0), rvy(0.0), rvz(0.0);
      rp = p0;
      rd = d0;
      rvx = amp;
      rvy = 0.0;
      rvz = 0.0;
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie; i++) {
            // Initialize (d, M, P)
            x1  = pcoord->x1v(i);
            phydro->u(IDN,k,j,i) = rd;
            phydro->u(IM1,k,j,i) = rd*rvx;
            phydro->u(IM2,k,j,i) = rd*rvy;
            if(!porb->orbital_advection_defined)
              phydro->u(IM2,k,j,i) -= rd*qshear*Omega_0*x1;
            phydro->u(IM3,k,j,i) = rd*rvz;
            if (NON_BAROTROPIC_EOS) {
              phydro->u(IEN,k,j,i) = rp/(gamma-1.0)
                                     + 0.5*(SQR(phydro->u(IM1,k,j,i))
                                            + SQR(phydro->u(IM2,k,j,i))
                                            + SQR(phydro->u(IM3,k,j,i)))/rd;
            }
          }
        }
      }
    } else if (ipert == 3) {
      std::stringstream msg;
      msg << "### FATAL ERROR in hgb.cpp ProblemGenerator" << std::endl
          << "ipert=3 (nonlinear density wave test of Fromang & Papaloizou)"
          << " not implemented yet!" << std::endl;
      ATHENA_ERROR(msg);
      // For FP density wave test, read data from file: not implemented yet.
    } else if (ipert == 4) {
      Real x1, x2;
      Real rd(0.0), rvx(0.0), rvy(0.0), rvz(0.0);
      if (NON_BAROTROPIC_EOS) {
        std::stringstream msg;
        msg << "### FATAL ERROR in hgb.cpp ProblemGenerator" << std::endl
            << "ipert=4 (nonlinear density wave test of Heinemann & Papaloizou)"
            << " requires the isothermal eos." << std::endl;
        ATHENA_ERROR(msg);
      } else { // isothermal
        Real kappa2 = 2.0*(2.0 - qshear)*Omega_0*Omega_0;
        Real aa = (kx*kx + ky*ky)*SQR(iso_cs) + kappa2;
        Real bb = 2.0*qshear*Omega_0*ky*iso_cs;
        Real denom = aa*aa + bb*bb;
        Real rd_hat = (ky*iso_cs*bb -2.0*Omega_0*aa)*amp/denom;
        Real px_hat = -d0*iso_cs*(ky*iso_cs*aa +2.0*Omega_0*bb)*amp/denom;
        Real py_hat = d0*(amp + ky*px_hat + (2.0-qshear)*Omega_0*rd_hat)/kx;
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              x1 = pcoord->x1v(i);
              x2 = pcoord->x2v(j);
              // Initialize (d, M, P)
              rd  = d0*(1.0 + rd_hat*std::cos(static_cast<Real>(kx*x1 + ky*x2)));
              rvx = px_hat*std::sin(static_cast<Real>(kx*x1 + ky*x2))/rd;
              rvy = py_hat*std::sin(static_cast<Real>(kx*x1 + ky*x2))/rd;
              rvz = 0.0;
              phydro->u(IDN,k,j,i) = rd;
              phydro->u(IM1,k,j,i) = rd*rvx;
              phydro->u(IM2,k,j,i) = rd*rvy;
              if(!porb->orbital_advection_defined)
                phydro->u(IM2,k,j,i) -= rd*qshear*Omega_0*x1;
              phydro->u(IM3,k,j,i) = rd*rvz;
            }
          }
        }
      }
    } // hydro

    // Initialize b.  For 3D shearing box B1=Bx, B2=By, B3=Bz
    // ifield = 1 - Bz=B0std::sin(ks*x1) field with zero-net-flux[default]
    // ifield = 2 - uniform Bz
    // ifield = 3 - B=(0,B0std::cos(kx*x1),B0std::sin(kx*x1))=zero-net flux w helicity
    // ifield = 4 - B=(0,B0/std::sqrt(2),B0/std::sqrt(2))= net toroidal+vertical field
    // ifield = 5 - uniform By
    if (ifield < 1 || ifield > 5) {
      std::stringstream msg;
      msg << "### FATAL ERROR in hb3.cpp ProblemGenerator" << std::endl
          << "ifield should be from 1 to 5 for ipert = 1-4." << std::endl;
      ATHENA_ERROR(msg);
    } else if (ifield != 3) {
      if (ifield == 1) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              x1  = pcoord->x1v(i);
              pfield->b.x1f(k,j,i) = 0.0;
              pfield->b.x2f(k,j,i) = 0.0;
              pfield->b.x3f(k,j,i) = B0*(std::sin(static_cast<Real>(kx)*x1));
              if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
              if (j==je) pfield->b.x2f(k,je+1,i) = 0.0;
              if (k==ke)
                pfield->b.x3f(ke+1,j,i) = B0*(std::sin(static_cast<Real>(kx)*x1));
            }
          }
        }
      } else if (ifield == 2) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              pfield->b.x1f(k,j,i) = 0.0;
              pfield->b.x2f(k,j,i) = 0.0;
              pfield->b.x3f(k,j,i) = B0*dir_sgn;
              if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
              if (j==je) pfield->b.x2f(k,je+1,i) = 0.0;
              if (k==ke) pfield->b.x3f(ke+1,j,i) = B0*dir_sgn;
            }
          }
        }
      } else if (ifield == 4) {
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              pfield->b.x1f(k,j,i) = 0.0;
              pfield->b.x2f(k,j,i) = B0/std::sqrt(2);
              pfield->b.x3f(k,j,i) = B0/std::sqrt(2);
              if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
              if (j==je) pfield->b.x2f(k,je+1,i) = B0/std::sqrt(2);
              if (k==ke) pfield->b.x3f(ke+1,j,i) = B0/std::sqrt(2);
            }
          }
        }
      } else { // ifield = 5
        for (int k=ks; k<=ke; k++) {
          for (int j=js; j<=je; j++) {
            for (int i=is; i<=ie; i++) {
              pfield->b.x1f(k,j,i) = 0.0;
              pfield->b.x2f(k,j,i) = B0;
              pfield->b.x3f(k,j,i) = 0.0;
              if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
              if (j==je) pfield->b.x2f(k,je+1,i) = B0;
              if (k==ke) pfield->b.x3f(ke+1,j,i) = 0.0;
            }
          }
        }
      }
    } else { // ifield = 4
      // nxN != ncellsN, in general. Allocate to extend through 2*ghost, regardless # dim
      int nx1 = block_size.nx1 + 2*NGHOST;
      int nx2 = block_size.nx2 + 2*NGHOST;
      int nx3 = block_size.nx3 + 2*NGHOST;

      // vector potential
      AthenaArray<Real> rax, ray, raz;
      rax.NewAthenaArray(nx3, nx2, nx1);
      ray.NewAthenaArray(nx3, nx2, nx1);
      raz.NewAthenaArray(nx3, nx2, nx1);

      // set rax
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie; i++) {
            rax(k,j,i) = 0.0;
          }
        }
      }

      // set ray
      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je; j++) {
          for (int i=is; i<=ie+1; i++) {
            x1 = pcoord->x1f(i);
            ray(k,j,i) = -B0*std::cos(kx*x1)/kx;
          }
        }
      }

      // set raz
      for (int k=ks; k<=ke; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie+1; i++) {
            x1 = pcoord->x1f(i);
            raz(k,j,i) = -B0*std::sin(kx*x1)/kx;
          }
        }
      }

      // calculate b from vector potential
      for (int k=ks; k<=ke  ; k++) {
        for (int j=js; j<=je  ; j++) {
          for (int i=is; i<=ie+1; i++) {
            pfield->b.x1f(k,j,i) = (raz(k,j+1,i)-raz(k,j,i))/pcoord->dx2f(j)
                                   -(ray(k+1,j,i)-ray(k,j,i))/pcoord->dx3f(k);
          }
        }
      }
      for (int k=ks; k<=ke  ; k++) {
        for (int j=js; j<=je+1; j++) {
          for (int i=is; i<=ie  ; i++) {
            pfield->b.x2f(k,j,i) = (rax(k+1,j,i)-rax(k,j,i))/pcoord->dx3f(k)
                                   -(raz(k,j,i+1)-raz(k,j,i))/pcoord->dx1f(i);
          }
        }
      }

      for (int k=ks; k<=ke+1; k++) {
        for (int j=js; j<=je  ; j++) {
          for (int i=is; i<=ie  ; i++) {
            pfield->b.x3f(k,j,i) = (ray(k,j,i+1)-ray(k,j,i))/pcoord->dx1f(i)
                                   -(rax(k,j+1,i)-rax(k,j,i))/pcoord->dx2f(j);
          }
        }
      }
      rax.DeleteAthenaArray();
      ray.DeleteAthenaArray();
      raz.DeleteAthenaArray();
    }
  } else { // ipert = 5
    Real x1, x2, x3;
    Real rd(0.0), rp(0.0), rvx(0.0), rvy(0.0), rvz(0.0);
    Real rbx(0.0), rby(0.0), rbz(0.0);
    // Following JG
    // amp (epsilon in JG) = 1.0e-6
    // Omega_0 = d0 = iso_cs = H = 1.0
    // Lx = Ly = Lz = 0.5H
    // nwx = -2, nwy = nwz = 1
    Real B02   =  0.5*SQR(B0); // the definition of beta is different in JG
    Real k2    =  SQR(kx)+SQR(ky)+SQR(kz);
    Real alpha =  std::sqrt(B02/static_cast<Real>(nwx*nwx+nwy*nwy));
    rbx   =  alpha*static_cast<Real>(nwy);
    rby   = -alpha*static_cast<Real>(nwx);
    rbz   =  0.0;

    Real sch   =  iso_cs/Omega_0; // H = 1
    Real cf1   =  std::sqrt(B02*(1.0+beta)); // va*sqrt(1+beta)
    Real cf2   =  amp*std::sqrt(sch*std::sqrt(k2*beta/(1.0+beta)));
    Real vd    =  cf1/std::sqrt(k2)*cf2;

    // hydro
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          x1 = pcoord->x1v(i);
          x2 = pcoord->x2v(j);
          x3 = pcoord->x3v(k);
          Real CS = std::cos(kx*x1+ky*x2+kz*x3);
          rd = d0*(1.0+cf2*CS);
          phydro->u(IDN,k,j,i) = rd;
          phydro->u(IM1,k,j,i) = rd*vd*kx*CS;
          phydro->u(IM2,k,j,i) = rd*vd*ky*CS;
          if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega_0*x1;
          phydro->u(IM3,k,j,i) = rd*vd*kz*CS;
        }
      }
    }

    // magnetic fields
    // Calculate magnetic fields using the vector potential.
    int nx1 = block_size.nx1 + 2*NGHOST;
    int nx2 = block_size.nx2 + 2*NGHOST;
    int nx3 = block_size.nx3 + 2*NGHOST;
    AthenaArray<Real> rax, ray, raz;
    rax.NewAthenaArray(nx3, nx2, nx1);
    ray.NewAthenaArray(nx3, nx2, nx1);
    raz.NewAthenaArray(nx3, nx2, nx1);
    // set rax
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          x1 = pcoord->x1v(i);
          x2 = pcoord->x2f(j);
          x3 = pcoord->x3f(k);
          Real temp = cf2/k2*std::sin(kx*x1+ky*x2+kz*x3);
          rax(k,j,i) = temp*(rby*kz-rbz*ky);
        }
      }
    }

    // set ray
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          x1 = pcoord->x1f(i);
          x2 = pcoord->x2v(j);
          x3 = pcoord->x3f(k);
          Real temp = cf2/k2*std::sin(kx*x1+ky*x2+kz*x3);
          ray(k,j,i) = temp*(rbz*kx-rbx*kz);
        }
      }
    }

    // set raz
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie+1; i++) {
          x1 = pcoord->x1f(i);
          x2 = pcoord->x2f(j);
          x3 = pcoord->x3v(k);
          Real temp = cf2/k2*std::sin(kx*x1+ky*x2+kz*x3);
          raz(k,j,i) = temp*(rbx*ky-rby*kx);
        }
      }
    }

    // initialize interface B
    // set bx
    for (int k=ks; k<=ke  ; k++) {
      for (int j=js; j<=je  ; j++) {
        for (int i=is; i<=ie+1; i++) {
          pfield->b.x1f(k,j,i) = rbx
                                 +(raz(k,j+1,i)-raz(k,j,i))/pcoord->dx2f(j)
                                 -(ray(k+1,j,i)-ray(k,j,i))/pcoord->dx3f(k);
        }
      }
    }

    // set by
    for (int k=ks; k<=ke  ; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie  ; i++) {
          pfield->b.x2f(k,j,i) = rby
                                 +(rax(k+1,j,i)-rax(k,j,i))/pcoord->dx3f(k)
                                 -(raz(k,j,i+1)-raz(k,j,i))/pcoord->dx1f(i);
        }
      }
    }

    // set bz
    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je  ; j++) {
        for (int i=is; i<=ie  ; i++) {
          pfield->b.x3f(k,j,i) = rbz
                                 +(ray(k,j,i+1)-ray(k,j,i))/pcoord->dx1f(i)
                                 -(rax(k,j+1,i)-rax(k,j,i))/pcoord->dx2f(j);
        }
      }
    }
    rax.DeleteAthenaArray();
    ray.DeleteAthenaArray();
    raz.DeleteAthenaArray();
  }

  return;
}

namespace {

Real HistoryBxBy(MeshBlock *pmb, int iout) {
  Real bxby = 0;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &b = pmb->pfield->bcc;
  AthenaArray<Real> volume; // 1D array of volumes
  // allocate 1D array for cell volume used in usr def history
  volume.NewAthenaArray(pmb->ncells1);

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        bxby-=volume(i)*b(IB1,k,j,i)*b(IB2,k,j,i);
      }
    }
  }
  return bxby;
}

Real HistorydVxVy(MeshBlock *pmb, int iout) {
  Real dvxvy = 0.0;
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &w = pmb->phydro->w;
  Real vshear = 0.0;
  AthenaArray<Real> volume; // 1D array of volumes
  // allocate 1D array for cell volume used in usr def history
  volume.NewAthenaArray(pmb->ncells1);

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        if(!pmb->porb->orbital_advection_defined) {
          vshear = -qshear*Omega_0*pmb->pcoord->x1v(i);
        } else {
          vshear = 0.0;
        }
        dvxvy += volume(i)*w(IDN,k,j,i)*w(IVX,k,j,i)*(w(IVY,k,j,i)+vshear);
      }
    }
  }
  return dvxvy;
}

Real HistorydBy(MeshBlock *pmb, int iout) {
  Real dby = 0;
  Real x1, x2, x3;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &b = pmb->pfield->bcc;

  Real kx = (TWO_PI/Lx)*(static_cast<Real>(nwx));
  kx += qshear*Omega_0*pmb->pmy_mesh->time*(TWO_PI/static_cast<Real>(Lx));
  Real ky = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/Lz)*(static_cast<Real>(nwz));
  if (pmb->pmy_mesh->mesh_size.nx3 == 1)
    kz = 0.0;
  Real by_analytic = -kx*std::sqrt(d0*SQR(iso_cs)/(beta*(SQR(kx)+SQR(ky))));
  Real total_volume = Lx*Ly*Lz;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        x1 = pmb->pcoord->x1v(i);
        x2 = pmb->pcoord->x2v(j);
        x3 = pmb->pcoord->x3v(k);
        Real CS = std::cos(kx*x1+ky*x2+kz*x3);
        dby += volume(i)*2.0*(b(IB2,k,j,i)-by_analytic)*CS/total_volume;
      }
    }
  }
  return dby;
}
} // namespace
