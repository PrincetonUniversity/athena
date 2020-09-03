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
//- ipert = 3 - J&G vortical shwave (hydro test)
//- ipert = 4 - nonlinear density wave test of Fromang & Papaloizou
//- ipert = 5 - 2nd MHD shwave test of JGG (2008) -- their figure 9
//- ipert = 6 - 3rd MHD shwave test of JGG (2008) -- their figure 11
//- ipert = 7 - nonlinear shearing wave test of Heinemann & Papaloizou (2008)
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
#include "../parameter_input.hpp"
#include "../utils/utils.hpp" // ran2()
#include "../scalars/scalars.hpp" //added scalars

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif
#if !SHEARING_BOX
#error "This problem generator requires shearing box"
#endif
#if NSCALARS!=2
#error "This problem generator requires 2 passive scalars"
#endif

namespace {
Real HistoryBxBy(MeshBlock *pmb, int iout);
Real HistorydVxVy(MeshBlock *pmb, int iout);
Real HistorydBy(MeshBlock *pmb, int iout);
// void ScalarSource(MeshBlock *pmb, const Real time, const Real dt, 
//               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
//               AthenaArray<Real> &cons);
Real Lx, Ly, Lz; // root grid size, global to share with output functions
Real Omega_0, qshear,tscal;
Real s_0_layer1,s_0_layer2,s_0_bkgd,s_1_layer1,s_1_layer2,s_1_bkgd;
Real layer_bound1,layer_bound2,layer_bound3;
} // namespace

// ===================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  AllocateUserHistoryOutput(3);
  EnrollUserHistoryOutput(0, HistoryBxBy, "-BxBy");
  EnrollUserHistoryOutput(1, HistorydVxVy, "dVxVy");
  EnrollUserHistoryOutput(2, HistorydBy, "dBy");
  // Read problem parameters
  Omega_0 = pin->GetOrAddReal("problem","Omega0",1.0);
  qshear  = pin->GetOrAddReal("problem","qshear",1.5);
  tscal = pin->GetOrAddReal("problem","tscal",0.0);
  static const std::string species_names[2] = {"C12","Mg24"};
  s_0_layer1 = pin->GetOrAddReal("problem","s_init_layer1_"+species_names[0], -1);
  s_0_layer2 = pin->GetOrAddReal("problem","s_init_layer2_"+species_names[0], -1);
  s_0_bkgd = pin->GetOrAddReal("problem","s_init_bkgd_"+species_names[0], -1);
  s_1_layer1 = pin->GetOrAddReal("problem","s_init_layer1_"+species_names[1], -1);
  s_1_layer2 = pin->GetOrAddReal("problem","s_init_layer2_"+species_names[1], -1);
  s_1_bkgd = pin->GetOrAddReal("problem","s_init_bkgd_"+species_names[1], -1);
  layer_bound1 = pin->GetOrAddReal("problem","layer_bound1", 0.6);
  layer_bound2 = pin->GetOrAddReal("problem","layer_bound2", 0.7);
  layer_bound3 = pin->GetOrAddReal("problem","layer_bound3", 0.8);
  // Enroll source funtion for scalars
  // EnrollUserExplicitSourceFunction(ScalarSource);
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
    std::cout << "[hgb.cpp]: HGB only works on a 2D or 3D grid" << std::endl;
  }

  // Read problem parameters for initial conditions
  Real amp = pin->GetReal("problem","amp");
  int ipert = pin->GetOrAddInteger("problem","ipert", 1);

  Real beta, dir_sgn;
  int ifield, Bdir;
  if (MAGNETIC_FIELDS_ENABLED) {
    beta = pin->GetReal("problem","beta");
    ifield = pin->GetOrAddInteger("problem","ifield", 1);
    // For net-flux calc, provide the direction of the B field
    Bdir = pin->GetOrAddInteger("problem","Bdir",1);
    if (Bdir > 0)
      dir_sgn = 1.0;
    else
      dir_sgn = -1.0;
  }

  // get initial scalar value
  const Real s_init = pin->GetOrAddReal("problem", "s_init", 0.);

  // Compute pressure based on the EOS.
  Real den = pin->GetReal("problem","den");
  Real iso_cs = peos->GetIsoSoundSpeed();
  Real pres =1.0, gamma=1.0;
  if (NON_BAROTROPIC_EOS) {
    gamma = peos->GetGamma();
    pres = pin->GetReal("problem","pres");
  } else {
    iso_cs =peos->GetIsoSoundSpeed();
    pres = den*SQR(iso_cs);
  }
  // Compute field strength based on beta.
  Real B0  = 0.0;
  if (MAGNETIC_FIELDS_ENABLED)
    B0 = std::sqrt(static_cast<Real>(2.0*pres/beta));

  // Ensure a different initial random seed for each meshblock.
  std::int64_t iseed = -1 - gid;

  // Initialize boxsize
  Lx = pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min;
  Ly = pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min;
  Lz = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;

  // initialize wavenumbers
  int nwx = pin->GetOrAddInteger("problem","nwx",1);
  int nwy = pin->GetOrAddInteger("problem","nwy",1);
  int nwz = pin->GetOrAddInteger("problem","nwz",1);
  Real kx = (TWO_PI/Lx)*(static_cast<Real>(nwx));// nxw=-ve for leading wave
  Real ky = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/Lz)*(static_cast<Real>(nwz));

  // For PF density wave test, read data from file: not implemented yet.


  // Rescale amp to sound speed for ipert 2,3
  if (NON_BAROTROPIC_EOS) {
    if (ipert == 2 || ipert == 3)
      amp *= std::sqrt(gamma*pres/den);
  } else {
    if (ipert == 2 || ipert == 3)
      amp *= iso_cs;
  }

  Real x1, x2, x3;  // xmin, xmax;
  Real x1f, x2f, x3f;
  Real rd(0.0), rp(0.0), rvx(0.0), rvy(0.0), rvz(0.0), rbx(0.0), rby(0.0), rbz(0.0);
  Real rval;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        x1 = pcoord->x1v(i);
        x2 = pcoord->x2v(j);
        x3 = pcoord->x3v(k);
        x1f = pcoord->x1f(i);
        x2f = pcoord->x2f(j);
        x3f = pcoord->x3f(k);

        //Initialize perturbations
        // ipert = 1 - random perturbations to P and V [default, used by HGB]
        // ipert = 2 - uniform Vx=amp (epicyclic wave test)
        // ipert = 3 - vortical shwave (hydro test)
        // ipert = 4 - Fromang & Papaloizou nonlinear density wave (hydro test)
        // ipert = 5 & 6 - JGG MHD shwave tests
        // ipert = 7 - Heinemann & Papaloizou (2008) nonlinear shwave (hydro test)
        if (ipert == 1) {
          rval = amp*(ran2(&iseed) - 0.5);
          if (NON_BAROTROPIC_EOS) {
            rp = pres*(1.0 + 2.0*rval);
            rd = den;
          } else {
            rd = den; //den*(1.0 + 2.0*rval);
          }
          // Follow HGB: the perturbations to V/Cs are
          // (1/5)amp/std::sqrt(gamma)
          rval = amp*(ran2(&iseed) - 0.5);
          rvx = (0.4/std::sqrt(3.0)) *rval*1e-3/std::sqrt(gamma);
          //rvx = 0.4*rval*std::sqrt(pres/den);
          SumRvx += rvx;

          rval = amp*(ran2(&iseed) - 0.5);
          rvy = (0.4/std::sqrt(3.0)) *rval*1e-3/std::sqrt(gamma);
          //rvy = 0.4*rval*std::sqrt(pres/den);
          SumRvy += rvy;

          rval = amp*(ran2(&iseed) - 0.5);
          rvz = (0.4/std::sqrt(3.0)) *rval*1e-3/std::sqrt(gamma);
          //rvz = 0.4*rval*std::sqrt(pres/den);
          SumRvz += rvz;
        }

        // Initialize (d, M, P)
        // for_the_future: if FARGO do not initialize the bg shear
        phydro->u(IDN,k,j,i) = rd;
        phydro->u(IM1,k,j,i) = rd*rvx;
        phydro->u(IM2,k,j,i) = rd*rvy;
        phydro->u(IM2,k,j,i) -= rd*(qshear*Omega_0*x1);
        phydro->u(IM3,k,j,i) = rd*rvz;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = rp/(gamma-1.0)
                                 + 0.5*(SQR(phydro->u(IM1,k,j,i))
                                        + SQR(phydro->u(IM2,k,j,i))
                                        + SQR(phydro->u(IM3,k,j,i)))/rd;
        } // Hydro

        // Initialize b.  For 3D shearing box B1=Bx, B2=By, B3=Bz
        // ifield = 0 - used with ipert=5 or 6
        // ifield = 1 - Bz=B0std::sin(x1) field with zero-net-flux[default]
        // ifield = 2 - uniform Bz
        // ifield = 3 - B=(0,B0std::cos(kx*x1),B0std::sin(kx*x1))=zero-net flux w helicity
        // ifield = 4 - B=(0,B0/std::sqrt(2),B0/std::sqrt(2))= net toroidal+vertical field
        if (MAGNETIC_FIELDS_ENABLED) {
          if (ifield == 1) {
            pfield->b.x1f(k,j,i) = 0.0;
            pfield->b.x2f(k,j,i) = 0.0;
            pfield->b.x3f(k,j,i) = B0*(std::sin(static_cast<Real>(kx)*x1));
            if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
            if (j==je) pfield->b.x2f(k,je+1,i) = 0.0;
            if (k==ke) pfield->b.x3f(ke+1,j,i) = B0*(std::sin(static_cast<Real>(kx)*x1));
          }
          if (ifield == 2) {
            pfield->b.x1f(k,j,i) = 0.0;
            pfield->b.x2f(k,j,i) = 0.0;
            pfield->b.x3f(k,j,i) = B0*dir_sgn;
            if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
            if (j==je) pfield->b.x2f(k,je+1,i) = 0.0;
            if (k==ke) pfield->b.x3f(ke+1,j,i) = B0*dir_sgn;
          }
          if (ifield == 5) {
            pfield->b.x1f(k,j,i) = 0.0;
            pfield->b.x2f(k,j,i) = B0;
            pfield->b.x3f(k,j,i) = 0.0;
            if (i==ie) pfield->b.x1f(k,j,ie+1) = 0.0;
            if (j==je) pfield->b.x2f(k,je+1,i) = B0;
            if (k==ke) pfield->b.x3f(ke+1,j,i) = 0.0;
          }
        } // MHD
        // initialize scalars
        if (NSCALARS > 0) {
          for (int ispec=0; ispec < NSCALARS; ++ispec) {
            pscalars->s(ispec, k, j, i) = s_init*rd;
          }
        }
      }
    }
  }


  // For random perturbations as in HGB, ensure net momentum is zero by
  // subtracting off mean of perturbations
  if (ipert == 1) {
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
  }


  return;
}

void MeshBlock::UserWorkInLoop() {
  if ((pmy_mesh->time >= tscal)&&(pmy_mesh->time < tscal+(pmy_mesh->dt))) {
//    printf("ACTIVATING SCALARS from %.2e, %.2e, %.2e, %.2e, %.2e, \n", layer_bound1, layer_bound3,
//      s_0_layer1, s_0_layer2, s_0_bkgd);
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real& u_d  = phydro->u(IDN,k,j,i);
        Real x1 = pcoord->x1v(i);
        // initialize scalars
          pscalars->s(0, k, j, i) = s_0_bkgd*u_d;
          pscalars->s(1, k, j, i) = s_1_bkgd*u_d;
          if ((x1 < layer_bound2*pmy_mesh->mesh_size.x1max) && (x1 > layer_bound1*pmy_mesh->mesh_size.x1max)){
            pscalars->s(0, k, j, i) = s_0_layer1*u_d;
            pscalars->s(1, k, j, i) = s_1_layer1*u_d;
          }
          else if ((x1 < layer_bound3*pmy_mesh->mesh_size.x1max) && (x1 > layer_bound2*pmy_mesh->mesh_size.x1max)){
              pscalars->s(0, k, j, i) = s_0_layer2*u_d;
              pscalars->s(1, k, j, i) = s_1_layer2*u_d;
          }
        }
      }
      }
  }
  return;
}

// void ScalarSource(MeshBlock *pmb, const Real time, const Real dt,const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
//               AthenaArray<Real> &cons)
// {
//   if (time > tscal && time < tscal + dt) {
//     int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
//     for (int k=ks; k<=ke; k++) {
//       for (int j=js; j<=je; j++) {
//         for (int i=is; i<=ie; i++) {
//           Real& u_d  = pmb->phydro->u(IDN,k,j,i);
//           Real x1 = pmb->pcoord->x1v(i);
//           static const std::string species_names[2] = {"C12","Mg24"};
//           pmb->pscalars->s(0, k, j, i) = s_0_bkgd*u_d;
//           pmb->pscalars->s(1, k, j, i) = s_1_bkgd*u_d;
//           if ((x1 < layer_bound2*pmb->pmy_mesh->mesh_size.x1max) && 
//             (x1 > layer_bound1*pmb->pmy_mesh->mesh_size.x1max)){
//             pmb->pscalars->s(0, k, j, i) = s_0_layer1*u_d;
//             pmb->pscalars->s(1, k, j, i) = s_1_layer1*u_d;
//           }
//           else if ((x1 < layer_bound3*pmb->pmy_mesh->mesh_size.x1max) && 
//             (x1 > layer_bound2*pmb->pmy_mesh->mesh_size.x1max)){
//               pmb->pscalars->s(0, k, j, i) = s_0_layer2*u_d;
//               pmb->pscalars->s(1, k, j, i) = s_1_layer2*u_d;
//           }
//         }
//       }
//     }
//   }
//   return;
// }

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
        vshear = qshear*Omega_0*pmb->pcoord->x1v(i);
        dvxvy += volume(i)*w(IDN,k,j,i)*w(IVX,k,j,i)*(w(IVY,k,j,i)+vshear);
      }
    }
  }
  return dvxvy;
}

Real HistorydBy(MeshBlock *pmb, int iout) {
  Real dby = 0;
  Real fkx, fky, fkz; // Fourier kx, ky
  Real x1, x2, x3;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> &b = pmb->pfield->bcc;

  fky = TWO_PI/Ly;
  fkx = -4.0*PI/Lx + qshear*Omega_0*fky*pmb->pmy_mesh->time;
  fkz = TWO_PI/Lz;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        x1 = pmb->pcoord->x1v(i);
        x2 = pmb->pcoord->x2v(j);
        x3 = pmb->pcoord->x3v(k);
        dby += (2.0
                * volume(i)
                * (b(IB2, k, j, i) - (0.2-0.15*Omega_0*pmb->pmy_mesh->time))
                * std::cos(fkx*x1 + fky*x2 + fkz*x3));
      }
    }
  }
  return dby;
}
} // namespace

