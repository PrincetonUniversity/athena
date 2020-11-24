//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file jgg.cpp
//
//  \brief Problem generator for linear MHD in shearing sheet.
//
// PURPOSE:  Problem generator for linear MHD in shearing sheet. Based on the initial
//   conditions described in "ORBITAL ADVECTION BY INTERPOLATION: A FAST AND ACCURATE
//   NUMERICAL SCHEME FOR SUPER-FAST MHD FLOWS" by Johnson, Guan, & Gammie, or JGG.
//
// Two kinds of perturbations are possible:
//
//- ipert = 1 - MHD simple shwave test of JGG -- similar to their figures 5 - 7
//- ipert = 2 - MHD compressive shwave test of JGG -- their figure 11
//
// REFERENCE: Johnson, Guan, & Gammie, ApJS, 177, 373 (2008)
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
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../outputs/outputs.hpp"
#include "../parameter_input.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

namespace {
Real iso_cs, gm1, d0, p0;
Real nwx, nwy, nwz; // Wavenumbers
Real Lx, Ly, Lz; // root grid size, global to share with output functions
Real amp, beta;
int ipert;
Real hst_dt, hst_next_time;
bool error_output;

Real HistoryBxc(MeshBlock *pmb, int iout);
Real HistoryBxs(MeshBlock *pmb, int iout);
Real HistorydBxc(MeshBlock *pmb, int iout);
Real HistorydByc(MeshBlock *pmb, int iout);
} // namespace

// ===================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (!shear_periodic) {
    std::stringstream msg;
    msg << "### FATAL ERROR in jgg.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box." << std::endl;
    ATHENA_ERROR(msg);
  }

  if (mesh_size.nx2 == 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in jgg.cpp ProblemGenerator" << std::endl
        << "This problem generator works only in 2D or 3D." << std::endl;
    ATHENA_ERROR(msg);
  }

  // Read parameters
  d0    = pin->GetOrAddReal("problem","d0", 1.0);
  if (NON_BAROTROPIC_EOS) {
    p0 = pin->GetReal("problem","p0");
    gm1 = (pin->GetReal("hydro","gamma")-1.0);
    iso_cs = std::sqrt((gm1+1.0)*p0/d0);
  } else {
    iso_cs = pin->GetReal("hydro","iso_sound_speed");
    p0 = d0*SQR(iso_cs);
    gm1 = 0.0;
  }
  ipert = pin->GetOrAddInteger("problem","ipert", 1);
  error_output = pin->GetOrAddBoolean("problem","error_output",false);
  amp = pin->GetReal("problem","amp");
  if (ipert == 1) {
    nwx = 0;
    nwy = pin->GetOrAddInteger("problem","nwy",1);
    nwz = 0;
  } else { // ipert = 2
    beta    = pin->GetReal("problem","beta");
    nwx = pin->GetOrAddInteger("problem","nwx",1);
    nwy = pin->GetOrAddInteger("problem","nwy",1);
    if (mesh_size.nx3 == 1) { // 2D
      nwz = 0;
    } else { // 3D
      nwz = pin->GetOrAddInteger("problem","nwz",1);
    }
  }
  // Initialize boxsize
  Lx = mesh_size.x1max - mesh_size.x1min;
  Ly = mesh_size.x2max - mesh_size.x2min;
  Lz = mesh_size.x3max - mesh_size.x3min;

  if (error_output) {
    if (ipert == 1) {
      // allocateDataField
      AllocateRealUserMeshDataField(1);
      ruser_mesh_data[0].NewAthenaArray(mesh_size.nx3,
                                        mesh_size.nx2, mesh_size.nx1);
      // allocate User-defined History Output
      AllocateUserHistoryOutput(2);
      EnrollUserHistoryOutput(0, HistoryBxc, "Bxc");
      EnrollUserHistoryOutput(1, HistoryBxs, "Bxs");

      // read history output timing
      InputBlock *pib = pin->pfirst_block;
      while (pib != nullptr) {
        if (pib->block_name.compare(0, 6, "output") == 0) {
          OutputParameters op;
          std::string outn = pib->block_name.substr(6);
          op.block_number = atoi(outn.c_str());
          op.block_name.assign(pib->block_name);
          op.next_time = pin->GetOrAddReal(op.block_name,"next_time", time);
          op.dt = pin->GetReal(op.block_name,"dt");
          op.file_type = pin->GetString(op.block_name,"file_type");
          if (op.file_type.compare("hst") == 0) {
            hst_dt = op.dt;
            hst_next_time = op.next_time;
          }
        }
        pib = pib->pnext;
      }
    } else if (ipert == 2) {
      // allocate User-defined History Output
      AllocateUserHistoryOutput(2);
      EnrollUserHistoryOutput(0, HistorydBxc, "dBxc");
      EnrollUserHistoryOutput(1, HistorydByc, "dByc");
    }
  }
  if ((ipert-1)*(ipert-2)!=0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in jgg.cpp ProblemGenerator" << std::endl
        << "Parameter problem/ipert should be chosen from 1  or 2." << std::endl;
    ATHENA_ERROR(msg);
  }
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief linear MHD waves in shearing box
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // shearing box parameter
  if (porb->shboxcoord != 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in jgg.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box in x-y plane." << std::endl;
    ATHENA_ERROR(msg);
  }
  Real qom = porb->Omega0*porb->qshear;

  Real kx = (TWO_PI/Lx)*(static_cast<Real>(nwx));
  Real ky = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/Lz)*(static_cast<Real>(nwz));

  // Initialize perturbations
  Real x1, x2, x3;
  Real rp(0.0);
  // ipert = 1 - MHD linear shwave test of JGG -- their figures 5 - 7
  // ipert = 2 - MHD compressive shwave test of JGG -- their figure 11
  if (ipert == 1) {
    // hydro
    Real rd = d0;
    rp = p0;
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IDN,k,j,i) = rd;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qom*pcoord->x1v(i);
          phydro->u(IM3,k,j,i) = 0.0;
        }
      }
    }
    // set b
    // initialize interface B
    // set bx
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          x2 = pcoord->x2v(j);
          pfield->b.x1f(k,j,i) = amp*std::cos(ky*x2);
        }
      }
    }
    // set by
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          pfield->b.x2f(k,j,i) = 0.0;
        }
      }
    }

    // set bz
     for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          pfield->b.x3f(k,j,i) = 0.0;
        }
      }
    }
  } else { // ipert == 2
    Real rd(0.0);
    Real rbx(0.0), rby(0.0), rbz(0.0);
    rp = p0;
    // In JGG
    // amp (epsilon in JGG)  = 1.0e-6
    // Omega_0 = d0 = iso_cs = H = 1.0
    // Lx = Ly = Lz = 0.5H
    // nwx = -2, nwy = nwz = 1
    Real B02  = static_cast<Real>(p0/beta);
    Real k2    =  SQR(kx)+SQR(ky)+SQR(kz);
    Real alpha =  std::sqrt(B02/(kx*kx+ky*ky));
    rbx   =  alpha*ky;
    rby   = -alpha*kx;
    rbz   =  0.0;

    Real sch   =  iso_cs/porb->Omega0; // H = 1
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
            phydro->u(IM2,k,j,i) -= rd*qom*x1;
          phydro->u(IM3,k,j,i) = rd*vd*kz*CS;
        }
      }
    }

    // Calculate magnetic fields using the vector potential.
    AthenaArray<Real> rax, ray, raz;
    int nz = (ncells3>1)? ncells3: 2;
    rax.NewAthenaArray(nz, ncells2, ncells1);
    ray.NewAthenaArray(nz, ncells2, ncells1);
    raz.NewAthenaArray(nz, ncells2, ncells1);
    // vector potential
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

  // initialize total energy
  if (NON_BAROTROPIC_EOS) {
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IEN,k,j,i) =
              rp/gm1 +
              0.5*(SQR(0.5*(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1))) +
                   SQR(0.5*(pfield->b.x2f(k,j,i) + pfield->b.x2f(k,j+1,i))) +
                   SQR(0.5*(pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i)))) + (0.5)*
              (SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i))
               + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkInLoop()
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void Mesh::UserWorkInLoop() {
  bool flag = false;
  // check output
  Real present_time = time + dt;
  if (ipert==1 && error_output) {
    if ((present_time < tlim) && (nlim < 0 || (ncycle + 1) < nlim)
        && (present_time > hst_next_time)) {
      flag = true;
      hst_next_time += hst_dt;
    }
    if ((present_time >= tlim) || (nlim >= 0 && (ncycle + 1) >= nlim)) {
      flag = true;
    }
  }
  // calculate dbxs, dbys
  if (flag) {
    AthenaArray<Real> &bs = ruser_mesh_data[0];
    Real kx0 = (TWO_PI/Lx)*(static_cast<Real>(nwx));
    Real ky  = (TWO_PI/Ly)*(static_cast<Real>(nwy));
    Real kz = (TWO_PI/Lz)*(static_cast<Real>(nwz));
    Real rbx(0.0), rby0(0.0), rby(0.0);

    // initialize bs
    for (int k=0; k<mesh_size.nx3; k++) {
      for (int j=0; j<mesh_size.nx2; j++) {
        for (int i=0; i<mesh_size.nx1; i++) {
          bs(0,k,j,i) = 0.0;
        }
      }
    }

    for (int bn=0; bn<nblocal; ++bn) {
      MeshBlock *pmb = my_blocks(bn);
      Real qom = pmb->porb->qshear*pmb->porb->Omega0;
      Real kx = kx0+qom*present_time*ky;
      LogicalLocation loc0;
      loc0.lx1 = pmb->loc.lx1;
      loc0.lx2 = pmb->loc.lx2;
      loc0.lx3 = pmb->loc.lx3;
      loc0.level = pmb->loc.level;
      if (loc0.level == root_level) { // root level
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int j=pmb->js; j<=pmb->je; j++) {
            for (int i=pmb->is; i<=pmb->ie; i++) {
              int ti = loc0.lx1*pmb->block_size.nx1+(i-pmb->is);
              int tj = loc0.lx2*pmb->block_size.nx2+(j-pmb->js);
              int tk = loc0.lx3*pmb->block_size.nx3+(k-pmb->ks);
              Real x1 = pmb->pcoord->x1v(i);
              Real x2 = pmb->pcoord->x2v(j);
              Real x3 = pmb->pcoord->x3v(k);
              Real vol = pmb->pcoord->dx3f(k)
                         *pmb->pcoord->dx2f(j)*pmb->pcoord->dx1f(i);
              Real SN = std::sin(kx*x1+ky*x2+kz*x3);

              Real bx = pmb->pfield->bcc(IB1,k,j,i);
              bs(0,tk,tj,ti) = 2.0*bx*vol*SN;
            }
          }
        }
      } else if (loc0.level-1 == root_level) { // level difference 1
        if (pmb->block_size.nx3==1) { // 2D
          int k = pmb->ks;
          for (int j=pmb->cjs; j<=pmb->cje; j++) {
            for (int i=pmb->cis; i<=pmb->cie; i++) {
              int ii = (i-pmb->cis)*2+pmb->is;
              int jj = (j-pmb->cjs)*2+pmb->js;
              int kk = k;
              int ti = (loc0.lx1>>1)*pmb->block_size.nx1
                       +(loc0.lx1%2)*(pmb->block_size.nx1/2)+(i-pmb->cis);
              int tj = (loc0.lx2>>1)*pmb->block_size.nx2
                       +(loc0.lx2%2)*(pmb->block_size.nx2/2)+(j-pmb->cjs);
              int tk = k;
              Real vol00 = pmb->pcoord->dx3f(kk)
                           *pmb->pcoord->dx2f(jj  )
                           *pmb->pcoord->dx1f(ii  );
              Real vol01 = pmb->pcoord->dx3f(kk)
                           *pmb->pcoord->dx2f(jj  )
                           *pmb->pcoord->dx1f(ii+1);
              Real vol10 = pmb->pcoord->dx3f(kk)
                           *pmb->pcoord->dx2f(jj+1)
                           *pmb->pcoord->dx1f(ii  );
              Real vol11 = pmb->pcoord->dx3f(kk)
                           *pmb->pcoord->dx2f(jj+1)
                           *pmb->pcoord->dx1f(ii+1);
              Real SN = std::sin(kx*pmb->pcoord->x1f(ii+1)
                                 +ky*pmb->pcoord->x2f(jj+1)
                                 +kz*pmb->pcoord->x3v(kk));

              Real bx00 = pmb->pfield->bcc(IB1,kk,jj  ,ii  );
              Real bx01 = pmb->pfield->bcc(IB1,kk,jj  ,ii+1);
              Real bx10 = pmb->pfield->bcc(IB1,kk,jj+1,ii  );
              Real bx11 = pmb->pfield->bcc(IB1,kk,jj  ,ii+1);
              Real bx_vol = bx00*vol00+bx01*vol01
                            +bx10*vol10+bx11*vol11;
              bs(0,tk,tj,ti) = 2.0*bx_vol*SN;
            }
          }
        } else { // 3D
          for (int k=pmb->cks; k<=pmb->cke; k++) {
            for (int j=pmb->cjs; j<=pmb->cje; j++) {
              for (int i=pmb->cis; i<=pmb->cie; i++) {
                int ii = (i-pmb->cis)*2+pmb->is;
                int jj = (j-pmb->cjs)*2+pmb->js;
                int kk = (k-pmb->cks)*2+pmb->ks;
                int ti = (loc0.lx1>>1)*pmb->block_size.nx1
                         +(loc0.lx1%2)*(pmb->block_size.nx1/2)+(i-pmb->cis);
                int tj = (loc0.lx2>>1)*pmb->block_size.nx2
                         +(loc0.lx2%2)*(pmb->block_size.nx2/2)+(j-pmb->cjs);
                int tk = (loc0.lx3>>1)*pmb->block_size.nx3
                         +(loc0.lx3%2)*(pmb->block_size.nx3/2)+(k-pmb->cks);
                Real vol000 = pmb->pcoord->dx3f(kk  )
                              *pmb->pcoord->dx2f(jj  )
                              *pmb->pcoord->dx1f(ii  );
                Real vol001 = pmb->pcoord->dx3f(kk  )
                              *pmb->pcoord->dx2f(jj  )
                              *pmb->pcoord->dx1f(ii+1);
                Real vol010 = pmb->pcoord->dx3f(kk  )
                              *pmb->pcoord->dx2f(jj+1)
                              *pmb->pcoord->dx1f(ii  );
                Real vol011 = pmb->pcoord->dx3f(kk  )
                              *pmb->pcoord->dx2f(jj+1)
                              *pmb->pcoord->dx1f(ii+1);
                Real vol100 = pmb->pcoord->dx3f(kk+1)
                              *pmb->pcoord->dx2f(jj  )
                              *pmb->pcoord->dx1f(ii  );
                Real vol101 = pmb->pcoord->dx3f(kk+1)
                              *pmb->pcoord->dx2f(jj  )
                              *pmb->pcoord->dx1f(ii+1);
                Real vol110 = pmb->pcoord->dx3f(kk+1)
                              *pmb->pcoord->dx2f(jj+1)
                              *pmb->pcoord->dx1f(ii  );
                Real vol111 = pmb->pcoord->dx3f(kk+1)
                              *pmb->pcoord->dx2f(jj+1)
                              *pmb->pcoord->dx1f(ii+1);
                Real SN = std::sin(kx*pmb->pcoord->x1f(ii+1)
                                   +ky*pmb->pcoord->x2f(jj+1)
                                   +kz*pmb->pcoord->x3f(kk+1));

                Real bx000 = pmb->pfield->bcc(IB1,kk  ,jj  ,ii  );
                Real bx001 = pmb->pfield->bcc(IB1,kk  ,jj  ,ii+1);
                Real bx010 = pmb->pfield->bcc(IB1,kk  ,jj+1,ii  );
                Real bx011 = pmb->pfield->bcc(IB1,kk  ,jj  ,ii+1);
                Real bx100 = pmb->pfield->bcc(IB1,kk+1,jj  ,ii  );
                Real bx101 = pmb->pfield->bcc(IB1,kk+1,jj  ,ii+1);
                Real bx110 = pmb->pfield->bcc(IB1,kk+1,jj+1,ii  );
                Real bx111 = pmb->pfield->bcc(IB1,kk+1,jj  ,ii+1);
                Real bx_vol = bx000*vol000+bx001*vol001
                              +bx010*vol010+bx011*vol011
                              +bx100*vol100+bx101*vol101
                              +bx110*vol110+bx111*vol111;
                bs(0,tk,tj,ti) = 2.0*bx_vol*SN;
              }
            }
          }
        }
      } else { // level difference 2
        std::stringstream msg;
        msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator"   << std::endl
            << "This problem prohibits level > 1 for ipert=3 "
            << "with error_output=true."  << std::endl;
          ATHENA_ERROR(msg);
      }
    } // pmb
#ifdef MPI_PARALLEL
    int ntot = mesh_size.nx3*mesh_size.nx2*mesh_size.nx1;
    if (Globals::my_rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, vs.data(), ntot, MPI_ATHENA_REAL,
                 MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
      MPI_Reduce(vs.data(), vs.data(), ntot, MPI_ATHENA_REAL,
                 MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif
  } // flag
  return;
}

namespace {

Real HistoryBxc(MeshBlock *pmb, int iout) {
  Real x1, x2, x3;
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);

  Real qom = pmb->porb->qshear*pmb->porb->Omega0;
  Real bxc = 0.0;

  AthenaArray<Real> &b = pmb->pfield->bcc;

  Real kx = (TWO_PI/Lx)
            *(static_cast<Real>(nwx)+qom*pmb->pmy_mesh->time);
  Real ky = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/Lz)*(static_cast<Real>(nwz));
  Real total_volume = Lx*Ly*Lz;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        x1 = pmb->pcoord->x1v(i);
        x2 = pmb->pcoord->x2v(j);
        x3 = pmb->pcoord->x3v(k);
        Real CS = std::cos(kx*x1+ky*x2+kz*x3);
        bxc += volume(i)*2.0*b(IB1,k,j,i)*CS;
      }
    }
  }
  return bxc/(amp*total_volume);
}

Real HistoryBxs(MeshBlock *pmb, int iout) {
  if (Globals::my_rank != 0) return 0.0;
  if (pmb->lid != 0) return 0.0;
  AthenaArray<Real> &bs = pmb->pmy_mesh->ruser_mesh_data[0];
  Real bxs = 0.0;
  Real bxs_temp;
  int nx1 = pmb->pmy_mesh->mesh_size.nx1;
  int nx2 = pmb->pmy_mesh->mesh_size.nx2;
  int nx3 = pmb->pmy_mesh->mesh_size.nx3;
  Real tvol = Lx*Ly*Ly;
  for (int k=0; k<nx3; k++) {
    for (int i=0; i<nx1; i++) {
      bxs_temp = 0.0;
      for (int j=0; j<nx2; j++) {
        bxs_temp += bs(0,k,j,i);
      }
      bxs += std::fabs(bxs_temp);
    }
  }
  return bxs/(amp*tvol);
}

Real HistorydBxc(MeshBlock *pmb, int iout) {
  Real x1, x2, x3;
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);

  Real qom = pmb->porb->qshear*pmb->porb->Omega0;
  Real sch = iso_cs/pmb->porb->Omega0;
  Real dbxc = 0.0;

  AthenaArray<Real> &b = pmb->pfield->bcc;

  Real kx0 = (TWO_PI/Lx)*(static_cast<Real>(nwx));
  Real ky  = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz  = (TWO_PI/Lz)*(static_cast<Real>(nwz));
  Real k20 = SQR(kx0)+SQR(ky)+SQR(kz);
  Real rbx  =  ky*iso_cs*std::sqrt(d0/(beta*(SQR(kx0)+SQR(ky))));
  Real kx  = kx0+qom*pmb->pmy_mesh->time*ky;
  Real total_volume = Lx*Ly*Lz;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        x1 = pmb->pcoord->x1v(i);
        x2 = pmb->pcoord->x2v(j);
        x3 = pmb->pcoord->x3v(k);
        Real CS = std::cos(kx*x1+ky*x2+kz*x3);
        dbxc += volume(i)*2.0*(b(IB1,k,j,i)-rbx)*CS;
      }
    }
  }
  Real dbxc0 = amp*rbx*std::sqrt(sch*std::sqrt(k20*beta/(1.0+beta)));
  return dbxc/(total_volume*dbxc0);
}

Real HistorydByc(MeshBlock *pmb, int iout) {
  Real x1, x2, x3;
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);

  Real qom = pmb->porb->qshear*pmb->porb->Omega0;
  Real sch = iso_cs/pmb->porb->Omega0;
  Real dbyc = 0.0;

  AthenaArray<Real> &b = pmb->pfield->bcc;

  Real kx0 = (TWO_PI/Lx)*(static_cast<Real>(nwx));
  Real ky  = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz  = (TWO_PI/Lz)*(static_cast<Real>(nwz));
  Real k20 = SQR(kx0)+SQR(ky)+SQR(kz);
  Real rbx  =  ky*iso_cs*std::sqrt(d0/(beta*(SQR(kx0)+SQR(ky))));
  Real rby0 = -kx0*iso_cs*std::sqrt(d0/(beta*(SQR(kx0)+SQR(ky))));
  Real kx  = kx0+qom*pmb->pmy_mesh->time*ky;
  Real rby = rby0-qom*pmb->pmy_mesh->time*rbx;
  Real total_volume = Lx*Ly*Lz;
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, volume);
      for (int i=is; i<=ie; i++) {
        x1 = pmb->pcoord->x1v(i);
        x2 = pmb->pcoord->x2v(j);
        x3 = pmb->pcoord->x3v(k);
        Real CS = std::cos(kx*x1+ky*x2+kz*x3);
        dbyc += volume(i)*2.0*(b(IB2,k,j,i)-rby)*CS;
      }
    }
  }
  Real dbyc0 = amp*rby0*std::sqrt(sch*std::sqrt(k20*beta/(1.0+beta)));
  return dbyc/(total_volume*dbyc0);
}
} // namespace
