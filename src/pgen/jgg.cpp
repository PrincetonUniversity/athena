//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file jgg.cpp
//! \brief Problem generator for linear MHD in shearing sheet.
//!
//! PURPOSE:  Problem generator for linear MHD in shearing sheet. Based on the initial
//!   conditions described in "ORBITAL ADVECTION BY INTERPOLATION: A FAST AND ACCURATE
//!   NUMERICAL SCHEME FOR SUPER-FAST MHD FLOWS" by Johnson, Guan, & Gammie, or JGG.
//!
//! Two kinds of perturbations are possible:
//!
//! - ipert = 1 - MHD simple shwave test
//! - ipert = 2 - MHD compressive shwave test of JGG -- their figure 11
//!
//! REFERENCE: Johnson, Guan, & Gammie, ApJS, 177, 373 (2008)
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

#ifdef FFT
#include "../fft/athena_fft.hpp"
#endif

#if !MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires magnetic fields"
#endif

namespace {
Real iso_cs, gm1, d0, p0;
Real nwx, nwy, nwz; // Wavenumbers
Real Lx, Ly, Lz; // root grid size, global to share with output functions
Real amp, beta;
Real qshear, Omega0;
int ipert;
Real hst_dt, hst_next_time;
bool error_output;

#ifdef FFT
  fftw_plan fplan;
  fftw_complex *fft_data;
#endif

Real HistoryBxc(MeshBlock *pmb, int iout);
Real HistoryBzc(MeshBlock *pmb, int iout);
Real HistoryBs(MeshBlock *pmb, int iout);
Real HistorydBxc(MeshBlock *pmb, int iout);
Real HistorydByc(MeshBlock *pmb, int iout);
Real HistoryxidBx(MeshBlock *pmb, int iout);
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

  // shearing box parameters
  qshear = pin->GetReal("orbital_advection","qshear");
  Omega0 = pin->GetReal("orbital_advection","Omega0");

  if (error_output) {
    if (ipert == 1) {
      // allocateDataField
      AllocateRealUserMeshDataField(1);
      ruser_mesh_data[0].NewAthenaArray(mesh_size.nx3, mesh_size.nx2, mesh_size.nx1);
      // allocate User-defined History Output
      AllocateUserHistoryOutput(2);
      if (mesh_size.nx3 > 1) { // 3D
        EnrollUserHistoryOutput(0, HistoryBzc, "Bzc");
        EnrollUserHistoryOutput(1, HistoryBs, "Bzs");
      } else { // 2D
        EnrollUserHistoryOutput(0, HistoryBxc, "Bxc");
        EnrollUserHistoryOutput(1, HistoryBs, "Bxs");
      }
    } else if (ipert == 2) {
      // preparation for fft
#ifndef FFT
        std::stringstream msg;
        msg << "### FATAL ERROR in jgg.cpp ProblemGenerator" << std::endl
            << "FFT option is required for error output."    << std::endl;
        ATHENA_ERROR(msg);
#else
      if (Globals::my_rank==0 || (Globals::nranks>1 && Globals::my_rank == 1)) {
        fft_data = new fftw_complex[mesh_size.nx2];
        fplan =  fftw_plan_dft_1d(mesh_size.nx2, fft_data, fft_data,
                                  FFTW_FORWARD, FFTW_ESTIMATE);
      }
#endif

      // allocateDataField
      AllocateRealUserMeshDataField(4);
      ruser_mesh_data[0].NewAthenaArray(mesh_size.nx3, mesh_size.nx2, mesh_size.nx1);
      ruser_mesh_data[1].NewAthenaArray(mesh_size.nx3, mesh_size.nx2, mesh_size.nx1);
      ruser_mesh_data[2].NewAthenaArray(2); // [dbx] 0: amplitude, 1: error rate
      ruser_mesh_data[2](0) = 1.0;
      ruser_mesh_data[2](1) = 0.0;
      ruser_mesh_data[3].NewAthenaArray(1); // [dby] 0: amplitude
      ruser_mesh_data[3](0) = 1.0;
      // allocate User-defined History Output
      AllocateUserHistoryOutput(3);
      EnrollUserHistoryOutput(0, HistorydBxc, "dBxc");
      EnrollUserHistoryOutput(1, HistorydByc, "dByc");
      EnrollUserHistoryOutput(2, HistoryxidBx, "xidBx");
    } // ipert == 2
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
  }
  if (ipert!=1 && ipert!=2) {
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

  Real kx = (TWO_PI/Lx)*(static_cast<Real>(nwx));
  Real ky = (TWO_PI/Ly)*(static_cast<Real>(nwy));
  Real kz = (TWO_PI/Lz)*(static_cast<Real>(nwz));

  // Initialize perturbations
  Real x1, x2, x3;
  Real rp(0.0);
  // ipert = 1 - MHD linear shwave test
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
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega0*pcoord->x1v(i);
          phydro->u(IM3,k,j,i) = 0.0;
        }
      }
    }
    // initialize interface B
    // set bx
    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          if (block_size.nx3>1) { // 3D
            pfield->b.x1f(k,j,i) = 0.0;
          } else { // 2D
            x2 = pcoord->x2v(j);
            pfield->b.x1f(k,j,i) = amp*std::cos(ky*x2);
          }
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
          if (block_size.nx3>1) { // 3D
            x2 = pcoord->x2v(j);
            pfield->b.x3f(k,j,i) = amp*std::cos(ky*x2);
          } else { // 2D
            pfield->b.x3f(k,j,i) = 0.0;
          }
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
    Real B02  = p0/beta;
    Real k2    =  SQR(kx)+SQR(ky)+SQR(kz);

    rbx =  ky*std::sqrt(B02/(SQR(kx)+SQR(ky)));
    rby = -kx*std::sqrt(B02/(SQR(kx)+SQR(ky)));
    rbz =  0.0;

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
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega0*x1;
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
  if (error_output) {
    // check flag
    if ((present_time < tlim) && (nlim < 0 || (ncycle + 1) < nlim)
        && (present_time > hst_next_time)) {
      flag = true;
      hst_next_time += hst_dt;
    }
    if ((present_time >= tlim) || (nlim >= 0 && (ncycle + 1) >= nlim)) {
      flag = true;
    }
  }

  if (flag) {
    Real rbx(0.0), rby0(0.0), rby(0.0);
    Real dbx0(0.0), dby0(0.0);

    Real kx0 = (TWO_PI/Lx)*(static_cast<Real>(nwx));
    Real ky  = (TWO_PI/Ly)*(static_cast<Real>(nwy));
    Real kz  = (TWO_PI/Lz)*(static_cast<Real>(nwz));
    Real kx  = kx0+qshear*Omega0*present_time*ky;
    if (ipert == 2) {
      Real B02 = p0/beta;
      Real sch = iso_cs/Omega0;
      Real k20 = SQR(kx0)+SQR(ky)+SQR(kz);
      rbx  =  ky*std::sqrt(B02/(SQR(kx0)+SQR(ky)));
      rby0 = -kx0*std::sqrt(B02/(SQR(kx0)+SQR(ky)));
      rby  =  rby0-qshear*Omega0*present_time*rbx;
      dbx0 = amp*rbx*std::sqrt(sch*std::sqrt(k20*beta/(1.0+beta)));
      dby0 = amp*rby0*std::sqrt(sch*std::sqrt(k20*beta/(1.0+beta)));
    }

    // initialize ruser_mesh_data[0] and ruser_mesh_data[1]
    for (int k=0; k<mesh_size.nx3; k++) {
      for (int j=0; j<mesh_size.nx2; j++) {
        for (int i=0; i<mesh_size.nx1; i++) {
          ruser_mesh_data[0](k,j,i) = 0.0;
          if (ipert==2) {
            ruser_mesh_data[1](k,j,i) = 0.0;
          }
        }
      }
    }

    for (int bn=0; bn<nblocal; ++bn) {
      MeshBlock *pmb = my_blocks(bn);
      LogicalLocation &loc = pmb->loc;
      if (loc.level == root_level) { // root level
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int j=pmb->js; j<=pmb->je; j++) {
            for (int i=pmb->is; i<=pmb->ie; i++) {
              int ti = static_cast<int>(loc.lx1)*pmb->block_size.nx1+(i-pmb->is);
              int tj = static_cast<int>(loc.lx2)*pmb->block_size.nx2+(j-pmb->js);
              int tk = static_cast<int>(loc.lx3)*pmb->block_size.nx3+(k-pmb->ks);
              if (ipert == 1) {
                Real x1  = pmb->pcoord->x1v(i);
                Real x2  = pmb->pcoord->x2v(j);
                Real x3  = pmb->pcoord->x3v(k);
                Real vol = pmb->pcoord->GetCellVolume(k,j,i);
                Real SN = std::sin(kx*x1+ky*x2+kz*x3);
                if (mesh_size.nx3>1) { // 3D
                  ruser_mesh_data[0](tk,tj,ti) = 2.0*pmb->pfield->bcc(IB3,k,j,i)*vol*SN;
                } else { // 2D
                  ruser_mesh_data[0](tk,tj,ti) = 2.0*pmb->pfield->bcc(IB1,k,j,i)*vol*SN;
                }
              } else { // ipert == 2
                ruser_mesh_data[0](tk,tj,ti) = (pmb->pfield->bcc(IB1,k,j,i)-rbx)/dbx0;
                ruser_mesh_data[1](tk,tj,ti) = (pmb->pfield->bcc(IB2,k,j,i)-rby)/dby0;
              }
            }
          }
        }
      } else if (loc.level-1 == root_level) { // level difference 1
        if (pmb->block_size.nx3==1) { // 2D
          int k = pmb->ks;
          for (int j=pmb->cjs; j<=pmb->cje; j++) {
            for (int i=pmb->cis; i<=pmb->cie; i++) {
              int ii = (i-pmb->cis)*2+pmb->is;
              int jj = (j-pmb->cjs)*2+pmb->js;
              int kk = k;
              int ti = static_cast<int>(loc.lx1>>1)*pmb->block_size.nx1
                       +static_cast<int>(loc.lx1%2)*(pmb->block_size.nx1/2)+(i-pmb->cis);
              int tj = static_cast<int>(loc.lx2>>1)*pmb->block_size.nx2
                       +static_cast<int>(loc.lx2%2)*(pmb->block_size.nx2/2)+(j-pmb->cjs);
              int tk = k;
              if (ipert == 1) {
                Real vol00 = pmb->pcoord->GetCellVolume(kk  ,jj  ,ii  );
                Real vol01 = pmb->pcoord->GetCellVolume(kk  ,jj  ,ii+1);
                Real vol10 = pmb->pcoord->GetCellVolume(kk  ,jj+1,ii  );
                Real vol11 = pmb->pcoord->GetCellVolume(kk  ,jj+1,ii+1);
                Real SN = std::sin(kx*pmb->pcoord->x1f(ii+1)
                                   +ky*pmb->pcoord->x2f(jj+1)+kz*pmb->pcoord->x3v(kk));
                Real bx00 = pmb->pfield->bcc(IB1,kk,jj  ,ii  );
                Real bx01 = pmb->pfield->bcc(IB1,kk,jj  ,ii+1);
                Real bx10 = pmb->pfield->bcc(IB1,kk,jj+1,ii  );
                Real bx11 = pmb->pfield->bcc(IB1,kk,jj+1,ii+1);
                Real bx_vol = bx00*vol00+bx01*vol01
                              +bx10*vol10+bx11*vol11;
                ruser_mesh_data[0](tk,tj,ti) = 2.0*bx_vol*SN;
              } else { // ipert == 2
                Real dbx00 = pmb->pfield->bcc(IB1,kk,jj  ,ii  )-rbx;
                Real dbx01 = pmb->pfield->bcc(IB1,kk,jj  ,ii+1)-rbx;
                Real dbx10 = pmb->pfield->bcc(IB1,kk,jj+1,ii  )-rbx;
                Real dbx11 = pmb->pfield->bcc(IB1,kk,jj+1,ii+1)-rbx;
                ruser_mesh_data[0](tk,tj,ti) = 0.25*(dbx00+dbx01+dbx10+dbx11)/dbx0;
                Real dby00 = pmb->pfield->bcc(IB2,kk,jj  ,ii  )-rby;
                Real dby01 = pmb->pfield->bcc(IB2,kk,jj  ,ii+1)-rby;
                Real dby10 = pmb->pfield->bcc(IB2,kk,jj+1,ii  )-rby;
                Real dby11 = pmb->pfield->bcc(IB2,kk,jj+1,ii+1)-rby;
                ruser_mesh_data[1](tk,tj,ti) = 0.25*(dby00+dby01+dby10+dby11)/dby0;
              }
            }
          }
        } else { // 3D
          for (int k=pmb->cks; k<=pmb->cke; k++) {
            for (int j=pmb->cjs; j<=pmb->cje; j++) {
              for (int i=pmb->cis; i<=pmb->cie; i++) {
                int ii = (i-pmb->cis)*2+pmb->is;
                int jj = (j-pmb->cjs)*2+pmb->js;
                int kk = (k-pmb->cks)*2+pmb->ks;
                int ti = static_cast<int>(loc.lx1>>1)*pmb->block_size.nx1
                         +static_cast<int>(loc.lx1%2)*(pmb->block_size.nx1/2)
                         +(i-pmb->cis);
                int tj = static_cast<int>(loc.lx2>>1)*pmb->block_size.nx2
                         +static_cast<int>(loc.lx2%2)*(pmb->block_size.nx2/2)
                         +(j-pmb->cjs);
                int tk = static_cast<int>(loc.lx3>>1)*pmb->block_size.nx3
                         +static_cast<int>(loc.lx3%2)*(pmb->block_size.nx3/2)
                         +(k-pmb->cks);
                if (ipert == 1) {
                  Real vol000 = pmb->pcoord->GetCellVolume(kk  ,jj  ,ii  );
                  Real vol001 = pmb->pcoord->GetCellVolume(kk  ,jj  ,ii+1);
                  Real vol010 = pmb->pcoord->GetCellVolume(kk  ,jj+1,ii  );
                  Real vol011 = pmb->pcoord->GetCellVolume(kk  ,jj+1,ii+1);
                  Real vol100 = pmb->pcoord->GetCellVolume(kk+1,jj  ,ii  );
                  Real vol101 = pmb->pcoord->GetCellVolume(kk+1,jj  ,ii+1);
                  Real vol110 = pmb->pcoord->GetCellVolume(kk+1,jj+1,ii  );
                  Real vol111 = pmb->pcoord->GetCellVolume(kk+1,jj+1,ii+1);
                  Real SN = std::sin(kx*pmb->pcoord->x1f(ii+1)+ky*pmb->pcoord->x2f(jj+1)
                                     +kz*pmb->pcoord->x3f(kk+1));

                  Real bz000 = pmb->pfield->bcc(IB3,kk  ,jj  ,ii  );
                  Real bz001 = pmb->pfield->bcc(IB3,kk  ,jj  ,ii+1);
                  Real bz010 = pmb->pfield->bcc(IB3,kk  ,jj+1,ii  );
                  Real bz011 = pmb->pfield->bcc(IB3,kk  ,jj+1,ii+1);
                  Real bz100 = pmb->pfield->bcc(IB3,kk+1,jj  ,ii  );
                  Real bz101 = pmb->pfield->bcc(IB3,kk+1,jj  ,ii+1);
                  Real bz110 = pmb->pfield->bcc(IB3,kk+1,jj+1,ii  );
                  Real bz111 = pmb->pfield->bcc(IB3,kk+1,jj+1,ii+1);
                  Real bz_vol = bz000*vol000+bz001*vol001
                                +bz010*vol010+bz011*vol011
                                +bz100*vol100+bz101*vol101
                                +bz110*vol110+bz111*vol111;
                  ruser_mesh_data[0](tk,tj,ti) = 2.0*bz_vol*SN;
                } else {
                  Real dbx000 = pmb->pfield->bcc(IB1,kk  ,jj  ,ii  )-rbx;
                  Real dbx001 = pmb->pfield->bcc(IB1,kk  ,jj  ,ii+1)-rbx;
                  Real dbx010 = pmb->pfield->bcc(IB1,kk  ,jj+1,ii  )-rbx;
                  Real dbx011 = pmb->pfield->bcc(IB1,kk  ,jj+1,ii+1)-rbx;
                  Real dbx100 = pmb->pfield->bcc(IB1,kk+1,jj  ,ii  )-rbx;
                  Real dbx101 = pmb->pfield->bcc(IB1,kk+1,jj  ,ii+1)-rbx;
                  Real dbx110 = pmb->pfield->bcc(IB1,kk+1,jj+1,ii  )-rbx;
                  Real dbx111 = pmb->pfield->bcc(IB1,kk+1,jj+1,ii+1)-rbx;
                  ruser_mesh_data[0](tk,tj,ti) = 0.125*(dbx000+dbx001+dbx010+dbx011
                                                       +dbx100+dbx101+dbx110+dbx111)/dbx0;
                  Real dby000 = pmb->pfield->bcc(IB2,kk  ,jj  ,ii  )-rby;
                  Real dby001 = pmb->pfield->bcc(IB2,kk  ,jj  ,ii+1)-rby;
                  Real dby010 = pmb->pfield->bcc(IB2,kk  ,jj+1,ii  )-rby;
                  Real dby011 = pmb->pfield->bcc(IB2,kk  ,jj+1,ii+1)-rby;
                  Real dby100 = pmb->pfield->bcc(IB2,kk+1,jj  ,ii  )-rby;
                  Real dby101 = pmb->pfield->bcc(IB2,kk+1,jj  ,ii+1)-rby;
                  Real dby110 = pmb->pfield->bcc(IB2,kk+1,jj+1,ii  )-rby;
                  Real dby111 = pmb->pfield->bcc(IB2,kk+1,jj+1,ii+1)-rby;
                  ruser_mesh_data[1](tk,tj,ti) = 0.125*(dby000+dby001+dby010+dby011
                                                       +dby100+dby101+dby110+dby111)/dby0;
                }
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
    if (Globals::nranks > 1) {
      int ntot = mesh_size.nx3*mesh_size.nx2*mesh_size.nx1;
      if (Globals::my_rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, ruser_mesh_data[0].data(), ntot,
                   MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
      } else {
        MPI_Reduce(ruser_mesh_data[0].data(), ruser_mesh_data[0].data(), ntot,
                   MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
      }
      if (ipert == 2) {
        if (Globals::my_rank == 1) {
          MPI_Reduce(MPI_IN_PLACE, ruser_mesh_data[1].data(), ntot,
                     MPI_ATHENA_REAL, MPI_SUM, 1, MPI_COMM_WORLD);
        } else {
          MPI_Reduce(ruser_mesh_data[1].data(), ruser_mesh_data[1].data(), ntot,
                     MPI_ATHENA_REAL, MPI_SUM, 1, MPI_COMM_WORLD);
        }
      }
    }
#endif
#ifdef FFT
    // fft execution
    if (ipert == 2) {
      int exe_rank_dby = (Globals::nranks>1)? 1 : 0;
      Real inv_nik = 1.0/static_cast<Real>(mesh_size.nx3*mesh_size.nx1);
      Real inv_half_nx2 = 1.0/static_cast<Real>(mesh_size.nx2/2);
      int abs_nwy = std::abs(nwy);
      // calculate FFT of dbx
      if (Globals::my_rank == 0) {
        ruser_mesh_data[2](0) = 0.0; // amplitude of main mode
        ruser_mesh_data[2](1) = 0.0; // error rate
        for (int k=0; k<mesh_size.nx3; k++) {
          for (int i=0; i<mesh_size.nx1; i++) {
            // set into buffer
            for (int j=0; j<mesh_size.nx2; j++) {
              fft_data[j][0] = ruser_mesh_data[0](k,j,i);
              fft_data[j][1] = 0.0;
            }
            fftw_execute(fplan);
            // calculation
            Real amp_wave = 0.0;
            Real amp_others = 0.0;
            for (int j=0; j<mesh_size.nx2/2; j++) {
              if (j == abs_nwy) {
                amp_wave = std::sqrt(SQR(fft_data[j][0])+SQR(fft_data[j][1]))
                           *inv_half_nx2;
              } else {
                amp_others += std::sqrt(SQR(fft_data[j][0])+SQR(fft_data[j][1]))
                              *inv_half_nx2;
              }
            }
            ruser_mesh_data[2](0) += amp_wave*inv_nik;
            ruser_mesh_data[2](1) += (amp_others/amp_wave)*inv_nik;
          }
        }
      }

      // calculate FFT of dby
      if (Globals::my_rank == exe_rank_dby) {
        ruser_mesh_data[3](0) = 0.0; // amplitude of main mode
        for (int k=0; k<mesh_size.nx3; k++) {
          for (int i=0; i<mesh_size.nx1; i++) {
            // set into buffer
            for (int j=0; j<mesh_size.nx2; j++) {
              fft_data[j][0] = ruser_mesh_data[1](k,j,i);
              fft_data[j][1] = 0.0;
            }
            fftw_execute(fplan);
            // calculation
            Real amp_wave = std::sqrt(SQR(fft_data[abs_nwy][0])
                                     +SQR(fft_data[abs_nwy][1]))*inv_half_nx2;
            ruser_mesh_data[3](0) += amp_wave*inv_nik;
          }
        }
      }
    }
#endif
  } // flag
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Function called after main loop is finished for user-defined work.
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
#ifdef FFT
  if(fplan) fftw_destroy_plan(fplan);
  if(fft_data) delete[] fft_data;
#endif
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

  Real qom = qshear*Omega0;
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

Real HistoryBzc(MeshBlock *pmb, int iout) {
  Real x1, x2, x3;
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);

  Real qom = qshear*Omega0;
  Real bzc = 0.0;

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
        bzc += volume(i)*2.0*b(IB3,k,j,i)*CS;
      }
    }
  }
  return bzc/(amp*total_volume);
}

Real HistoryBs(MeshBlock *pmb, int iout) {
  if (Globals::my_rank != 0) return 0.0;
  if (pmb->lid != 0) return 0.0;
  AthenaArray<Real> &bb = pmb->pmy_mesh->ruser_mesh_data[0];
  Real bs = 0.0;
  Real bs_temp;
  int nx1 = pmb->pmy_mesh->mesh_size.nx1;
  int nx2 = pmb->pmy_mesh->mesh_size.nx2;
  int nx3 = pmb->pmy_mesh->mesh_size.nx3;
  Real total_volume = Lx*Ly*Ly;
  for (int k=0; k<nx3; k++) {
    for (int i=0; i<nx1; i++) {
      bs_temp = 0.0;
      for (int j=0; j<nx2; j++) {
        bs_temp += bb(0,k,j,i);
      }
      bs += std::abs(bs_temp);
    }
  }
  return bs/(amp*total_volume);
}

Real HistorydBxc(MeshBlock *pmb, int iout) {
  if (Globals::my_rank != 0) return 0.0;
  if (pmb->lid != 0) return 0.0;
  return pmb->pmy_mesh->ruser_mesh_data[2](0);
}

Real HistorydByc(MeshBlock *pmb, int iout) {
  int exe_rank_dby = (Globals::nranks>1)? 1 : 0;
  if (Globals::my_rank != exe_rank_dby) return 0.0;
  if (pmb->lid != 0) return 0.0;
  return pmb->pmy_mesh->ruser_mesh_data[3](0);
}

Real HistoryxidBx(MeshBlock *pmb, int iout) {
  if (Globals::my_rank != 0) return 0.0;
  if (pmb->lid != 0) return 0.0;
  return pmb->pmy_mesh->ruser_mesh_data[2](1);
}
} // namespace
