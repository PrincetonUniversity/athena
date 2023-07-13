//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file ssheet.cpp
//! \brief Shearing wave problem generator.
//! Several different initial conditions:
//!  - ipert = 1  pure shearing background flow
//!  - ipert = 2  epicycle motion (0.1 c_s initial kick in radial)
//!  - ipert = 3  shwave test for compressible flow
//!               JG: Johnson & Gammie 2005, ApJ, 626, 978
//======================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <fstream>    // ofstream
#include <iomanip>    // setprecision
#include <iostream>   // cout, endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../outputs/outputs.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires does not work with MHD."
#endif

namespace {
Real iso_cs, gm1, d0, p0;
Real amp; // amplitude
int nwx, nwy; // Wavenumbers
Real x1size,x2size,x3size;
int ipert; // initial pattern
Real qshear, Omega0;
Real hst_dt, hst_next_time;
bool error_output;

Real Historydvyc(MeshBlock *pmb, int iout);
Real Historyvxs(MeshBlock *pmb, int iout);
Real Historydvys(MeshBlock *pmb, int iout);
} // namespace

//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (!shear_periodic) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box."   << std::endl;
    ATHENA_ERROR(msg);
  }

  if (mesh_size.nx2 == 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator" << std::endl
        << "This problem does NOT work on a 1D grid." << std::endl;
    ATHENA_ERROR(msg);
  }

  // read ipert parameter
  d0 = 1.0;
  p0 = 1e-6;
  if (NON_BAROTROPIC_EOS) {
    gm1 = (pin->GetReal("hydro","gamma") - 1.0);
    iso_cs = std::sqrt((gm1+1.0)*p0/d0);
  } else {
    iso_cs = pin->GetReal("hydro","iso_sound_speed");
    p0 = d0*SQR(iso_cs);
  }
  ipert = pin->GetInteger("problem","ipert");
  x1size = mesh_size.x1max - mesh_size.x1min;
  x2size = mesh_size.x2max - mesh_size.x2min;
  x3size = mesh_size.x3max - mesh_size.x3min;

  // shearing box parameters
  qshear = pin->GetReal("orbital_advection","qshear");
  Omega0 = pin->GetReal("orbital_advection","Omega0");

  if (ipert == 3) {
    amp = pin->GetReal("problem","amp");
    nwx = pin->GetInteger("problem","nwx");
    nwy = pin->GetInteger("problem","nwy");
    if (nwx == 0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator"   << std::endl
          << "Parameter nwx must be non-zero." << std::endl;
      ATHENA_ERROR(msg);
    }
    error_output = pin->GetOrAddBoolean("problem","error_output",false);
    if (error_output) {
      // allocateDataField
      AllocateRealUserMeshDataField(2);
      ruser_mesh_data[0].NewAthenaArray(mesh_size.nx3, mesh_size.nx2, mesh_size.nx1);
      ruser_mesh_data[1].NewAthenaArray(mesh_size.nx3, mesh_size.nx2, mesh_size.nx1);

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

      // allocate User-defined History Output
      AllocateUserHistoryOutput(3);
      EnrollUserHistoryOutput(0, Historydvyc, "dvyc",
                              UserHistoryOperation::sum);
      EnrollUserHistoryOutput(1, Historyvxs,  "vxs",
                              UserHistoryOperation::sum);
      EnrollUserHistoryOutput(2, Historydvys, "dvys",
                              UserHistoryOperation::sum);
    }
  } else if (ipert == 1 || ipert == 2) {
    amp = 0.0;
    nwx = 0;
    nwy = 0;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator"   << std::endl
        << "This problem requires that ipert is from 1 to 3." << std::endl;
    ATHENA_ERROR(msg);
  }

  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  int shboxcoord = porb->shboxcoord;
  int il = is - NGHOST; int iu = ie + NGHOST;
  int jl = js - NGHOST; int ju = je + NGHOST;
  int kl = ks;          int ku = ke;
  if (block_size.nx3 > 1) {
    kl = ks - NGHOST;
    ku = ke + NGHOST;
  }

  if (gid == 0) {
    std::cout << "iso_cs = " << iso_cs << std::endl;
    std::cout << "d0 = " << d0 << std::endl;
    std::cout << "p0 = " << p0 << std::endl;
    std::cout << "ipert  = " << ipert  << std::endl;
    std::cout << "[ssheet.cpp]: [Lx,Ly,Lz] = [" <<x1size <<","<<x2size
              <<","<<x3size<<"]"<<std::endl;
  }

  // calculate wave number just for ipert = 3
  Real kx(0.0), ky(0.0);
  if (ipert==3) {
    kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
    ky = (TWO_PI/x2size)*(static_cast<Real>(nwy));
  }

  Real x1, x2, rd, rp, rvx, rvy;
  // update the physical variables as initial conditions
  for (int k=kl; k<=ku; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=il; i<=iu; i++) {
        x1 = pcoord->x1v(i);
        x2 = pcoord->x2v(j);
        rd = d0;
        rp = p0;
        if (ipert == 1) {
          // 1) pure shear bg flow:
          phydro->u(IDN,k,j,i) = rd;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega0*x1;
          phydro->u(IM3,k,j,i) = 0.0;
        } else if (ipert == 2) {
          // 2) epicyclic oscillation
          if (shboxcoord == 1) { // x-y shear
            rvx = 0.1*iso_cs;
            rvy = 0.0;
            phydro->u(IDN,k,j,i) = rd;
            phydro->u(IM1,k,j,i) = rd*rvx;
            phydro->u(IM2,k,j,i) = -rd*rvy;
            if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega0*x1;
            phydro->u(IM3,k,j,i) = 0.0;
          } else { // x-z plane
            rvx = 0.1*iso_cs;
            rvy = 0.0;
            phydro->u(IDN,k,j,i) = rd;
            phydro->u(IM1,k,j,i) = rd*rvx;
            phydro->u(IM2,k,j,i) = 0.0;
            phydro->u(IM3,k,j,i) = -rd*(rvy+qshear*Omega0*x1);
          }
        } else if (ipert == 3) {
          // 3) JG HD shwave test
          rvx = amp*iso_cs*std::cos(kx*x1 + ky*x2);
          rvy = amp*iso_cs*(ky/kx)*std::cos(kx*x1 + ky*x2);
          phydro->u(IDN,k,j,i) = rd;
          phydro->u(IM1,k,j,i) = -rd*rvx;
          phydro->u(IM2,k,j,i) = -rd*rvy;
          if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega0*x1;
          phydro->u(IM3,k,j,i) = 0.0;
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator" << std::endl
              << "Shearing wave sheet ipert=" << ipert << " is unrecognized" << std::endl;
          ATHENA_ERROR(msg);
        }
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = rp/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                                               SQR(phydro->u(IM2,k,j,i)) +
                                               SQR(phydro->u(IM3,k,j,i))
                                              ) / phydro->u(IDN,k,j,i);
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
  int cncycle = ncycle + 1;
  if (error_output) {
    if ((present_time < tlim) && (nlim < 0 || cncycle < nlim)
        && (present_time > hst_next_time)) {
      flag = true;
      hst_next_time += hst_dt;
    }
    if ((present_time >= tlim) || (nlim >= 0 && cncycle >= nlim)) {
      flag = true;
    }
  }
  // calculate vxs, dvys
  if (flag) {
    AthenaArray<Real> &vxs = ruser_mesh_data[0];
    AthenaArray<Real> &dvys = ruser_mesh_data[1];
    Real qom = qshear*Omega0;
    Real kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
    Real ky = (TWO_PI/x2size)*(static_cast<Real>(nwy));
    kx += qom*present_time*ky;

    // initialize vs
    for (int k=0; k<mesh_size.nx3; k++) {
      for (int j=0; j<mesh_size.nx2; j++) {
        for (int i=0; i<mesh_size.nx1; i++) {
          vxs(k,j,i) = 0.0;
          dvys(k,j,i) = 0.0;
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
              Real x1 = pmb->pcoord->x1v(i);
              Real x2 = pmb->pcoord->x2v(j);
              Real vx = pmb->phydro->w(IVX,k,j,i);
              Real dvy = pmb->phydro->w(IVY,k,j,i);
              if(!pmb->porb->orbital_advection_defined)
                dvy += qom*x1;
              Real vol = pmb->pcoord->GetCellVolume(k,j,i);
              Real SN = std::sin(kx*x1+ky*x2);
              vxs(tk,tj,ti) = 2.0*vx*vol*SN;
              dvys(tk,tj,ti) = 2.0*dvy*vol*SN;
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
              Real vol00 = pmb->pcoord->GetCellVolume(kk  ,jj  ,ii  );
              Real vol01 = pmb->pcoord->GetCellVolume(kk  ,jj  ,ii+1);
              Real vol10 = pmb->pcoord->GetCellVolume(kk  ,jj+1,ii  );
              Real vol11 = pmb->pcoord->GetCellVolume(kk  ,jj+1,ii+1);
              Real vx00 = pmb->phydro->w(IVX,kk  ,jj  ,ii  );
              Real vx01 = pmb->phydro->w(IVX,kk  ,jj  ,ii+1);
              Real vx10 = pmb->phydro->w(IVX,kk  ,jj+1,ii  );
              Real vx11 = pmb->phydro->w(IVX,kk  ,jj+1,ii+1);
              Real dvy00 = pmb->phydro->w(IVY,kk  ,jj  ,ii  );
              Real dvy01 = pmb->phydro->w(IVY,kk  ,jj  ,ii+1);
              Real dvy10 = pmb->phydro->w(IVY,kk  ,jj+1,ii  );
              Real dvy11 = pmb->phydro->w(IVY,kk  ,jj+1,ii+1);
              if(!pmb->porb->orbital_advection_defined) {
                dvy00 += qom*pmb->pcoord->x1v(ii  );
                dvy01 += qom*pmb->pcoord->x1v(ii+1);
                dvy10 += qom*pmb->pcoord->x1v(ii  );
                dvy11 += qom*pmb->pcoord->x1v(ii+1);
              }
              Real SN = std::sin(kx*pmb->pcoord->x1f(ii+1)
                                 +ky*pmb->pcoord->x2f(jj+1));
              Real vx_vol  = vx00*vol00+vx01*vol01
                             +vx10*vol10+vx11*vol11;
              Real dvy_vol = dvy00*vol00+dvy01*vol01
                             +dvy10*vol10+dvy11*vol11;
              vxs(tk,tj,ti) = 2.0*SN*vx_vol;
              dvys(tk,tj,ti) = 2.0*SN*dvy_vol;
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
                Real vol000 = pmb->pcoord->GetCellVolume(kk  ,jj  ,ii  );
                Real vol001 = pmb->pcoord->GetCellVolume(kk  ,jj  ,ii+1);
                Real vol010 = pmb->pcoord->GetCellVolume(kk  ,jj+1,ii  );
                Real vol011 = pmb->pcoord->GetCellVolume(kk  ,jj+1,ii+1);
                Real vol100 = pmb->pcoord->GetCellVolume(kk+1,jj  ,ii  );
                Real vol101 = pmb->pcoord->GetCellVolume(kk+1,jj  ,ii+1);
                Real vol110 = pmb->pcoord->GetCellVolume(kk+1,jj+1,ii  );
                Real vol111 = pmb->pcoord->GetCellVolume(kk+1,jj+1,ii+1);
                Real vx000 = pmb->phydro->w(IVX,kk  ,jj  ,ii  );
                Real vx001 = pmb->phydro->w(IVX,kk  ,jj  ,ii+1);
                Real vx010 = pmb->phydro->w(IVX,kk  ,jj+1,ii  );
                Real vx011 = pmb->phydro->w(IVX,kk  ,jj+1,ii+1);
                Real vx100 = pmb->phydro->w(IVX,kk+1,jj  ,ii  );
                Real vx101 = pmb->phydro->w(IVX,kk+1,jj  ,ii+1);
                Real vx110 = pmb->phydro->w(IVX,kk+1,jj+1,ii  );
                Real vx111 = pmb->phydro->w(IVX,kk+1,jj+1,ii+1);
                Real dvy000 = pmb->phydro->w(IVY,kk  ,jj  ,ii  );
                Real dvy001 = pmb->phydro->w(IVY,kk  ,jj  ,ii+1);
                Real dvy010 = pmb->phydro->w(IVY,kk  ,jj+1,ii  );
                Real dvy011 = pmb->phydro->w(IVY,kk  ,jj+1,ii+1);
                Real dvy100 = pmb->phydro->w(IVY,kk+1,jj  ,ii  );
                Real dvy101 = pmb->phydro->w(IVY,kk+1,jj  ,ii+1);
                Real dvy110 = pmb->phydro->w(IVY,kk+1,jj+1,ii  );
                Real dvy111 = pmb->phydro->w(IVY,kk+1,jj+1,ii+1);
                if(!pmb->porb->orbital_advection_defined) {
                  dvy000 += qom*pmb->pcoord->x1v(ii  );
                  dvy001 += qom*pmb->pcoord->x1v(ii+1);
                  dvy010 += qom*pmb->pcoord->x1v(ii  );
                  dvy011 += qom*pmb->pcoord->x1v(ii+1);
                  dvy100 += qom*pmb->pcoord->x1v(ii  );
                  dvy101 += qom*pmb->pcoord->x1v(ii+1);
                  dvy110 += qom*pmb->pcoord->x1v(ii  );
                  dvy111 += qom*pmb->pcoord->x1v(ii+1);
                }
                Real SN = std::sin(kx*pmb->pcoord->x1f(ii+1)
                                   +ky*pmb->pcoord->x2f(jj+1));
                Real vx_vol = vx000*vol000+vx001*vol001
                              +vx010*vol010+vx011*vol011
                              +vx100*vol100+vx101*vol101
                              +vx110*vol110+vx111*vol111;
                Real dvy_vol = dvy000*vol000+dvy001*vol001
                               +dvy010*vol010+dvy011*vol011
                               +dvy100*vol100+dvy101*vol101
                               +dvy110*vol110+dvy111*vol111;
                vxs(tk,tj,ti) = 2.0*SN*vx_vol;
                dvys(tk,tj,ti) = 2.0*SN*dvy_vol;
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
      // vxs
      if (Globals::my_rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, vxs.data(), ntot, MPI_ATHENA_REAL,
                   MPI_SUM, 0, MPI_COMM_WORLD);
      } else {
        MPI_Reduce(vxs.data(), vxs.data(), ntot, MPI_ATHENA_REAL,
                   MPI_SUM, 0, MPI_COMM_WORLD);
      }
      // dvys
      if (Globals::my_rank == 1) {
        MPI_Reduce(MPI_IN_PLACE, dvys.data(), ntot, MPI_ATHENA_REAL,
                   MPI_SUM, 1, MPI_COMM_WORLD);
      } else {
        MPI_Reduce(dvys.data(), dvys.data(), ntot, MPI_ATHENA_REAL,
                   MPI_SUM, 1, MPI_COMM_WORLD);
      }
    }
#endif
  } // flag
  return;
}

namespace {

Real Historydvyc(MeshBlock *pmb, int iout) {
  Real qom = qshear*Omega0;
  Real kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
  Real ky = (TWO_PI/x2size)*(static_cast<Real>(nwy));
  kx += qom*pmb->pmy_mesh->time*ky;
  Real dvyc = 0.0;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);
  Real tvol = x1size*x2size*x3size;
  for (int k=pmb->ks; k<=pmb->ke  ; k++) {
    for (int j=pmb->js; j<=pmb->je  ; j++) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, volume);
      for (int i=pmb->is; i<=pmb->ie  ; i++) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real CS = std::cos(kx*x1+ky*x2);
        Real dvy = pmb->phydro->w(IVY,k,j,i);
        if(!pmb->porb->orbital_advection_defined)
          dvy += qom*x1;
        dvyc += volume(i)*2.0*dvy*CS;
      }
    }
  }
  Real dvy0 = iso_cs*amp
              *std::abs(static_cast<Real>(nwy)/static_cast<Real>(nwx));
  return dvyc/(dvy0*tvol);
}

Real Historyvxs(MeshBlock *pmb, int iout) {
  if (Globals::my_rank != 0) return 0.0;
  if (pmb->lid != 0) return 0.0;
  AthenaArray<Real> &vs = pmb->pmy_mesh->ruser_mesh_data[0];
  Real vxs = 0.0;
  Real vxs_temp;
  int nx1 = pmb->pmy_mesh->mesh_size.nx1;
  int nx2 = pmb->pmy_mesh->mesh_size.nx2;
  int nx3 = pmb->pmy_mesh->mesh_size.nx3;
  Real tvol = x1size*x2size*x3size;
  for (int k=0; k<nx3; k++) {
    for (int i=0; i<nx1; i++) {
      vxs_temp = 0.0;
      for (int j=0; j<nx2; j++) {
        vxs_temp += vs(k,j,i);
      }
      vxs += std::abs(vxs_temp);
    }
  }
  Real vx0 = iso_cs*amp;
  return vxs/(vx0*tvol);
}

Real Historydvys(MeshBlock *pmb, int iout) {
  int exe_rank_dvy = (Globals::nranks>1)? 1 : 0;
  if (Globals::my_rank != exe_rank_dvy) return 0.0;
  if (pmb->lid != 0) return 0.0;
  AthenaArray<Real> &vs = pmb->pmy_mesh->ruser_mesh_data[1];
  Real dvys = 0.0;
  Real dvys_temp;
  int nx1 = pmb->pmy_mesh->mesh_size.nx1;
  int nx2 = pmb->pmy_mesh->mesh_size.nx2;
  int nx3 = pmb->pmy_mesh->mesh_size.nx3;
  Real tvol = x1size*x2size*x3size;
  for (int k=0; k<nx3; k++) {
    for (int i=0; i<nx1; i++) {
      dvys_temp = 0.0;
      for (int j=0; j<nx2; j++) {
        dvys_temp += vs(k,j,i);
      }
      dvys += std::abs(dvys_temp);
    }
  }
  Real dvy0 = iso_cs*amp
              *std::abs(static_cast<Real>(nwy)/static_cast<Real>(nwx));
  return dvys/(dvy0*tvol);
}
} // namespace
