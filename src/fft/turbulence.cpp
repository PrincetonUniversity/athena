//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turbulence.cpp
//  \brief implementation of functions in class Turbulence

// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <random>     // mt19937, normal_distribution, uniform_real_distribution
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "athena_fft.hpp"
#include "turbulence.hpp"

//----------------------------------------------------------------------------------------
//! \fn TurbulenceDriver::TurbulenceDriver(Mesh *pm, ParameterInput *pin)
//! \brief TurbulenceDriver constructor
//!
//! \note
//! turb_flag is initialzed in the Mesh constructor to 0 by default;
//! - turb_flag = 1 for decaying turbulence
//! - turb_flag = 2 for impulsively driven turbulence
//! - turb_flag = 3 for continuously driven turbulence

TurbulenceDriver::TurbulenceDriver(Mesh *pm, ParameterInput *pin) :
    FFTDriver(pm, pin),
    rseed(pin->GetOrAddInteger("turbulence", "rseed", -1)), // seed for RNG
    // cut-off wavenumbers, low and high:
    nlow(pin->GetOrAddInteger("turbulence", "nlow", 0)),
    nhigh(pin->GetOrAddInteger("turbulence", "nhigh", pm->mesh_size.nx1/2)),
    tdrive(pm->time),
    // driving interval must be set manually:
    dtdrive(pm->turb_flag == 2 ? pin->GetReal("turbulence", "dtdrive") : 0.0),
    // correlation time scales for OU smoothing:
    tcorr(pm->turb_flag > 1 ? pin->GetReal("turbulence", "tcorr") : 0.0),
    f_shear(pin->GetOrAddReal("turbulence", "f_shear", -1)), // ratio of shear component
    expo(pin->GetOrAddReal("turbulence", "expo", 2)), // power-law exponent
    dedt(pin->GetReal("turbulence", "dedt")), // turbulence amplitude
    // TODO(changgoo): this assumes 3D and should not work with 1D, 2D. Add check.
    vel{ {nmb, pm->my_blocks(0)->ncells3,
               pm->my_blocks(0)->ncells2, pm->my_blocks(0)->ncells1},
         {nmb, pm->my_blocks(0)->ncells3,
               pm->my_blocks(0)->ncells2, pm->my_blocks(0)->ncells1},
         {nmb, pm->my_blocks(0)->ncells3,
               pm->my_blocks(0)->ncells2, pm->my_blocks(0)->ncells1} },
    fv_new_(nullptr) {
  if (f_shear > 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "The ratio between shear and compressible components should be less than one"
        << std::endl;
    ATHENA_ERROR(msg);
    return;
  }

  if (pm->turb_flag == 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "Turbulence flag is set to zero! Shouldn't reach here!" << std::endl;
    ATHENA_ERROR(msg);
    return;
  } else {
#ifndef FFT
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "non zero Turbulence flag is set without FFT!" << std::endl;
    ATHENA_ERROR(msg);
    return;
#endif
  }

  InitializeFFTBlock(true);
  // note, pmy_fb won't be defined until InitializeFFTBlock is called:
  dvol = pmy_fb->dx1*pmy_fb->dx2*pmy_fb->dx3;
  QuickCreatePlan();

  fv_ = new std::complex<Real>*[3];
  fv_sh_ = new std::complex<Real>*[3];
  fv_co_ = new std::complex<Real>*[3];
  if (pm->turb_flag > 1) fv_new_ = new std::complex<Real>*[3];
  for (int nv=0; nv<3; nv++) {
    fv_[nv] = new std::complex<Real>[pmy_fb->cnt_];
    fv_sh_[nv] = new std::complex<Real>[pmy_fb->cnt_];
    fv_co_[nv] = new std::complex<Real>[pmy_fb->cnt_];
    if (pm->turb_flag > 1) fv_new_[nv] = new std::complex<Real>[pmy_fb->cnt_];
  }

  // initialize MT19937 random number generator
  if (rseed < 0) {
    std::random_device device;
    rseed = static_cast<std::int64_t>(device());
  } else {
    // If rseed is specified with a non-negative value,
    // PS is generated with a global random number sequence.
    // This would make perturbation identical irrespective of number of MPI ranks,
    // but the cost of the PowerSpectrum() function call is huge.
    // Not recommended with turb_flag = 3 or turb_flag = with small dtdrive
    global_ps_ = true;
    if ((pm->turb_flag == 3) & (Globals::my_rank == 0)) {
      std::cout << "### Warning: continuous turbulence driving (turb_flag == 3)"
                << std::endl << " with a specific rseed (rseed >= 0) will be slow"
                << std::endl;
    }
  }
  rng_generator.seed(rseed);
}

// destructor
TurbulenceDriver::~TurbulenceDriver() {
  for (int nv=0; nv<3; nv++) {
    delete [] fv_[nv];
    delete [] fv_sh_[nv];
    delete [] fv_co_[nv];
    if (fv_new_ != nullptr) delete [] fv_new_[nv];
  }
  delete [] fv_;
  delete [] fv_sh_;
  delete [] fv_co_;
  if (fv_new_ != nullptr) delete [] fv_new_;
}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::Driving()
//! \brief Generate and Perturb the velocity field

void TurbulenceDriver::Driving() {
  Mesh *pm = pmy_mesh_;

  // check driving time interval to generate new perturbation
  switch(pm->turb_flag) {
    case 1: // turb_flag == 1 : decaying turbulence
      Generate();
      Perturb(0);
      break;
    case 2: // turb_flag == 2 : impulsively driven turbulence with OU smoothing
      if (pm->time >= tdrive) {
        tdrive = pm->time + dtdrive;
        Generate();
        Perturb(dtdrive);
      }
      break;
    case 3: // turb_flag == 3 : continuously driven turbulence with OU smoothing
      Generate();
      Perturb(pm->dt);
      break;
    default:
      std::stringstream msg;
      msg << "### FATAL ERROR in TurbulenceDriver::Driving" << std::endl
          << "Turbulence flag " << pm->turb_flag << " is not supported!" << std::endl;
      ATHENA_ERROR(msg);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::Generate()
//! \brief Generate velocity pertubation.

void TurbulenceDriver::Generate() {
  Mesh *pm = pmy_mesh_;
  FFTBlock *pfb = pmy_fb;
  AthenaFFTPlan *plan = pfb->bplan_;

  int nbs = nslist_[Globals::my_rank];
  int nbe = nbs + nblist_[Globals::my_rank] - 1;


  // For driven turbulence (turb_flag == 2 or 3),
  // Ornstein-Uhlenbeck (OU) process is implemented.
  // fv_ are set initially (or in restart) and kept
  // unless tcorr == 0
  if (!initialized_) {
    for (int nv=0; nv<3; nv++) {
      std::complex<Real> *fv = fv_[nv];
      PowerSpectrum(fv);
    }
    if (f_shear >= 0) Project(fv_, f_shear);
    if (tcorr > 0.) initialized_ = true;
  } else {
    Real OUdt = pm->dt;
    if (pm->turb_flag == 2) OUdt=dtdrive;
    OUProcess(OUdt);
  }

  for (int nv=0; nv<3; nv++) {
    AthenaArray<Real> &dv = vel[nv], dv_mb;
    for (int kidx=0; kidx<pfb->cnt_; kidx++) pfb->in_[kidx] = fv_[nv][kidx];
    pfb->Execute(plan);
    for (int nb=0; nb<pm->nblocal; ++nb) {
      MeshBlock *pmb = pm->my_blocks(nb);
      dv_mb.InitWithShallowSlice(dv, 4, nb, 1);
      pfb->RetrieveResult(dv_mb, 0, NGHOST, pmb->loc, pmb->block_size);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::OUProcess(Real dt)
//! \brief Generate velocity pertubation.
//!
//! Ornstein–Uhlenbeck (OU) process based on Eq. 26 in Lynn et al. (2012)
//! original formalism for
//! \f$ f=exp(-dt/tcorr) \f$
//! \f[ dv_k(t+dt) = f*dv_k(t) + sqrt(1-f^2)*dv_k' \f]

void TurbulenceDriver::OUProcess(Real dt) {
  FFTBlock *pfb = pmy_fb;
  Real factor = std::exp(-dt/tcorr);
  Real sqrt_factor = std::sqrt(1 - factor*factor);

  for (int nv=0; nv<3; nv++) PowerSpectrum(fv_new_[nv]);

  if (f_shear >= 0) Project(fv_new_, f_shear);

  for (int nv=0; nv<3; nv++) {
    for (int k=0; k<pfb->cnt_; k++) {
      fv_[nv][k] = factor * fv_[nv][k] + sqrt_factor * fv_new_[nv][k];
    }
  }
}


//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::PowerSpectrum(std::complex<Real> *amp)
//! \brief Generate Power spectrum in Fourier space with power-law

void TurbulenceDriver::PowerSpectrum(std::complex<Real> *amp) {
  Real pcoeff;
  FFTBlock *pfb = pmy_fb;
  AthenaFFTIndex *idx = pfb->b_in_;
  int kNx1 = pfb->kNx[0], kNx2 = pfb->kNx[1], kNx3 = pfb->kNx[2];
  int knx1 = pfb->knx[0], knx2 = pfb->knx[1], knx3 = pfb->knx[2];
  int kdisp1 = pfb->kdisp[0], kdisp2 = pfb->kdisp[1], kdisp3 = pfb->kdisp[2];

  std::normal_distribution<Real> ndist(0.0,1.0); // standard normal distribution
  std::uniform_real_distribution<Real> udist(0.0,1.0); // uniform in [0,1)

  // set random amplitudes with gaussian deviation
  // loop over entire Mesh
  if (global_ps_) {
    for (int gk=0; gk<kNx3; gk++) {
      for (int gj=0; gj<kNx2; gj++) {
        for (int gi=0; gi<kNx1; gi++) {
          int k = gk - kdisp3;
          int j = gj - kdisp2;
          int i = gi - kdisp1;
          if ((k >= 0) && (k < knx3) &&
              (j >= 0) && (j < knx2) &&
              (i >= 0) && (i < knx1)) {
            // Box-Muller Method
            //Real q1 = ran2(&rseed);
            //Real q2 = ran2(&rseed);
            //Real A = std::sqrt(-2.0*std::log(q1 + 1.e-20))*std::cos(TWO_PI*q2);
            //Real ph = ran2(&rseed)*TWO_PI;
            std::int64_t kidx = pfb->GetIndex(i,j,k,idx);
            Real A = ndist(rng_generator);
            Real ph = udist(rng_generator)*TWO_PI;
            amp[kidx] = A*std::complex<Real>(std::cos(ph), std::sin(ph));
          } else { // if it is not in FFTBlock, just burn unused random numbers
            Real A = ndist(rng_generator);
            Real ph = udist(rng_generator)*TWO_PI;
            //ran2(&rseed);
            //ran2(&rseed);
            //ran2(&rseed);
          }
        }
      }
    }
  }

  // set power spectrum: only power-law

  // find the reference 2PI/L along the longest axis
  Real dkx = pfb->dkx[0];
  if (knx2>1) dkx = dkx < pfb->dkx[1] ? dkx : pfb->dkx[1];
  if (knx3>2) dkx = dkx < pfb->dkx[2] ? dkx : pfb->dkx[2];

  for (int k=0; k<knx3; k++) {
    for (int j=0; j<knx2; j++) {
      for (int i=0; i<knx1; i++) {
        std::int64_t nx = GetKcomp(i,pfb->kdisp[0],pfb->kNx[0]);
        std::int64_t ny = GetKcomp(j,pfb->kdisp[1],pfb->kNx[1]);
        std::int64_t nz = GetKcomp(k,pfb->kdisp[2],pfb->kNx[2]);
        Real nmag = std::sqrt(nx*nx+ny*ny+nz*nz);
        Real kx = nx*pfb->dkx[0];
        Real ky = ny*pfb->dkx[1];
        Real kz = nz*pfb->dkx[2];
        Real kmag = std::sqrt(kx*kx+ky*ky+kz*kz);

        std::int64_t gidx = pfb->GetGlobalIndex(i,j,k);

        if (gidx == 0) {
          pcoeff = 0.0;
        } else {
          if ((kmag/dkx > nlow) && (kmag/dkx < nhigh)) {
            pcoeff = 1.0/std::pow(kmag,(expo+2.0)/2.0);
          } else {
            pcoeff = 0.0;
          }
        }
        std::int64_t kidx=pfb->GetIndex(i,j,k,idx);

        if (global_ps_) {
          amp[kidx] *= pcoeff;
        } else {
          Real A = ndist(rng_generator);
          Real ph = udist(rng_generator)*TWO_PI;
          amp[kidx] = pcoeff*A*std::complex<Real>(std::cos(ph), std::sin(ph));
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::Perturb(Real dt)
//! \brief Add velocity perturbation to the hydro variables

void TurbulenceDriver::Perturb(Real dt) {
  Mesh *pm = pmy_mesh_;
  int nbs = nslist_[Globals::my_rank];
  int nbe = nbs+nblist_[Globals::my_rank]-1;

  int il = pm->my_blocks(0)->is, iu = pm->my_blocks(0)->ie;
  int jl = pm->my_blocks(0)->js, ju = pm->my_blocks(0)->je;
  int kl = pm->my_blocks(0)->ks, ku = pm->my_blocks(0)->ke;

  Real aa, b, c, s, de, v1, v2, v3, den, M1, M2, M3;
  Real m[4] = {0};
  AthenaArray<Real> &dv1 = vel[0], &dv2 = vel[1], &dv3 = vel[2];

  for (int nb=0; nb<pm->nblocal; ++nb) {
    MeshBlock *pmb = pm->my_blocks(nb);
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=il; i<=iu; i++) {
          den = pmb->phydro->u(IDN,k,j,i);
          m[0] += den;
          m[1] += den*dv1(nb,k,j,i);
          m[2] += den*dv2(nb,k,j,i);
          m[3] += den*dv3(nb,k,j,i);
        }
      }
    }
  }

#ifdef MPI_PARALLEL
  int mpierr;
  // Sum the perturbations over all processors
  mpierr = MPI_Allreduce(MPI_IN_PLACE, m, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) {
    std::stringstream msg;
    msg << "[normalize]: MPI_Allreduce error = " << mpierr << std::endl;
    ATHENA_ERROR(msg);
  }
#endif // MPI_PARALLEL

  for (int nb=0; nb<nmb; nb++) {
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=il; i<=iu; i++) {
          dv1(nb,k,j,i) -= m[1]/m[0];
          dv2(nb,k,j,i) -= m[2]/m[0];
          dv3(nb,k,j,i) -= m[3]/m[0];
        }
      }
    }
  }

  // Calculate unscaled energy of perturbations
  m[0] = 0.0;
  m[1] = 0.0;
  for (int nb=0; nb<pm->nblocal; ++nb) {
    MeshBlock *pmb = pm->my_blocks(nb);
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=il; i<=iu; i++) {
          v1 = dv1(nb,k,j,i);
          v2 = dv2(nb,k,j,i);
          v3 = dv3(nb,k,j,i);
          den = pmb->phydro->u(IDN,k,j,i);
          M1 = pmb->phydro->u(IM1,k,j,i);
          M2 = pmb->phydro->u(IM2,k,j,i);
          M3 = pmb->phydro->u(IM3,k,j,i);
          m[0] += den*(SQR(v1) + SQR(v2) + SQR(v3));
          m[1] += M1*v1 + M2*v2 + M3*v3;
        }
      }
    }
  }

#ifdef MPI_PARALLEL
  // Sum the perturbations over all processors
  mpierr = MPI_Allreduce(MPI_IN_PLACE, m, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) {
    std::stringstream msg;
    msg << "[normalize]: MPI_Allreduce error = "
        << mpierr << std::endl;
    ATHENA_ERROR(msg);
  }
#endif // MPI_PARALLEL

  // Rescale to give the correct energy injection rate
  if (pm->turb_flag > 1) {
    // driven turbulence
    de = dedt*dt;
  } else {
    // decaying turbulence (all in one shot)
    de = dedt;
  }
  aa = 0.5*m[0];
  aa = std::max(aa,static_cast<Real>(1.0e-20));
  b = m[1];
  c = -de/dvol;
  if (b >= 0.0)
    s = (-2.0*c)/(b + std::sqrt(b*b - 4.0*aa*c));
  else
    s = (-b + std::sqrt(b*b - 4.0*aa*c))/(2.0*aa);

  if (std::isnan(s)) std::cout << "[perturb]: s is NaN!" << std::endl;

  // Apply momentum pertubations
  for (int nb=0; nb<pm->nblocal; ++nb) {
    MeshBlock *pmb = pm->my_blocks(nb);
    for (int k=kl; k<=ku; k++) {
      for (int j=jl; j<=ju; j++) {
        for (int i=il; i<=iu; i++) {
          v1 = dv1(nb,k,j,i);
          v2 = dv2(nb,k,j,i);
          v3 = dv3(nb,k,j,i);
          den = pmb->phydro->u(IDN,k,j,i);
          M1 = pmb->phydro->u(IM1,k,j,i);
          M2 = pmb->phydro->u(IM2,k,j,i);
          M3 = pmb->phydro->u(IM3,k,j,i);

          if (NON_BAROTROPIC_EOS) {
            pmb->phydro->u(IEN,k,j,i) += s*(M1*v1 + M2*v2+M3*v3)
                                         + 0.5*s*s*den*(SQR(v1) + SQR(v2) + SQR(v3));
          }
          pmb->phydro->u(IM1,k,j,i) += s*den*v1;
          pmb->phydro->u(IM2,k,j,i) += s*den*v2;
          pmb->phydro->u(IM3,k,j,i) += s*den*v3;
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::Project(std::complex<Real> **fv, Real f_shear)
//! \brief calculate velocity field with a given ratio of shear to comp.

void TurbulenceDriver::Project(std::complex<Real> **fv, Real f_shear) {
  FFTBlock *pfb = pmy_fb;
  Project(fv, fv_sh_, fv_co_);
  for (int nv=0; nv<3; nv++) {
    for (int kidx=0; kidx<pfb->cnt_; kidx++) {
      fv[nv][kidx] = (1-f_shear)*fv_co_[nv][kidx] + f_shear*fv_sh_[nv][kidx];
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::Project(std::complex<Real> **fv,
//!                                    std::complex<Real> **fv_sh,
//!                                    std::complex<Real> **fv_co)
//! \brief calculates shear and compressible components
void TurbulenceDriver::Project(std::complex<Real> **fv, std::complex<Real> **fv_sh,
                               std::complex<Real> **fv_co) {
  FFTBlock *pfb = pmy_fb;
  AthenaFFTIndex *idx = pfb->b_in_;
  int knx1 = pfb->knx[0], knx2 = pfb->knx[1], knx3 = pfb->knx[2];

  for (int k=0; k<knx3; k++) {
    for (int j=0; j<knx2; j++) {
      for (int i=0; i<knx1; i++) {
        // Get khat
        std::int64_t nx = GetKcomp(i, pfb->kdisp[0], pfb->kNx[0]);
        std::int64_t ny = GetKcomp(j, pfb->kdisp[1], pfb->kNx[1]);
        std::int64_t nz = GetKcomp(k, pfb->kdisp[2], pfb->kNx[2]);
        Real kx = nx*pfb->dkx[0];
        Real ky = ny*pfb->dkx[1];
        Real kz = nz*pfb->dkx[2];
        Real kmag = std::sqrt(kx*kx+ky*ky+kz*kz);

        std::int64_t kidx = pfb->GetIndex(i, j, k, idx);
        std::int64_t gidx = pfb->GetGlobalIndex(i,j,k);
        if (gidx == 0.0) {
          fv_co[0][kidx] = std::complex<Real>(0,0);
          fv_co[1][kidx] = std::complex<Real>(0,0);
          fv_co[2][kidx] = std::complex<Real>(0,0);

          fv_sh[0][kidx] = std::complex<Real>(0,0);
          fv_sh[1][kidx] = std::complex<Real>(0,0);
          fv_sh[2][kidx] = std::complex<Real>(0,0);
        } else {
          kx /= kmag;
          ky /= kmag;
          kz /= kmag;
          // Form (khat.f)
          std::complex<Real> kdotf = kx*fv[0][kidx] + ky*fv[1][kidx] + kz*fv[2][kidx];

          fv_co[0][kidx] = kdotf * kx;
          fv_co[1][kidx] = kdotf * ky;
          fv_co[2][kidx] = kdotf * kz;

          fv_sh[0][kidx] = fv[0][kidx] - fv_co[0][kidx];
          fv_sh[1][kidx] = fv[1][kidx] - fv_co[1][kidx];
          fv_sh[2][kidx] = fv[2][kidx] - fv_co[2][kidx];
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void TurbulenceDriver::GetKcomp(int idx, int disp, int Nx)
//! \brief Get k index, which runs from 0, 1, ... Nx/2-1, -Nx/2, -Nx/2+1, ..., -1.

std::int64_t TurbulenceDriver::GetKcomp(int idx, int disp, int Nx) {
  return ((idx+disp) - static_cast<std::int64_t>(2*(idx+disp)/Nx)*Nx);
}
