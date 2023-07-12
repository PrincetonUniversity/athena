//=======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file forward_euler.cpp
//! \brief implementation of the forward Euler solver for chemistry module

// C header

// C++ header
#include <array>
#include <ctime>
#include <iostream>   // endl, ostream
#include <limits>
#include <sstream>    // stringstream
#include <stdexcept>
#include <string>

// Athena++ classes headers
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

// this class header
#include "ode_wrapper.hpp" // NOLINT

namespace {
  Real cfl_cool_sub; //cfl number for subcycling
  Real yfloor; //species abundance floor for calculating the cooling time
  int nsub_max; //maximum number of substeps
  Real GetChemTime(const Real *y, const Real *ydot,
                   const Real E, const Real Edot);
  void IntegrateOneSubstep(Real tsub, Real *y, const Real *ydot,
                           Real &E, const Real Edot);
} // namespace

//----------------------------------------------------------------------------------------
//! \brief ODEWrapper constructor
ODEWrapper::ODEWrapper(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block_ = pmb;
  if (NON_BAROTROPIC_EOS) {
    dim_ = NSPECIES + 1;
  } else {
    dim_ = NSPECIES;
  }
  output_zone_sec_ = pin->GetOrAddBoolean("chemistry", "output_zone_sec", false);
  cfl_cool_sub = pin->GetOrAddReal("chemistry","cfl_cool_sub",0.1);
  yfloor = pin->GetOrAddReal("chemistry","yfloor",1e-3);
  nsub_max = pin->GetOrAddInteger("chemistry","nsub_max",1e5);
}

//----------------------------------------------------------------------------------------
//! \brief ODEWrapper destructor
ODEWrapper::~ODEWrapper() {
}

//----------------------------------------------------------------------------------------
//! \fn void ODEWrapper::Initialize(ParameterInput *pin)
//! \brief Initialize ODE solver parameters
void ODEWrapper::Initialize(ParameterInput *pin) {
  //Note: this cannot be in the constructor, since it needs the PassiveScalars
  //class, and the ODEWrapper class is constructed in the PassiveScalars constructor.
  pmy_spec_ = pmy_block_->pscalars;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ODEWrapper::Integrate(const Real tinit, const Real dt)
//! \brief Integrate the ODE forward for time dt
void ODEWrapper::Integrate(const Real tinit, const Real dt) {
  int is = pmy_block_->is;
  int js = pmy_block_->js;
  int ks = pmy_block_->ks;
  int ie = pmy_block_->ie;
  int je = pmy_block_->je;
  int ke = pmy_block_->ke;
  int ncycle = pmy_block_->pmy_mesh->ncycle;
  clock_t tstart=0.0;
  clock_t tstop;
  if (output_zone_sec_) {
    tstart = std::clock();
  }
  const Real scalar_floor = pmy_block_->peos->GetScalarFloor();
  // primitive conserved variables
  AthenaArray<Real> &u = pmy_block_->phydro->u;
  AthenaArray<Real> &bcc = pmy_block_->pfield->bcc;

  // chemical species and rates

  // TODO(KGF): alternative 1
  //std::array<Real, NSPECIES> y, ydot;

  // TODO(KGF): alternative 2
  Real *y = new Real[NSPECIES];
  Real *ydot = new Real[NSPECIES];

  // internal energy and rates
  Real E = 0.;
  Real Edot = 0.;
  Real time = pmy_block_->pmy_mesh->time;// code time
  Real dt_mhd = pmy_block_->pmy_mesh->dt; // MHD timestep
  // subcycling variables
  Real tend, tsub, tnow, tleft;
  Real icount;
  // loop over all cells
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        pmy_spec_->chemnet.InitializeNextStep(k, j, i);
        // copy species abundance
        for (int ispec=0; ispec<NSPECIES; ispec++) {
          y[ispec] = pmy_spec_->s(ispec,k,j,i)/u(IDN,k,j,i);
        }
        // assign internal energy, if not isothermal eos
        if (NON_BAROTROPIC_EOS) {
          E = u(IEN,k,j,i)
            - 0.5*( SQR(u(IM1,k,j,i)) + SQR(u(IM2,k,j,i)) + SQR(u(IM3,k,j,i))
                   )/u(IDN,k,j,i);
          if (MAGNETIC_FIELDS_ENABLED) {
            E -= 0.5*(
                SQR(bcc(IB1,k,j,i)) + SQR(bcc(IB2,k,j,i)) + SQR(bcc(IB3,k,j,i)) );
          }
        }
        // subcycling
        icount = 0;
        tnow = time;
        tend = time + dt_mhd;
        tleft = dt_mhd;
        while (tnow < tend) {
          // calculate reaction rates
          pmy_spec_->chemnet.RHS(time, y, E, ydot);
          // calculate heating and cooling rats
          if (NON_BAROTROPIC_EOS) {
            Edot = pmy_spec_->chemnet.Edot(time, y, E);
          }
          // advance one substep
          tsub = cfl_cool_sub * GetChemTime(y, ydot, E, Edot);
          tsub = std::min(tsub, tleft);
          IntegrateOneSubstep(tsub, y, ydot, E, Edot);
          // update timing
          tnow += tsub;
          tleft = tend - tnow;
          icount++;
          // check maximum number of steps
          if (icount > nsub_max) {
            std::stringstream msg;
            msg << "### FATAL ERROR in function ODEWrapper::Integrate: "
              << "Maximum number of substeps = " << nsub_max
              << " exceeded for forward Euler solver." << std::endl;
            ATHENA_ERROR(msg);
          }
        }
        // copy species abundance back to s
        for (int ispec=0; ispec<NSPECIES; ispec++) {
          // apply floor to passive scalar concentrations
          y[ispec] = (y[ispec] < scalar_floor) ?  scalar_floor : y[ispec];
          pmy_spec_->s(ispec,k,j,i) = y[ispec]*u(IDN,k,j,i);
        }
        // assign internal energy, if not isothermal eos
        if (NON_BAROTROPIC_EOS) {
          u(IEN,k,j,i) = E
            + 0.5*( SQR(u(IM1,k,j,i)) + SQR(u(IM2,k,j,i)) + SQR(u(IM3,k,j,i))
                )/u(IDN,k,j,i);
          if (MAGNETIC_FIELDS_ENABLED) {
            u(IEN,k,j,i) += 0.5*(
                SQR(bcc(IB1,k,j,i)) + SQR(bcc(IB2,k,j,i)) + SQR(bcc(IB3,k,j,i)) );
          }
        }
      }
    }
  }
  if (output_zone_sec_) {
    tstop = std::clock();
    double cpu_time = (tstop>tstart ? static_cast<double> (tstop-tstart) :
                       1.0)/static_cast<double> (CLOCKS_PER_SEC);
    std::uint64_t nzones =
      static_cast<std::uint64_t> (pmy_block_->GetNumberOfMeshBlockCells());
    // double zone_sec = static_cast<double> (nzones) / cpu_time;
    printf("chemistry ODE integration: ");
    printf("ncycle = %d, total time in sec = %.2e, zone/sec=%.2e\n",
        ncycle, cpu_time, Real(nzones)/cpu_time);
  }
  delete[] y;
  delete[] ydot;
  return;
}

namespace {

//----------------------------------------------------------------------------------------
//! \fn Real GetChemTime(const Real *y, const Real *ydot,
//                       const Real E, const Real Edot)
//! \brief calculate chemistry timescale
Real GetChemTime(const Real *y, const Real *ydot, const Real E, const Real Edot) {
  const Real small_ = 1024 * std::numeric_limits<float>::min();
  // put floor in species abundance
  std::array<Real, NSPECIES> yf;
  for (int ispec=0; ispec<NSPECIES; ispec++) {
    yf[ispec] = std::max(y[ispec], yfloor);
  }
  // calculate chemistry timescale
  Real tchem = std::abs( yf[0]/(ydot[0] + small_) );
  for (int ispec=1; ispec<NSPECIES; ispec++) {
    tchem = std::min( tchem, std::abs(yf[ispec]/(ydot[ispec]+small_)) );
  }
  if (NON_BAROTROPIC_EOS) {
    tchem = std::min( tchem, std::abs(E/(Edot+small_)) );
  }
  return tchem;
}

//----------------------------------------------------------------------------------------
//! \fn void IntegrateOneSubstep(Real tsub, Real *y, const Real *ydot,
//!                              Real &E, const Real Edot)
//! \brief advance chemical abundance and energy for tsub
void IntegrateOneSubstep(Real tsub, Real *y, const Real *ydot,
                         Real &E, const Real Edot) {
  for (int ispec=0; ispec<NSPECIES; ispec++) {
    y[ispec] += ydot[ispec] * tsub;
  }
  if (NON_BAROTROPIC_EOS) {
    E += Edot * tsub;
  }
  return;
}

} // namespace
