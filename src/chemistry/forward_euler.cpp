//=======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file forward_euler.cpp
//! \brief implementation of the forward Euler solver

// this class header
#include "ode_wrapper.hpp"

//c header
#include <stdio.h> //c style io

//c++ header
#include <ctime> //time
#include <stdexcept> //throw exceptions
#include <string>

// Athena++ classes headers
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

//----------------------------------------------------------------------------------------
//! \brief ODEWrapper constructor
ODEWrapper::ODEWrapper(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block_ = pmb;
  if (NON_BAROTROPIC_EOS) {
    dim_ = NSCALARS + 1;
  } else {
    dim_ = NSCALARS;
  }
  output_zone_sec_ = pin->GetOrAddBoolean("chemistry", "output_zone_sec", false);
}

//----------------------------------------------------------------------------------------
//! \brief ODEWrapper destructor
ODEWrapper::~ODEWrapper() {}

//----------------------------------------------------------------------------------------
//! \fn void ODEWrapper::Initialize(ParameterInput *pin)
//! \brief Initialize ODE solver parameters
void ODEWrapper::Initialize(ParameterInput *pin) {
  //TODO: initialization
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ODEWrapper::Integrate(const Real tinit, const Real dt)
//! \brief Integrate the ODE forward for time dt
void ODEWrapper::Integrate(const Real tinit, const Real dt) {
  clock_t tstart, tstop;
  tstart = std::clock();
  //TODO: integration
  tstop = std::clock();
  if (output_zone_sec_) {
    double cpu_time = (tstop>tstart ? static_cast<double> (tstop-tstart) :
                       1.0)/static_cast<double> (CLOCKS_PER_SEC);
    std::uint64_t nzones = static_cast<std::uint64_t> (pmy_block_->GetNumberOfMeshBlockCells());
    double zone_sec = static_cast<double> (nzones) / cpu_time;
    printf("chemistry ODE integration: ");
    printf("ncycle = %d, total time in sec = %.2e, zone/sec=%.2e\n",
        ncycle, cpu_time, Real(nzones)/cpu_time);
  }
  return;
}
