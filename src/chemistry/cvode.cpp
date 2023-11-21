//=======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cvode.cpp
//! \brief implementation of the CVODE solver

// C header

// C++ header
#include <ctime>
#include <stdexcept> // throw exceptions
#include <string>

// Athena++ classes headers
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

// this class header
#include "ode_wrapper.hpp"

namespace {
  // solver tolerance
  Real reltol_; // relative tolerance
  AthenaArray<Real> abstol_;
  // solver stepsize
  Real h_init_;
  bool use_previous_h_;
  // factor of the max timestep in solver relative to the hydrostep
  Real fac_dtmax_;
  void CheckFlag(const void *flagvalue, const char *funcname,
                 const int opt);
} // namespace

//----------------------------------------------------------------------------------------
//! \brief ODEWrapper constructor
ODEWrapper::ODEWrapper(MeshBlock *pmb, ParameterInput *pin) {
  int flag;
  pmy_block_ = pmb;
  dense_matrix_ = NULL;
  dense_ls_ = NULL;
  if (NON_BAROTROPIC_EOS) {
    dim_ = NSPECIES + 1;
  } else {
    dim_ = NSPECIES;
  }
  abstol_.NewAthenaArray(dim_);
  // Create the SUNDIALS context
  flag = SUNContext_Create(NULL, &sunctx_);
  CheckFlag(&flag, "SUNContext_Create", 1);
  // allocate y_
  y_ = N_VNew_Serial(dim_, sunctx_);
  CheckFlag(static_cast<void *>(y_), "N_VNew_Serial", 0);
  ydata_ = NV_DATA_S(y_);
  reltol_ = pin->GetOrAddReal("chemistry", "reltol", 1.0e-2);
  output_zone_sec_ = pin->GetOrAddBoolean("chemistry", "output_zone_sec", false);
  fac_dtmax_ = pin->GetOrAddReal("chemistry", "fac_dtmax", 10.);
}

//----------------------------------------------------------------------------------------
//! \brief ODEWrapper destructor
ODEWrapper::~ODEWrapper() {
  NV_DATA_S(y_) = ydata_;
  // Free y_ vector
  N_VDestroy(y_);
  // Free integrator memory
  CVodeFree(&cvode_mem_);
  // Free linear solver memory
  SUNLinSolFree(dense_ls_);
  // Free matrix memory
  SUNMatDestroy(dense_matrix_);
  // Free SUNDIALS context
  SUNContext_Free(&sunctx_);
}

//----------------------------------------------------------------------------------------
//! \fn void ODEWrapper::Initialize(ParameterInput *pin)
//! \brief Initialize ODE solver parameters
void ODEWrapper::Initialize(ParameterInput *pin) {
  // Note: this cannot be in the constructor, since it needs the PassiveScalars
  // class, and the ODEWrapper class is constructed in the PassiveScalars constructor.
  int flag;
  pmy_spec_ = pmy_block_->pscalars;

  // tolerance
  Real abstol_all = pin->GetOrAddReal("chemistry", "abstol", 1.0e-12);
  for (int i=0; i<NSPECIES; i++) {
    abstol_(i) = pin->GetOrAddReal("chemistry",
        "abstol_"+pmy_spec_->chemnet.species_names[i], -1);
    if (abstol_(i) < 0) {
      abstol_(i) = abstol_all;
    }
  }
  if (NON_BAROTROPIC_EOS) {
    abstol_(dim_-1) = pin->GetOrAddReal("chemistry", "abstol_E", -1);
    if (abstol_(dim_-1) < 0) {
      abstol_(dim_-1) = abstol_all;
    }
  }
  // read initial step, default -1 uses CVODE estimation
  h_init_ = pin->GetOrAddReal("chemistry", "h_init", -1);
  // choice of using the previous stepsize or default estimation from CVODE
  use_previous_h_ = pin->GetOrAddBoolean("chemistry", "use_previous_h", true);
  // user input Jacobian flag
  int user_jac = pin->GetOrAddBoolean("chemistry", "user_jac", false);
  // maximum number of steps
  int64_t maxsteps = pin->GetOrAddInteger("chemistry", "maxsteps", 1000);
  // maximum order
  int maxorder = pin->GetOrAddInteger("chemistry", "maxorder", 2);
  // maximum number of error test fails
  int maxerrtest = pin->GetOrAddInteger("chemistry", "maxerrtest", 7);
  // maximum number of t = h+t warnings
  int maxhnil = pin->GetOrAddInteger("chemistry", "maxhnil", 10);
  // stability limit detection
  int stldet = pin->GetOrAddInteger("chemistry", "stldet", 0);

  // -----------Initialize absolute value vector----------
  N_Vector abstol_vec = N_VNew_Serial(dim_, sunctx_);
  CheckFlag(static_cast<void *>(abstol_vec), "N_VNew_Serial", 0);
  for (int i=0; i<dim_; i++) {
    NV_Ith_S(abstol_vec, i) = abstol_(i);
  }

  //-------------initialize CVODE------------------
  // Call CVodeCreate to create the solver memory and specify the
  // Backward Differentiation Formula and the use of a Newton iteration
  cvode_mem_ = CVodeCreate(CV_BDF, sunctx_);
  CheckFlag(static_cast<void *>(cvode_mem_), "CVodeCreate", 0);

  // Set the user data pointer to NetworkWrapper
  flag = CVodeSetUserData(cvode_mem_, &(pmy_spec_->chemnet));
  CheckFlag(&flag, "CVodeSetUserData", 1);

  // Call CVodeInit to initialize the integrator memory and specify the
  // user's right hand side function in y'=f(t,y), the inital time T0, and
  // the initial dependent variable vector y.
  flag = CVodeInit(cvode_mem_,  pmy_spec_->chemnet.WrapRHS,
                   pmy_block_->pmy_mesh->time, y_);
  CheckFlag(&flag, "CVodeInit", 1);

  // Call CVodeSVtolerances to specify the scalar relative tolerance
  // and vector absolute tolerances
  flag = CVodeSVtolerances(cvode_mem_, reltol_, abstol_vec);
  CheckFlag(&flag, "CVodeSVtolerances", 1);

  // Create dense SUNMatrix for use in linear solves
  dense_matrix_ = SUNDenseMatrix(dim_, dim_, sunctx_);
  CheckFlag(static_cast<void *>(dense_matrix_), "SUNDenseMatrix", 0);

  // Create dense SUNLinearSolver object for use by CVode
  dense_ls_ = SUNLinSol_Dense(y_, dense_matrix_, sunctx_);
  CheckFlag(static_cast<void *>(dense_ls_), "SUNDenseLinearSolver", 0);

  // Call CVodeSetLinearSolver to attach the matrix and linear solver
  flag = CVodeSetLinearSolver(cvode_mem_, dense_ls_, dense_matrix_);
  CheckFlag(&flag, "CVodeSetLinearSolver", 1);

  // set maximum number of steps
  flag = CVodeSetMaxNumSteps(cvode_mem_, maxsteps);
  CheckFlag(&flag, "CVodeSetMaxNumSteps", 1);

  // set maximum number of convergence failure
  flag = CVodeSetMaxConvFails(cvode_mem_, 100000);
  CheckFlag(&flag, "CVodeSetMaxNumSteps", 1);

  // set maximum order
  flag = CVodeSetMaxOrd(cvode_mem_, maxorder);
  CheckFlag(&flag, "CVodeSetMaxOrd", 1);

  // set maximum number of error test fails
  flag = CVodeSetMaxErrTestFails(cvode_mem_, maxerrtest);
  CheckFlag(&flag, "CVodeSetMaxErrTestFails", 1);

  // set maximum number of t+h=t warnings
  flag = CVodeSetMaxHnilWarns(cvode_mem_, maxhnil);
  CheckFlag(&flag, "CVodeSetMaxHnilWarns", 1);

  // set whether enable stability limit detection
  // Only used to reduce the order of the order is larger than 3
  flag = CVodeSetStabLimDet(cvode_mem_, stldet);
  CheckFlag(&flag, "CVodeSetStabLimDet", 1);

  // Set the Jacobian routine to Jac (user-supplied)
  if (user_jac) {
    flag = CVodeSetJacFn(cvode_mem_, pmy_spec_->chemnet.WrapJacobian);
    CheckFlag(&flag, "CVDlsSetDenseJacFn", 1);
  }

  // Free abstol_ vector
  N_VDestroy(abstol_vec);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ODEWrapper::Integrate(const Real tinit, const Real dt)
//! \brief Integrate the ODE forward for time dt
//!
//! Update abundance in PassiveScalars over time dt.
//! For each cell:
//! Step 1: Set the radiation field strength in ChemNetwork.
//! Depends on the data structure of radiation field, this can be copying
//! the value from ChemRadiation class to ChemNetwork class, or just pass a pointer.
//!
//! Step 2: re-initialize CVODE with starting time t, and starting abundance
//! y. If x(k, j, i, ispec), we can just pass a pointer to CVODE, otherwise,
//! we need to copy the abundance of PassiveScalars to an array.
//!
//! Step 3: Integration. Update the array of PassiveScalars abundance in that
//! cell over time dt.
//!
//! Note that this will be not vectorizable(?).

void ODEWrapper::Integrate(const Real tinit, const Real dt) {
  int is = pmy_block_->is;
  int js = pmy_block_->js;
  int ks = pmy_block_->ks;
  int ie = pmy_block_->ie;
  int je = pmy_block_->je;
  int ke = pmy_block_->ke;
  Real *pdata_r_copy = pmy_spec_->r_copy.data();
  int ncycle = pmy_block_->pmy_mesh->ncycle;
  Real tfinal = tinit + dt;
  Real treturn = 0;
  int flag;
  // primitive conserved variables
  AthenaArray<Real> &u = pmy_block_->phydro->u;
  AthenaArray<Real> &bcc = pmy_block_->pfield->bcc;
  // constants
  //const Real gm1 = pmy_block_->peos->GetGamma() - 1.0;
  const Real scalar_floor = pmy_block_->peos->GetScalarFloor();
  // timing of the chemistry in each cycle
  clock_t tstart=0.0;
  clock_t tstop;
  if (output_zone_sec_) {
    tstart = std::clock();
  }
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      // copy s to r_copy
      for (int ispec=0; ispec<NSPECIES; ispec++) {
        for (int i=is; i<=ie; ++i) {
          pmy_spec_->r_copy(i, ispec) = pmy_spec_->s(ispec,k,j,i)/u(IDN,k,j,i);
        }
      }
      // assign internal energy, if not isothermal eos
      if (NON_BAROTROPIC_EOS) {
        for (int i=is; i<=ie; ++i) {
          pmy_spec_->r_copy(i, NSPECIES) = u(IEN,k,j,i)
            - 0.5*( SQR(u(IM1,k,j,i)) + SQR(u(IM2,k,j,i)) + SQR(u(IM3,k,j,i))
                   )/u(IDN,k,j,i);
          if (MAGNETIC_FIELDS_ENABLED) {
            pmy_spec_->r_copy(i, NSPECIES) -= 0.5*(
                SQR(bcc(IB1,k,j,i)) + SQR(bcc(IB2,k,j,i)) + SQR(bcc(IB3,k,j,i)) );
          }
        }
      }
      // loop over each cell
      for (int i=is; i<=ie; ++i) {
        // step 1: initialize chemistry network, eg: density, radiation
        NV_DATA_S(y_) = pdata_r_copy + i*dim_;
        pmy_spec_->chemnet.InitializeNextStep(k, j, i);
        // step 2: re-initialize CVODE with starting time t, and vector y
        // allocate r_copy(i, *) to y_.
        // TODO(Gong): make sure Real and realtype are the same.
        flag = CVodeReInit(cvode_mem_, tinit, y_);
        CheckFlag(&flag, "CVodeReInit", 1);
        // set initial step
        if (ncycle == 0 && h_init_ > 0.) {
          SetInitStep(h_init_);
        }
        if (ncycle != 0 && use_previous_h_) {
          SetInitStep(pmy_spec_->h(k, j, i));
        }
        // set maximum step to be a factor times dt
        flag = CVodeSetMaxStep(cvode_mem_, dt*fac_dtmax_);
        CheckFlag(&flag, "CVodeSetMaxStep", 1);
        // step 3: integration. update array abundance over time dt
        // in CV_NORMAL model, treturn=tfinal (the time of output)
        flag = CVode(cvode_mem_, tfinal, y_, &treturn, CV_NORMAL);
        CheckFlag(&flag, "CVode", 3);

        // update next step size
        if (ncycle != 0) {
          pmy_spec_->h(k, j, i) = GetNextStep();
        }
      }

      // copy r_copy back to s
      for (int ispec=0; ispec<NSPECIES; ispec++) {
        for (int i=is; i<=ie; ++i) {
          Real& r_copy_i  = pmy_spec_->r_copy(i,ispec);
          // apply floor to passive scalar concentrations
          r_copy_i = (r_copy_i < scalar_floor) ?  scalar_floor : r_copy_i;
          pmy_spec_->s(ispec,k,j,i) = r_copy_i*u(IDN,k,j,i);
        }
      }

      // assign internal energy, if not isothermal eos
      if (NON_BAROTROPIC_EOS) {
        for (int i=is; i<=ie; ++i) {
          u(IEN,k,j,i) = pmy_spec_->r_copy(i, NSPECIES)
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
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ODEWrapper::SetInitStep(const Real h_init) const
//! \brief Set initial stepsize for ODE solver

void ODEWrapper::SetInitStep(const Real h_init) const {
  int flag = CVodeSetInitStep(cvode_mem_, h_init);
  CheckFlag(&flag, "CVodeSetInitStep", 1);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Real ODEWrapper::GetLastStep() const
//! \brief Get last stepsize

Real ODEWrapper::GetLastStep() const {
  Real hlast;
  int flag;
  flag = CVodeGetLastStep(cvode_mem_, &hlast);
  CheckFlag(&flag, "CVodeGetLastStep", 1);
  return hlast;
}

//----------------------------------------------------------------------------------------
//! \fn Real ODEWrapper::GetNextStep() const
//! \brief Get next stepsize
Real ODEWrapper::GetNextStep() const {
  Real hlast;
  int flag;
  flag = CVodeGetCurrentStep(cvode_mem_, &hlast);
  CheckFlag(&flag, "CVodeGetLastStep", 1);
  return hlast;
}

//----------------------------------------------------------------------------------------
//! \fn long int ODEWrapper::GetNsteps() const
//! \brief Get the number of steps between two reinits
long int ODEWrapper::GetNsteps() const { // NOLINT (runtime/int)
  long int nst; // NOLINT (runtime/int)
  int flag;
  flag = CVodeGetNumSteps(cvode_mem_, &nst);
  CheckFlag(&flag, "CVodeGetNumSteps", 1);
  return nst;
}

namespace {
//----------------------------------------------------------------------------------------
//! \fn void CheckFlag(const void *flagvalue, const char *funcname,
//!                    const int opt)
//! \brief CVODE flag check
void CheckFlag(const void *flagvalue, const char *funcname,
               const int opt) {
  const int *errflag;

  // Check if SUNDIALS function returned NULL pointer - no memory allocated
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
    throw std::runtime_error("SUNDIALS:Sundials error.");
    return;
  } else if (opt == 1) { // Check if flag < 0
    errflag = static_cast<const int *>(flagvalue);
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
          funcname, *errflag);
      throw std::runtime_error("SUNDIALS:Sundials error.");
      return;
    }
  } else if (opt == 2 && flagvalue == NULL) {
    // Check if function returned NULL pointer - no memory allocated
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
    throw std::runtime_error("SUNDIALS:Memory error.");
    return;
  } else if (opt == 3) { // Check if CV_SUCCESS for integration.
    errflag = static_cast<const int *>(flagvalue);
    if (*errflag != CV_SUCCESS) {
      fprintf(stderr, "\nCV_SUCCESS error: %s() failed with flag = %d\n\n",
          funcname, *errflag);
      throw std::runtime_error("SUNDIALS:CV_SUCCESS error");
      return;
    }
  }
}

} // namespace
