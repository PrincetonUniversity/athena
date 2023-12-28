//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file time_integrator.cpp
//! \brief derived class for time integrator task list. Can create task lists for one
//! of many different time integrators (e.g. van Leer, RK2, RK3, etc.)

// C headers

// C++ headers
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../bvals/bvals.hpp"
#include "../chemistry/network/network.hpp"
#include "../cr/cr.hpp"
#include "../cr/integrators/cr_integrators.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../scalars/scalars.hpp"
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
//! TimeIntegratorTaskList constructor

TimeIntegratorTaskList::TimeIntegratorTaskList(ParameterInput *pin, Mesh *pm) {
  //! \note
  //! First, define each time-integrator by setting weights for each step of
  //! the algorithm and the CFL number stability limit when coupled to the single-stage
  //! spatial operator.
  //! Currently, the explicit, multistage time-integrators must be expressed as 2S-type
  //! algorithms as in Ketcheson (2010) Algorithm 3, which incudes 2N (Williamson) and 2R
  //! (van der Houwen) popular 2-register low-storage RK methods. The 2S-type integrators
  //! depend on a bidiagonally sparse Shu-Osher representation; at each stage l:
  //! \f[
  //!   U^{l} = a_{l,l-2}*U^{l-2} + a_{l-1}*U^{l-1}
  //!         + b_{l,l-2}*dt*Div(F_{l-2}) + b_{l,l-1}*dt*Div(F_{l-1}),
  //! \f]
  //! where \f$U^{l-1}\f$ and \f$U^{l-2}\f$ are previous stages and
  //! \f$a_{l,l-2}\f$, \f$a_{l,l-1}=(1-a_{l,l-2})\f$,
  //! and \f$b_{l,l-2}\f$, \f$b_{l,l-1}\f$
  //! are weights that are different for each stage and
  //! integrator. Previous timestep \f$U^{0} = U^n\f$ is given, and the integrator solves
  //! for \f$U^{l}\f$ for 1 <= l <= nstages.
  //!
  //! \note
  //! The 2x RHS evaluations of Div(F) and source terms per stage is avoided by adding
  //! another weighted average / caching of these terms each stage. The API and framework
  //! is extensible to three register 3S* methods,
  //! although none are currently implemented.
  //!
  //! \note
  //! Notation: exclusively using "stage", equivalent in lit. to "substage" or "substep"
  //! (infrequently "step"), to refer to the intermediate values of U^{l} between each
  //! "timestep" = "cycle" in explicit, multistage methods. This is to disambiguate the
  //! temporal integration from other iterative sequences;  generic
  //! "Step" is often used for sequences in code, e.g. main.cpp: "Step 1: MPI"
  //!
  //! \note
  //! main.cpp invokes the tasklist in a for () loop from stage=1 to stage=ptlist->nstages
  //!
  //! \todo (felker):
  //! - validate Field and Hydro diffusion with RK3, RK4, SSPRK(5,4)
  integrator = pin->GetOrAddString("time", "integrator", "vl2");

  // nr_radiation enabled but not implicit_radiation
  bool radiation_flag = (NR_RADIATION_ENABLED && (!IM_RADIATION_ENABLED));

  // Read a flag for orbital advection
  ORBITAL_ADVECTION = (pm->orbital_advection != 0)? true : false;

  // Read a flag for shear periodic
  SHEAR_PERIODIC = pm->shear_periodic;

  if (integrator == "rk4" || integrator == "ssprk5_4") {
    // shear periodic not work with rk4 or ssprk5_4
    if (SHEAR_PERIODIC) {
      std::stringstream msg;
      msg << "### FATAL ERROR in TimeIntegratorTaskList constructor" << std::endl
          << "integrator=" << integrator << " does not work with shear periodic boundary."
          << std::endl;
      ATHENA_ERROR(msg);
    }
    // orbital advection + mesh refinement not work with rk4 or ssprk5_4
    if (ORBITAL_ADVECTION && pm->multilevel) {
      std::stringstream msg;
      msg << "### FATAL ERROR in TimeIntegratorTaskList constructor" << std::endl
          << "integrator=" << integrator << " does not work with orbital advection and "
          << "mesh refinement"
          << std::endl;
      ATHENA_ERROR(msg);
    }
  }

  if (integrator == "vl2") {
    //! \note `integrator == "vl2"`
    //! - VL: second-order van Leer integrator (Stone & Gardiner, NewA 14, 139 2009)
    //! - Simple predictor-corrector scheme similar to MUSCL-Hancock
    //! - Expressed in 2S or 3S* algorithm form

    // set number of stages and time coeff.
    nstages_main = 2;
    if (ORBITAL_ADVECTION) {
      // w/ orbital advection
      if (SHEAR_PERIODIC || pm->multilevel) {
        // w/ shear_periodic or refinements
        if (pm->orbital_advection==1) { // first order splitting
          nstages = nstages_main+1;
          for (int l=0; l<nstages; l++) {
            if (l == nstages-1) { // last stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].ebeta = 0.0;
            }
            stage_wghts[l].sbeta = 0.0;
          }
          stage_wghts[0].beta = 0.5;
          stage_wghts[1].beta = 1.0;
          stage_wghts[2].beta = 0.0;
        } else { // second order splitting
          nstages = nstages_main+2;
          for (int l=0; l<nstages; l++) {
            if (l == 0) { // first stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.0;
            } else if (l == nstages-1) { // last stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
            }
          }
          stage_wghts[0].beta = 0.0;
          stage_wghts[1].beta = 0.5;
          stage_wghts[2].beta = 1.0;
          stage_wghts[3].beta = 0.0;
        }
      } else { // w/o shear periodic and refinements
        if (pm->orbital_advection==1) { // first order splitting
          nstages = nstages_main;
          for (int l=0; l<nstages; l++) {
            stage_wghts[l].main_stage = true;
            if (l == nstages-1) { // last stage
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].ebeta = 0.0;
            }
            stage_wghts[l].sbeta = 0.0;
          }
          stage_wghts[0].beta = 0.5;
          stage_wghts[1].beta = 1.0;
        } else { // second order splitting
          nstages = nstages_main+1;
          for (int l=0; l<nstages; l++) {
            if (l == 0) { // first stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.0;
            } else if (l == nstages-1) { // last stage
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
            }
          }
          stage_wghts[0].beta = 0.0;
          stage_wghts[1].beta = 0.5;
          stage_wghts[2].beta = 1.0;
        }
      }
    } else { // w/o orbital advection
      nstages = nstages_main;
      for (int l=0; l<nstages; l++) {
        stage_wghts[l].main_stage = true;
        stage_wghts[l].orbital_stage = false;
      }
      stage_wghts[0].sbeta = 0.0;
      stage_wghts[0].ebeta = 0.5;
      stage_wghts[1].sbeta = 0.5;
      stage_wghts[1].ebeta = 1.0;
      stage_wghts[0].beta = 0.5;
      stage_wghts[1].beta = 1.0;
    }
    cfl_limit = 1.0;
    // Modify VL2 stability limit in 2D, 3D
    if (pm->ndim == 2) cfl_limit = 0.5;
    if (pm->ndim == 3) cfl_limit = 0.5;

    // set delta and gamma at each stage
    int n_main = 0;
    for (int n=0; n<nstages; n++) {
      if (stage_wghts[n].main_stage) {
        if (n_main == 0) {
          stage_wghts[n].delta = 1.0; // required for consistency
          stage_wghts[n].gamma_1 = 0.0;
          stage_wghts[n].gamma_2 = 1.0;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        } else if (n_main == 1) {
          stage_wghts[n].delta = 0.0;
          stage_wghts[n].gamma_1 = 0.0;
          stage_wghts[n].gamma_2 = 1.0;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        }
      }
    }
  } else if (integrator == "rk1") {
    //! \note `integrator == "rk1"`
    //! - RK1: first-order Runge-Kutta / the forward Euler (FE) method

    // set number of stages and time coeff.
    nstages_main = 1;
    if (ORBITAL_ADVECTION) {
      // w/ orbital advection
      if (SHEAR_PERIODIC || pm->multilevel) {
        // w/ shear_periodic or refinements
        if (pm->orbital_advection==1) { // first order splitting
          nstages = nstages_main+1;
          for (int l=0; l<nstages; l++) {
            if (l == nstages-1) { // last stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].ebeta = 0.0;
            }
            stage_wghts[l].sbeta = 0.0;
          }
          stage_wghts[0].beta = 1.0;
          stage_wghts[1].beta = 0.0;
        } else { // second order splitting
          nstages = nstages_main+2;
          for (int l=0; l<nstages; l++) {
            if (l == 0) { // first stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.0;
            } else if (l == nstages-1) { // last stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
            }
          }
          stage_wghts[0].beta = 0.0;
          stage_wghts[1].beta = 1.0;
          stage_wghts[2].beta = 0.0;
        }
      } else { // w/o shear periodic and refinements
        if (pm->orbital_advection==1) { // first order splitting
          nstages = nstages_main;
          stage_wghts[0].main_stage = true;
          stage_wghts[0].orbital_stage = true;
          stage_wghts[0].sbeta = 0.0;
          stage_wghts[0].ebeta = 1.0;
          stage_wghts[0].beta = 1.0;
        } else { // second order splitting
          nstages = nstages_main+1;
          stage_wghts[0].main_stage = false;
          stage_wghts[1].main_stage = true;
          stage_wghts[0].orbital_stage = true;
          stage_wghts[1].orbital_stage = true;
          stage_wghts[0].sbeta = 0.0;
          stage_wghts[0].ebeta = 0.5;
          stage_wghts[1].sbeta = 0.5;
          stage_wghts[1].ebeta = 1.0;
          stage_wghts[0].beta = 0.0;
          stage_wghts[1].beta = 1.0;
        }
      }
    } else { // w/o orbital advection
      nstages = nstages_main;
      stage_wghts[0].main_stage = true;
      stage_wghts[0].orbital_stage = false;
      stage_wghts[0].sbeta = 0.0;
      stage_wghts[0].ebeta = 1.0;
      stage_wghts[0].beta = 1.0;
    }
    cfl_limit = 1.0;

    // set delta and gamma at each stage
    int n_main = 0;
    for (int n=0; n<nstages; n++) {
      if (stage_wghts[n].main_stage) {
        if (n_main == 0) {
          stage_wghts[n].delta = 1.0;
          stage_wghts[n].gamma_1 = 0.0;
          stage_wghts[n].gamma_2 = 1.0;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        }
      }
    }
  } else if (integrator == "rk2") {
    //! \note `integrator == "rk2"`
    //! - Heun's method / SSPRK (2,2): Gottlieb (2009) equation 3.1
    //! - Optimal (in error bounds) explicit two-stage, second-order SSPRK

    // set number of stages and time coeff.
    nstages_main = 2;
    if (ORBITAL_ADVECTION) {
      // w/ orbital advection
      if (SHEAR_PERIODIC || pm->multilevel) {
        // w/ shear_periodic or refinements
        if (pm->orbital_advection==1) { // first order splitting
          nstages = nstages_main+1;
          for (int l=0; l<nstages; l++) {
            if (l == nstages-1) { // last stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].ebeta = 0.0;
            }
            stage_wghts[l].sbeta = 0.0;
          }
          stage_wghts[0].beta = 1.0;
          stage_wghts[1].beta = 0.5;
          stage_wghts[2].beta = 0.0;
        } else { // second order splitting
          nstages = nstages_main+2;
          for (int l=0; l<nstages; l++) {
            if (l == 0) { // first stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.0;
            } else if (l == nstages-1) { // last stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
            }
          }
          stage_wghts[0].beta = 0.0;
          stage_wghts[1].beta = 1.0;
          stage_wghts[2].beta = 0.5;
          stage_wghts[3].beta = 0.0;
        }
      } else { // w/o shear periodic and refinements
        if (pm->orbital_advection==1) { // first order splitting
          nstages = nstages_main;
          for (int l=0; l<nstages; l++) {
            stage_wghts[l].main_stage = true;
            if (l == nstages-1) { // last stage
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].ebeta = 0.0;
            }
            stage_wghts[l].sbeta = 0.0;
          }
          stage_wghts[0].beta = 1.0;
          stage_wghts[1].beta = 0.5;
        } else { // second order splitting
          nstages = nstages_main+1;
          for (int l=0; l<nstages; l++) {
            if (l == 0) { // first stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.0;
            } else if (l == nstages-1) { // last stage
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
            }
          }
          stage_wghts[0].beta = 0.0;
          stage_wghts[1].beta = 1.0;
          stage_wghts[2].beta = 0.5;
        }
      }
    } else { // w/o orbital advection
      nstages = nstages_main;
      for (int l=0; l<nstages; l++) {
        stage_wghts[l].main_stage = true;
        stage_wghts[l].orbital_stage = false;
      }
      stage_wghts[0].sbeta = 0.0;
      stage_wghts[0].ebeta = 1.0;
      stage_wghts[1].sbeta = 1.0;
      stage_wghts[1].ebeta = 1.0;
      stage_wghts[0].beta = 1.0;
      stage_wghts[1].beta = 0.5;
    }
    cfl_limit = 1.0;  // c_eff = c/nstages = 1/2 (Gottlieb (2009), pg 271)

    // set delta and gamma at each stage
    int n_main = 0;
    for (int n=0; n<nstages; n++) {
      if (stage_wghts[n].main_stage) {
        if (n_main == 0) {
          stage_wghts[n].delta = 1.0;
          stage_wghts[n].gamma_1 = 0.0;
          stage_wghts[n].gamma_2 = 1.0;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        } else if (n_main == 1) {
          stage_wghts[n].delta = 0.0;
          stage_wghts[n].gamma_1 = 0.5;
          stage_wghts[n].gamma_2 = 0.5;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        }
      }
    }
  } else if (integrator == "rk3") {
    //! \note `integrator == "rk3"`
    //! - SSPRK (3,3): Gottlieb (2009) equation 3.2
    //! - Optimal (in error bounds) explicit three-stage, third-order SSPRK

    // set number of stages and time coeff.
    nstages_main = 3;
    if (ORBITAL_ADVECTION) {
      // w/ orbital advection
      if (SHEAR_PERIODIC || pm->multilevel) {
        // w/ shear_periodic or refinements
        if (pm->orbital_advection==1) { // first order splitting
          nstages = nstages_main+1;
          for (int l=0; l<nstages; l++) {
            if (l == nstages-1) { // last stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].ebeta = 0.0;
            }
            stage_wghts[l].sbeta = 0.0;
          }
          stage_wghts[0].beta = 1.0;
          stage_wghts[1].beta = 0.25;
          stage_wghts[2].beta = TWO_3RD;
          stage_wghts[3].beta = 0.0;
        } else { // second order splitting
          nstages = nstages_main+2;
          for (int l=0; l<nstages; l++) {
            if (l == 0) { // first stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.0;
            } else if (l == nstages-1) { // last stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
            }
          }
          stage_wghts[0].beta = 0.0;
          stage_wghts[1].beta = 1.0;
          stage_wghts[2].beta = 0.25;
          stage_wghts[3].beta = TWO_3RD;
          stage_wghts[4].beta = 0.0;
        }
      } else { // w/o shear periodic and refinements
        if (pm->orbital_advection==1) { // first order splitting
          nstages = nstages_main;
          for (int l=0; l<nstages; l++) {
            stage_wghts[l].main_stage = true;
            if (l == nstages-1) { // last stage
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].ebeta = 0.0;
            }
            stage_wghts[l].sbeta = 0.0;
          }
          stage_wghts[0].beta = 1.0;
          stage_wghts[1].beta = 0.25;
          stage_wghts[2].beta = TWO_3RD;
        } else { // second order splitting
          nstages = nstages_main+1;
          for (int l=0; l<nstages; l++) {
            if (l == 0) { // first stage
              stage_wghts[l].main_stage = false;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.0;
            } else if (l == nstages-1) { // last stage
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = true;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
              stage_wghts[l].ebeta = 1.0;
            } else {
              stage_wghts[l].main_stage = true;
              stage_wghts[l].orbital_stage = false;
              stage_wghts[l].sbeta = 0.5;
              stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
            }
          }
          stage_wghts[0].beta = 0.0;
          stage_wghts[1].beta = 1.0;
          stage_wghts[2].beta = 0.25;
          stage_wghts[3].beta = TWO_3RD;
        }
      }
    } else { // w/o orbital advection
      nstages = nstages_main;
      for (int l=0; l<nstages; l++) {
        stage_wghts[l].main_stage = true;
        stage_wghts[l].orbital_stage = false;
      }
      stage_wghts[0].sbeta = 0.0;
      stage_wghts[0].ebeta = 1.0;
      stage_wghts[1].sbeta = 1.0;
      stage_wghts[1].ebeta = 0.5;
      stage_wghts[2].sbeta = 0.5;
      stage_wghts[2].ebeta = 1.0;
      stage_wghts[0].beta = 1.0;
      stage_wghts[1].beta = 0.25;
      stage_wghts[2].beta = TWO_3RD;
    }
    cfl_limit = 1.0;  // c_eff = c/nstages = 1/3 (Gottlieb (2009), pg 271)

    // set delta and gamma at each stage
    int n_main = 0;
    for (int n=0; n<nstages; n++) {
      if (stage_wghts[n].main_stage) {
        if (n_main == 0) {
          stage_wghts[n].delta = 1.0;
          stage_wghts[n].gamma_1 = 0.0;
          stage_wghts[n].gamma_2 = 1.0;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        } else if (n_main == 1) {
          stage_wghts[n].delta = 0.0;
          stage_wghts[n].gamma_1 = 0.25;
          stage_wghts[n].gamma_2 = 0.75;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        } else if (n_main == 2) {
          stage_wghts[n].delta = 0.0;
          stage_wghts[n].gamma_1 = TWO_3RD;
          stage_wghts[n].gamma_2 = ONE_3RD;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        }
      }
    }
  } else if (integrator == "rk4") {
    //! \note `integorator == "rk4"`
    //! - RK4()4[2S] from Table 2 of Ketcheson (2010)
    //! - Non-SSP, explicit four-stage, fourth-order RK
    //! - Stability properties are similar to classical (non-SSP) RK4
    //!   (but ~2x L2 principal error norm).
    //! - Refer to Colella (2011) for linear stability analysis of constant
    //!   coeff. advection of classical RK4 + 4th or 1st order (limiter engaged) fluxes
    nstages_main = 4;
    cfl_limit = 1.3925; // Colella (2011) eq 101; 1st order flux is most severe constraint

    if (ORBITAL_ADVECTION) {
      if (pm->orbital_advection==1) { // first order splitting
        nstages = nstages_main;
        for (int l=0; l<nstages; l++) {
          stage_wghts[l].main_stage = true;
          if (l == nstages-1) { // last stage
            stage_wghts[l].orbital_stage = true;
          } else {
            stage_wghts[l].orbital_stage = false;
          }
        }
        stage_wghts[0].beta = 1.193743905974738;
        stage_wghts[1].beta = 0.099279895495783;
        stage_wghts[2].beta = 1.131678018054042;
        stage_wghts[3].beta = 0.310665766509336;
      } else { // second order splitting
        nstages = nstages_main+1;
        for (int l=0; l<nstages; l++) {
          if (l == 0) { // first stage
            stage_wghts[l].main_stage = false;
            stage_wghts[l].orbital_stage = true;
          } else if (l == nstages-1) { // last stage
            stage_wghts[l].main_stage = true;
            stage_wghts[l].orbital_stage = true;
          } else {
            stage_wghts[l].main_stage = true;
            stage_wghts[l].orbital_stage = false;
          }
        }
        stage_wghts[0].beta = 0.0;
        stage_wghts[1].beta = 1.193743905974738;
        stage_wghts[2].beta = 0.099279895495783;
        stage_wghts[3].beta = 1.131678018054042;
        stage_wghts[4].beta = 0.310665766509336;
      }
    } else { // w/o orbital advection
      nstages = nstages_main;
      for (int l=0; l<nstages; l++) {
        stage_wghts[l].main_stage = true;
        stage_wghts[l].orbital_stage = false;
      }
      stage_wghts[0].beta = 1.193743905974738;
      stage_wghts[1].beta = 0.099279895495783;
      stage_wghts[2].beta = 1.131678018054042;
      stage_wghts[3].beta = 0.310665766509336;
    }

    // set delta and gamma at each stage
    int n_main = 0;
    for (int n=0; n<nstages; n++) {
      if (stage_wghts[n].main_stage) {
        if (n_main == 0) {
          stage_wghts[n].delta = 1.0;
          stage_wghts[n].gamma_1 = 0.0;
          stage_wghts[n].gamma_2 = 1.0;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        } else if (n_main == 1) {
          stage_wghts[n].delta = 0.217683334308543;
          stage_wghts[n].gamma_1 = 0.121098479554482;
          stage_wghts[n].gamma_2 = 0.721781678111411;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        } else if (n_main == 2) {
          stage_wghts[n].delta = 1.065841341361089;
          stage_wghts[n].gamma_1 = -3.843833699660025;
          stage_wghts[n].gamma_2 = 2.121209265338722;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        } else if (n_main == 3) {
          stage_wghts[n].delta = 0.0;
          stage_wghts[n].gamma_1 = 0.546370891121863;
          stage_wghts[n].gamma_2 = 0.198653035682705;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        }
      }
    }

    // set sbeta & ebeta
    if (ORBITAL_ADVECTION) {
      if (pm->orbital_advection==1) { // first order splitting
        for (int l=0; l<nstages; l++) {
          if (l == nstages-1) { // last stage
            stage_wghts[l].ebeta = 1.0;
          } else {
            stage_wghts[l].ebeta = 0.0;
          }
          stage_wghts[l].sbeta = 0.0;
        }
      } else { // second order splitting
        for (int l=0; l<nstages; l++) {
          if (l == 0) { // first stage
            stage_wghts[l].sbeta = 0.0;
          } else if (l == nstages-1) { // last stage
            stage_wghts[l].sbeta = 0.5;
            stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
            stage_wghts[l].ebeta = 1.0;
          } else {
            stage_wghts[l].sbeta = 0.5;
            stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
          }
        }
      }
    } else { // w/o orbital advection ///
      Real temp = 0.0;
      Real temp_prev = 0.0;
      stage_wghts[0].sbeta = 0.0;
      for (int l=0; l<nstages-1; l++) {
        temp_prev = temp;
        temp = temp_prev + stage_wghts[l].delta*stage_wghts[l].sbeta;
        stage_wghts[l].ebeta = stage_wghts[l].gamma_1*temp_prev
                               + stage_wghts[l].gamma_2*temp
                               + stage_wghts[l].gamma_3*0.0
                               + stage_wghts[l].beta;
        stage_wghts[l+1].sbeta = stage_wghts[l].ebeta;
      }
      stage_wghts[nstages-1].ebeta = 1.0;
    }
  } else if (integrator == "ssprk5_4") {
    //! \note `integrator == "ssprk5_4"`
    //! - SSPRK (5,4): Gottlieb (2009) section 3.1; between eq 3.3 and 3.4
    //! - Optimal (in error bounds) explicit five-stage, fourth-order SSPRK
    //!   3N method, but there is no 3S* formulation due to irregular sparsity
    //!   of Shu-Osher form matrix, alpha.
    //! - Because it is an SSP method, we can use the SSP coefficient c=1.508 to to
    //!   trivially relate the CFL constraint to the RK1 CFL=1 (for first-order fluxes).
    //! - There is no need to perform stability analysis from scratch
    //!   (unlike e.g. the linear stability analysis for classical/non-SSP RK4 in
    //!   Colella (2011)).
    //! - However, PLM and PPM w/o the limiter engaged are unconditionally unstable
    //!   under RK1 integration, so the SSP guarantees do not hold for
    //!   the Athena++ spatial discretizations.
    nstages_main = 5;
    cfl_limit = 1.508;         //  (effective SSP coeff = 0.302) Gottlieb (2009) pg 272

    if (ORBITAL_ADVECTION) {
      if (pm->orbital_advection==1) { // first order splitting
        nstages = nstages_main;
        for (int l=0; l<nstages; l++) {
          stage_wghts[l].main_stage = true;
          if (l == nstages-1) { // last stage
            stage_wghts[l].orbital_stage = true;
          } else {
            stage_wghts[l].orbital_stage = false;
          }
        }
        stage_wghts[0].beta = 0.391752226571890;
        stage_wghts[1].beta = 0.368410593050371;
        stage_wghts[2].beta = 0.251891774271694;
        stage_wghts[3].beta = 0.544974750228521;
        stage_wghts[4].beta = 0.226007483236906; // F(u^(4)) coeff.
      } else { // second order splitting
        nstages = nstages_main+1;
        for (int l=0; l<nstages; l++) {
          if (l == 0) { // first stage
            stage_wghts[l].main_stage = false;
            stage_wghts[l].orbital_stage = true;
          } else if (l == nstages-1) { // last stage
            stage_wghts[l].main_stage = true;
            stage_wghts[l].orbital_stage = true;
          } else {
            stage_wghts[l].main_stage = true;
            stage_wghts[l].orbital_stage = false;
          }
        }
        stage_wghts[0].beta = 0.0;
        stage_wghts[1].beta = 0.391752226571890;
        stage_wghts[2].beta = 0.368410593050371;
        stage_wghts[3].beta = 0.251891774271694;
        stage_wghts[4].beta = 0.544974750228521;
        stage_wghts[5].beta = 0.226007483236906; // F(u^(4)) coeff.
      }
    } else { // w/o orbital advection
      nstages = nstages_main;
      for (int l=0; l<nstages; l++) {
        stage_wghts[l].main_stage = true;
        stage_wghts[l].orbital_stage = false;
      }
      stage_wghts[0].beta = 0.391752226571890;
      stage_wghts[1].beta = 0.368410593050371;
      stage_wghts[2].beta = 0.251891774271694;
      stage_wghts[3].beta = 0.544974750228521;
      stage_wghts[4].beta = 0.226007483236906; // F(u^(4)) coeff.
    }

    // set delta and gamma at each stage
    int n_main = 0;
    for (int n=0; n<nstages; n++) {
      if (stage_wghts[n].main_stage) {
        if (n_main == 0) {
          stage_wghts[n].delta = 1.0; // u1 = u^n
          stage_wghts[n].gamma_1 = 0.0;
          stage_wghts[n].gamma_2 = 1.0;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        } else if (n_main == 1) {
          stage_wghts[n].delta = 0.0; // u1 = u^n
          stage_wghts[n].gamma_1 = 0.555629506348765;
          stage_wghts[n].gamma_2 = 0.444370493651235;
          stage_wghts[n].gamma_3 = 0.0;
          n_main++;
        } else if (n_main == 2) {
          stage_wghts[n].delta = 0.517231671970585; // u1 <- (u^n + d*u^(2))
          stage_wghts[n].gamma_1 = 0.379898148511597;
          stage_wghts[n].gamma_2 = 0.0;
          stage_wghts[n].gamma_3 = 0.620101851488403; // u^(n) coeff =  u2
          n_main++;
        } else if (n_main == 3) {
          stage_wghts[n].delta = 0.096059710526147; // u1 <- (u^n + d*u^(2) + d'*u^(3))
          stage_wghts[n].gamma_1 = 0.821920045606868;
          stage_wghts[n].gamma_2 = 0.0;
          stage_wghts[n].gamma_3 = 0.178079954393132; // u^(n) coeff =  u2
          n_main++;
        } else if (n_main == 4) {
          stage_wghts[n].delta = 0.0;
          // 1 ulp lower than Gottlieb u^(4) coeff
          stage_wghts[n].gamma_1 = 0.386708617503268;
          // u1 <- (u^n + d*u^(2) + d'*u^(3))
          stage_wghts[n].gamma_2 = 1.0;
          // partial sum from hardcoded extra stage=4
          stage_wghts[n].gamma_3 = 1.0;
          n_main++;
        }
      }
    }

    // set sbeta & ebeta
    if (ORBITAL_ADVECTION) {
      if (pm->orbital_advection==1) { // first order splitting
        for (int l=0; l<nstages; l++) {
          if (l == nstages-1) { // last stage
            stage_wghts[l].ebeta = 1.0;
          } else {
            stage_wghts[l].ebeta = 0.0;
          }
          stage_wghts[l].sbeta = 0.0;
        }
      } else { // second order splitting
        for (int l=0; l<nstages; l++) {
          if (l == 0) { // first stage
            stage_wghts[l].sbeta = 0.0;
          } else if (l == nstages-1) { // last stage
            stage_wghts[l].sbeta = 0.5;
            stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
            stage_wghts[l].ebeta = 1.0;
          } else {
            stage_wghts[l].sbeta = 0.5;
            stage_wghts[l-1].ebeta = stage_wghts[l].sbeta;
          }
        }
      }
    } else { // w/o orbital advection ///
      Real temp = 0.0;
      Real temp_prev = 0.0;
      stage_wghts[0].sbeta = 0.0;
      for (int l=0; l<nstages-1; l++) {
        temp_prev = temp;
        temp = temp_prev + stage_wghts[l].delta*stage_wghts[l].sbeta;
        stage_wghts[l].ebeta = stage_wghts[l].gamma_1*temp_prev
                               + stage_wghts[l].gamma_2*temp
                               + stage_wghts[l].gamma_3*0.0
                               + stage_wghts[l].beta;
        stage_wghts[l+1].sbeta = stage_wghts[l].ebeta;
      }
      stage_wghts[nstages-1].ebeta = 1.0;
    }
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in TimeIntegratorTaskList constructor" << std::endl
        << "integrator=" << integrator << " not valid time integrator" << std::endl;
    ATHENA_ERROR(msg);
  }

  // Set cfl_number based on user input and time integrator CFL limit
  Real cfl_number = pin->GetReal("time", "cfl_number");
  if (cfl_number > cfl_limit
      && pm->fluid_setup == FluidFormulation::evolve) {
    std::cout << "### Warning in TimeIntegratorTaskList constructor" << std::endl
              << "User CFL number " << cfl_number << " must be smaller than " << cfl_limit
              << " for integrator=" << integrator << " in " << pm->ndim
              << "D simulation" << std::endl << "Setting to limit" << std::endl;
    cfl_number = cfl_limit;
  }
  // Save to Mesh class
  pm->cfl_number = cfl_number;

  // Now assemble list of tasks for each stage of time integrator
  {using namespace HydroIntegratorTaskNames; // NOLINT (build/namespace)
    // calculate hydro/field diffusive fluxes
    if (!STS_ENABLED) {
      AddTask(DIFFUSE_HYD,NONE);
      if (MAGNETIC_FIELDS_ENABLED) {
        AddTask(DIFFUSE_FLD,NONE);
        // compute hydro fluxes, integrate hydro variables
        AddTask(CALC_HYDFLX,(DIFFUSE_HYD|DIFFUSE_FLD));
      } else { // Hydro
        AddTask(CALC_HYDFLX,DIFFUSE_HYD);
      }
      if (NSCALARS > 0) {
        AddTask(DIFFUSE_SCLR,NONE);
        AddTask(CALC_SCLRFLX,(CALC_HYDFLX|DIFFUSE_SCLR));
      }
    } else { // STS enabled:
      AddTask(CALC_HYDFLX,NONE);
      if (NSCALARS > 0)
        AddTask(CALC_SCLRFLX,CALC_HYDFLX);
    }
    if (pm->multilevel || SHEAR_PERIODIC) { // SMR or AMR or shear periodic
      AddTask(SEND_HYDFLX,CALC_HYDFLX);
      AddTask(RECV_HYDFLX,CALC_HYDFLX);
      if (SHEAR_PERIODIC) {
        AddTask(SEND_HYDFLXSH,RECV_HYDFLX);
        AddTask(RECV_HYDFLXSH,(SEND_HYDFLX|RECV_HYDFLX));
        AddTask(INT_HYD,RECV_HYDFLXSH);
      } else {
        AddTask(INT_HYD,RECV_HYDFLX);
      }
    } else {
      AddTask(INT_HYD, CALC_HYDFLX);
    }

    if (radiation_flag) {
      AddTask(CALC_RADFLX,NONE);
      if (pm->multilevel || SHEAR_PERIODIC) { // SMR or AMR or shear periodic
        AddTask(SEND_RADFLX,CALC_RADFLX);
        AddTask(RECV_RADFLX,CALC_RADFLX);
        if (SHEAR_PERIODIC) {
          AddTask(SEND_RADFLXSH,RECV_RADFLX);
          AddTask(RECV_RADFLXSH,(SEND_RADFLX|RECV_RADFLX));
          AddTask(INT_RAD,RECV_RADFLXSH);
        } else {
          AddTask(INT_RAD,RECV_RADFLX);
        }
      } else {
        AddTask(INT_RAD, CALC_RADFLX);
      }
      AddTask(SRCTERM_RAD,INT_RAD);
      AddTask(SEND_RAD,SRCTERM_RAD);
      AddTask(RECV_RAD,NONE);
      AddTask(SETB_RAD,(RECV_RAD|SRCTERM_RAD));
    }

    if (CR_ENABLED) {
      AddTask(CALC_CRTCFLX,NONE);
      if (pm->multilevel) { // SMR or AMR
        AddTask(SEND_CRTCFLX,CALC_CRTCFLX);
        AddTask(RECV_CRTCFLX,CALC_CRTCFLX);
        AddTask(INT_CRTC,RECV_CRTCFLX);
      } else {
        AddTask(INT_CRTC, CALC_CRTCFLX);
      }
      AddTask(SRCTERM_CRTC,INT_CRTC);
      AddTask(SEND_CRTC,SRCTERM_CRTC);
      AddTask(RECV_CRTC,NONE);
      AddTask(SETB_CRTC,(RECV_CRTC|SRCTERM_CRTC));
    }

    if (NSCALARS > 0) {
      AddTask(SRC_TERM,(INT_HYD|INT_SCLR|INT_CHM));
    } else {
      AddTask(SRC_TERM,INT_HYD);
    }

    // Hydro will also be updated with radiation source term
    TaskID src_aterm = SRC_TERM;
    if (radiation_flag)
      src_aterm = (src_aterm | SRCTERM_RAD);

    if (CR_ENABLED)
      src_aterm = (src_aterm | SRCTERM_CRTC);

    if (ORBITAL_ADVECTION) {
      AddTask(SEND_HYDORB,src_aterm);
      AddTask(RECV_HYDORB,NONE);
      AddTask(CALC_HYDORB,(SEND_HYDORB|RECV_HYDORB));
      AddTask(SEND_HYD,CALC_HYDORB);
      AddTask(RECV_HYD,NONE);
      AddTask(SETB_HYD,(RECV_HYD|CALC_HYDORB));
    } else {
      AddTask(SEND_HYD,src_aterm);
      AddTask(RECV_HYD,NONE);
      AddTask(SETB_HYD,(RECV_HYD|SRC_TERM));
    }

    if (SHEAR_PERIODIC) {
      AddTask(SEND_HYDSH,SETB_HYD);
      AddTask(RECV_HYDSH,SEND_HYDSH);

      if (radiation_flag) {
        AddTask(SEND_RADSH,SETB_RAD);
        // shearing periodic boundary of radiation requires
        // shearing velocity to be set first
        AddTask(RECV_RADSH,SEND_RADSH|RECV_HYDSH);
      }
    }

    if (NSCALARS > 0) {
      if (pm->multilevel || SHEAR_PERIODIC) {
        AddTask(SEND_SCLRFLX,CALC_SCLRFLX);
        AddTask(RECV_SCLRFLX,CALC_SCLRFLX);
        if (SHEAR_PERIODIC) {
          AddTask(SEND_SCLRFLXSH,RECV_SCLRFLX);
          AddTask(RECV_SCLRFLXSH,(SEND_SCLRFLX|RECV_SCLRFLX));
          AddTask(INT_SCLR,RECV_SCLRFLXSH);
        } else {
          AddTask(INT_SCLR,RECV_SCLRFLX);
        }
      } else {
        AddTask(INT_SCLR,CALC_SCLRFLX);
      }
      AddTask(INT_CHM,INT_SCLR);
      if (ORBITAL_ADVECTION) {
        AddTask(SEND_SCLR,CALC_HYDORB);
        AddTask(RECV_SCLR,NONE);
        AddTask(SETB_SCLR,(RECV_SCLR|CALC_HYDORB));
      } else {
        AddTask(SEND_SCLR,SRC_TERM);
        AddTask(RECV_SCLR,NONE);
        AddTask(SETB_SCLR,(RECV_SCLR|SRC_TERM));
      }
      if (SHEAR_PERIODIC) {
        AddTask(SEND_SCLRSH,SETB_SCLR);
        AddTask(RECV_SCLRSH,SEND_SCLRSH);
      }
    }

    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      // compute MHD fluxes, integrate field
      AddTask(CALC_FLDFLX,CALC_HYDFLX);
      AddTask(SEND_FLDFLX,CALC_FLDFLX);
      AddTask(RECV_FLDFLX,SEND_FLDFLX);
      if (SHEAR_PERIODIC) {
        AddTask(SEND_EMFSH,RECV_FLDFLX);
        AddTask(RECV_EMFSH,RECV_FLDFLX);
        AddTask(INT_FLD,RECV_EMFSH);
      } else {
        AddTask(INT_FLD,RECV_FLDFLX);
      }

      if (ORBITAL_ADVECTION) {
        AddTask(SEND_FLDORB,INT_FLD);
        AddTask(RECV_FLDORB,NONE);
        AddTask(CALC_FLDORB,(SEND_FLDORB|RECV_FLDORB));
        AddTask(SEND_FLD,CALC_FLDORB);
        AddTask(RECV_FLD,NONE);
        AddTask(SETB_FLD,(RECV_FLD|CALC_FLDORB));
      } else {
        AddTask(SEND_FLD,INT_FLD);
        AddTask(RECV_FLD,NONE);
        AddTask(SETB_FLD,(RECV_FLD|INT_FLD));
      }
      if (SHEAR_PERIODIC) {
        AddTask(SEND_FLDSH,SETB_FLD);
        AddTask(RECV_FLDSH,SEND_FLDSH);
      }

      // TODO(felker): these nested conditionals are horrible now. Add option to AddTask
      // for "wait for all previously added tasks"?

      // prolongate, compute new primitives
      if (pm->multilevel) { // SMR or AMR
       TaskID setb=(SEND_HYD|SEND_FLD);
        if (SHEAR_PERIODIC) {
          if (NSCALARS > 0) {
            setb=(setb|RECV_HYDSH|SEND_SCLR|RECV_SCLRSH|RECV_FLDSH);
          } else {
            setb=(setb|RECV_HYDSH|RECV_FLDSH);
          }
        } else {
          if (NSCALARS > 0) {
            setb=(setb|SETB_HYD|SETB_FLD|SEND_SCLR|SETB_SCLR);
          } else {
            setb=(setb|SETB_HYD|SETB_FLD);
          }
        }

        if (radiation_flag) {
          setb=(setb|SEND_RAD|SETB_RAD);
          if (SHEAR_PERIODIC)
            setb=(setb|RECV_RADSH);
        }

        if (CR_ENABLED) {
          setb=(setb|SEND_CRTC|SETB_CRTC);
        }

        AddTask(PROLONG,setb);
        AddTask(CONS2PRIM,PROLONG);
      } else {
        if (SHEAR_PERIODIC) {
          if (NSCALARS > 0) {
            AddTask(CONS2PRIM,(RECV_HYDSH|RECV_FLDSH|RECV_SCLRSH));
          } else {
            AddTask(CONS2PRIM,(RECV_HYDSH|RECV_FLDSH));
          }
        } else {
          if (NSCALARS > 0) {
            AddTask(CONS2PRIM,(SETB_HYD|SETB_FLD|SETB_SCLR));
          } else {
            AddTask(CONS2PRIM,(SETB_HYD|SETB_FLD));
          }
        }
      }
    } else {  // HYDRO
      // prolongate, compute new primitives
      if (pm->multilevel) { // SMR or AMR
        TaskID setb=(SEND_HYD);
        if (SHEAR_PERIODIC) {
          if (NSCALARS > 0) {
            setb=(setb|RECV_HYDSH|SEND_SCLR|RECV_SCLRSH);
          } else {
            setb=(setb|RECV_HYDSH);
          }
        } else {
          if (NSCALARS > 0) {
            setb=(setb|SETB_HYD|SEND_SCLR|SETB_SCLR);
          } else {
            setb=(setb|SETB_HYD);
          }
        }

        if (radiation_flag) {
          setb=(setb|SEND_RAD|SETB_RAD);
          if (SHEAR_PERIODIC)
            setb=(setb|RECV_RADSH);
        }

        if (CR_ENABLED)
          setb=(setb|SEND_CRTC|SETB_CRTC);

        AddTask(PROLONG,setb);
        AddTask(CONS2PRIM,PROLONG);
      } else {
        if (SHEAR_PERIODIC) {
          if (NSCALARS > 0) {
            AddTask(CONS2PRIM,(RECV_HYDSH|RECV_SCLRSH));
          } else {
            AddTask(CONS2PRIM,RECV_HYDSH);
          }
        } else {
          if (NSCALARS > 0) {
            AddTask(CONS2PRIM,(SETB_HYD|SETB_SCLR));
          } else {
            AddTask(CONS2PRIM,SETB_HYD);
          }
        }
      }
    }

    // everything else

    TaskID before_bval = CONS2PRIM;
    TaskID before_userwork = PHY_BVAL;
    if (radiation_flag) {
      before_bval = (before_bval|SETB_RAD|SEND_RAD);
      if (SHEAR_PERIODIC)
        before_bval = (before_bval|RECV_RADSH);

      before_userwork = (before_userwork|RAD_MOMOPACITY);
    }

    if (CR_ENABLED) {
      before_bval = (before_bval|SETB_CRTC|SEND_CRTC);
      before_userwork = (before_userwork|CRTC_OPACITY);
    }

    AddTask(PHY_BVAL,before_bval);

    if (radiation_flag)
      AddTask(RAD_MOMOPACITY,PHY_BVAL);

    if (CR_ENABLED)
      AddTask(CRTC_OPACITY,PHY_BVAL);

    if (!STS_ENABLED || pm->sts_integrator == "rkl1") {
      AddTask(USERWORK,before_userwork);
      AddTask(NEW_DT,USERWORK);
      if (pm->adaptive) {
        AddTask(FLAG_AMR,USERWORK);
        AddTask(CLEAR_ALLBND,FLAG_AMR);
      } else {
        AddTask(CLEAR_ALLBND,NEW_DT);
      }
    } else {
      AddTask(CLEAR_ALLBND,PHY_BVAL);
    }
  } // end of using namespace block
}

//----------------------------------------------------------------------------------------
//!  Sets id and dependency for "ntask" member of task_list_ array, then iterates value of
//!  ntask.

void TimeIntegratorTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id = id;
  task_list_[ntasks].dependency = dep;
  //! \todo (felker):
  //! - change naming convention of either/both of TASK_NAME and TaskFunc
  //! - There are some issues with the current names:
  //!   VERB_OBJECT is confusing with ObjectVerb(). E.g. seeing SEND_HYD in the task list
  //!   assembly would lead the user to believe the corresponding function is
  //!   SendHydro(), when it is actually HydroSend() ---
  //!   Probaby change function names to active voice
  //!   VerbObject() since "HydroFluxCalculate()" doesn't sound quite right.
  //! - There are exceptions to the "verb+object" convention in some TASK_NAMES and
  //!   TaskFunc, e.g. NEW_DT + NewBlockTimeStep() and AMR_FLAG + CheckRefinement(),
  //!   SRC_TERM and SourceTerms(), USERWORK, PHY_BVAL, PROLONG, CONS2PRIM,
  //!   ... Although, AMR_FLAG = "flag blocks for AMR" should be FLAG_AMR in VERB_OBJECT
  using namespace HydroIntegratorTaskNames; // NOLINT (build/namespace)
  if (id == CLEAR_ALLBND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ClearAllBoundary);
    task_list_[ntasks].lb_time = false;
  } else if (id == CALC_HYDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateHydroFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == CALC_FLDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateEMF);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_HYDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendHydroFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_FLDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendEMF);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_HYDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveAndCorrectHydroFlux);
    task_list_[ntasks].lb_time = false;
  } else if (id == RECV_FLDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveAndCorrectEMF);
    task_list_[ntasks].lb_time = false;
  } else if (id == INT_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::IntegrateHydro);
    task_list_[ntasks].lb_time = true;
  } else if (id == INT_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::IntegrateField);
    task_list_[ntasks].lb_time = true;
  } else if (id == SRC_TERM) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::AddSourceTerms);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendHydro);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendField);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveHydro);
    task_list_[ntasks].lb_time = false;
  } else if (id == RECV_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveField);
    task_list_[ntasks].lb_time = false;
  } else if (id == SETB_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SetBoundariesHydro);
    task_list_[ntasks].lb_time = true;
  } else if (id == SETB_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SetBoundariesField);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_HYDFLXSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendHydroFluxShear);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_HYDFLXSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveHydroFluxShear);
    task_list_[ntasks].lb_time = false;
  } else if (id == SEND_HYDSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendHydroShear);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_HYDSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveHydroShear);
    task_list_[ntasks].lb_time = false;
  } else if (id == SEND_FLDSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendFieldShear);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_FLDSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveFieldShear);
    task_list_[ntasks].lb_time = false;
  } else if (id == SEND_EMFSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendEMFShear);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_EMFSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveEMFShear);
    task_list_[ntasks].lb_time = false;
  } else if (id == PROLONG) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Prolongation);
    task_list_[ntasks].lb_time = true;
  } else if (id == CONS2PRIM) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::Primitives);
    task_list_[ntasks].lb_time = true;
  } else if (id == PHY_BVAL) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::PhysicalBoundary);
    task_list_[ntasks].lb_time = true;
  } else if (id == USERWORK) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::UserWork);
    task_list_[ntasks].lb_time = true;
  } else if (id == NEW_DT) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::NewBlockTimeStep);
    task_list_[ntasks].lb_time = true;
  } else if (id == FLAG_AMR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CheckRefinement);
    task_list_[ntasks].lb_time = true;
  } else if (id == DIFFUSE_HYD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::DiffuseHydro);
    task_list_[ntasks].lb_time = true;
  } else if (id == DIFFUSE_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::DiffuseField);
    task_list_[ntasks].lb_time = true;
  } else if (id == CALC_SCLRFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateScalarFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_SCLRFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendScalarFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_SCLRFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveScalarFlux);
    task_list_[ntasks].lb_time = false;
  } else if (id == INT_SCLR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::IntegrateScalars);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_SCLR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendScalars);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_SCLR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveScalars);
    task_list_[ntasks].lb_time = false;
  } else if (id == SETB_SCLR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SetBoundariesScalars);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_SCLRSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendScalarsShear);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_SCLRSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveScalarsShear);
    task_list_[ntasks].lb_time = false;
  } else if (id == SEND_SCLRFLXSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendScalarsFluxShear);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_SCLRFLXSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveScalarsFluxShear);
    task_list_[ntasks].lb_time = false;
  } else if (id == DIFFUSE_SCLR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::DiffuseScalars);
    task_list_[ntasks].lb_time = true;
  } else if (id == INT_CHM) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::IntegrateChemistry);
  } else if (id == SEND_HYDORB) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendHydroOrbital);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_HYDORB) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveHydroOrbital);
    task_list_[ntasks].lb_time = false;
  } else if (id == CALC_HYDORB) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateHydroOrbital);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_FLDORB) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendFieldOrbital);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_FLDORB) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveFieldOrbital);
    task_list_[ntasks].lb_time = false;
  } else if (id == CALC_FLDORB) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateFieldOrbital);
    task_list_[ntasks].lb_time = true;
  } else if (id == CALC_RADFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateRadFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == INT_RAD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::IntegrateRad);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_RADFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendRadFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_RADFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveAndCorrectRadFlux);
    task_list_[ntasks].lb_time = false;
  } else if (id == SRCTERM_RAD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::AddSourceTermsRad);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_RAD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendRad);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_RAD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveRad);
    task_list_[ntasks].lb_time = true;
  } else if (id == SETB_RAD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SetBoundariesRad);
    task_list_[ntasks].lb_time = true;
  } else if (id == RAD_MOMOPACITY ) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::RadMomOpacity);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_RADFLXSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendRadFluxShear);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_RADFLXSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveRadFluxShear);
    task_list_[ntasks].lb_time = false;
  } else if (id == SEND_RADSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendRadShear);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_RADSH) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveRadShear);
    task_list_[ntasks].lb_time = false;
  } else if (id == CALC_CRTCFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CalculateCRTCFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == INT_CRTC) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::IntegrateCRTC);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_CRTCFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendCRTCFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_CRTCFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveAndCorrectCRTCFlux);
    task_list_[ntasks].lb_time = true;
  } else if (id == SEND_CRTC) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SendCRTC);
    task_list_[ntasks].lb_time = true;
  } else if (id == RECV_CRTC) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::ReceiveCRTC);
    task_list_[ntasks].lb_time = true;
  } else if (id == SRCTERM_CRTC) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::AddSourceTermsCRTC);
    task_list_[ntasks].lb_time = true;
  } else if (id == SETB_CRTC) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::SetBoundariesCRTC);
    task_list_[ntasks].lb_time = true;
  } else if (id == CRTC_OPACITY) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::CRTCOpacity);
    task_list_[ntasks].lb_time = true;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in AddTask" << std::endl
        << "Invalid Task is specified" << std::endl;
    ATHENA_ERROR(msg);
  }
  ntasks++;
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void TimeIntegratorTaskList::StartupTaskList(MeshBlock *pmb, int stage)
//! \brief Initialize time abscissae

void TimeIntegratorTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  bool radiation_flag = (NR_RADIATION_ENABLED && (!IM_RADIATION_ENABLED));
  if (stage == 1) {
    // Initialize storage registers
    Hydro *ph = pmb->phydro;
    ph->u1.ZeroClear();
    if (integrator == "ssprk5_4")
      ph->u2 = ph->u;

    if (MAGNETIC_FIELDS_ENABLED) { // MHD
      Field *pf = pmb->pfield;
      pf->b1.x1f.ZeroClear();
      pf->b1.x2f.ZeroClear();
      pf->b1.x3f.ZeroClear();
      if (integrator == "ssprk5_4") {
        std::stringstream msg;
        msg << "### FATAL ERROR in TimeIntegratorTaskList::StartupTaskList\n"
            << "integrator=" << integrator << " is currently incompatible with MHD."
            << std::endl;
        ATHENA_ERROR(msg);
      }
    }
    if (NSCALARS > 0) {
      PassiveScalars *ps = pmb->pscalars;
      ps->s1.ZeroClear();
      if (integrator == "ssprk5_4")
        ps->s2 = ps->s;
    }

    if (radiation_flag) {
      pmb->pnrrad->ir1.ZeroClear();
      if (integrator == "ssprk5_4")
        pmb->pnrrad->ir2 = pmb->pnrrad->ir;
    }

    if (CR_ENABLED) {
      pmb->pcr->u_cr1.ZeroClear();
      if (integrator == "ssprk5_4")
        pmb->pcr->u_cr2 = pmb->pcr->u_cr;
    }
  }

  if (SHEAR_PERIODIC) {
    Real dt_fc   = pmb->pmy_mesh->dt*(stage_wghts[stage-1].sbeta);
    Real dt_int  = pmb->pmy_mesh->dt*(stage_wghts[stage-1].ebeta);
    Real time = pmb->pmy_mesh->time;
    if (stage==1 ||
        ((stage_wghts[stage-1].sbeta != stage_wghts[stage-2].sbeta)
        || (stage_wghts[stage-1].ebeta != stage_wghts[stage-2].ebeta)))
      pmb->pbval->ComputeShear(time+dt_fc, time+dt_int);
  }

  if (stage_wghts[stage-1].main_stage) {
    pmb->pbval->StartReceivingSubset(BoundaryCommSubset::all, pmb->pbval->bvars_main_int);
  } else {
    pmb->pbval->StartReceivingSubset(BoundaryCommSubset::orbital,
                                     pmb->pbval->bvars_main_int);
  }
  if (stage_wghts[stage-1].orbital_stage && pmb->porb->orbital_advection_active) {
    Real dt = (stage_wghts[(stage-1)].ebeta-stage_wghts[(stage-1)].sbeta)
              *pmb->pmy_mesh->dt;
    pmb->porb->orb_bc->ComputeOrbit(dt);
    pmb->porb->orb_bc->StartReceiving(BoundaryCommSubset::all);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! Functions to end MPI communication

TaskStatus TimeIntegratorTaskList::ClearAllBoundary(MeshBlock *pmb, int stage) {
  if (stage_wghts[stage-1].main_stage) {
    pmb->pbval->ClearBoundarySubset(BoundaryCommSubset::all,
                                    pmb->pbval->bvars_main_int);
  } else {
    pmb->pbval->ClearBoundarySubset(BoundaryCommSubset::orbital,
                                    pmb->pbval->bvars_main_int);
  }
  if (stage_wghts[stage-1].orbital_stage && pmb->porb->orbital_advection_active) {
    pmb->porb->orb_bc->ClearBoundary(BoundaryCommSubset::all);
  }

  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
// Functions to calculates Hydro fluxes

TaskStatus TimeIntegratorTaskList::CalculateHydroFlux(MeshBlock *pmb, int stage) {
  Hydro *phydro = pmb->phydro;
  Field *pfield = pmb->pfield;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      if ((integrator == "vl2") && (stage-stage_wghts[0].orbital_stage == 1)) {
        phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, 1);
      } else {
        phydro->CalculateFluxes(phydro->w,  pfield->b,  pfield->bcc, pmb->precon->xorder);
      }
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
// Functions to calculates EMFs

TaskStatus TimeIntegratorTaskList::CalculateEMF(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      pmb->pfield->ComputeCornerE(pmb->phydro->w,  pmb->pfield->bcc);
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
// Functions to communicate fluxes between MeshBlocks for flux correction with AMR

TaskStatus TimeIntegratorTaskList::SendHydroFlux(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      pmb->phydro->hbvar.SendFluxCorrection();
    }
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
// Functions to communicate emf between MeshBlocks for flux correction with AMR

TaskStatus TimeIntegratorTaskList::SendEMF(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      pmb->pfield->fbvar.SendFluxCorrection();
    }
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
// Functions to receive fluxes between MeshBlocks

TaskStatus TimeIntegratorTaskList::ReceiveAndCorrectHydroFlux(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      if (pmb->phydro->hbvar.ReceiveFluxCorrection()) {
        return TaskStatus::next;
      } else {
        return TaskStatus::fail;
      }
    } else {
      return TaskStatus::next;
    }
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
// Functions to receive emf between MeshBlocks

TaskStatus TimeIntegratorTaskList::ReceiveAndCorrectEMF(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      if (pmb->pfield->fbvar.ReceiveFluxCorrection()) {
        return TaskStatus::next;
      } else {
        return TaskStatus::fail;
      }
    } else {
      return TaskStatus::next;
    }
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
// Functions to integrate Hydro variables

TaskStatus TimeIntegratorTaskList::IntegrateHydro(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Field *pf = pmb->pfield;

  if (pmb->pmy_mesh->fluid_setup != FluidFormulation::evolve) return TaskStatus::next;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      // This time-integrator-specific averaging operation logic is identical to FieldInt
      Real ave_wghts[5];
      ave_wghts[0] = 1.0;
      ave_wghts[1] = stage_wghts[stage-1].delta;
      ave_wghts[2] = 0.0;
      ave_wghts[3] = 0.0;
      ave_wghts[4] = 0.0;
      pmb->WeightedAve(ph->u1, ph->u, ph->u2, ph->u0, ph->fl_div, ave_wghts);

      ave_wghts[0] = stage_wghts[stage-1].gamma_1;
      ave_wghts[1] = stage_wghts[stage-1].gamma_2;
      ave_wghts[2] = stage_wghts[stage-1].gamma_3;
      if (ave_wghts[0] == 0.0 && ave_wghts[1] == 1.0 && ave_wghts[2] == 0.0) {
        ph->u.SwapAthenaArray(ph->u1);
      } else {
        pmb->WeightedAve(ph->u, ph->u1, ph->u2, ph->u0, ph->fl_div, ave_wghts);
      }

      const Real wght = stage_wghts[stage-1].beta*pmb->pmy_mesh->dt;
      ph->AddFluxDivergence(wght, ph->u);
      // add coordinate (geometric) source terms
      pmb->pcoord->AddCoordTermsDivergence(wght, ph->flux, ph->w, pf->bcc, ph->u);

      // Hardcode an additional flux divergence weighted average for the penultimate
      // stage of SSPRK(5,4) since it cannot be expressed in a 3S* framework
      if (stage == 4 && integrator == "ssprk5_4") {
        // From Gottlieb (2009), u^(n+1) partial calculation
        ave_wghts[0] = -1.0; // -u^(n) coeff.
        ave_wghts[1] = 0.0;
        ave_wghts[2] = 0.0;
        const Real beta = 0.063692468666290; // F(u^(3)) coeff.
        const Real wght_ssp = beta*pmb->pmy_mesh->dt;
        // writing out to u2 register
        pmb->WeightedAve(ph->u2, ph->u1, ph->u2, ph->u0, ph->fl_div, ave_wghts);
        ph->AddFluxDivergence(wght_ssp, ph->u2);
        // add coordinate (geometric) source terms
        pmb->pcoord->AddCoordTermsDivergence(wght_ssp, ph->flux, ph->w, pf->bcc, ph->u2);
      }
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
// Functions to integrate Field variables


TaskStatus TimeIntegratorTaskList::IntegrateField(MeshBlock *pmb, int stage) {
  Field *pf = pmb->pfield;

  if (pmb->pmy_mesh->fluid_setup != FluidFormulation::evolve) return TaskStatus::next;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      // This time-integrator-specific averaging operation logic is identical to HydroInt
      Real ave_wghts[5];
      ave_wghts[0] = 1.0;
      ave_wghts[1] = stage_wghts[stage-1].delta;
      ave_wghts[2] = 0.0;
      ave_wghts[3] = 0.0;
      ave_wghts[4] = 0.0;
      pmb->WeightedAve(pf->b1, pf->b, pf->b2, pf->b0, pf->ct_update, ave_wghts);

      ave_wghts[0] = stage_wghts[stage-1].gamma_1;
      ave_wghts[1] = stage_wghts[stage-1].gamma_2;
      ave_wghts[2] = stage_wghts[stage-1].gamma_3;
      if (ave_wghts[0] == 0.0 && ave_wghts[1] == 1.0 && ave_wghts[2] == 0.0) {
        pf->b.x1f.SwapAthenaArray(pf->b1.x1f);
        pf->b.x2f.SwapAthenaArray(pf->b1.x2f);
        pf->b.x3f.SwapAthenaArray(pf->b1.x3f);
      } else {
        pmb->WeightedAve(pf->b, pf->b1, pf->b2, pf->b0, pf->ct_update, ave_wghts);
      }

      pf->CT(stage_wghts[stage-1].beta*pmb->pmy_mesh->dt, pf->b);
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
//! Functions to add source terms

TaskStatus TimeIntegratorTaskList::AddSourceTerms(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Field *pf = pmb->pfield;
  PassiveScalars *ps = pmb->pscalars;

  // return if there are no source terms to be added
  if (!(ph->hsrc.hydro_sourceterms_defined)
      || pmb->pmy_mesh->fluid_setup != FluidFormulation::evolve) return TaskStatus::next;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      // Time at beginning of stage for u()
      Real t_start_stage = pmb->pmy_mesh->time
                           + stage_wghts[(stage-1)].sbeta*pmb->pmy_mesh->dt;
      // Scaled coefficient for RHS update
      Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
      // Evaluate the source terms at the time at the beginning of the stage
      ph->hsrc.AddSourceTerms(t_start_stage, dt, ph->flux, ph->w, ps->r, pf->bcc,
                                   ph->u, ps->s);
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
//! Functions to calculate hydro diffusion fluxes (stored in HydroDiffusion::visflx[],
//! cndflx[], added at the end of Hydro::CalculateFluxes()

TaskStatus TimeIntegratorTaskList::DiffuseHydro(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Field *pf = pmb->pfield;

  // return if there are no diffusion to be added
  if (!(ph->hdif.hydro_diffusion_defined)
      || pmb->pmy_mesh->fluid_setup != FluidFormulation::evolve) return TaskStatus::next;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      // if using orbital advection, put modified conservative into the function
      if (pmb->porb->orbital_advection_defined) {
        pmb->porb->ConvertOrbitalSystem(ph->w, ph->u, OrbitalTransform::prim);
        ph->hdif.CalcDiffusionFlux(ph->w, pmb->porb->w_orb, pf->bcc);
      } else {
        ph->hdif.CalcDiffusionFlux(ph->w, ph->w, pf->bcc);
      }
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
//! Functions to calculate diffusion EMF

TaskStatus TimeIntegratorTaskList::DiffuseField(MeshBlock *pmb, int stage) {
  Field *pf = pmb->pfield;

  // return if there are no diffusion to be added
  if (!(pf->fdif.field_diffusion_defined)) return TaskStatus::next;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      // TODO(pdmullen): DiffuseField is also called in SuperTimeStepTaskLsit.
      // It must skip Hall effect (once implemented) diffusion process in STS
      // and always calculate those terms in the main integrator.
      pf->fdif.CalcDiffusionEMF(pf->b, pf->bcc, pf->e);
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
//! Functions to communicate Hydro variables between MeshBlocks

TaskStatus TimeIntegratorTaskList::SendHydro(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    // Swap Hydro quantity in BoundaryVariable interface back to conserved var formulation
    // (also needed in SetBoundariesHydro(), since the tasks are independent)
    pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->u, HydroBoundaryQuantity::cons);
    pmb->phydro->hbvar.SendBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
//! Functions to communicate Field variables between MeshBlocks

TaskStatus TimeIntegratorTaskList::SendField(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pfield->fbvar.SendBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
//! Functions to receive Hydro variables between MeshBlocks

TaskStatus TimeIntegratorTaskList::ReceiveHydro(MeshBlock *pmb, int stage) {
  bool ret;
  if (stage <= nstages) {
    ret = pmb->phydro->hbvar.ReceiveBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}

//----------------------------------------------------------------------------------------
//! Functions to receive Field variables between MeshBlocks

TaskStatus TimeIntegratorTaskList::ReceiveField(MeshBlock *pmb, int stage) {
  bool ret;
  if (stage <= nstages) {
    ret = pmb->pfield->fbvar.ReceiveBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}

//----------------------------------------------------------------------------------------
//! Functions to set Hydro boundaries

TaskStatus TimeIntegratorTaskList::SetBoundariesHydro(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->u, HydroBoundaryQuantity::cons);
    pmb->phydro->hbvar.SetBoundaries();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
//! Functions to set Field boundaries

TaskStatus TimeIntegratorTaskList::SetBoundariesField(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pfield->fbvar.SetBoundaries();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
//! Functions to communicate Hydro variables between MeshBlocks with shear

TaskStatus TimeIntegratorTaskList::SendHydroShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->phydro->hbvar.SendShearingBoxBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
//! Functions to communicate Hydro variables between MeshBlocks with shear

TaskStatus TimeIntegratorTaskList::ReceiveHydroShear(MeshBlock *pmb, int stage) {
  bool ret;
  ret = false;
  if (stage <= nstages) {
    ret = pmb->phydro->hbvar.ReceiveShearingBoxBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    pmb->phydro->hbvar.SetShearingBoxBoundaryBuffers();
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}

//----------------------------------------------------------------------------------------
//! Functions to communicate Field variables between MeshBlocks with shear

TaskStatus TimeIntegratorTaskList::SendHydroFluxShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      pmb->phydro->hbvar.SendFluxShearingBoxBoundaryBuffers();
    }
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::ReceiveHydroFluxShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      if (pmb->phydro->hbvar.ReceiveFluxShearingBoxBoundaryBuffers()) {
        pmb->phydro->hbvar.SetFluxShearingBoxBoundaryBuffers();
        return TaskStatus::success;
      } else {
        return TaskStatus::fail;
      }
    } else {
      return TaskStatus::success;
    }
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendFieldShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pfield->fbvar.SendShearingBoxBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
//! Functions to communicate Field variables between MeshBlocks with shear

TaskStatus TimeIntegratorTaskList::ReceiveFieldShear(MeshBlock *pmb, int stage) {
  bool ret;
  ret = false;
  if (stage <= nstages) {
    ret = pmb->pfield->fbvar.ReceiveShearingBoxBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    pmb->pfield->fbvar.SetShearingBoxBoundaryBuffers();
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}

//----------------------------------------------------------------------------------------
//! Functions to communicate EMFs between MeshBlocks with shear

TaskStatus TimeIntegratorTaskList::SendEMFShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      pmb->pfield->fbvar.SendEMFShearingBoxBoundaryCorrection();
    }
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
//! Functions to communicate EMFs between MeshBlocks with shear

TaskStatus TimeIntegratorTaskList::ReceiveEMFShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      if (pmb->pfield->fbvar.ReceiveEMFShearingBoxBoundaryCorrection()) {
        pmb->pfield->fbvar.SetEMFShearingBoxBoundaryCorrection();
        return TaskStatus::success;
      } else {
        return TaskStatus::fail;
      }
    } else {
      return TaskStatus::success;
    }
  }
  return TaskStatus::fail;
}

//--------------------------------------------------------------------------------------
// Functions for everything else

TaskStatus TimeIntegratorTaskList::Prolongation(MeshBlock *pmb, int stage) {
  BoundaryValues *pbval = pmb->pbval;

  if (stage <= nstages) {
    // Time at the end of stage for (u, b) register pair
    Real t_end_stage = pmb->pmy_mesh->time
                       + stage_wghts[(stage-1)].ebeta*pmb->pmy_mesh->dt;
    // Scaled coefficient for RHS time-advance within stage
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
    pbval->ProlongateBoundaries(t_end_stage, dt, pmb->pbval->bvars_main_int);
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::Primitives(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Field *pf = pmb->pfield;
  PassiveScalars *ps = pmb->pscalars;
  BoundaryValues *pbval = pmb->pbval;

  int il = pmb->is, iu = pmb->ie, jl = pmb->js, ju = pmb->je, kl = pmb->ks, ku = pmb->ke;
  if (pbval->nblevel[1][1][0] != -1) il -= NGHOST;
  if (pbval->nblevel[1][1][2] != -1) iu += NGHOST;
  if (pbval->nblevel[1][0][1] != -1) jl -= NGHOST;
  if (pbval->nblevel[1][2][1] != -1) ju += NGHOST;
  if (pbval->nblevel[0][1][1] != -1) kl -= NGHOST;
  if (pbval->nblevel[2][1][1] != -1) ku += NGHOST;

  if (stage <= nstages) {
    // At beginning of this task, ph->w contains previous stage's W(U) output
    // and ph->w1 is used as a register to store the current stage's output.
    // For the second order integrators VL2 and RK2, the prim_old initial guess for the
    // Newton-Raphson solver in GR EOS uses the following abscissae:
    // stage=1: W at t^n and
    // stage=2: W at t^{n+1/2} (VL2) or t^{n+1} (RK2)
    pmb->peos->ConservedToPrimitive(ph->u, ph->w, pf->b,
                                    ph->w1, pf->bcc, pmb->pcoord,
                                    il, iu, jl, ju, kl, ku);
    if (pmb->porb->orbital_advection_defined) {
      pmb->porb->ResetOrbitalSystemConversionFlag();
    }
    if (NSCALARS > 0) {
      // r1/r_old for GR is currently unused:
      pmb->peos->PassiveScalarConservedToPrimitive(ps->s, ph->u, ps->r, ps->r,
                                                   pmb->pcoord, il, iu, jl, ju, kl, ku);
    }
    // fourth-order EOS:
    if (pmb->precon->xorder == 4) {
      // for hydro, shrink buffer by 1 on all sides
      if (pbval->nblevel[1][1][0] != -1) il += 1;
      if (pbval->nblevel[1][1][2] != -1) iu -= 1;
      if (pbval->nblevel[1][0][1] != -1) jl += 1;
      if (pbval->nblevel[1][2][1] != -1) ju -= 1;
      if (pbval->nblevel[0][1][1] != -1) kl += 1;
      if (pbval->nblevel[2][1][1] != -1) ku -= 1;
      // for MHD, shrink buffer by 3
      // TODO(felker): add MHD loop limit calculation for 4th order W(U)
      // Apply physical boundaries prior to 4th order W(U)
      // Time at the end of stage for (u, b) register pair
      Real t_end_stage = pmb->pmy_mesh->time
                         + stage_wghts[(stage-1)].ebeta*pmb->pmy_mesh->dt;
      // Scaled coefficient for RHS time-advance within stage
      Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
      // Swap Hydro and (possibly) passive scalar quantities in BoundaryVariable interface
      // from conserved to primitive formulations:
      ph->hbvar.SwapHydroQuantity(ph->w1, HydroBoundaryQuantity::prim);
      if (NSCALARS > 0) {
        ps->sbvar.var_cc = &(ps->r);
        if (pmb->pmy_mesh->multilevel) {
          ps->sbvar.coarse_buf = &(ps->coarse_r_);
        }
      }
      pbval->ApplyPhysicalBoundaries(t_end_stage, dt, pmb->pbval->bvars_main_int);
      // Perform 4th order W(U)
      pmb->peos->ConservedToPrimitiveCellAverage(ph->u, ph->w, pf->b,
                                                 ph->w1, pf->bcc, pmb->pcoord,
                                                 il, iu, jl, ju, kl, ku);
      if (NSCALARS > 0) {
        pmb->peos->PassiveScalarConservedToPrimitiveCellAverage(
            ps->s, ps->r, ps->r, pmb->pcoord, il, iu, jl, ju, kl, ku);
      }
    }
    // swap AthenaArray data pointers so that w now contains the updated w_out
    ph->w.SwapAthenaArray(ph->w1);
    // r1/r_old for GR is currently unused:
    // ps->r.SwapAthenaArray(ps->r1);
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::PhysicalBoundary(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  PassiveScalars *ps = pmb->pscalars;
  BoundaryValues *pbval = pmb->pbval;

  if (stage <= nstages) {
    // Time at the end of stage for (u, b) register pair
    Real t_end_stage = pmb->pmy_mesh->time
                       + stage_wghts[(stage-1)].ebeta*pmb->pmy_mesh->dt;
    // Scaled coefficient for RHS time-advance within stage
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
    // Swap Hydro and (possibly) passive scalar quantities in BoundaryVariable interface
    // from conserved to primitive formulations:
    ph->hbvar.SwapHydroQuantity(ph->w, HydroBoundaryQuantity::prim);
    if (NSCALARS > 0) {
      ps->sbvar.var_cc = &(ps->r);
      if (pmb->pmy_mesh->multilevel) {
        ps->sbvar.coarse_buf = &(ps->coarse_r_);
      }
    }
    pbval->ApplyPhysicalBoundaries(t_end_stage, dt, pmb->pbval->bvars_main_int);
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::UserWork(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  pmb->UserWorkInLoop();
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::NewBlockTimeStep(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  pmb->phydro->NewBlockTimeStep();
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::CheckRefinement(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  pmb->pmr->CheckRefinementCondition();
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::CalculateScalarFlux(MeshBlock *pmb, int stage) {
  PassiveScalars *ps = pmb->pscalars;
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      if ((integrator == "vl2") && (stage-stage_wghts[0].orbital_stage == 1)) {
        ps->CalculateFluxes(ps->r, 1);
      } else {
        ps->CalculateFluxes(ps->r, pmb->precon->xorder);
      }
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendScalarFlux(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      pmb->pscalars->sbvar.SendFluxCorrection();
    }
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::ReceiveScalarFlux(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      if (pmb->pscalars->sbvar.ReceiveFluxCorrection()) {
        return TaskStatus::next;
      } else {
        return TaskStatus::fail;
      }
    } else {
      return TaskStatus::next;
    }
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::IntegrateScalars(MeshBlock *pmb, int stage) {
  if (pmb->pmy_mesh->fluid_setup == FluidFormulation::fixed) return TaskStatus::next;

  PassiveScalars *ps = pmb->pscalars;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      // This time-integrator-specific averaging operation logic is identical to
      // IntegrateHydro, IntegrateField
      Real ave_wghts[5];
      ave_wghts[0] = 1.0;
      ave_wghts[1] = stage_wghts[stage-1].delta;
      ave_wghts[2] = 0.0;
      ave_wghts[3] = 0.0;
      ave_wghts[4] = 0.0;
      pmb->WeightedAve(ps->s1, ps->s, ps->s2, ps->s0, ps->s_fl_div, ave_wghts);

      ave_wghts[0] = stage_wghts[stage-1].gamma_1;
      ave_wghts[1] = stage_wghts[stage-1].gamma_2;
      ave_wghts[2] = stage_wghts[stage-1].gamma_3;
      if (ave_wghts[0] == 0.0 && ave_wghts[1] == 1.0 && ave_wghts[2] == 0.0) {
        ps->s.SwapAthenaArray(ps->s1);
      } else {
        pmb->WeightedAve(ps->s, ps->s1, ps->s2, ps->s0, ps->s_fl_div, ave_wghts);
      }

      const Real wght = stage_wghts[stage-1].beta*pmb->pmy_mesh->dt;
      ps->AddFluxDivergence(wght, ps->s);

      // Hardcode an additional flux divergence weighted average for the penultimate
      // stage of SSPRK(5,4) since it cannot be expressed in a 3S* framework
      if (stage == 4 && integrator == "ssprk5_4") {
        // From Gottlieb (2009), u^(n+1) partial calculation
        ave_wghts[0] = -1.0; // -u^(n) coeff.
        ave_wghts[1] = 0.0;
        ave_wghts[2] = 0.0;
        const Real beta = 0.063692468666290; // F(u^(3)) coeff.
        const Real wght_ssp = beta*pmb->pmy_mesh->dt;
        // writing out to s2 register
        pmb->WeightedAve(ps->s2, ps->s1, ps->s2, ps->s0, ps->s_fl_div, ave_wghts);
        ps->AddFluxDivergence(wght_ssp, ps->s2);
      }
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

TaskStatus TimeIntegratorTaskList::IntegrateChemistry(MeshBlock *pmb, int stage) {
  // integrate chemistry reactions
  if (CHEMISTRY_ENABLED) {
    if (stage != nstages) return TaskStatus::success; // only do on last stage

    pmb->pscalars->odew.Integrate(pmb->pmy_mesh->time, pmb->pmy_mesh->dt);
  }
  return TaskStatus::next;
}

TaskStatus TimeIntegratorTaskList::SendScalars(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    // Swap PassiveScalars quantity in BoundaryVariable interface back to conserved var
    // formulation (also needed in SetBoundariesScalars() since the tasks are independent)
    pmb->pscalars->sbvar.var_cc = &(pmb->pscalars->s);
    if (pmb->pmy_mesh->multilevel) {
      pmb->pscalars->sbvar.coarse_buf = &(pmb->pscalars->coarse_s_);
    }
    pmb->pscalars->sbvar.SendBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::ReceiveScalars(MeshBlock *pmb, int stage) {
  bool ret;
  if (stage <= nstages) {
    ret = pmb->pscalars->sbvar.ReceiveBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::SetBoundariesScalars(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    // Set PassiveScalars quantity in BoundaryVariable interface to cons var formulation
    pmb->pscalars->sbvar.var_cc = &(pmb->pscalars->s);
    if (pmb->pmy_mesh->multilevel) {
      pmb->pscalars->sbvar.coarse_buf = &(pmb->pscalars->coarse_s_);
    }
    pmb->pscalars->sbvar.SetBoundaries();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::DiffuseScalars(MeshBlock *pmb, int stage) {
  if (pmb->pmy_mesh->fluid_setup == FluidFormulation::fixed) return TaskStatus::next;

  PassiveScalars *ps = pmb->pscalars;
  Hydro *ph = pmb->phydro;
  // return if there are no diffusion to be added
  if (!(ps->scalar_diffusion_defined))
    return TaskStatus::next;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      // TODO(felker): adapted directly from HydroDiffusion::ClearFlux. Deduplicate
      ps->diffusion_flx[X1DIR].ZeroClear();
      ps->diffusion_flx[X2DIR].ZeroClear();
      ps->diffusion_flx[X3DIR].ZeroClear();

      // unlike HydroDiffusion, only 1x passive scalar diffusive process is allowed, so
      // there is no need for counterpart to wrapper fn HydroDiffusion::CalcDiffusionFlux
      ps->DiffusiveFluxIso(ps->r, ph->w, ps->diffusion_flx);
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendScalarsShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pscalars->sbvar.SendShearingBoxBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::ReceiveScalarsShear(MeshBlock *pmb, int stage) {
  bool ret;
  ret = false;
  if (stage <= nstages) {
    ret = pmb->pscalars->sbvar.ReceiveShearingBoxBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    pmb->pscalars->sbvar.SetShearingBoxBoundaryBuffers();
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}


TaskStatus TimeIntegratorTaskList::SendScalarsFluxShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      pmb->pscalars->sbvar.SendFluxShearingBoxBoundaryBuffers();
    }
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::ReceiveScalarsFluxShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_before ||
        pmb->pmy_mesh->sts_loc == TaskType::op_split_after) {
      if (pmb->pscalars->sbvar.ReceiveFluxShearingBoxBoundaryBuffers()) {
        pmb->pscalars->sbvar.SetFluxShearingBoxBoundaryBuffers();
        return TaskStatus::success;
      } else {
        return TaskStatus::fail;
      }
    } else {
      return TaskStatus::success;
    }
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendHydroOrbital(MeshBlock *pmb, int stage) {
  if (!stage_wghts[stage-1].orbital_stage) {
    return TaskStatus::success;
  } else {
    OrbitalAdvection *porb = pmb->porb;
    if (!porb->orbital_advection_active) return TaskStatus::success;
    Hydro *ph = pmb->phydro;
    PassiveScalars *ps = pmb->pscalars;
    porb->SetOrbitalAdvectionCC(ph->u, ps->s);
    porb->orb_bc->SendBoundaryBuffersCC();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendFieldOrbital(MeshBlock *pmb, int stage) {
  if (!stage_wghts[stage-1].orbital_stage) {
    return TaskStatus::success;
  } else {
    OrbitalAdvection *porb = pmb->porb;
    if (!porb->orbital_advection_active) return TaskStatus::success;
    Field *pf = pmb->pfield;
    porb->SetOrbitalAdvectionFC(pf->b);
    porb->orb_bc->SendBoundaryBuffersFC();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::ReceiveHydroOrbital(MeshBlock *pmb, int stage) {
  if (!stage_wghts[stage-1].orbital_stage) {
    return TaskStatus::success;
  } else {
    OrbitalAdvection *porb = pmb->porb;
    if (!porb->orbital_advection_active) return TaskStatus::success;
    if (porb->orb_bc->ReceiveBoundaryBuffersCC()) {
      return TaskStatus::success;
    }
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::ReceiveFieldOrbital(MeshBlock *pmb, int stage) {
  if (!stage_wghts[stage-1].orbital_stage) {
    return TaskStatus::success;
  } else {
    OrbitalAdvection *porb = pmb->porb;
    if (!porb->orbital_advection_active) return TaskStatus::success;
    if (porb->orb_bc->ReceiveBoundaryBuffersFC()) {
      return TaskStatus::success;
    }
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::CalculateHydroOrbital(MeshBlock *pmb, int stage) {
  if (!stage_wghts[stage-1].orbital_stage) {
    return TaskStatus::success;
  } else {
    OrbitalAdvection *porb = pmb->porb;
    Hydro *ph = pmb->phydro;
    PassiveScalars *ps = pmb->pscalars;
    Real dt = pmb->pmy_mesh->dt
              *(stage_wghts[(stage-1)].ebeta-stage_wghts[(stage-1)].sbeta);
    porb->CalculateOrbitalAdvectionCC(dt, ph->u, ps->s);
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::CalculateFieldOrbital(MeshBlock *pmb, int stage) {
  if (!stage_wghts[stage-1].orbital_stage) {
    return TaskStatus::success;
  } else {
    OrbitalAdvection *porb = pmb->porb;
    if (!porb->orbital_advection_active) return TaskStatus::success;
    Field *pf = pmb->pfield;
    Real dt = pmb->pmy_mesh->dt
              *(stage_wghts[(stage-1)].ebeta-stage_wghts[(stage-1)].sbeta);
    porb->CalculateOrbitalAdvectionFC(dt, pf->e);
    pf->CT(1.0, pf->b);
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


//----------------------------------------------------------------------------------------
// New steps necessary for radiation/cosmic ray modules

TaskStatus TimeIntegratorTaskList::IntegrateRad(MeshBlock *pmb, int stage) {
  NRRadiation *prad = pmb->pnrrad;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      Real ave_wghts[3];
      ave_wghts[0] = 1.0;
      ave_wghts[1] = stage_wghts[stage-1].delta;
      ave_wghts[2] = 0.0;
      pmb->WeightedAve(prad->ir1, prad->ir, prad->ir2, ave_wghts,1);

      ave_wghts[0] = stage_wghts[stage-1].gamma_1;
      ave_wghts[1] = stage_wghts[stage-1].gamma_2;
      ave_wghts[2] = stage_wghts[stage-1].gamma_3;
      if (ave_wghts[0] == 0.0 && ave_wghts[1] == 1.0 && ave_wghts[2] == 0.0) {
        prad->ir.SwapAthenaArray(prad->ir1);
      } else {
        pmb->WeightedAve(prad->ir, prad->ir1, prad->ir2, ave_wghts,1);
      }
      const Real wght = stage_wghts[stage-1].beta*pmb->pmy_mesh->dt;
      // ir is already partially updated
      prad->pradintegrator->FluxDivergence(wght, prad->ir, prad->ir);

      // no geometric source term for radiation

      // Hardcode an additional flux divergence weighted average for the penultimate
      // stage of SSPRK(5,4) since it cannot be expressed in a 3S* framework
      if (stage == 4 && integrator == "ssprk5_4") {
        // From Gottlieb (2009), u^(n+1) partial calculation
        ave_wghts[0] = -1.0; // -u^(n) coeff.
        ave_wghts[1] = 0.0;
        ave_wghts[2] = 0.0;
        const Real beta = 0.063692468666290; // F(u^(3)) coeff.
        const Real wght_ssp = beta*pmb->pmy_mesh->dt;
        // writing out to u2 register
        pmb->WeightedAve(prad->ir2, prad->ir1, prad->ir2, ave_wghts,1);
        prad->pradintegrator->FluxDivergence(wght_ssp, prad->ir2, prad->ir2);
      }
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

// tasks for radiation transport
TaskStatus TimeIntegratorTaskList::CalculateRadFlux(MeshBlock *pmb, int stage) {
  Hydro *phydro = pmb->phydro;
  NRRadiation *prad = pmb->pnrrad;
  if (stage_wghts[stage-1].main_stage) {
    Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);
    prad->pradintegrator->GetTgasVel(pmb,dt,phydro->u,phydro->w,
                                   pmb->pfield->bcc,prad->ir);
  }

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      if ((stage == 1) && (integrator == "vl2")) {
        prad->pradintegrator->CalculateFluxes(phydro->w,  prad->ir, 1);
      } else {
        prad->pradintegrator->CalculateFluxes(phydro->w,  prad->ir,
                                     prad->pradintegrator->rad_xorder);
      }
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}


//----------------------------------------------------------------------------------------
// Functions to add source terms

TaskStatus TimeIntegratorTaskList::AddSourceTermsRad(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  NRRadiation *prad = pmb->pnrrad;

  int is=pmb->is, ie=pmb->ie;
  int js=pmb->js, je=pmb->je;
  int ks=pmb->ks, ke=pmb->ke;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      // Real t_start_stage = pmb->pmy_mesh->time
      //                      + stage_wghts[(stage-1)].sbeta*pmb->pmy_mesh->dt;

      // Scaled coefficient for RHS update
      Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);

      // Evaluate the time-dependent source terms
      // Both u and ir are partially updated, only w is from the beginning of the step
      prad->ir_old = prad->ir;

      for(int k=ks; k<=ke; ++k)
        for(int j=js; j<=je; ++j)
          for(int i=is; i<=ie; ++i) {
            prad->pradintegrator->CalSourceTerms(pmb, dt, k, j, i, ph->u,
                                                       prad->ir, prad->ir);
          }

      if (prad->set_source_flag > 0) {
        prad->pradintegrator->GetHydroSourceTerms(pmb, prad->ir_old, prad->ir);
        prad->pradintegrator->AddSourceTerms(pmb, ph->u);
      }
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}


//----------------------------------------------------------------------------------------
// Functions to communicate fluxes between MeshBlocks for flux correction with AMR

TaskStatus TimeIntegratorTaskList::SendRadFlux(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      pmb->pnrrad->rad_bvar.SendFluxCorrection();
    }
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::ReceiveAndCorrectRadFlux(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      if (pmb->pnrrad->rad_bvar.ReceiveFluxCorrection()) {
        return TaskStatus::next;
      } else {
        return TaskStatus::fail;
      }
    } else {
      return TaskStatus::next;
    }
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendRad(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pnrrad->rad_bvar.SendBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::ReceiveRad(MeshBlock *pmb, int stage) {
  bool ret;
  if (stage <= nstages) {
    ret = pmb->pnrrad->rad_bvar.ReceiveBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}

TaskStatus TimeIntegratorTaskList::SetBoundariesRad(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pnrrad->rad_bvar.SetBoundaries();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}

TaskStatus TimeIntegratorTaskList::RadMomOpacity(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    // pmb->pnrrad->CalculateMoment(pmb->pnrrad->ir);
    pmb->pnrrad->UpdateOpacity(pmb,pmb->phydro->w);
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


//----------------------------------------------------------------------------------------
//! Functions to communicate Field variables between MeshBlocks with shear

TaskStatus TimeIntegratorTaskList::SendRadFluxShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      pmb->pnrrad->rad_bvar.SendFluxShearingBoxBoundaryBuffers();
    }
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::ReceiveRadFluxShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      if (pmb->pnrrad->rad_bvar.ReceiveFluxShearingBoxBoundaryBuffers()) {
        pmb->pnrrad->rad_bvar.SetFluxShearingBoxBoundaryBuffers();
        return TaskStatus::success;
      } else {
        return TaskStatus::fail;
      }
    } else {
      return TaskStatus::success;
    }
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendRadShear(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    pmb->pnrrad->rad_bvar.SendShearingBoxBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}


TaskStatus TimeIntegratorTaskList::ReceiveRadShear(MeshBlock *pmb, int stage) {
  bool ret;
  ret = false;
  if (stage <= nstages) {
    ret = pmb->pnrrad->rad_bvar.ReceiveShearingBoxBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret) {
    pmb->pnrrad->rad_bvar.SetShearingBoxBoundaryBuffers();
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}

//----------------------------------------------------------------------------------
//function for cosmic ray transport

TaskStatus TimeIntegratorTaskList::CalculateCRTCFlux(MeshBlock *pmb, int stage) {
  Hydro *phydro = pmb->phydro;
  CosmicRay *pcr = pmb->pcr;
  Field *pf = pmb->pfield;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      if ((stage == 1) && (integrator == "vl2")) {
        if (CR_ENABLED)
          pcr->pcrintegrator->CalculateFluxes(phydro->w, pf->bcc, pcr->u_cr, 1);
      } else {
        if (CR_ENABLED)
          pcr->pcrintegrator->CalculateFluxes(phydro->w, pf->bcc, pcr->u_cr,
                                           pcr->pcrintegrator->cr_xorder);
      }
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}



//----------------------------------------------------------------------------------------
// Functions to integrate conserved variables

TaskStatus TimeIntegratorTaskList::IntegrateCRTC(MeshBlock *pmb, int stage) {
  CosmicRay *pcr = pmb->pcr;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      Real ave_wghts[5];
      ave_wghts[0] = 1.0;
      ave_wghts[1] = stage_wghts[stage-1].delta;
      ave_wghts[2] = 0.0;
      ave_wghts[3] = 0.0;
      ave_wghts[4] = 0.0;
      if (CR_ENABLED) {
        pmb->WeightedAve(pcr->u_cr1, pcr->u_cr, pcr->u_cr2, pcr->u_cr2, pcr->u_cr2,
                         ave_wghts);
      }
      ave_wghts[0] = stage_wghts[stage-1].gamma_1;
      ave_wghts[1] = stage_wghts[stage-1].gamma_2;
      ave_wghts[2] = stage_wghts[stage-1].gamma_3;
      if (ave_wghts[0] == 0.0 && ave_wghts[1] == 1.0 && ave_wghts[2] == 0.0) {
        if (CR_ENABLED)
          pcr->u_cr1.SwapAthenaArray(pcr->u_cr);
      } else {
        // ave_wghts[3] and ave_wght[4] = 0
        if (CR_ENABLED)
          pmb->WeightedAve(pcr->u_cr, pcr->u_cr1, pcr->u_cr2, pcr->u_cr2, pcr->u_cr2,
                           ave_wghts);
      }
      const Real wght = stage_wghts[stage-1].beta*pmb->pmy_mesh->dt;
      if (CR_ENABLED) {
        pcr->pcrintegrator->FluxDivergence(wght, pcr->u_cr);
      }
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

TaskStatus TimeIntegratorTaskList::SendCRTCFlux(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      if (CR_ENABLED)
        pmb->pcr->cr_bvar.SendFluxCorrection();
    }
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
// Functions to receive fluxes between MeshBlocks

TaskStatus TimeIntegratorTaskList::ReceiveAndCorrectCRTCFlux(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      bool flag_cr = true;
      if (CR_ENABLED) {
        flag_cr = pmb->pcr->cr_bvar.ReceiveFluxCorrection();
      }

      if (flag_cr) {
        return TaskStatus::next;
      } else {
        return TaskStatus::fail;
      }
    } else {
      return TaskStatus::next;
    }
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::SendCRTC(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    // Swap Hydro quantity in BoundaryVariable interface back to conserved var formulation
    // (also needed in SetBoundariesHydro(), since the tasks are independent)
    if (CR_ENABLED)
      pmb->pcr->cr_bvar.SendBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
// Functions to receive conserved variables between MeshBlocks

TaskStatus TimeIntegratorTaskList::ReceiveCRTC(MeshBlock *pmb, int stage) {
  bool ret_cr = true;
  if (stage <= nstages) {
    if (CR_ENABLED)
      ret_cr = pmb->pcr->cr_bvar.ReceiveBoundaryBuffers();
  } else {
    return TaskStatus::fail;
  }
  if (ret_cr) {
    return TaskStatus::success;
  } else {
    return TaskStatus::fail;
  }
}




TaskStatus TimeIntegratorTaskList::SetBoundariesCRTC(MeshBlock *pmb, int stage) {
  if (stage <= nstages) {
    if (CR_ENABLED)
      pmb->pcr->cr_bvar.SetBoundaries();
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::CRTCOpacity(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  CosmicRay *pcr = pmb->pcr;
  Field *pf = pmb->pfield;
  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      if (CR_ENABLED)
        pcr->UpdateOpacity(pmb, pcr->u_cr, ph->w, pf->bcc);
    }
    return TaskStatus::success;
  }
  return TaskStatus::fail;
}


TaskStatus TimeIntegratorTaskList::AddSourceTermsCRTC(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  CosmicRay *pcr = pmb->pcr;
  Field *pf = pmb->pfield;

  if (stage <= nstages) {
    if (stage_wghts[stage-1].main_stage) {
      // Real t_start_stage = pmb->pmy_mesh->time
      //                      + stage_wghts[(stage-1)].sbeta*pmb->pmy_mesh->dt;
      // Scaled coefficient for RHS update
      Real dt = (stage_wghts[(stage-1)].beta)*(pmb->pmy_mesh->dt);

      // Evaluate the time-dependent source terms
      // Both u and ir are partially updated, only w is from the beginning of the step
      if (CR_ENABLED)
        pcr->pcrintegrator->AddSourceTerms(pmb, dt, ph->u, ph->w, pf->bcc, pcr->u_cr);
    }
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}
