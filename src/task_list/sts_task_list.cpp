//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file sts_task_list.cpp
//! \brief RKL Super-Time-Stepping
//!
//! REFERENCE:
//! Meyer, C. D., Balsara, D. S., & Aslam, T. D. 2014, J. Comput. Phys., 257A, 594-626

// C headers

// C++ headers
#include <algorithm>  // std::binary_search
#include <cstring>    // strcmp()
#include <iomanip>    // std::setprecision()
#include <iostream>   // endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <vector>     // std::vector

// Athena++ headers
#include "../athena.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "../scalars/scalars.hpp"
#include "task_list.hpp"

//----------------------------------------------------------------------------------------
//! SuperTimeStepTaskList constructor

SuperTimeStepTaskList::SuperTimeStepTaskList(
    ParameterInput *pin, Mesh *pm, TimeIntegratorTaskList *ptlist) :
    sts_max_dt_ratio(pin->GetOrAddReal("time", "sts_max_dt_ratio", -1.0)),
    ptlist_(ptlist) {
  // Read a flag for shear periodic
  SHEAR_PERIODIC = pm->shear_periodic;

  // Check for valid STS integrator:
  if (pm->sts_integrator != "rkl1" && pm->sts_integrator != "rkl2") {
    std::stringstream msg;
    msg << "### FATAL ERROR in SuperTimeStepTaskList constructor" << std::endl
        << "sts_integrator=" << pm->sts_integrator
        << " not valid STS integrator" << std::endl;
    ATHENA_ERROR(msg);
  }

  MeshBlock *pmb = pm->my_blocks(0);
  // Check for STS incompatiblities:
  if (!(pmb->phydro->hdif.hydro_diffusion_defined)
      // short-circuit evaluation makes these safe (won't dereference pscalars=nullptr):
      && !(MAGNETIC_FIELDS_ENABLED && pmb->pfield->fdif.field_diffusion_defined)
      && !(NSCALARS > 0 && pmb->pscalars->scalar_diffusion_defined)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in SuperTimeStepTaskList" << std::endl
        << "Super-time-stepping requires setting parameters for "
        << "at least one diffusive process in the input file." << std::endl;
    ATHENA_ERROR(msg);
  }

  // identify what diffusion processes will be integrated with STS and
  // assemble a vector containing the subset of NHYDRO indices updated inside STS
  do_sts_hydro = false;
  do_sts_field = false;
  do_sts_scalar = false;
  if (pmb->phydro->hdif.hydro_diffusion_defined) {
    do_sts_hydro = true;
    if (pmb->phydro->hdif.nu_iso > 0.0
        || pmb->phydro->hdif.nu_aniso > 0.0) {
      sts_idx_subset.push_back(IM1);
      sts_idx_subset.push_back(IM2);
      sts_idx_subset.push_back(IM3);
      if (NON_BAROTROPIC_EOS) {
        sts_idx_subset.push_back(IEN);
      }
    }
    if (pmb->phydro->hdif.kappa_iso > 0.0
        || pmb->phydro->hdif.kappa_aniso > 0.0) {
      if (!std::binary_search(sts_idx_subset.begin(), sts_idx_subset.end(), IEN)) {
        sts_idx_subset.push_back(IEN);
      }
    }
  }
  if (MAGNETIC_FIELDS_ENABLED) {
    if (pmb->pfield->fdif.field_diffusion_defined) {
      do_sts_field = true;
      if (NON_BAROTROPIC_EOS) {
        // Below designations are needed to account for Poynting flux
        do_sts_hydro = true;
        if (!std::binary_search(sts_idx_subset.begin(), sts_idx_subset.end(), IEN)) {
          sts_idx_subset.push_back(IEN);
        }
      }
    }
  }
  if (NSCALARS > 0) {
    if (pmb->pscalars->scalar_diffusion_defined) {
      do_sts_scalar = true;
    }
  }

  // Now assemble list of tasks for each stage of SuperTimeStep integrator
  {using namespace HydroIntegratorTaskNames; // NOLINT (build/namespace)
    // calculate hydro/field diffusive fluxes
    if (do_sts_hydro) {
      AddTask(DIFFUSE_HYD,NONE);
    }
    if (do_sts_field) {
      AddTask(DIFFUSE_FLD,NONE);
      if (do_sts_hydro) {
        // compute hydro fluxes
        AddTask(CALC_HYDFLX,(DIFFUSE_HYD|DIFFUSE_FLD));
      }
    } else if (do_sts_hydro) {
      AddTask(CALC_HYDFLX,DIFFUSE_HYD);
    }
    if (do_sts_scalar) {
      AddTask(DIFFUSE_SCLR,NONE);
      if (do_sts_hydro) {
        AddTask(CALC_SCLRFLX,(CALC_HYDFLX|DIFFUSE_SCLR));
      } else {
        AddTask(CALC_SCLRFLX,DIFFUSE_SCLR);
      }
    }

    if (do_sts_hydro) {
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
      AddTask(SEND_HYD,INT_HYD);
      AddTask(RECV_HYD,NONE);
      AddTask(SETB_HYD,(RECV_HYD|INT_HYD));
      if (SHEAR_PERIODIC) {
        AddTask(SEND_HYDSH,SETB_HYD);
        AddTask(RECV_HYDSH,SEND_HYDSH);
      }
    }

    if (do_sts_scalar) {
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
      AddTask(SEND_SCLR,INT_SCLR);
      AddTask(RECV_SCLR,NONE);
      AddTask(SETB_SCLR,(RECV_SCLR|INT_SCLR));
      if (SHEAR_PERIODIC) {
        AddTask(SEND_SCLRSH,SETB_SCLR);
        AddTask(RECV_SCLRSH,SEND_SCLRSH);
      }
    }

    if (do_sts_field) { // field diffusion
      // compute MHD fluxes, integrate field
      if (do_sts_hydro) {
        AddTask(CALC_FLDFLX,CALC_HYDFLX);
      } else {
        AddTask(CALC_FLDFLX,DIFFUSE_FLD);
      }
      AddTask(SEND_FLDFLX,CALC_FLDFLX);
      AddTask(RECV_FLDFLX,SEND_FLDFLX);
      if (SHEAR_PERIODIC) {
        AddTask(SEND_EMFSH,RECV_FLDFLX);
        AddTask(RECV_EMFSH,RECV_FLDFLX);
        AddTask(INT_FLD,RECV_EMFSH);
      } else {
        AddTask(INT_FLD,RECV_FLDFLX);
      }
      AddTask(SEND_FLD,INT_FLD);
      AddTask(RECV_FLD,NONE);
      AddTask(SETB_FLD,(RECV_FLD|INT_FLD));
      if (SHEAR_PERIODIC) {
        AddTask(SEND_FLDSH,SETB_FLD);
        AddTask(RECV_FLDSH,SEND_FLDSH);
      }

      // prolongate, compute new primitives
      if (pm->multilevel) { // SMR or AMR
        if (do_sts_scalar) {
          if (do_sts_hydro) {
            if (SHEAR_PERIODIC) {
              AddTask(PROLONG,(SEND_HYD|RECV_HYDSH|SEND_FLD|RECV_FLDSH
                               |SEND_SCLR|RECV_SCLRSH));
            } else {
              AddTask(PROLONG,(SEND_HYD|SETB_HYD|SEND_FLD|SETB_FLD|SEND_SCLR|SETB_SCLR));
            }
          } else {
            if (SHEAR_PERIODIC) {
              AddTask(PROLONG,(SEND_FLD|RECV_FLDSH|SEND_SCLR|RECV_SCLRSH));
            } else {
              AddTask(PROLONG,(SEND_FLD|SETB_FLD|SEND_SCLR|SETB_SCLR));
            }
          }
        } else {
          if (do_sts_hydro) {
            if (SHEAR_PERIODIC) {
              AddTask(PROLONG,(SEND_HYD|RECV_HYDSH|SEND_FLD|RECV_FLDSH));
            } else {
              AddTask(PROLONG,(SEND_HYD|SETB_HYD|SEND_FLD|SETB_FLD));
            }
          } else {
            if (SHEAR_PERIODIC) {
              AddTask(PROLONG,(SEND_FLD|RECV_FLDSH));
            } else {
              AddTask(PROLONG,(SEND_FLD|SETB_FLD));
            }
          }
        }
        AddTask(CONS2PRIM,PROLONG);
      } else {
        if (do_sts_scalar) {
          if (do_sts_hydro) {
            if (SHEAR_PERIODIC) {
              AddTask(CONS2PRIM,(RECV_HYDSH|RECV_FLDSH|RECV_SCLRSH));
            } else {
              AddTask(CONS2PRIM,(SETB_HYD|SETB_FLD|SETB_SCLR));
            }
          } else {
            if (SHEAR_PERIODIC) {
              AddTask(CONS2PRIM,(RECV_FLDSH|RECV_SCLRSH));
            } else {
              AddTask(CONS2PRIM,(SETB_FLD|SETB_SCLR));
            }
          }
        } else {
          if (do_sts_hydro) {
            if (SHEAR_PERIODIC) {
              AddTask(CONS2PRIM,(RECV_HYDSH|RECV_FLDSH));
            } else {
              AddTask(CONS2PRIM,(SETB_HYD|SETB_FLD));
            }
          } else {
            if (SHEAR_PERIODIC) {
              AddTask(CONS2PRIM,RECV_FLDSH);
            } else {
              AddTask(CONS2PRIM,SETB_FLD);
            }
          }
        }
      }
    } else if (do_sts_hydro) {  // no field diffusion
      // prolongate, compute new primitives
      if (pm->multilevel) { // SMR or AMR
        if (do_sts_scalar) {
          if (SHEAR_PERIODIC) {
            AddTask(PROLONG,(SEND_HYD|RECV_HYDSH|SEND_SCLR|RECV_SCLRSH));
          } else {
            AddTask(PROLONG,(SEND_HYD|SETB_HYD|SEND_SCLR|SETB_SCLR));
          }
        } else {
          if (SHEAR_PERIODIC) {
            AddTask(PROLONG,(SEND_HYD|RECV_HYDSH));
          } else {
            AddTask(PROLONG,(SEND_HYD|SETB_HYD));
          }
        }
        AddTask(CONS2PRIM,PROLONG);
      } else {
        if (do_sts_scalar) {
          if (SHEAR_PERIODIC) {
            AddTask(CONS2PRIM,(RECV_HYDSH|RECV_SCLRSH));
          } else {
            AddTask(CONS2PRIM,(SETB_HYD|SETB_SCLR));
          }
        } else {
          if (SHEAR_PERIODIC) {
            AddTask(CONS2PRIM,RECV_HYDSH);
          } else {
            AddTask(CONS2PRIM,SETB_HYD);
          }
        }
      }
    } else {  // only scalar diffusion
      // prolongate, compute new primitives
      if (pm->multilevel) { // SMR or AMR
        if (SHEAR_PERIODIC) {
          AddTask(PROLONG,(SEND_SCLR|RECV_SCLRSH));
        } else {
          AddTask(PROLONG,(SEND_SCLR|SETB_SCLR));
        }
        AddTask(CONS2PRIM,PROLONG);
      } else {
        if (SHEAR_PERIODIC) {
          AddTask(CONS2PRIM,RECV_SCLRSH);
        } else {
          AddTask(CONS2PRIM,SETB_SCLR);
        }
      }
    }

    // everything else
    AddTask(PHY_BVAL,CONS2PRIM);
    if (pm->sts_integrator == "rkl2") {
      // TaskType::op_split_after tasks:
      AddTask(USERWORK,PHY_BVAL);
      AddTask(NEW_DT,USERWORK);
      if (pm->adaptive) {
        AddTask(FLAG_AMR,USERWORK);
        AddTask(CLEAR_ALLBND,FLAG_AMR);
      } else {
        AddTask(CLEAR_ALLBND,NEW_DT);
      }
    } else { // rkl1, no TaskType::op_split_after tasks
      AddTask(CLEAR_ALLBND,PHY_BVAL);
    }
  } // end of using namespace block
}

//---------------------------------------------------------------------------------------
//! Sets id and dependency for "ntask" member of task_list_ array, then iterates value of
//! ntask.

void SuperTimeStepTaskList::AddTask(const TaskID& id, const TaskID& dep) {
  task_list_[ntasks].task_id = id;
  task_list_[ntasks].dependency = dep;

  using namespace HydroIntegratorTaskNames; // NOLINT (build/namespace)

  if (id == CLEAR_ALLBND) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&SuperTimeStepTaskList::ClearAllBoundary_STS);
    task_list_[ntasks].lb_time = false;
  } else if (id == CALC_HYDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&SuperTimeStepTaskList::CalculateHydroFlux_STS);
    task_list_[ntasks].lb_time = true;
  } else if (id == CALC_FLDFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&SuperTimeStepTaskList::CalculateEMF_STS);
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
        (&SuperTimeStepTaskList::IntegrateHydro_STS);
    task_list_[ntasks].lb_time = true;
  } else if (id == INT_FLD) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&SuperTimeStepTaskList::IntegrateField_STS);
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
  } else if (id == CALC_SCLRFLX) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&SuperTimeStepTaskList::CalculateScalarFlux_STS);
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
        (&SuperTimeStepTaskList::IntegrateScalars_STS);
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
  } else if (id == PROLONG) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&SuperTimeStepTaskList::Prolongation_STS);
    task_list_[ntasks].lb_time = true;
  } else if (id == CONS2PRIM) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&SuperTimeStepTaskList::Primitives_STS);
    task_list_[ntasks].lb_time = true;
  } else if (id == USERWORK) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&SuperTimeStepTaskList::UserWork_STS);
    task_list_[ntasks].lb_time = true;
  } else if (id == NEW_DT) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&SuperTimeStepTaskList::NewBlockTimeStep_STS);
    task_list_[ntasks].lb_time = true;
  } else if (id == FLAG_AMR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&SuperTimeStepTaskList::CheckRefinement_STS);
    task_list_[ntasks].lb_time = true;
  } else if (id == PHY_BVAL) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&SuperTimeStepTaskList::PhysicalBoundary_STS);
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
  } else if (id == DIFFUSE_SCLR) {
    task_list_[ntasks].TaskFunc=
        static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
        (&TimeIntegratorTaskList::DiffuseScalars);
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

//---------------------------------------------------------------------------------------
//! \brief

void SuperTimeStepTaskList::StartupTaskList(MeshBlock *pmb, int stage) {
  Mesh *pm =  pmb->pmy_mesh;
#pragma omp single
  {
    if (pm->sts_integrator == "rkl2") { // Set RKL2 params
      Real bj    = (std::pow((stage   ), 2.) + (stage   ) - 2.)
                   / (2.*(stage   )*((stage   ) + 1.));
      Real bj_m1 = (std::pow((stage-1.), 2.) + (stage-1.) - 2.)
                   / (2.*(stage-1.)*((stage-1.) + 1.));
      Real bj_m2 = (std::pow((stage-2.), 2.) + (stage-2.) - 2.)
                   / (2.*(stage-2.)*((stage-2.) + 1.));
      if (stage == 1 || stage == 2) {
        bj = bj_m1 = bj_m2 = 1./3.;
      } else if (stage == 3) {
        bj_m1 = bj_m2 = 1./3.;
      } else if (stage == 4) {
        bj_m2 = 1./3.;
      }
      pm->muj = (2.*stage - 1.)/stage*bj/bj_m1;
      pm->nuj = -1.*(stage - 1.)/stage*bj/bj_m2;
      if (stage == 1) {
        pm->muj_tilde = bj*4./(std::pow(nstages, 2.) + nstages - 2.);
        pm->gammaj_tilde = 0.;
      } else {
        pm->muj_tilde = pm->muj*4./(std::pow(nstages, 2.) + nstages - 2.);
        pm->gammaj_tilde = -1.*(1. - bj_m1)*pm->muj_tilde;
      }
    } else { // Set RKL1 params
      pm->muj = (2.*stage - 1.)/stage;
      pm->nuj = (1. - stage)/stage;
      pm->muj_tilde = pm->muj*2./(std::pow(nstages, 2.) + nstages);
    }

    Real dt_ratio = pm->dt / pm->dt_parabolic;
    Real nstages_time_int = ptlist_->nstages;
    Real tot_nstages_sts = nstages;
    if (pm->sts_integrator == "rkl2") {
      tot_nstages_sts *= 2.;
    }
    Real stage_ratio = nstages_time_int*dt_ratio/(tot_nstages_sts + nstages_time_int);

    if (Globals::my_rank == 0) {
      // output additional diagnostics indiciating progress through STS stages:
      if (pm->dt_diagnostics != -1 && pm->ncycle_out != 0
          && pm->ncycle % pm->ncycle_out == 0) {
        const int ratio_precision = 3;
        const int dt_precision = std::numeric_limits<Real>::max_digits10 - 1;
        if (pm->dt_diagnostics == 0) {
          if (stage == nstages) {
            std::cout << "stage=" << stage << "/" << nstages
                      << " dt_parabolic=" << pm->dt_parabolic
                      << " ratio=" << std::setprecision(ratio_precision) <<  dt_ratio
                      << " stage_ratio=" << stage_ratio
                      << std::setprecision(dt_precision)
                      << std::endl;
          }
        } else {
          if (stage % pm->dt_diagnostics == 0) {
            std::cout << "stage=" << stage << "/" << nstages
                      << " dt_parabolic=" << pm->dt_parabolic
                      << " ratio=" << std::setprecision(ratio_precision) << dt_ratio
                      << " stage_ratio=" << stage_ratio
                      << std::setprecision(dt_precision)
                      << std::endl;
          }
        }
      }
    }
  }

  // Compute Shear for Shearing Box at stage = 1
  if (SHEAR_PERIODIC) {
    if (stage == 1) {
      if (pmb->pmy_mesh->sts_loc == TaskType::op_split_before) {
        pmb->pbval->ComputeShear(pm->time, pm->time);
      } else { // op_split_after
        pmb->pbval->ComputeShear(pm->time+pm->dt, pm->time+pm->dt);
      }
    }
  }

  // Clear flux arrays from previous stage
  pmb->phydro->hdif.ClearFlux(pmb->phydro->flux);
  if (MAGNETIC_FIELDS_ENABLED)
    pmb->pfield->fdif.ClearEMF(pmb->pfield->e);
  if (NSCALARS > 0) {
    PassiveScalars *ps = pmb->pscalars;
    ps->s_flux[X1DIR].ZeroClear();
    ps->s_flux[X2DIR].ZeroClear();
    ps->s_flux[X3DIR].ZeroClear();
  }

  pmb->pbval->StartReceivingSubset(BoundaryCommSubset::all,
                                   pmb->pbval->bvars_sts);

  return;
}

//----------------------------------------------------------------------------------------
//! Functions to end MPI communication

TaskStatus SuperTimeStepTaskList::ClearAllBoundary_STS(MeshBlock *pmb, int stage) {
  pmb->pbval->ClearBoundarySubset(BoundaryCommSubset::all,
                                  pmb->pbval->bvars_sts);
  return TaskStatus::success;
}

//----------------------------------------------------------------------------------------
//! Functions to calculates hydro flux

TaskStatus SuperTimeStepTaskList::CalculateHydroFlux_STS(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  if (stage <= nstages) {
    ph->CalculateFluxes_STS();
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
//! Functions to calculates scalar flux

TaskStatus SuperTimeStepTaskList::CalculateScalarFlux_STS(MeshBlock *pmb, int stage) {
  PassiveScalars *ps = pmb->pscalars;
  if (stage <= nstages) {
    ps->CalculateFluxes_STS();
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
//! Functions to calculates emf

TaskStatus SuperTimeStepTaskList::CalculateEMF_STS(MeshBlock *pmb, int stage) {
  Field *pf = pmb->pfield;
  if (stage <= nstages) {
    pf->ComputeCornerE_STS();
    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
//! Functions to integrate hydro variables

TaskStatus SuperTimeStepTaskList::IntegrateHydro_STS(MeshBlock *pmb, int stage) {
  Hydro *ph = pmb->phydro;
  Field *pf = pmb->pfield;

  if (pmb->pmy_mesh->fluid_setup != FluidFormulation::evolve) return TaskStatus::next;

  // set registers
  if (pmb->pmy_mesh->sts_integrator == "rkl2" && stage == 1) {
    ph->u0 = ph->u;
  }
  ph->u2.SwapAthenaArray(ph->u1);
  ph->u1.SwapAthenaArray(ph->u);

  // update u
  if (stage <= nstages) {
    Real ave_wghts[5];
    ave_wghts[0] = 0.;
    ave_wghts[1] = pmb->pmy_mesh->muj;
    ave_wghts[2] = pmb->pmy_mesh->nuj;
    ave_wghts[3] = 0.;
    ave_wghts[4] = 0.;
    if (pmb->pmy_mesh->sts_integrator == "rkl2") {
      ave_wghts[3] = 1. - pmb->pmy_mesh->muj - pmb->pmy_mesh->nuj;
      ave_wghts[4] = pmb->pmy_mesh->gammaj_tilde;
    }
    pmb->WeightedAve(ph->u, ph->u1, ph->u2, ph->u0, ph->fl_div, ave_wghts);

    Real wght = pmb->pmy_mesh->muj_tilde*pmb->pmy_mesh->dt;
    if (pmb->pmy_mesh->sts_integrator == "rkl2") {
      wght *= 0.5;
    }
    ph->AddFluxDivergence_STS(wght, stage, ph->u, ph->fl_div, sts_idx_subset);
    pmb->pcoord->AddCoordTermsDivergence_STS(wght, stage, ph->flux,
                                             ph->u, ph->fl_div);

    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
//! Functions to integrate scalars

TaskStatus SuperTimeStepTaskList::IntegrateScalars_STS(MeshBlock *pmb, int stage) {
  PassiveScalars *ps = pmb->pscalars;
  // set registers
  if (pmb->pmy_mesh->sts_integrator == "rkl2" && stage == 1) {
    ps->s0 = ps->s;
  }
  ps->s2.SwapAthenaArray(ps->s1);
  ps->s1.SwapAthenaArray(ps->s);

  // update s
  if (stage <= nstages) {
    Real ave_wghts[5];
    ave_wghts[0] = 0.;
    ave_wghts[1] = pmb->pmy_mesh->muj;
    ave_wghts[2] = pmb->pmy_mesh->nuj;
    ave_wghts[3] = 0.;
    ave_wghts[4] = 0.;
    if (pmb->pmy_mesh->sts_integrator == "rkl2") {
      ave_wghts[3] = 1. - pmb->pmy_mesh->muj - pmb->pmy_mesh->nuj;
      ave_wghts[4] = pmb->pmy_mesh->gammaj_tilde;
    }

    pmb->WeightedAve(ps->s, ps->s1, ps->s2, ps->s0, ps->s_fl_div, ave_wghts);

    Real wght = pmb->pmy_mesh->muj_tilde*pmb->pmy_mesh->dt;
    if (pmb->pmy_mesh->sts_integrator == "rkl2") {
      wght *= 0.5;
    }
    ps->AddFluxDivergence_STS(wght, stage, ps->s, ps->s_fl_div);

    return TaskStatus::next;
  }

  return TaskStatus::fail;
}

//----------------------------------------------------------------------------------------
//! Functions to integrate field variables

TaskStatus SuperTimeStepTaskList::IntegrateField_STS(MeshBlock *pmb, int stage) {
  Field *pf = pmb->pfield;

  if (pmb->pmy_mesh->fluid_setup != FluidFormulation::evolve) return TaskStatus::next;

  // set reigsters
  if (pmb->pmy_mesh->sts_integrator == "rkl2" && stage == 1) {
    pf->b0.x1f = pf->b.x1f;
    pf->b0.x2f = pf->b.x2f;
    pf->b0.x3f = pf->b.x3f;
  }
  pf->b2.x1f.SwapAthenaArray(pf->b1.x1f);
  pf->b2.x2f.SwapAthenaArray(pf->b1.x2f);
  pf->b2.x3f.SwapAthenaArray(pf->b1.x3f);

  pf->b1.x1f.SwapAthenaArray(pf->b.x1f);
  pf->b1.x2f.SwapAthenaArray(pf->b.x2f);
  pf->b1.x3f.SwapAthenaArray(pf->b.x3f);


  // update b
  if (stage <= nstages) {
    Real ave_wghts[5];
    ave_wghts[0] = 0.;
    ave_wghts[1] = pmb->pmy_mesh->muj;
    ave_wghts[2] = pmb->pmy_mesh->nuj;
    ave_wghts[3] = 0.;
    ave_wghts[4] = 0.;
    if (pmb->pmy_mesh->sts_integrator == "rkl2") {
      ave_wghts[3] = 1. - pmb->pmy_mesh->muj - pmb->pmy_mesh->nuj;
      ave_wghts[4] = pmb->pmy_mesh->gammaj_tilde;
    }
    pmb->WeightedAve(pf->b, pf->b1, pf->b2, pf->b0, pf->ct_update, ave_wghts);

    Real wght = pmb->pmy_mesh->muj_tilde*pmb->pmy_mesh->dt;
    if (pmb->pmy_mesh->sts_integrator == "rkl2") {
      wght *= 0.5;
    }
    pf->CT_STS(wght, stage, pf->b, pf->ct_update);

    return TaskStatus::next;
  }
  return TaskStatus::fail;
}

//--------------------------------------------------------------------------------------
// Functions for everything else

TaskStatus SuperTimeStepTaskList::Prolongation_STS(MeshBlock *pmb,
                                                   int stage) {
  BoundaryValues *pbval = pmb->pbval;
  if (stage <= nstages) {
    Real time = pmb->pmy_mesh->time;
    if (pmb->pmy_mesh->sts_loc == TaskType::op_split_after) time += pmb->pmy_mesh->dt;
    pbval->ProlongateBoundaries(time, 0.0, pbval->bvars_sts);
  } else {
    return TaskStatus::fail;
  }

  return TaskStatus::success;
}

TaskStatus SuperTimeStepTaskList::Primitives_STS(MeshBlock *pmb, int stage) {
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
    //! \note
    //! At beginning of this task, ph->w and pf->bcc contain the previous stage's output.
    //! Now, we wish to update ph->w and pf->bcc with the present ph->u and pf->b.
    //! In TimeIntegratorTaskList::Primitives, ph->w1 is used to contain the present
    //! stage's output, however, after performing w.SwapAthenaArray(ph->w1), the
    //! resulting ph->w does not contain correct values for the ghost zones.
    //! This is resolved by calling TimeIntegratorTaskList::ApplyPhysicalBoundaries,
    //! however, in STS, ApplyPhysicalBoundaries only cycles through boundary types
    //! integrated in SuperTimeStepTaskList, **not all** boundary types.  We still
    //! need the correct values in ghost zones for all other non-STS-integrated
    //! boundaries, thus we choose to directly update ph->w, instead of storing in an
    //! intermediate ph->w1. ph->w1 is only used in RELATAVISTIC DYNAMICS,
    //! which is not yet compatible with STS, thus making this a safe choice.
    if (do_sts_hydro || do_sts_field) {
      pmb->peos->ConservedToPrimitive(ph->u, ph->w, pf->b,
                                      ph->w, pf->bcc, pmb->pcoord,
                                      il, iu, jl, ju, kl, ku);
      if (pmb->porb->orbital_advection_defined) {
        pmb->porb->ResetOrbitalSystemConversionFlag();
      }
    }

    if (do_sts_scalar) {
      pmb->peos->PassiveScalarConservedToPrimitive(ps->s, ph->w,
                                                   ps->r, ps->r,
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
      if (do_sts_hydro || do_sts_field) {
        pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->w, HydroBoundaryQuantity::prim);
      }
      if (do_sts_scalar) {
        pmb->pscalars->sbvar.var_cc = &(pmb->pscalars->r);
        if (pmb->pmy_mesh->multilevel) {
          pmb->pscalars->sbvar.coarse_buf = &(pmb->pscalars->coarse_r_);
        }
      }
      Real time = pmb->pmy_mesh->time;
      if (pmb->pmy_mesh->sts_loc == TaskType::op_split_after) time += pmb->pmy_mesh->dt;
      pbval->ApplyPhysicalBoundaries(time, 0.0, pbval->bvars_sts);
      // Perform 4th order W(U)
      if (do_sts_hydro || do_sts_field) {
        pmb->peos->ConservedToPrimitiveCellAverage(ph->u, ph->w, pf->b,
                                                   ph->w, pf->bcc, pmb->pcoord,
                                                   il, iu, jl, ju, kl, ku);
      }
      if (do_sts_scalar) {
        pmb->peos->PassiveScalarConservedToPrimitiveCellAverage(
            ps->s, ps->r, ps->r, pmb->pcoord, il, iu, jl, ju, kl, ku);
      }
    }
  } else {
    return TaskStatus::fail;
  }

  return TaskStatus::success;
}


TaskStatus SuperTimeStepTaskList::PhysicalBoundary_STS(MeshBlock *pmb, int stage) {
  BoundaryValues *pbval=pmb->pbval;
  if (stage <= nstages) {
    if (do_sts_hydro || do_sts_field) {
      pmb->phydro->hbvar.SwapHydroQuantity(pmb->phydro->w, HydroBoundaryQuantity::prim);
    }
    if (do_sts_scalar) {
      pmb->pscalars->sbvar.var_cc = &(pmb->pscalars->r);
      if (pmb->pmy_mesh->multilevel) {
        pmb->pscalars->sbvar.coarse_buf = &(pmb->pscalars->coarse_r_);
      }
    }
    Real time = pmb->pmy_mesh->time;
    if (pmb->pmy_mesh->sts_loc == TaskType::op_split_after) time += pmb->pmy_mesh->dt;
    pbval->ApplyPhysicalBoundaries(time, 0.0, pbval->bvars_sts);
  } else {
    return TaskStatus::fail;
  }
  return TaskStatus::success;
}

TaskStatus SuperTimeStepTaskList::UserWork_STS(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  if (pmb->pmy_mesh->sts_loc == TaskType::op_split_after)
    pmb->UserWorkInLoop();
  return TaskStatus::success;
}


TaskStatus SuperTimeStepTaskList::NewBlockTimeStep_STS(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  if (pmb->pmy_mesh->sts_loc == TaskType::op_split_after)
    pmb->phydro->NewBlockTimeStep();
  return TaskStatus::success;
}


TaskStatus SuperTimeStepTaskList::CheckRefinement_STS(MeshBlock *pmb, int stage) {
  if (stage != nstages) return TaskStatus::success; // only do on last stage

  if (pmb->pmy_mesh->sts_loc == TaskType::op_split_after)
    pmb->pmr->CheckRefinementCondition();
  return TaskStatus::success;
}
