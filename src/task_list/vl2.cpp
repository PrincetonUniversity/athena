//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file vl2.hpp
//  \brief implement steps for [second-order] Van Leer integrator
//======================================================================================

// Athena++ classes headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../bvals/bvals.hpp"
#include "../eos/eos.hpp"
#include "../hydro/srcterms/srcterms.hpp"

// this class header
#include "task_list.hpp"

//--------------------------------------------------------------------------------------
//! \fn
//  \brief

void TaskList::CreateVL2Integrator(Mesh *pm)
{
  if (MAGNETIC_FIELDS_ENABLED) { // MHD
    // MHD predict
    AddTask(1,CALC_FLX,1,NONE);
    if(pm->multilevel==true) { // SMR or AMR
      AddTask(1,FLX_SEND,1,CALC_FLX);
      AddTask(1,FLX_RECV,1,CALC_FLX);
      AddTask(1,HYD_INT, 1,FLX_RECV);
      AddTask(1,HYD_SEND,1,HYD_INT);
      AddTask(1,CALC_EMF,1,CALC_FLX);
      AddTask(1,EMF_SEND,1,CALC_EMF);
      AddTask(1,EMF_RECV,1,EMF_SEND);
      AddTask(1,FLD_INT, 1,EMF_RECV);
    } else {
      AddTask(1,HYD_INT, 1,CALC_FLX);
      AddTask(1,HYD_SEND,1,HYD_INT);
      AddTask(1,CALC_EMF,1,CALC_FLX);
      AddTask(1,EMF_SEND,1,CALC_EMF);
      AddTask(1,EMF_RECV,1,EMF_SEND);
      AddTask(1,FLD_INT, 1,EMF_RECV);
    }
    AddTask(1,FLD_SEND,1,FLD_INT);
    AddTask(1,HYD_RECV,1,NONE);
    AddTask(1,FLD_RECV,1,NONE);
    if(pm->multilevel==true) { // SMR or AMR
      AddTask(1,PROLONG,1,(HYD_SEND|HYD_RECV|FLD_SEND|FLD_RECV));
      AddTask(1,CON2PRIM,1,PROLONG);
    }
    else {
      AddTask(1,CON2PRIM,1,(HYD_INT|HYD_RECV|FLD_INT|FLD_RECV));
    }
    AddTask(1,PHY_BVAL,1,CON2PRIM);

    // MHD correct
    AddTask(2,CALC_FLX,1,PHY_BVAL);
    if(pm->multilevel==true) { // SMR or AMR
      AddTask(2,FLX_SEND,2,CALC_FLX);
      AddTask(2,FLX_RECV,2,CALC_FLX);
      AddTask(2,HYD_INT, 2,FLX_RECV);
      AddTask(2,HYD_SEND,2,HYD_INT);
      AddTask(2,CALC_EMF,2,CALC_FLX);
      AddTask(2,EMF_SEND,2,CALC_EMF);
      AddTask(2,EMF_RECV,2,EMF_SEND);
      AddTask(2,FLD_INT, 2,EMF_RECV);
    } else {
      AddTask(2,HYD_INT, 2,CALC_FLX);
      AddTask(2,HYD_SEND,2,HYD_INT);
      AddTask(2,CALC_EMF,2,CALC_FLX);
      AddTask(2,EMF_SEND,2,CALC_EMF);
      AddTask(2,EMF_RECV,2,EMF_SEND);
      AddTask(2,FLD_INT, 2,EMF_RECV);
    }
    AddTask(2,FLD_SEND,2,FLD_INT);
    AddTask(2,HYD_RECV,1,PHY_BVAL);
    AddTask(2,FLD_RECV,1,PHY_BVAL);
    if(pm->multilevel==true) { // SMR or AMR
      AddTask(2,PROLONG,2,(HYD_SEND|HYD_RECV|FLD_SEND|FLD_RECV));
      AddTask(2,CON2PRIM,2,PROLONG);
    }
    else {
      AddTask(2,CON2PRIM,2,(HYD_INT|HYD_RECV|FLD_INT|FLD_RECV));
    }
    AddTask(2,PHY_BVAL,2,CON2PRIM);
  }
  else {
    // Hydro predict
    AddTask(1,CALC_FLX,1,NONE);
    if(pm->multilevel==true) { // SMR or AMR
      AddTask(1,FLX_SEND,1,CALC_FLX);
      AddTask(1,FLX_RECV,1,CALC_FLX);
      AddTask(1,HYD_INT, 1,FLX_RECV);
    }
    else {
      AddTask(1,HYD_INT, 1,CALC_FLX);
    }
    AddTask(1,HYD_SEND,1,HYD_INT);
    AddTask(1,HYD_RECV,1,NONE);
    if(pm->multilevel==true) { // SMR or AMR
      AddTask(1,PROLONG,1,(HYD_SEND|HYD_RECV));
      AddTask(1,CON2PRIM,1,PROLONG);
    } else {
      AddTask(1,CON2PRIM,1,(HYD_INT|HYD_RECV));
    }
    AddTask(1,PHY_BVAL,1,CON2PRIM);

    // Hydro correct
    AddTask(2,CALC_FLX,1,PHY_BVAL);
    if(pm->multilevel==true) { // SMR or AMR
      AddTask(2,FLX_SEND,2,CALC_FLX);
      AddTask(2,FLX_RECV,2,CALC_FLX);
      AddTask(2,HYD_INT, 2,FLX_RECV);
    }
    else {
      AddTask(2,HYD_INT, 2,CALC_FLX);
    }
    AddTask(2,HYD_SEND,2,HYD_INT);
    AddTask(2,HYD_RECV,1,PHY_BVAL);
    if(pm->multilevel==true) { // SMR or AMR
      AddTask(2,PROLONG,2,(HYD_SEND|HYD_RECV));
      AddTask(2,CON2PRIM,2,PROLONG);
    } else {
      AddTask(2,CON2PRIM,2,(HYD_INT|HYD_RECV));
    }
    AddTask(2,PHY_BVAL,2,CON2PRIM);
  }

  AddTask(2,USERWORK,2,PHY_BVAL);

  // New timestep on mesh block
  AddTask(2,NEW_DT,2,USERWORK);
  if(pm->adaptive==true) {
    AddTask(2,AMR_FLAG,2,USERWORK);
  }
}
