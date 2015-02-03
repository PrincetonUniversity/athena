//======================================================================================
//! \file tasklist.hpp
//  \brief task functions
//======================================================================================
#include <iostream>   // cout, endl
#include <sstream>    // sstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

#include "athena.hpp"
#include "tasklist.hpp"
#include "mesh.hpp"
#include "fluid/fluid.hpp"
#include "field/field.hpp"
#include "bvals/bvals.hpp"
#include "fluid/eos/eos.hpp"
#include "fluid/integrators/fluid_integrator.hpp"
#include "field/integrators/field_integrator.hpp"

namespace taskfunc {

bool Primitives(MeshBlock *pmb, int task_arg)
{
  Fluid *pfluid=pmb->pfluid;
  Field *pfield=pmb->pfield;
  if(task_arg==0)
    pfluid->pf_eos->ConservedToPrimitive(pfluid->u, pfluid->w1, pfield->b,
                                         pfluid->w, pfield->bcc);
  else if(task_arg==1)
    pfluid->pf_eos->ConservedToPrimitive(pfluid->u1, pfluid->w, pfield->b1,
                                         pfluid->w1, pfield->bcc1);

  return true;
}

bool FluidIntegrateSendX1(MeshBlock *pmb, int task_arg)
{
  Fluid *pfluid=pmb->pfluid;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  if(task_arg==0) {
    pfluid->pf_integrator->OneStep(pmb, pfluid->u, pfluid->w1, pfield->b1,
                                   pfield->bcc1, 2);
    if (MAGNETIC_FIELDS_ENABLED)
      pbval->LoadAndSendEFluxBoundaryBuffer(pfield->ei,pfield->wght,0);
    pbval->LoadAndSendFluidBoundaryBuffer(outer_x1,pfluid->u,0);
    pbval->LoadAndSendFluidBoundaryBuffer(inner_x1,pfluid->u,0);
  }
  else if(task_arg==1) {
    pfluid->u1 = pfluid->u;
    pfluid->pf_integrator->OneStep(pmb, pfluid->u1, pfluid->w, pfield->b,
                                   pfield->bcc, 1);
    if (MAGNETIC_FIELDS_ENABLED)
      pbval->LoadAndSendEFluxBoundaryBuffer(pfield->ei,pfield->wght,1);
    pbval->LoadAndSendFluidBoundaryBuffer(outer_x1,pfluid->u1,1);
    pbval->LoadAndSendFluidBoundaryBuffer(inner_x1,pfluid->u1,1);
  }

  return true;
}


bool EFluxReceive(MeshBlock *pmb, int task_arg)
{
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  return pbval->ReceiveAndSetEFluxBoundary(pfield->ei,pfield->wght,task_arg);
}


bool FieldIntegrateSendX1(MeshBlock *pmb, int task_arg)
{
  Fluid *pfluid=pmb->pfluid;
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  if(task_arg==0) {
    pfield->pint->CT(pmb, pfield->b, pfluid->w1, pfield->bcc1, 2);
    pbval->LoadAndSendFieldBoundaryBuffer(outer_x1,pfield->b,0);
    pbval->LoadAndSendFieldBoundaryBuffer(inner_x1,pfield->b,0);
  }
  else if(task_arg==1) {
    pfield->b1.x1f = pfield->b.x1f;
    pfield->b1.x2f = pfield->b.x2f;
    pfield->b1.x3f = pfield->b.x3f;
    pfield->pint->CT(pmb, pfield->b1, pfluid->w, pfield->bcc, 1);
    pbval->LoadAndSendFieldBoundaryBuffer(outer_x1,pfield->b1,1);
    pbval->LoadAndSendFieldBoundaryBuffer(inner_x1,pfield->b1,1);
  }
  return true;
}

bool FluidRecvX1SendX2(MeshBlock *pmb, int task_arg)
{
  Fluid *pfluid=pmb->pfluid;
  BoundaryValues *pbval=pmb->pbval;
  bool lf, rf;
  if(task_arg==0) {
    lf=pbval->ReceiveAndSetFluidBoundary(inner_x1,pfluid->u,0);
    rf=pbval->ReceiveAndSetFluidBoundary(outer_x1,pfluid->u,0);
    if(lf==true && rf==true) {
      if(pmb->pmy_mesh->mesh_size.nx2>1) { // 2D or 3D
        pbval->LoadAndSendFluidBoundaryBuffer(outer_x2,pfluid->u,0);
        pbval->LoadAndSendFluidBoundaryBuffer(inner_x2,pfluid->u,0);
      }
      return true;
    }
  }
  else if(task_arg==1) {
    lf=pbval->ReceiveAndSetFluidBoundary(inner_x1,pfluid->u1,1);
    rf=pbval->ReceiveAndSetFluidBoundary(outer_x1,pfluid->u1,1);
    if(lf==true && rf==true) {
      if(pmb->pmy_mesh->mesh_size.nx2>1) { // 2D or 3D
        pbval->LoadAndSendFluidBoundaryBuffer(outer_x2,pfluid->u1,1);
        pbval->LoadAndSendFluidBoundaryBuffer(inner_x2,pfluid->u1,1);
      }
      return true;
    }
  }
  return false;
}

bool FieldRecvX1SendX2(MeshBlock *pmb, int task_arg)
{
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  bool lf, rf;
  if(task_arg==0) {
    lf=pbval->ReceiveAndSetFieldBoundary(inner_x1,pfield->b,0);
    rf=pbval->ReceiveAndSetFieldBoundary(outer_x1,pfield->b,0);
    if(lf==true && rf==true) {
      if(pmb->pmy_mesh->mesh_size.nx2>1) { // 2D or 3D
        pbval->LoadAndSendFieldBoundaryBuffer(outer_x2,pfield->b,0);
        pbval->LoadAndSendFieldBoundaryBuffer(inner_x2,pfield->b,0);
      }
      return true;
    }
  }
  else if(task_arg==1) {
    lf=pbval->ReceiveAndSetFieldBoundary(inner_x1,pfield->b1,1);
    rf=pbval->ReceiveAndSetFieldBoundary(outer_x1,pfield->b1,1);
    if(lf==true && rf==true) {
      if(pmb->pmy_mesh->mesh_size.nx2>1) { // 2D or 3D
        pbval->LoadAndSendFieldBoundaryBuffer(outer_x2,pfield->b1,1);
        pbval->LoadAndSendFieldBoundaryBuffer(inner_x2,pfield->b1,1);
      }
      return true;
    }
  }
  return false;
}

bool FluidRecvX2SendX3(MeshBlock *pmb, int task_arg)
{
  Fluid *pfluid=pmb->pfluid;
  BoundaryValues *pbval=pmb->pbval;
  bool lf, rf;
  if(task_arg==0) {
    lf=pbval->ReceiveAndSetFluidBoundary(inner_x2,pfluid->u,0);
    rf=pbval->ReceiveAndSetFluidBoundary(outer_x2,pfluid->u,0);
    if(lf==true && rf==true) {
      if(pmb->pmy_mesh->mesh_size.nx3>1) { // 3D
        pbval->LoadAndSendFluidBoundaryBuffer(outer_x3,pfluid->u,0);
        pbval->LoadAndSendFluidBoundaryBuffer(inner_x3,pfluid->u,0);
      }
      return true;
    }
  }
  else if(task_arg==1) {
    lf=pbval->ReceiveAndSetFluidBoundary(inner_x2,pfluid->u1,1);
    rf=pbval->ReceiveAndSetFluidBoundary(outer_x2,pfluid->u1,1);
    if(lf==true && rf==true) {
      if(pmb->pmy_mesh->mesh_size.nx3>1) { // 3D
        pbval->LoadAndSendFluidBoundaryBuffer(outer_x3,pfluid->u1,1);
        pbval->LoadAndSendFluidBoundaryBuffer(inner_x3,pfluid->u1,1);
      }
      return true;
    }
  }
  return false;
}

bool FieldRecvX2SendX3(MeshBlock *pmb, int task_arg)
{
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  bool lf, rf;
  if(task_arg==0) {
    lf=pbval->ReceiveAndSetFieldBoundary(inner_x2,pfield->b,0);
    rf=pbval->ReceiveAndSetFieldBoundary(outer_x2,pfield->b,0);
    if(lf==true && rf==true) {
      if(pmb->pmy_mesh->mesh_size.nx3>1) { // 2D or 3D
        pbval->LoadAndSendFieldBoundaryBuffer(outer_x3,pfield->b,0);
        pbval->LoadAndSendFieldBoundaryBuffer(inner_x3,pfield->b,0);
      }
      return true;
    }
  }
  else if(task_arg==1) {
    lf=pbval->ReceiveAndSetFieldBoundary(inner_x2,pfield->b1,1);
    rf=pbval->ReceiveAndSetFieldBoundary(outer_x2,pfield->b1,1);
    if(lf==true && rf==true) {
      if(pmb->pmy_mesh->mesh_size.nx3>1) { // 3D
        pbval->LoadAndSendFieldBoundaryBuffer(outer_x3,pfield->b1,1);
        pbval->LoadAndSendFieldBoundaryBuffer(inner_x3,pfield->b1,1);
      }
      return true;
    }
  }
  return false;
}

bool FluidRecvX3(MeshBlock *pmb, int task_arg)
{
  Fluid *pfluid=pmb->pfluid;
  BoundaryValues *pbval=pmb->pbval;
  bool lf, rf;
  if(task_arg==0) {
    lf=pbval->ReceiveAndSetFluidBoundary(inner_x3,pfluid->u,0);
    rf=pbval->ReceiveAndSetFluidBoundary(outer_x3,pfluid->u,0);
  }
  else if(task_arg==1) {
    lf=pbval->ReceiveAndSetFluidBoundary(inner_x3,pfluid->u1,1);
    rf=pbval->ReceiveAndSetFluidBoundary(outer_x3,pfluid->u1,1);
  }
  if(lf==true && rf==true)
    return true;
  return false;
}

bool FieldRecvX3(MeshBlock *pmb, int task_arg)
{
  Field *pfield=pmb->pfield;
  BoundaryValues *pbval=pmb->pbval;
  bool lf, rf;
  if(task_arg==0) {
    lf=pbval->ReceiveAndSetFieldBoundary(inner_x3,pfield->b,0);
    rf=pbval->ReceiveAndSetFieldBoundary(outer_x3,pfield->b,0);
  }
  else if(task_arg==1) {
    lf=pbval->ReceiveAndSetFieldBoundary(inner_x3,pfield->b1,1);
    rf=pbval->ReceiveAndSetFieldBoundary(outer_x3,pfield->b1,1);
  }
  if(lf==true && rf==true)
    return true;
  return false;
}

bool NewBlockTimeStep(MeshBlock *pmb, int task_arg)
{
  pmb->pfluid->NewBlockTimeStep(pmb);
  return true;
}

} // namespace task

unsigned long DependFlag(enum task dep)
{
  return 1L<<dep;
}

void TaskList::AddTask(enum task t, unsigned long int dependence)
{
  std::stringstream msg;
  task[ntask].taskid=t;
  task[ntask].depend=dependence;
  switch(t)
  {
  case primitives_0:
    task[ntask].TaskFunc=taskfunc::Primitives;
    task[ntask].task_arg=0;
    break;

  case fluid_integrate_sendx1_0: // fluid correction
    task[ntask].TaskFunc=taskfunc::FluidIntegrateSendX1;
    task[ntask].task_arg=0;
    break;

  case eflux_recv_0:
    task[ntask].TaskFunc=taskfunc::EFluxReceive;
    task[ntask].task_arg=0;
    break;

  case field_integrate_sendx1_0: // field correction
    task[ntask].TaskFunc=taskfunc::FieldIntegrateSendX1;
    task[ntask].task_arg=0;
    break;

  case fluid_recvx1_sendx2_0: // same as fluid_recvx1_0
    task[ntask].TaskFunc=taskfunc::FluidRecvX1SendX2;
    task[ntask].task_arg=0;
    break;

  case field_recvx1_sendx2_0: // same as field_recvx1_0
    task[ntask].TaskFunc=taskfunc::FieldRecvX1SendX2;
    task[ntask].task_arg=0;
    break;

  case fluid_recvx2_sendx3_0:
    task[ntask].TaskFunc=taskfunc::FluidRecvX2SendX3;
    task[ntask].task_arg=0;
    break;

  case field_recvx2_sendx3_0:
    task[ntask].TaskFunc=taskfunc::FieldRecvX2SendX3;
    task[ntask].task_arg=0;
    break;

  case fluid_recvx3_0:
    task[ntask].TaskFunc=taskfunc::FluidRecvX3;
    task[ntask].task_arg=0;
    break;

  case field_recvx3_0:
    task[ntask].TaskFunc=taskfunc::FieldRecvX3;
    task[ntask].task_arg=0;
    break;


  case primitives_1:
    task[ntask].TaskFunc=taskfunc::Primitives;
    task[ntask].task_arg=1;
    break;

  case fluid_integrate_sendx1_1: // fluid prediction
    task[ntask].TaskFunc=taskfunc::FluidIntegrateSendX1;
    task[ntask].task_arg=1;
    break;

  case eflux_recv_1:
    task[ntask].TaskFunc=taskfunc::EFluxReceive;
    task[ntask].task_arg=1;
    break;
  
  case field_integrate_sendx1_1: // field prediction
    task[ntask].TaskFunc=taskfunc::FieldIntegrateSendX1;
    task[ntask].task_arg=1;
    break;

  case fluid_recvx1_sendx2_1:
    task[ntask].TaskFunc=taskfunc::FluidRecvX1SendX2;
    task[ntask].task_arg=1;
    break;

  case field_recvx1_sendx2_1:
    task[ntask].TaskFunc=taskfunc::FieldRecvX1SendX2;
    task[ntask].task_arg=1;
    break;

  case fluid_recvx2_sendx3_1:
    task[ntask].TaskFunc=taskfunc::FluidRecvX2SendX3;
    task[ntask].task_arg=1;
    break;

  case field_recvx2_sendx3_1:
    task[ntask].TaskFunc=taskfunc::FieldRecvX2SendX3;
    task[ntask].task_arg=1;
    break;

  case fluid_recvx3_1:
    task[ntask].TaskFunc=taskfunc::FluidRecvX3;
    task[ntask].task_arg=1;
    break;

  case field_recvx3_1:
    task[ntask].TaskFunc=taskfunc::FieldRecvX3;
    task[ntask].task_arg=1;
    break;

  case new_blocktimestep:
    task[ntask].TaskFunc=taskfunc::NewBlockTimeStep;
    task[ntask].task_arg=0;
    break;

  default:
    msg << "### FATAL ERROR in DoOneTask" << std::endl
        << "Invalid Task "<< t << " is specified" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  
  }
  ntask++;
  return;
}

