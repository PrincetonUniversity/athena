//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file poisson.cpp
//  \brief Problem generator to test Poisson's solver
//

// C++ headers
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <ctime>

// Athena++ headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../gravity/gravity.hpp"
#include "../gravity/fftgravity.hpp"
#include "../mesh/mesh.hpp"

#ifdef OPENMP_PARALLEL
#include "omp.h"
#endif

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  Real four_pi_G = pin->GetReal("problem","four_pi_G");
  SetFourPiG(four_pi_G);
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real x0=0.0, y0=0.0, z0=0.0;
  Real x, y, z;
  Real r, r2;
  Real den, phia;
  RegionSize &mesh_size = pmy_mesh->mesh_size;

  Real x1size = mesh_size.x1max - mesh_size.x1min;
  Real x2size = mesh_size.x2max - mesh_size.x2min;
  Real x3size = mesh_size.x3max - mesh_size.x3min;

  Real four_pi_G = pin->GetReal("problem","four_pi_G");
  Real gconst = four_pi_G / (4.0*PI);
  Real grav_mean_rho = 0.0;
  
  int iprob = pin->GetOrAddInteger("problem","iprob",1);
  int nlim = pin->GetInteger("time","nlim");

  int dim = 1;
  if(mesh_size.nx2 > 1) dim=2;
  if(mesh_size.nx3 > 1) dim=3;
  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
    x = pcoord->x1v(i);
    y = pcoord->x2v(j);
    z = pcoord->x3v(k);
    r2 = SQR(x - x0);
    if(dim > 1) r2 += SQR(y - y0);
    if(dim > 2) r2 += SQR(z - z0);
    if(iprob == 1){
      den = std::sin(2*PI*x/x1size);
      if(dim > 1) den *= std::sin(2*PI*y/x2size);
      if(dim > 2) den *= std::sin(2*PI*z/x3size); // dkx=2*PI/Nx
      phia = SQR(2*PI/x1size);
      if(dim > 1) phia += SQR(2*PI/x2size);
      if(dim > 2) phia += SQR(2*PI/x3size);
      phia = -den*four_pi_G/phia;
    } else if (iprob == 2){
      Real M = pin->GetOrAddReal("problem","M",1.0);
      Real a0 = pin->GetOrAddReal("problem","a0",1.0);
      Real da = 3.0*M/(4*PI*a0*a0*a0);
      den = da*std::pow((1.0+r2/SQR(a0)),-2.5);
      phia = - gconst*M/sqrt(r2+SQR(a0));
    } else if (iprob == 3){
      Real a0 = pin->GetOrAddReal("problem","a0",1.0);
      den = (4.0*SQR(a0)*r2-6.0*a0)*exp(-a0*r2);
      phia = four_pi_G*exp(-a0*r2);
    }

    if(nlim > 0){
      phydro->u(IDN,k,j,i) = den;
      phydro->u(IM1,k,j,i) = 0.0;
      phydro->u(IM2,k,j,i) = 0.0;
      phydro->u(IM3,k,j,i) = 0.0;
    } else {
      phydro->u(IDN,k,j,i) = den;
      phydro->u(IM1,k,j,i) = den;
      phydro->u(IM2,k,j,i) = phia;
      phydro->u(IM3,k,j,i) = 0.0;
    }
  }}}
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void MeshBlock::UserWorkInLoop(void)
{
  // do nothing
  return;
}
//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  Hydro *phydro = pblock->phydro;
  Coordinates *pcoord = pblock->pcoord;
  Real x0=0.0, y0=0.0, z0=0.0;
  int is=pblock->is, ie=pblock->ie;
  int js=pblock->js, je=pblock->je;
  int ks=pblock->ks, ke=pblock->ke;
  int cnt = (ke-ks+1)*(je-js+1)*(ie-is+1);

  int nlim = pin->GetInteger("time","nlim");

  if(SELF_GRAVITY_ENABLED){
    Gravity *pgrav = pblock->pgrav;
    if(nlim == 0){
      Real err1=0.0,err2=0.0;
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        err1 += std::abs(pgrav->phi(k,j,i) - phydro->u(IM2,k,j,i));
        err2 += pgrav->phi(k,j,i)/phydro->u(IM2,k,j,i);
      }}} // for-loop
 
      err1 = err1/cnt;
      err2 = err2/cnt;
 
      if(Globals::my_rank == 0){
        std::cout << "=====================================================" << std::endl;
        std::cout << "L1 : " << err1 <<" L2: " << err2 << std::endl;
        std::cout << "=====================================================" << std::endl;
      }
// timing measure after loop
      int ncycle = pin->GetOrAddInteger("problem","ncycle",100);
      if(Globals::my_rank == 0){
        std::cout << "=====================================================" << std::endl;
        std::cout << "Call Poisson Solver  " << ncycle << " times          " << std::endl;
        std::cout << "=====================================================" << std::endl;
      }
      clock_t tstart = clock();
#ifdef OPENMP_PARALLEL
      double omp_start_time = omp_get_wtime();
#endif
      for (int n=0; n <= ncycle; n++) pgrav->Solver(phydro->u);
#ifdef OPENMP_PARALLEL
      double omp_time = omp_get_wtime() - omp_start_time;;
#endif
      clock_t tstop = clock();
      float cpu_time = (tstop>tstart ? (float)(tstop-tstart) : 1.0)/(float)CLOCKS_PER_SEC;
      int64_t zones = GetTotalCells();
      int64_t mb_zones = GetTotalCells()/nbtotal;
      float zc_cpus = (float)(mb_zones*ncycle)/cpu_time;
      float zc_cpus2 = (float)(mb_zones*log2(mb_zones)*ncycle)/cpu_time;

      if(Globals::my_rank == 0){
        std::cout << "Timing Possison Solver                               " << std::endl;
        std::cout << "number of zones in Mesh = " << zones << std::endl;
        std::cout << "Mesh configuration = " << mesh_size.nx1 << "x" 
                  << mesh_size.nx2 << "x" << mesh_size.nx3 << std::endl;
        std::cout << "number of zones in MeshBlock = " << mb_zones << std::endl;
        std::cout << "MeshBlock configuration = " << pblock->block_size.nx1 << "x" 
                  << pblock->block_size.nx2 << "x" << pblock->block_size.nx3 << std::endl;
        std::cout << "number of processors = " << Globals::nranks << std::endl;
        std::cout << "processor configuration = " 
                  << nrbx1 << "x" << nrbx2 << "x" << nrbx3 << std::endl;
        std::cout << "cpu time used  = " << cpu_time << std::endl;
        std::cout << "cpu time used/cycle  = " << cpu_time/(float)(ncycle) << std::endl;
        std::cout << "zone-cycles/cpu_second = " << zc_cpus << std::endl;
        std::cout << "zone-cycles(NlogN)/cpu_second = " << zc_cpus2 << std::endl;
#ifdef OPENMP_PARALLEL
        std::cout << "=====================================================" << std::endl;
        float zc_omps = (float)(zones*ncycle)/omp_time;
        std::cout << "omp number of threads = " << GetNumMeshThreads() << std::endl;
        std::cout << "omp wtime used = " << omp_time << std::endl;
        std::cout << "zone-cycles/omp_wsecond = " << zc_omps << std::endl;
#endif
        std::cout << "=====================================================" << std::endl;
      }
    }
  }

  return;
}
