//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file trackers.cpp
//  \brief implementation of functions in the Tracker classes

#include <cstdio>
#include <stdexcept>
#include <sstream>
#include <iostream>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#include "trackers.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "z4c.hpp"
#include "z4c_macro.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../lagrange_interp.hpp"

Tracker::Tracker(Mesh * pmesh, ParameterInput * pin):
    pmesh(pmesh), pofile(NULL) {
  int punct_idx;
  std::cout<<"Tracker is being constructed\n";
  ofname = pin->GetOrAddString("z4c", "tracker_filename", "punctures_position.txt");
  npunct = pin->GetOrAddInteger("z4c", "n_punct", 2);

  Initialize(pmesh, pin);
  std::cout<<pos_body[0].pos[0]<<"\n";
//  root = pin->GetOrAddInteger("wave", "mpi_root", 0);
  root = 0;
#ifdef MPI_PARALLEL
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ioproc = (rank == 0);
#else
  ioproc = true;
#endif
  if (ioproc) {
    pofile = fopen(ofname.c_str(), "w");
    if (NULL == pofile) {
      std::stringstream msg;
      msg << "### FATAL ERROR in Tracker constructor" << std::endl;
      msg << "Could not open file '" << ofname << "' for writing!";
      throw std::runtime_error(msg.str().c_str());
    }
    fprintf(pofile, "#%-12s%-13s", "1:iter", "2:time");
    for (int i_file = 1; i_file <= npunct; ++i_file) {
      punct_idx = 3*(i_file);   
      fprintf(pofile, "%d:P-x-%-7d%d:P-y-%-7d%d:P-z-%-7d", punct_idx, i_file, punct_idx+1, i_file, punct_idx+2, i_file);
    }
    fprintf(pofile, "\n");
  }
}

void Tracker::Initialize(Mesh * pmesh, ParameterInput * pin){
    
  // Initialize position for twopunctures, TODO generalize
  Real par_b = pin->GetOrAddReal("problem", "par_b", 2.);
  Real of1 = pin->GetOrAddReal("problem", "center_offset1", 0.);
  Real of2 = pin->GetOrAddReal("problem", "center_offset2", 0.);
  Real of3 = pin->GetOrAddReal("problem", "center_offset3", 0.);

  //TODO: to be generalized for many punctures, this is written for twopunctures
  for (int i_punc = 0; i_punc < npunct; ++i_punc) {    
    pos_body[i_punc].pos[0] = pow(-1, i_punc)*par_b + of1;
    pos_body[i_punc].pos[1] = of2;
    pos_body[i_punc].pos[2] = of3;

    pos_body[i_punc].betap[0] = 0.;
    pos_body[i_punc].betap[1] = 0.;
    pos_body[i_punc].betap[2] = 0.;
  }
}

void Tracker::ReduceTracker() {
  std::cout<<"In EvolveTracker\n";
  //Reduce to get betap
  MeshBlock const * pmb = pmesh->pblock;
  Tracker * ptracker = pmesh->pz4c_tracker;
  //Real times_in_block[npunct];
  Real total_betap;

  // Initialize counter and collective betap
  for (int i_punc = 0; i_punc < npunct; ++i_punc) {
    times_in_block[i_punc] = 0;
    for (int i_dim = 0; i_dim < NDIM; ++i_dim) {
      pos_body[i_punc].betap[i_dim] = 0.;
    }
  }
  
  // Loop over all meshblocks to sum all the contributions to betap
  while (pmb != NULL) {
    for (int i_punc = 0; i_punc < npunct; ++i_punc) {
      //total_betap = 0;
      if (pmb->pz4c_tracker_loc->betap[i_punc].inblock) {  
        for (int i_dim = 0; i_dim < NDIM; ++i_dim) {
          pos_body[i_punc].betap[i_dim] += pmb->pz4c_tracker_loc->betap[i_punc].betap[i_dim];
          //total_betap += pos_body[i_punc].betap[i_dim];
        }
        times_in_block[i_punc] += 1.;
      }
    }
    pmb = pmb->next;
  }
  
  //Collect contributions from all MPI ranks  
#ifdef MPI_PARALLEL
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout<<"\n<================================>\nI am rank "<<rank<<"\n<=============================>\n\n";
  for (int i_punc = 0; i_punc < npunct; ++i_punc) {
    for (int i_dim = 0; i_dim < NDIM; ++i_dim) {
      if (root == rank) {
	if (i_dim == 0) MPI_Reduce(MPI_IN_PLACE, &times_in_block[i_punc], 1, MPI_ATHENA_REAL, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, &pos_body[i_punc].betap[i_dim], 1, MPI_ATHENA_REAL, MPI_SUM, root, MPI_COMM_WORLD);
      }
      else {
	if (i_dim == 0) MPI_Reduce(&times_in_block[i_punc], &times_in_block[i_punc], 1, MPI_ATHENA_REAL, MPI_SUM, root, MPI_COMM_WORLD);
        MPI_Reduce(&pos_body[i_punc].betap[i_dim], &pos_body[i_punc].betap[i_dim], 1, MPI_ATHENA_REAL, MPI_SUM, root, MPI_COMM_WORLD);
      }
	
      
    }
  }
#endif
}

void Tracker::EvolveTracker() 
{
  if (ioproc) {
    for (int i_punc = 0; i_punc < npunct; ++i_punc) {
      for (int i_dim = 0; i_dim < NDIM; ++i_dim) {
        times_in_block[i_punc] = (times_in_block[i_punc]==0) ? 1 : times_in_block[i_punc]; 
        pos_body[i_punc].betap[i_dim] /= times_in_block[i_punc]; 
      }
      std::cout<<'\n'<<"===== puncture "<<i_punc<<" ====="<<'\n';
      std::cout<<"I was in a block "<<times_in_block[i_punc]<<" times.\n";
      std::cout<<"beta0 = "<<pos_body[i_punc].betap[0]<<'\n';
      std::cout<<"====================\n";
    }
  EvolveTrackerIntegrateEuler();
  }
}

//----------------------------------------------------------------------------------------
// \!fn void Z4c::Z4cEvolveTrackerIntegrate(AthenaArray<Real> & u_pos, AthenaArray<Real> & u)
// \brief compute the RHS given the state vector and matter state
//
// This function implements puncture tracker with Euler Timestep
void Tracker::EvolveTrackerIntegrateEuler()
{
  std::cout<<"In EvolveTrackerIntegrate Euler\n";
  
  Real dt = pmesh->dt;
  Tracker * ptracker = pmesh->pz4c_tracker;
  for(int i_punct = 0; i_punct < NPUNCT; ++i_punct) {
    for(int i_dim = 0; i_dim < NDIM; ++i_dim) {
        ptracker->pos_body[i_punct].pos[i_dim] += dt * (-ptracker->pos_body[i_punct].betap[i_dim]);
    }
    std::cout<<"Pos 0, punct "<<i_punct<<'\n'<<"Val: ";
    std::cout<<ptracker->pos_body[i_punct].pos[0];
    std::cout<<"; Betap: "<<pos_body[i_punct].betap[0]<<'\n';
  }

}

void Tracker::WriteTracker(int iter, Real time) const {
  if (ioproc) {
    fprintf(pofile, "%-13d%-13.5e", iter, time);
    for (int i_file = 0; i_file < npunct; ++i_file) {
      fprintf(pofile, "%-13.5e%-13.5e%-13.5e", pos_body[i_file].pos[0], pos_body[i_file].pos[1], pos_body[i_file].pos[2]);
    }
    fprintf(pofile, "\n");
  }
}

Tracker::~Tracker() {
  if (ioproc) {
     fclose(pofile);
  }
  std::cout<<"Tracker is being deleted\n";
}


//TrackerLocal class 

TrackerLocal::TrackerLocal(MeshBlock * pmb, ParameterInput * pin) {
  pmy_block = pmb;
  Coordinates * pco = pmb->pcoord;
  std::cout<<"Tracker Local is being constructed\n";
}

//----------------------------------------------------------------------------------------
// \!fn Z4c::StoreBetaPrev(AthenaArray<Real> & u, AthenaArray<Real> & u_pos) 
// \brief Store Beta at previous timestep, needed for integration.
//
void TrackerLocal::StoreBetaPrev(Betap_vars betap[2], AthenaArray<Real> & u, int body)
{
  Real dt;
  Real origin[3];
  Real delta[3];
  int size[3];
  Real coord[3];
  AthenaArray<Real> mySrc;
  int const nvars = u.GetDim4();
  LagrangeInterpND<2*NGHOST-1, 3> * pinterp = nullptr;
   
  MeshBlock * pmb = pmy_block;
  Tracker * ptracker = pmb->pmy_mesh->pz4c_tracker;
  Coordinates const * pmc = pmy_block->pcoord;
  
  std::cout<<"In StoreBetaPrev\n";
  if (!((pmy_block->block_size.x1min <= ptracker->pos_body[body].pos[0]) && (ptracker->pos_body[body].pos[0] <= pmy_block->block_size.x1max) &&
      (pmy_block->block_size.x2min <= ptracker->pos_body[body].pos[1]) && (ptracker->pos_body[body].pos[1] <= pmy_block->block_size.x2max) &&
      (pmy_block->block_size.x3min <= ptracker->pos_body[body].pos[2]) && (ptracker->pos_body[body].pos[2] <= pmy_block->block_size.x3max))) {   
    for (int i_dim = 0; i_dim < NDIM; ++i_dim) {
      betap[body].betap[i_dim] = std::nan("1");
      betap[body].inblock = 0;
    }
  }
  else {
    betap[body].inblock = 1;
//#if PREFER_VC
//    origin[0] = pmc->x1f(0);
//    origin[1] = pmc->x2f(0);
//    origin[2] = pmc->x3f(0);
//    size[0] = pmb->block_size.nx1 + 2*(NGHOST) + 1;
//    size[1] = pmb->block_size.nx2 + 2*(NGHOST) + 1;
//    size[2] = pmb->block_size.nx3 + 2*(NGHOST) + 1;
//    delta[0] = pmc->dx1f(0);
//    delta[1] = pmc->dx2f(0);
//    delta[2] = pmc->dx3f(0);
//#else
//    origin[0] = pmc->x1v(0);
//    origin[1] = pmc->x2v(0);
//    origin[2] = pmc->x3v(0);
//    size[0] = pmb->block_size.nx1 + 2*(NGHOST);
//    size[1] = pmb->block_size.nx2 + 2*(NGHOST);
//    size[2] = pmb->block_size.nx3 + 2*(NGHOST);
//    delta[0] = pmc->dx1v(0);
//    delta[1] = pmc->dx2v(0);
//    delta[2] = pmc->dx3v(0);
//#endif
    origin[0] = pmb->pz4c->mbi.x1(0);
    origin[1] = pmb->pz4c->mbi.x2(0);
    origin[2] = pmb->pz4c->mbi.x3(0);
    size[0] = pmb->pz4c->mbi.nn1 + 2*(NGHOST);
    size[1] = pmb->pz4c->mbi.nn2 + 2*(NGHOST);
    size[2] = pmb->pz4c->mbi.nn3 + 2*(NGHOST);
#if PREFER_VC
    delta[0] = pmc->dx1f(0);
    delta[1] = pmc->dx2f(0);
    delta[2] = pmc->dx3f(0);
#else
    delta[0] = pmc->dx1v(0);
    delta[1] = pmc->dx2v(0);
    delta[2] = pmc->dx3v(0);
#endif
    
    for (int i_dim = 0; i_dim < NDIM; ++i_dim) {
      coord[i_dim] = ptracker->pos_body[body].pos[i_dim];
    }

    std::cout<<"\nPunc position: x ="<<coord[0]<<", y = "<<coord[1]<<", z = "<<coord[2]<<'\n';
    std::cout<<"Edges:\n"<<pmy_block->block_size.x1min<<"<=x<="<<pmy_block->block_size.x1max<<'\n';
    std::cout<<pmy_block->block_size.x2min<<"<=y<="<<pmy_block->block_size.x2max<<'\n';
    std::cout<<pmy_block->block_size.x3min<<"<=z<="<<pmy_block->block_size.x3max<<'\n';
    pinterp = new LagrangeInterpND<2*NGHOST-1, 3>(origin, delta, size, coord);
    
    for (int iv = 19; iv < nvars; ++iv) {
      mySrc.InitWithShallowSlice(const_cast<AthenaArray<Real>&>(u), iv, 1);
      betap[body].betap[iv-19] = pinterp->eval(mySrc.data());
    }


    delete pinterp;
  }

}
                                                                                      
TrackerLocal::~TrackerLocal() {
  std::cout<<"Tracker Local is being deleted\n";
}
