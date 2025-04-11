//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file amr_loadbalance.cpp
//! \brief implementation of Mesh::AdaptiveMeshRefinement() and related utilities

// C headers

// C++ headers
#include <algorithm>  // std::sort()
#include <cstdint>
#include <iostream>
#include <limits>
#include <sstream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../utils/buffer_utils.hpp"
#include "mesh.hpp"
#include "mesh_refinement.hpp"
#include "meshblock_tree.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


//----------------------------------------------------------------------------------------
//! \fn void Mesh::LoadBalancingAndAdaptiveMeshRefinement(ParameterInput *pin)
//! \brief Main function for adaptive mesh refinement

void Mesh::LoadBalancingAndAdaptiveMeshRefinement(ParameterInput *pin) {
  int nnew = 0, ndel = 0;
  amr_updated = false;

  if (adaptive) {
    UpdateMeshBlockTree(nnew, ndel);
    nbnew += nnew; nbdel += ndel;
  }

  lb_flag_ |= lb_automatic_;

  UpdateCostList();

  if (nnew != 0 || ndel != 0) { // at least one (de)refinement happened
    amr_updated = true;
    GatherCostListAndCheckBalance();
    RedistributeAndRefineMeshBlocks(pin, nbtotal + nnew - ndel);
  } else if (lb_flag_ && step_since_lb >= lb_interval_) {
    if (!GatherCostListAndCheckBalance()) { // load imbalance detected
      amr_updated = true;
      RedistributeAndRefineMeshBlocks(pin, nbtotal);
    }
    lb_flag_ = false;
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Mesh::CalculateLoadBalance(double *clist, int *rlist, int *slist,
//!                                     int *nlist, int nb)
//! \brief Calculate distribution of MeshBlocks based on the cost list

void Mesh::CalculateLoadBalance(double *clist, int *rlist, int *slist, int *nlist,
                                int nb) {
  std::stringstream msg;
  double real_max  =  std::numeric_limits<double>::max();
  double totalcost = 0, maxcost = 0.0, mincost = (real_max);

  for (int i=0; i<nb; i++) {
    totalcost += clist[i];
    mincost = std::min(mincost,clist[i]);
    maxcost = std::max(maxcost,clist[i]);
  }

  int j = (Globals::nranks) - 1;
  double targetcost = totalcost/Globals::nranks;
  double mycost = 0.0;
  // create rank list from the end: the master MPI rank should have less load
  for (int i=nb-1; i>=0; i--) {
    if (targetcost == 0.0) {
      msg << "### FATAL ERROR in CalculateLoadBalance" << std::endl
          << "There is at least one process which has no MeshBlock" << std::endl
          << "Decrease the number of processes or use smaller MeshBlocks." << std::endl;
      ATHENA_ERROR(msg);
    }
    mycost += clist[i];
    rlist[i] = j;
    if (mycost >= targetcost && j>0) {
      j--;
      totalcost -= mycost;
      mycost = 0.0;
      targetcost = totalcost/(j+1);
    }
  }
  slist[0] = 0;
  j = 0;
  for (int i=1; i<nb; i++) { // make the list of nbstart and nblocks
    if (rlist[i] != rlist[i-1]) {
      nlist[j] = i-slist[j];
      slist[++j] = i;
    }
  }
  nlist[j] = nb-slist[j];

#ifdef MPI_PARALLEL
  if (nb % (Globals::nranks * num_mesh_threads_) != 0
      && !adaptive && !lb_flag_ && maxcost == mincost && Globals::my_rank == 0) {
    std::cout << "### Warning in CalculateLoadBalance" << std::endl
              << "The number of MeshBlocks cannot be divided evenly. "
              << "This will result in poor load balancing." << std::endl;
  }
#endif
  if ((Globals::nranks)*(num_mesh_threads_) > nb) {
    msg << "### FATAL ERROR in CalculateLoadBalance" << std::endl
        << "There are fewer MeshBlocks than OpenMP threads on each MPI rank" << std::endl
        << "Decrease the number of threads or use more MeshBlocks." << std::endl;
    ATHENA_ERROR(msg);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::ResetLoadBalanceVariables()
//! \brief reset counters and flags for load balancing

void Mesh::ResetLoadBalanceVariables() {
  if (lb_automatic_) {
    for (int i=0; i<nblocal; ++i) {
      MeshBlock *pmb = my_blocks(i);
      costlist[pmb->gid] = TINY_NUMBER;
      pmb->ResetTimeMeasurement();
    }
  }
  lb_flag_ = false;
  step_since_lb = 0;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::UpdateCostList()
//! \brief update the cost list

void Mesh::UpdateCostList() {
  if (lb_automatic_) {
    double w = static_cast<double>(lb_interval_-1)/static_cast<double>(lb_interval_);
    for (int i=0; i<nblocal; ++i) {
      MeshBlock *pmb = my_blocks(i);
      costlist[pmb->gid] = costlist[pmb->gid]*w+pmb->cost_;
    }
  } else if (lb_flag_) {
    for (int i=0; i<nblocal; ++i) {
      MeshBlock *pmb = my_blocks(i);
      costlist[pmb->gid] = pmb->cost_;
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::UpdateMeshBlockTree(int &nnew, int &ndel)
//! \brief collect refinement flags and manipulate the MeshBlockTree

void Mesh::UpdateMeshBlockTree(int &nnew, int &ndel) {
  // compute nleaf= number of leaf MeshBlocks per refined block
  int nleaf = 2;
  if (mesh_size.nx2 > 1) nleaf = 4;
  if (mesh_size.nx3 > 1) nleaf = 8;

  // collect refinement flags from all the meshblocks
  // count the number of the blocks to be (de)refined
  nref[Globals::my_rank] = 0;
  nderef[Globals::my_rank] = 0;
  for (int i=0; i<nblocal; ++i) {
    MeshBlock *pmb = my_blocks(i);
    if (pmb->pmr->refine_flag_ ==  1) nref[Globals::my_rank]++;
    if (pmb->pmr->refine_flag_ == -1) nderef[Globals::my_rank]++;
  }
#ifdef MPI_PARALLEL
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, nref,   1, MPI_INT, MPI_COMM_WORLD);
  MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, nderef, 1, MPI_INT, MPI_COMM_WORLD);
#endif

  // count the number of the blocks to be (de)refined and displacement
  int tnref = 0, tnderef = 0;
  for (int n=0; n<Globals::nranks; n++) {
    tnref  += nref[n];
    tnderef += nderef[n];
  }
  if (tnref == 0 && tnderef < nleaf) // nothing to do
    return;

  int rd = 0, dd = 0;
  for (int n=0; n<Globals::nranks; n++) {
    rdisp[n] = rd;
    ddisp[n] = dd;
    // technically could overflow, since sizeof() operator returns
    // std::size_t = long unsigned int > int
    // on many platforms (LP64). However, these are used below in MPI calls for
    // integer arguments (recvcounts, displs). MPI does not support > 64-bit count ranges
    bnref[n] = static_cast<int>(nref[n]*sizeof(LogicalLocation));
    bnderef[n] = static_cast<int>(nderef[n]*sizeof(LogicalLocation));
    brdisp[n] = static_cast<int>(rd*sizeof(LogicalLocation));
    bddisp[n] = static_cast<int>(dd*sizeof(LogicalLocation));
    rd += nref[n];
    dd += nderef[n];
  }

  // allocate memory for the location arrays
  LogicalLocation *lref{}, *lderef{}, *clderef{};
  if (tnref > 0)
    lref = new LogicalLocation[tnref];
  if (tnderef >= nleaf) {
    lderef = new LogicalLocation[tnderef];
    clderef = new LogicalLocation[tnderef/nleaf];
  }

  // collect the locations and costs
  int iref = rdisp[Globals::my_rank], ideref = ddisp[Globals::my_rank];
  for (int i=0; i<nblocal; ++i) {
    MeshBlock *pmb = my_blocks(i);
    if (pmb->pmr->refine_flag_ ==  1)
      lref[iref++] = pmb->loc;
    if (pmb->pmr->refine_flag_ == -1 && tnderef >= nleaf)
      lderef[ideref++] = pmb->loc;
  }
#ifdef MPI_PARALLEL
  if (tnref > 0) {
    MPI_Allgatherv(MPI_IN_PLACE, bnref[Globals::my_rank],   MPI_BYTE,
                   lref,   bnref,   brdisp, MPI_BYTE, MPI_COMM_WORLD);
  }
  if (tnderef >= nleaf) {
    MPI_Allgatherv(MPI_IN_PLACE, bnderef[Globals::my_rank], MPI_BYTE,
                   lderef, bnderef, bddisp, MPI_BYTE, MPI_COMM_WORLD);
  }
#endif

  // calculate the list of the newly derefined blocks
  int ctnd = 0;
  if (tnderef >= nleaf) {
    int lk = 0, lj = 0;
    if (mesh_size.nx2 > 1) lj = 1;
    if (mesh_size.nx3 > 1) lk = 1;
    for (int n=0; n<tnderef; n++) {
      if ((lderef[n].lx1 & 1LL) == 0LL &&
          (lderef[n].lx2 & 1LL) == 0LL &&
          (lderef[n].lx3 & 1LL) == 0LL) {
        int r = n, rr = 0;
        for (std::int64_t k=0; k<=lk; k++) {
          for (std::int64_t j=0; j<=lj; j++) {
            for (std::int64_t i=0; i<=1; i++) {
              if (r < tnderef) {
                if ((lderef[n].lx1+i) == lderef[r].lx1
                    && (lderef[n].lx2+j) == lderef[r].lx2
                    && (lderef[n].lx3+k) == lderef[r].lx3
                    &&  lderef[n].level  == lderef[r].level)
                  rr++;
                r++;
              }
            }
          }
        }
        if (rr == nleaf) {
          clderef[ctnd].lx1   = lderef[n].lx1>>1;
          clderef[ctnd].lx2   = lderef[n].lx2>>1;
          clderef[ctnd].lx3   = lderef[n].lx3>>1;
          clderef[ctnd].level = lderef[n].level-1;
          ctnd++;
        }
      }
    }
  }
  // sort the lists by level
  if (ctnd > 1)
    std::sort(clderef, &(clderef[ctnd-1]), LogicalLocation::Greater);

  if (tnderef >= nleaf)
    delete [] lderef;

  // Now the lists of the blocks to be refined and derefined are completed
  // Start tree manipulation
  // Step 1. perform refinement
  for (int n=0; n<tnref; n++) {
    MeshBlockTree *bt=tree.FindMeshBlock(lref[n]);
    bt->Refine(nnew);
  }
  if (tnref != 0)
    delete [] lref;

  // Step 2. perform derefinement
  for (int n=0; n<ctnd; n++) {
    MeshBlockTree *bt = tree.FindMeshBlock(clderef[n]);
    bt->Derefine(ndel);
  }
  if (tnderef >= nleaf)
    delete [] clderef;

  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool Mesh::GatherCostListAndCheckBalance()
//! \brief collect the cost from MeshBlocks and check the load balance

bool Mesh::GatherCostListAndCheckBalance() {
  if (lb_manual_ || lb_automatic_) {
#ifdef MPI_PARALLEL
    MPI_Allgatherv(MPI_IN_PLACE, nblist[Globals::my_rank], MPI_DOUBLE, costlist, nblist,
                   nslist, MPI_DOUBLE, MPI_COMM_WORLD);
#endif
    double maxcost = 0.0, avecost = 0.0;
    for (int rank=0; rank<Globals::nranks; rank++) {
      double rcost = 0.0;
      int ns = nslist[rank];
      int ne = ns + nblist[rank];
      for (int n=ns; n<ne; ++n)
        rcost += costlist[n];
      maxcost = std::max(maxcost,rcost);
      avecost += rcost;
    }
    avecost /= Globals::nranks;

    if (adaptive) lb_tolerance_ = 2.0*static_cast<double>(Globals::nranks)
                                     /static_cast<double>(nbtotal);

    if (maxcost > (1.0 + lb_tolerance_)*avecost)
      return false;
  }
  return true;
}


//----------------------------------------------------------------------------------------
//! \fn void Mesh::RedistributeAndRefineMeshBlocks(ParameterInput *pin, int ntot)
//! \brief redistribute MeshBlocks according to the new load balance

void Mesh::RedistributeAndRefineMeshBlocks(ParameterInput *pin, int ntot) {
  // compute nleaf= number of leaf MeshBlocks per refined block
  int nleaf = 2;
  if (mesh_size.nx2 > 1) nleaf = 4;
  if (mesh_size.nx3 > 1) nleaf = 8;

  // Step 1. construct new lists
  LogicalLocation *newloc = new LogicalLocation[ntot];
  int *newrank = new int[ntot];
  double *newcost = new double[ntot];
  int *newtoold = new int[ntot];
  int *oldtonew = new int[nbtotal];
  int nbtold = nbtotal;
  tree.GetMeshBlockList(newloc, newtoold, nbtotal);

  // create a list mapping the previous gid to the current one
  oldtonew[0] = 0;
  int mb_idx = 1;
  for (int n=1; n<ntot; n++) {
    if (newtoold[n] == newtoold[n-1] + 1) { // normal
      oldtonew[mb_idx++] = n;
    } else if (newtoold[n] == newtoold[n-1] + nleaf) { // derefined
      for (int j=0; j<nleaf-1; j++)
        oldtonew[mb_idx++] = n-1;
      oldtonew[mb_idx++] = n;
    }
  }
  // fill the last block
  for ( ; mb_idx<nbtold; mb_idx++)
    oldtonew[mb_idx] = ntot-1;

  current_level = 0;
  for (int n=0; n<ntot; n++) {
    // "on" = "old n" = "old gid" = "old global MeshBlock ID"
    int on = newtoold[n];
    if (newloc[n].level>current_level) // set the current max level
      current_level = newloc[n].level;
    if (newloc[n].level > loclist[on].level) // refined
      refined_.push_back(n);
    if (newloc[n].level >= loclist[on].level) { // same or refined
      newcost[n] = costlist[on];
    } else {
      double acost = 0.0;
      for (int l=0; l<nleaf; l++)
        acost += costlist[on+l];
      newcost[n] = acost/nleaf;
    }
  }

  // Step 2. Calculate new load balance
  CalculateLoadBalance(newcost, newrank, nslist, nblist, ntot);

  int nbs = nslist[Globals::my_rank];
  int nbe = nbs + nblist[Globals::my_rank] - 1;

#ifdef MPI_PARALLEL
  // Step 3. count the number of the blocks to be sent / received
  int nsend = 0, nrecv = 0;
  for (int n=nbs; n<=nbe; n++) {
    int on = newtoold[n];
    if (loclist[on].level > newloc[n].level) { // f2c
      for (int k=0; k<nleaf; k++) {
        if (ranklist[on+k] != Globals::my_rank)
          nrecv++;
      }
    } else {
      if (ranklist[on] != Globals::my_rank)
        nrecv++;
    }
  }
  for (int n=gids_; n<=gide_; n++) {
    int nn = oldtonew[n];
    if (loclist[n].level < newloc[nn].level) { // c2f
      for (int k=0; k<nleaf; k++) {
        if (newrank[nn+k] != Globals::my_rank)
          nsend++;
      }
    } else {
      if (newrank[nn] != Globals::my_rank)
        nsend++;
    }
  }

  // Step 4. allocate and start receiving buffers
  Real **sendbuf, **recvbuf;
  MPI_Request *req_send, *req_recv;
  if (nrecv != 0) {
    recvbuf = new Real*[nrecv];
    req_recv = new MPI_Request[nrecv];
    int rb_idx = 0;     // recv buffer index
    for (int n=nbs; n<=nbe; n++) {
      int on = newtoold[n];
      LogicalLocation const &oloc = loclist[on];
      LogicalLocation const &nloc = newloc[n];
      if (oloc.level > nloc.level) { // f2c
        for (int l=0; l<nleaf; l++) {
          if (ranklist[on+l] == Globals::my_rank) continue;
          LogicalLocation &lloc = loclist[on+l];
          int ox1 = ((lloc.lx1 & 1LL) == 1LL), ox2 = ((lloc.lx2 & 1LL) == 1LL),
              ox3 = ((lloc.lx3 & 1LL) == 1LL);
          recvbuf[rb_idx] = new Real[bsf2c];
          int tag = CreateAMRMPITag(n-nbs, ox1, ox2, ox3);
          MPI_Irecv(recvbuf[rb_idx], bsf2c, MPI_ATHENA_REAL, ranklist[on+l],
                    tag, MPI_COMM_WORLD, &(req_recv[rb_idx]));
          rb_idx++;
        }
      } else { // same level or c2f
        if (ranklist[on] == Globals::my_rank) continue;
        int size;
        if (oloc.level == nloc.level) {
          size = bssame;
        } else {
          size = bsc2f;
        }
        recvbuf[rb_idx] = new Real[size];
        int tag = CreateAMRMPITag(n-nbs, 0, 0, 0);
        MPI_Irecv(recvbuf[rb_idx], size, MPI_ATHENA_REAL, ranklist[on],
                  tag, MPI_COMM_WORLD, &(req_recv[rb_idx]));
        rb_idx++;
      }
    }
  }
  // Step 5. allocate, pack and start sending buffers
  if (nsend != 0) {
    sendbuf = new Real*[nsend];
    req_send = new MPI_Request[nsend];
    int sb_idx = 0;      // send buffer index
    for (int n=gids_; n<=gide_; n++) {
      int nn = oldtonew[n];
      LogicalLocation &oloc = loclist[n];
      LogicalLocation const &nloc = newloc[nn];
      MeshBlock* pb = FindMeshBlock(n);
      if (nloc.level == oloc.level) { // same level
        if (newrank[nn] == Globals::my_rank) continue;
        sendbuf[sb_idx] = new Real[bssame];
        PrepareSendSameLevel(pb, sendbuf[sb_idx]);
        int tag = CreateAMRMPITag(nn-nslist[newrank[nn]], 0, 0, 0);
        MPI_Isend(sendbuf[sb_idx], bssame, MPI_ATHENA_REAL, newrank[nn],
                  tag, MPI_COMM_WORLD, &(req_send[sb_idx]));
        sb_idx++;
      } else if (nloc.level > oloc.level) { // c2f
        // c2f must communicate to multiple leaf blocks (unlike f2c, same2same)
        for (int l=0; l<nleaf; l++) {
          if (newrank[nn+l] == Globals::my_rank) continue;
          sendbuf[sb_idx] = new Real[bsc2f];
          PrepareSendCoarseToFineAMR(pb, sendbuf[sb_idx], newloc[nn+l]);
          int tag = CreateAMRMPITag(nn+l-nslist[newrank[nn+l]], 0, 0, 0);
          MPI_Isend(sendbuf[sb_idx], bsc2f, MPI_ATHENA_REAL, newrank[nn+l],
                    tag, MPI_COMM_WORLD, &(req_send[sb_idx]));
          sb_idx++;
        } // end loop over nleaf (unique to c2f branch in this step 6)
      } else { // f2c: restrict + pack + send
        if (newrank[nn] == Globals::my_rank) continue;
        sendbuf[sb_idx] = new Real[bsf2c];
        PrepareSendFineToCoarseAMR(pb, sendbuf[sb_idx]);
        int ox1 = ((oloc.lx1 & 1LL) == 1LL), ox2 = ((oloc.lx2 & 1LL) == 1LL),
            ox3 = ((oloc.lx3 & 1LL) == 1LL);
        int tag = CreateAMRMPITag(nn-nslist[newrank[nn]], ox1, ox2, ox3);
        MPI_Isend(sendbuf[sb_idx], bsf2c, MPI_ATHENA_REAL, newrank[nn],
                  tag, MPI_COMM_WORLD, &(req_send[sb_idx]));
        sb_idx++;
      }
    }
  } // if (nsend !=0)
#endif // MPI_PARALLEL

  if (my_blocks(0)->pmr->pvars_fc_.size() > 0)
    PrepareAndSendFaceFieldCorrection(newloc, ranklist, newrank, nslist, nbtold);

  // Step 6. construct a new MeshBlock list (moving the data within the MPI rank)
  AthenaArray<MeshBlock*> newlist;
  newlist.NewAthenaArray(nblist[Globals::my_rank]);
  RegionSize block_size = my_blocks(0)->block_size;

  for (int n=nbs; n<=nbe; n++) {
    int on = newtoold[n];
    if ((ranklist[on] == Globals::my_rank) && (loclist[on].level == newloc[n].level)) {
      // on the same MPI rank and same level -> just move it
      MeshBlock* pob = FindMeshBlock(on);
      pob->gid = n;
      pob->lid = n - nbs;
      newlist(n-nbs) = pob;
      my_blocks(on-gids_) = nullptr;
    } else {
      // on a different refinement level or MPI rank - create a new block
      BoundaryFlag block_bcs[6];
      SetBlockSizeAndBoundaries(newloc[n], block_size, block_bcs);
      newlist(n-nbs) = new MeshBlock(n, n-nbs, newloc[n], block_size, block_bcs, this,
                                     pin, true);
      // fill the conservative variables
      if ((loclist[on].level > newloc[n].level)) { // fine to coarse (f2c)
        for (int ll=0; ll<nleaf; ll++) {
          if (ranklist[on+ll] != Globals::my_rank) continue;
          // fine to coarse on the same MPI rank (different AMR level) - restriction
          MeshBlock* pob = FindMeshBlock(on+ll);
          FillSameRankFineToCoarseAMR(pob, newlist(n-nbs), loclist[on+ll]);
        }
      } else if ((loclist[on].level < newloc[n].level) && // coarse to fine (c2f)
                 (ranklist[on] == Globals::my_rank)) {
        // coarse to fine on the same MPI rank (different AMR level) - prolongation
        MeshBlock* pob = FindMeshBlock(on);
        FillSameRankCoarseToFineAMR(pob, newlist(n-nbs), newloc[n]);
      }
    }
  }

  // discard remaining MeshBlocks
  // they could be reused, but for the moment, just throw them away for simplicity
  for (int n = 0; n<nblocal; n++) {
    delete my_blocks(n); // OK to delete even if it is nullptr
    my_blocks(n) = nullptr;
  }

  // Replace the MeshBlock list
  my_blocks.ExchangeAthenaArray(newlist);
  nblocal = nblist[Globals::my_rank];
  gids_ = nbs;
  gide_ = nbe;

  // Step 7. Receive the data and load into MeshBlocks
  // This is a test: try MPI_Waitall later.
#ifdef MPI_PARALLEL
  if (nrecv != 0) {
    int rb_idx = 0;     // recv buffer index
    for (int n=nbs; n<=nbe; n++) {
      int on = newtoold[n];
      LogicalLocation const &oloc = loclist[on];
      LogicalLocation const &nloc = newloc[n];
      MeshBlock *pb = FindMeshBlock(n);
      if (oloc.level == nloc.level) { // same
        if (ranklist[on] == Globals::my_rank) continue;
        MPI_Wait(&(req_recv[rb_idx]), MPI_STATUS_IGNORE);
        FinishRecvSameLevel(pb, recvbuf[rb_idx]);
        rb_idx++;
      } else if (oloc.level > nloc.level) { // f2c
        for (int l=0; l<nleaf; l++) {
          if (ranklist[on+l] == Globals::my_rank) continue;
          MPI_Wait(&(req_recv[rb_idx]), MPI_STATUS_IGNORE);
          FinishRecvFineToCoarseAMR(pb, recvbuf[rb_idx], loclist[on+l]);
          rb_idx++;
        }
      } else { // c2f
        if (ranklist[on] == Globals::my_rank) continue;
        MPI_Wait(&(req_recv[rb_idx]), MPI_STATUS_IGNORE);
        ReceiveCoarseToFineAMR(pb, recvbuf[rb_idx]);
        rb_idx++;
      }
    }
  }
#endif

  if (my_blocks(0)->pmr->pvars_fc_.size() > 0)
    ReceiveAndSetFaceFieldCorrection(newrank);

  // Step 7b. Prolongate MeshBlocks
  for (int n=nbs; n<=nbe; n++) {
    int on = newtoold[n];
    LogicalLocation const &oloc = loclist[on];
    LogicalLocation const &nloc = newloc[n];
    if (oloc.level < nloc.level) { // c2f
      MeshBlock *pb = FindMeshBlock(n);
      ProlongateMeshBlock(pb);
    }
  }

  // deallocate arrays
  delete [] loclist;
  delete [] ranklist;
  delete [] costlist;
  delete [] newtoold;
  delete [] oldtonew;
  if (my_blocks(0)->pmr->pvars_fc_.size() > 0) {
    refined_.clear();
    ffc_send_.clear();
    ffc_recv_.clear();
    for (int l = 0; l < max_level - root_level; ++l)
      locmap_[l].clear();
  }

#ifdef MPI_PARALLEL
  if (nsend != 0) {
    MPI_Waitall(nsend, req_send, MPI_STATUSES_IGNORE);
    for (int n=0; n<nsend; n++)
      delete [] sendbuf[n];
    delete [] sendbuf;
    delete [] req_send;
  }
  if (nrecv != 0) {
    for (int n=0; n<nrecv; n++)
      delete [] recvbuf[n];
    delete [] recvbuf;
    delete [] req_recv;
  }
#endif

  // update the lists
  loclist = newloc;
  ranklist = newrank;
  costlist = newcost;

  // re-initialize the MeshBlocks
  for (int i=0; i<nblocal; ++i)
    my_blocks(i)->pbval->SearchAndSetNeighbors(tree, ranklist, nslist);
  Initialize(2, pin);

  ResetLoadBalanceVariables();

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::PrepareSendSameLevel(MeshBlock* pb, Real *sendbuf)
//! \brief AMR: step 6, branch 1 (same2same: just pack+send)
//!
//! This helper fn is used for AMR and non-refinement load balancing of MeshBlocks.
//! Therefore, unlike PrepareSendCoarseToFineAMR(), etc., it loops over
//! MeshBlock::vars_cc/fc_ containers, not MeshRefinement::pvars_cc/fc_ containers

void Mesh::PrepareSendSameLevel(MeshBlock* pb, Real *sendbuf) {
  // pack
  int p = 0;

  //! \todo (felker):
  //! * add explicit check to ensure that elements of pb->vars_cc/fc_ and
  //!   pb->pmr->pvars_cc/fc_ v point to the same objects, if adaptive

  // (C++11) range-based for loop: (automatic type deduction fails when iterating over
  // container with std::reference_wrapper; could use auto var_cc_r = var_cc.get())
  for (AthenaArray<Real> &var_cc : pb->vars_cc_) {
    int nu = var_cc.GetDim4() - 1;
    BufferUtility::PackData(var_cc, sendbuf, 0, nu,
                            pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke, p);
  }
  if ((NR_RADIATION_ENABLED|| IM_RADIATION_ENABLED)) {
    int nu = pb->pnrrad->ir.GetDim1() - 1;
    BufferUtility::PackData(pb->pnrrad->ir, sendbuf, pb->ks, pb->ke,
                            0,  nu, pb->is, pb->ie, pb->js, pb->je, p);
  }

  for (FaceField &var_fc : pb->vars_fc_) {
    BufferUtility::PackData(var_fc.x1f, sendbuf,
                            pb->is, pb->ie+1, pb->js, pb->je, pb->ks, pb->ke, p);
    BufferUtility::PackData(var_fc.x2f, sendbuf,
                            pb->is, pb->ie, pb->js, pb->je+f2, pb->ks, pb->ke, p);
    BufferUtility::PackData(var_fc.x3f, sendbuf,
                            pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke+f3, p);
  }
  //! \warning (felker):
  //! * casting from "Real *" to "int *" in order to append single integer
  //!   to send buffer is slightly unsafe (especially if sizeof(int) > sizeof(Real))
  if (adaptive) {
    int *dcp = reinterpret_cast<int *>(&(sendbuf[p]));
    *dcp = pb->pmr->deref_count_;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::PrepareSendCoarseToFineAMR(MeshBlock* pb, Real *sendbuf,
//!                                           LogicalLocation &lloc)
//! \brief step 6, branch 2 (c2f: just pack+send)

void Mesh::PrepareSendCoarseToFineAMR(MeshBlock* pb, Real *sendbuf,
                                      LogicalLocation &lloc) {
  int ox1 = ((lloc.lx1 & 1LL) == 1LL), ox2 = ((lloc.lx2 & 1LL) == 1LL),
      ox3 = ((lloc.lx3 & 1LL) == 1LL);
  // pack
  int il, iu, jl, ju, kl, ku;
  if (ox1 == 0) il = pb->is - 1,               iu = pb->is + pb->block_size.nx1/2;
  else        il = pb->is + pb->block_size.nx1/2-1,  iu = pb->ie + 1;
  if (ox2 == 0) jl = pb->js - f2,              ju = pb->js + pb->block_size.nx2/2;
  else        jl = pb->js + pb->block_size.nx2/2 - f2, ju = pb->je + f2;
  if (ox3 == 0) kl = pb->ks - f3,              ku = pb->ks + pb->block_size.nx3/2;
  else        kl = pb->ks + pb->block_size.nx3/2 - f3, ku = pb->ke + f3;
  int p = 0;
  for (auto cc_pair : pb->pmr->pvars_cc_) {
    AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
    int nu = var_cc->GetDim4() - 1;
    BufferUtility::PackData(*var_cc, sendbuf, 0, nu,
                            il, iu, jl, ju, kl, ku, p);
  }
  if ((NR_RADIATION_ENABLED|| IM_RADIATION_ENABLED)) {
    int nu = pb->pnrrad->ir.GetDim1() - 1;
    BufferUtility::PackData(pb->pnrrad->ir, sendbuf, kl, ku,
                            0,  nu, il, iu, jl, ju, p);
  }


  for (auto fc_pair : pb->pmr->pvars_fc_) {
    FaceField *var_fc = std::get<0>(fc_pair);
    BufferUtility::PackData((*var_fc).x1f, sendbuf,
                            il, iu+1, jl, ju, kl, ku, p);
    BufferUtility::PackData((*var_fc).x2f, sendbuf,
                            il, iu, jl, ju+f2, kl, ku, p);
    BufferUtility::PackData((*var_fc).x3f, sendbuf,
                            il, iu, jl, ju, kl, ku+f3, p);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::PrepareSendFineToCoarseAMR(MeshBlock* pb, Real *sendbuf)
//! \brief step 6, branch 3 (f2c: restrict, pack, send)

void Mesh::PrepareSendFineToCoarseAMR(MeshBlock* pb, Real *sendbuf) {
  // restrict and pack
  MeshRefinement *pmr = pb->pmr;
  int p = 0;
  for (auto cc_pair : pmr->pvars_cc_) {
    AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
    AthenaArray<Real> *coarse_cc = std::get<1>(cc_pair);
    int nu = var_cc->GetDim4() - 1;
    pmr->RestrictCellCenteredValues(*var_cc, *coarse_cc,
                                    0, nu,
                                    pb->cis, pb->cie,
                                    pb->cjs, pb->cje,
                                    pb->cks, pb->cke);
    BufferUtility::PackData(*coarse_cc, sendbuf, 0, nu,
                            pb->cis, pb->cie,
                            pb->cjs, pb->cje,
                            pb->cks, pb->cke, p);
  }

  if ((NR_RADIATION_ENABLED|| IM_RADIATION_ENABLED)) {
    AthenaArray<Real> &var_cc = pb->pnrrad->ir;
    AthenaArray<Real> &coarse_cc = pb->pnrrad->coarse_ir_;
    int nu = var_cc.GetDim1() - 1;

    pmr->RestrictCellCenteredValues(var_cc, coarse_cc,-1,
                                    0, nu,
                                    pb->cis, pb->cie,
                                    pb->cjs, pb->cje,
                                    pb->cks, pb->cke);
    BufferUtility::PackData(coarse_cc, sendbuf,
                            pb->cks, pb->cke,
                            0, nu,
                            pb->cis, pb->cie,
                            pb->cjs, pb->cje, p);
  }

  for (auto fc_pair : pb->pmr->pvars_fc_) {
    FaceField *var_fc = std::get<0>(fc_pair);
    FaceField *coarse_fc = std::get<1>(fc_pair);
    pmr->RestrictFieldX1((*var_fc).x1f, (*coarse_fc).x1f,
                         pb->cis, pb->cie+1,
                         pb->cjs, pb->cje,
                         pb->cks, pb->cke);
    BufferUtility::PackData((*coarse_fc).x1f, sendbuf,
                            pb->cis, pb->cie+1,
                            pb->cjs, pb->cje,
                            pb->cks, pb->cke, p);
    pmr->RestrictFieldX2((*var_fc).x2f, (*coarse_fc).x2f,
                         pb->cis, pb->cie,
                         pb->cjs, pb->cje+f2,
                         pb->cks, pb->cke);
    BufferUtility::PackData((*coarse_fc).x2f, sendbuf,
                            pb->cis, pb->cie,
                            pb->cjs, pb->cje+f2,
                            pb->cks, pb->cke, p);
    pmr->RestrictFieldX3((*var_fc).x3f, (*coarse_fc).x3f,
                         pb->cis, pb->cie,
                         pb->cjs, pb->cje,
                         pb->cks, pb->cke+f3);
    BufferUtility::PackData((*coarse_fc).x3f, sendbuf,
                            pb->cis, pb->cie,
                            pb->cjs, pb->cje,
                            pb->cks, pb->cke+f3, p);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::FillSameRankFineToCoarseAMR(MeshBlock* pob, MeshBlock* pmb,
//!                                            LogicalLocation &loc)
//! \brief step 7: f2c, same MPI rank, different level (just restrict+copy, no pack/send)

void Mesh::FillSameRankFineToCoarseAMR(MeshBlock* pob, MeshBlock* pmb,
                                       LogicalLocation &loc) {
  MeshRefinement *pmr = pob->pmr;
  int il = pmb->is + ((loc.lx1 & 1LL) == 1LL)*pmb->block_size.nx1/2;
  int jl = pmb->js + ((loc.lx2 & 1LL) == 1LL)*pmb->block_size.nx2/2;
  int kl = pmb->ks + ((loc.lx3 & 1LL) == 1LL)*pmb->block_size.nx3/2;

  // absent a zip() feature for range-based for loops, manually advance the
  // iterator over "SMR/AMR-enrolled" cell-centered quantities on the new
  // MeshBlock in lock-step with pob
  auto pmb_cc_it = pmb->pmr->pvars_cc_.begin();
  // iterate MeshRefinement std::vectors on pob
  for (auto cc_pair : pmr->pvars_cc_) {
    AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
    AthenaArray<Real> *coarse_cc = std::get<1>(cc_pair);
    int nu = var_cc->GetDim4() - 1;
    pmr->RestrictCellCenteredValues(*var_cc, *coarse_cc,
                                    0, nu,
                                    pob->cis, pob->cie,
                                    pob->cjs, pob->cje,
                                    pob->cks, pob->cke);
    // copy from old/original/other MeshBlock (pob) to newly created block (pmb)
    AthenaArray<Real> const &src = *coarse_cc;
    AthenaArray<Real> &dst = *std::get<0>(*pmb_cc_it); // pmb->phydro->u;
    for (int nv=0; nv<=nu; nv++) {
      for (int k=kl, fk=pob->cks; fk<=pob->cke; k++, fk++) {
        for (int j=jl, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
          for (int i=il, fi=pob->cis; fi<=pob->cie; i++, fi++)
            dst(nv, k, j, i) = src(nv, fk, fj, fi);
        }
      }
    }
    pmb_cc_it++;
  }


  if ((NR_RADIATION_ENABLED|| IM_RADIATION_ENABLED)) {
    // restrict from pob block
    AthenaArray<Real> &var_cc = pob->pnrrad->ir;
    AthenaArray<Real> &coarse_cc = pob->pnrrad->coarse_ir_;
    int nu = var_cc.GetDim1() - 1;

    pmr->RestrictCellCenteredValues(var_cc, coarse_cc,-1,
                                    0, nu,
                                    pob->cis, pob->cie,
                                    pob->cjs, pob->cje,
                                    pob->cks, pob->cke);

    // now copy from pob to new pmb
    AthenaArray<Real> const &src = coarse_cc;
    AthenaArray<Real> &dst = pmb->pnrrad->ir;

    for (int k=kl, fk=pob->cks; fk<=pob->cke; k++, fk++) {
      for (int j=jl, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
        for (int i=il, fi=pob->cis; fi<=pob->cie; i++, fi++) {
          for (int nv=0; nv<=nu; nv++)
            dst(k, j, i, nv) = src(fk, fj, fi, nv);
        }
      }
    }
  }

  auto pmb_fc_it = pmb->pmr->pvars_fc_.begin();
  for (auto fc_pair : pmr->pvars_fc_) {
    FaceField *var_fc = std::get<0>(fc_pair);
    FaceField *coarse_fc = std::get<1>(fc_pair);
    pmr->RestrictFieldX1((*var_fc).x1f, (*coarse_fc).x1f,
                         pob->cis, pob->cie+1,
                         pob->cjs, pob->cje,
                         pob->cks, pob->cke);
    pmr->RestrictFieldX2((*var_fc).x2f, (*coarse_fc).x2f,
                         pob->cis, pob->cie,
                         pob->cjs, pob->cje+f2,
                         pob->cks, pob->cke);
    pmr->RestrictFieldX3((*var_fc).x3f, (*coarse_fc).x3f,
                         pob->cis, pob->cie,
                         pob->cjs, pob->cje,
                         pob->cks, pob->cke+f3);
    FaceField &src_b = *coarse_fc;
    FaceField &dst_b = *std::get<0>(*pmb_fc_it); // pmb->pfield->b;
    for (int k=kl, fk=pob->cks; fk<=pob->cke; k++, fk++) {
      for (int j=jl, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
        for (int i=il, fi=pob->cis; fi<=pob->cie+1; i++, fi++)
          dst_b.x1f(k, j, i) = src_b.x1f(fk, fj, fi);
      }
    }
    for (int k=kl, fk=pob->cks; fk<=pob->cke; k++, fk++) {
      for (int j=jl, fj=pob->cjs; fj<=pob->cje+f2; j++, fj++) {
        for (int i=il, fi=pob->cis; fi<=pob->cie; i++, fi++)
          dst_b.x2f(k, j, i) = src_b.x2f(fk, fj, fi);
      }
    }
    if (pmb->block_size.nx2 == 1) {
      int iu = il + pmb->block_size.nx1/2 - 1;
      for (int i=il; i<=iu; i++)
        dst_b.x2f(pmb->ks, pmb->js+1, i) = dst_b.x2f(pmb->ks, pmb->js, i);
    }
    for (int k=kl, fk=pob->cks; fk<=pob->cke+f3; k++, fk++) {
      for (int j=jl, fj=pob->cjs; fj<=pob->cje; j++, fj++) {
        for (int i=il, fi=pob->cis; fi<=pob->cie; i++, fi++)
          dst_b.x3f(k, j, i) = src_b.x3f(fk, fj, fi);
      }
    }
    if (pmb->block_size.nx3 == 1) {
      int iu = il + pmb->block_size.nx1/2 - 1, ju = jl + pmb->block_size.nx2/2 - 1;
      if (pmb->block_size.nx2 == 1) ju = jl;
      for (int j=jl; j<=ju; j++) {
        for (int i=il; i<=iu; i++)
          dst_b.x3f(pmb->ks+1, j, i) = dst_b.x3f(pmb->ks, j, i);
      }
    }
    pmb_fc_it++;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::FillSameRankCoarseToFineAMR(MeshBlock* pob, MeshBlock* pmb,
//!                                            LogicalLocation &newloc)
//! \brief step 7: c2f, same MPI rank, different level (just copy, no pack/send)

void Mesh::FillSameRankCoarseToFineAMR(MeshBlock* pob, MeshBlock* pmb,
                                       LogicalLocation &newloc) {
  MeshRefinement *pmr = pmb->pmr;
  int il = pob->cis - 1, iu = pob->cie + 1, jl = pob->cjs - f2,
      ju = pob->cje + f2, kl = pob->cks - f3, ku = pob->cke + f3;
  int cis = ((newloc.lx1 & 1LL) == 1LL)*pob->block_size.nx1/2 + pob->is - 1;
  int cjs = ((newloc.lx2 & 1LL) == 1LL)*pob->block_size.nx2/2 + pob->js - f2;
  int cks = ((newloc.lx3 & 1LL) == 1LL)*pob->block_size.nx3/2 + pob->ks - f3;

  auto pob_cc_it = pob->pmr->pvars_cc_.begin();
  // iterate MeshRefinement std::vectors on new pmb
  for (auto cc_pair : pmr->pvars_cc_) {
    AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
    AthenaArray<Real> *coarse_cc = std::get<1>(cc_pair);
    int nu = var_cc->GetDim4() - 1;

    AthenaArray<Real> const &src = *std::get<0>(*pob_cc_it);
    AthenaArray<Real> &dst = *coarse_cc;
    // fill the coarse buffer
    for (int nv=0; nv<=nu; nv++) {
      for (int k=kl, ck=cks; k<=ku; k++, ck++) {
        for (int j=jl, cj=cjs; j<=ju; j++, cj++) {
          for (int i=il, ci=cis; i<=iu; i++, ci++)
            dst(nv, k, j, i) = src(nv, ck, cj, ci);
        }
      }
    }
    pob_cc_it++;
  }

  if ((NR_RADIATION_ENABLED|| IM_RADIATION_ENABLED)) {
    // copy from pmb block
    AthenaArray<Real> &var_cc = pmb->pnrrad->ir;
    AthenaArray<Real> &coarse_cc = pmb->pnrrad->coarse_ir_;
    int nu = var_cc.GetDim1() - 1;

    // fill the coarse buffer
    AthenaArray<Real> &src = pob->pnrrad->ir;
    AthenaArray<Real> &dst = coarse_cc;

    for (int k=kl, ck=cks; k<=ku; k++, ck++) {
      for (int j=jl, cj=cjs; j<=ju; j++, cj++) {
        for (int i=il, ci=cis; i<=iu; i++, ci++) {
          for (int nv=0; nv<=nu; nv++)
            dst(k, j, i, nv) = src(ck, cj, ci, nv);
        }
      }
    }
  }

  auto pob_fc_it = pob->pmr->pvars_fc_.begin();
  // iterate MeshRefinement std::vectors on new pmb
  for (auto fc_pair : pmr->pvars_fc_) {
    FaceField *var_fc = std::get<0>(fc_pair);
    FaceField *coarse_fc = std::get<1>(fc_pair);

    FaceField &src_b = *std::get<0>(*pob_fc_it);
    FaceField &dst_b = *coarse_fc;
    for (int k=kl, ck=cks; k<=ku; k++, ck++) {
      for (int j=jl, cj=cjs; j<=ju; j++, cj++) {
        for (int i=il, ci=cis; i<=iu+1; i++, ci++)
          dst_b.x1f(k, j, i) = src_b.x1f(ck, cj, ci);
      }
    }
    for (int k=kl, ck=cks; k<=ku; k++, ck++) {
      for (int j=jl, cj=cjs; j<=ju+f2; j++, cj++) {
        for (int i=il, ci=cis; i<=iu; i++, ci++)
          dst_b.x2f(k, j, i) = src_b.x2f(ck, cj, ci);
      }
    }
    for (int k=kl, ck=cks; k<=ku+f3; k++, ck++) {
      for (int j=jl, cj=cjs; j<=ju; j++, cj++) {
        for (int i=il, ci=cis; i<=iu; i++, ci++)
          dst_b.x3f(k, j, i) = src_b.x3f(ck, cj, ci);
      }
    }
    pob_fc_it++;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::FinishRecvSameLevel(MeshBlock *pb, Real *recvbuf)
//! \brief step 8 (receive and load), branch 1 (same2same: unpack)
void Mesh::FinishRecvSameLevel(MeshBlock *pb, Real *recvbuf) {
  MeshRefinement *pmr = pb->pmr;
  int p = 0;
  for (AthenaArray<Real> &var_cc : pb->vars_cc_) {
    int nu = var_cc.GetDim4() - 1;
    BufferUtility::UnpackData(recvbuf, var_cc, 0, nu,
                              pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke, p);
  }

  if ((NR_RADIATION_ENABLED|| IM_RADIATION_ENABLED)) {
    int nu = pb->pnrrad->ir.GetDim1() - 1;
    BufferUtility::UnpackData(recvbuf, pb->pnrrad->ir, pb->ks, pb->ke, 0, nu,
                              pb->is, pb->ie, pb->js, pb->je, p);
  }

  for (FaceField &var_fc : pb->vars_fc_) {
    BufferUtility::UnpackData(recvbuf, var_fc.x1f,
                              pb->is, pb->ie+1, pb->js, pb->je, pb->ks, pb->ke, p);
    BufferUtility::UnpackData(recvbuf, var_fc.x2f,
                              pb->is, pb->ie, pb->js, pb->je+f2, pb->ks, pb->ke, p);
    BufferUtility::UnpackData(recvbuf, var_fc.x3f,
                              pb->is, pb->ie, pb->js, pb->je, pb->ks, pb->ke+f3, p);
    if (pb->block_size.nx2 == 1) {
      for (int i=pb->is; i<=pb->ie; i++)
        var_fc.x2f(pb->ks, pb->js+1, i) = var_fc.x2f(pb->ks, pb->js, i);
    }
    if (pb->block_size.nx3 == 1) {
      for (int j=pb->js; j<=pb->je; j++) {
        for (int i=pb->is; i<=pb->ie; i++)
          var_fc.x3f(pb->ks+1, j, i) = var_fc.x3f(pb->ks, j, i);
      }
    }
  }
  //! \warning (felker):
  //! * casting from "Real *" to "int *" in order to read single
  //!   appended integer from received buffer is slightly unsafe
  if (adaptive) {
    int *dcp = reinterpret_cast<int *>(&(recvbuf[p]));
    pmr->deref_count_ = *dcp;
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::FinishRecvFineToCoarseAMR(MeshBlock *pb, Real *recvbuf,
//!                                          LogicalLocation &lloc)
//! \brief step 8 (receive and load), branch 2 (f2c: unpack)

void Mesh::FinishRecvFineToCoarseAMR(MeshBlock *pb, Real *recvbuf,
                                     LogicalLocation &lloc) {
  MeshRefinement *pmr = pb->pmr;
  int ox1 = ((lloc.lx1 & 1LL) == 1LL), ox2 = ((lloc.lx2 & 1LL) == 1LL),
      ox3 = ((lloc.lx3 & 1LL) == 1LL);
  int p = 0, il, iu, jl, ju, kl, ku;
  if (ox1 == 0) il = pb->is,            iu = pb->is + pb->block_size.nx1/2 - 1;
  else        il = pb->is + pb->block_size.nx1/2, iu = pb->ie;
  if (ox2 == 0) jl = pb->js,            ju = pb->js + pb->block_size.nx2/2 - f2;
  else        jl = pb->js + pb->block_size.nx2/2, ju = pb->je;
  if (ox3 == 0) kl = pb->ks,            ku = pb->ks + pb->block_size.nx3/2 - f3;
  else        kl = pb->ks + pb->block_size.nx3/2, ku = pb->ke;

  for (auto cc_pair : pmr->pvars_cc_) {
    AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
    int nu = var_cc->GetDim4() - 1;
    BufferUtility::UnpackData(recvbuf, *var_cc, 0, nu,
                              il, iu, jl, ju, kl, ku, p);
  }

  if ((NR_RADIATION_ENABLED|| IM_RADIATION_ENABLED)) {
    int nu = pb->pnrrad->ir.GetDim1() - 1;
    BufferUtility::UnpackData(recvbuf, pb->pnrrad->ir, kl, ku, 0, nu,
                              il, iu, jl, ju, p);
  }

  for (auto fc_pair : pmr->pvars_fc_) {
    FaceField *var_fc = std::get<0>(fc_pair);
    FaceField &dst_b = *var_fc;
    BufferUtility::UnpackData(recvbuf, dst_b.x1f,
                              il, iu+1, jl, ju, kl, ku, p);
    BufferUtility::UnpackData(recvbuf, dst_b.x2f,
                              il, iu, jl, ju+f2, kl, ku, p);
    BufferUtility::UnpackData(recvbuf, dst_b.x3f,
                              il, iu, jl, ju, kl, ku+f3, p);
    if (pb->block_size.nx2 == 1) {
      for (int i=il; i<=iu; i++)
        dst_b.x2f(pb->ks, pb->js+1, i) = dst_b.x2f(pb->ks, pb->js, i);
    }
    if (pb->block_size.nx3 == 1) {
      for (int j=jl; j<=ju; j++) {
        for (int i=il; i<=iu; i++)
          dst_b.x3f(pb->ks+1, j, i) = dst_b.x3f(pb->ks, j, i);
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::ReceiveCoarseToFineAMR(MeshBlock *pb, Real *recvbuf)
//! \brief step 8, branch 2 (c2f: receive and unpack)

void Mesh::ReceiveCoarseToFineAMR(MeshBlock *pb, Real *recvbuf) {
  MeshRefinement *pmr = pb->pmr;
  int p = 0;
  int il = pb->cis - 1, iu = pb->cie+1, jl = pb->cjs - f2,
      ju = pb->cje + f2, kl = pb->cks - f3, ku = pb->cke + f3;
  for (auto cc_pair : pmr->pvars_cc_) {
    AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
    AthenaArray<Real> *coarse_cc = std::get<1>(cc_pair);
    int nu = var_cc->GetDim4() - 1;
    BufferUtility::UnpackData(recvbuf, *coarse_cc,
                              0, nu, il, iu, jl, ju, kl, ku, p);
  }

  if ((NR_RADIATION_ENABLED|| IM_RADIATION_ENABLED)) {
    AthenaArray<Real> &var_cc = pb->pnrrad->ir;
    AthenaArray<Real> &coarse_cc = pb->pnrrad->coarse_ir_;
    int nu = var_cc.GetDim1() - 1;
    BufferUtility::UnpackData(recvbuf, coarse_cc, kl, ku,
                              0, nu, il, iu, jl, ju, p);
  }

  for (auto fc_pair : pmr->pvars_fc_) {
    FaceField *var_fc = std::get<0>(fc_pair);
    FaceField *coarse_fc = std::get<1>(fc_pair);

    BufferUtility::UnpackData(recvbuf, (*coarse_fc).x1f,
                              il, iu+1, jl, ju, kl, ku, p);
    BufferUtility::UnpackData(recvbuf, (*coarse_fc).x2f,
                              il, iu, jl, ju+f2, kl, ku, p);
    BufferUtility::UnpackData(recvbuf, (*coarse_fc).x3f,
                              il, iu, jl, ju, kl, ku+f3, p);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::ProlongateMeshBlock(MeshBlock *pb)
//! \brief step 8, branch 2 (c2f: prolongate)

void Mesh::ProlongateMeshBlock(MeshBlock *pb) {
  MeshRefinement *pmr = pb->pmr;
  for (auto cc_pair : pmr->pvars_cc_) {
    AthenaArray<Real> *var_cc = std::get<0>(cc_pair);
    AthenaArray<Real> *coarse_cc = std::get<1>(cc_pair);
    int nu = var_cc->GetDim4() - 1;
    pmr->ProlongateCellCenteredValues(*coarse_cc, *var_cc, 0, nu,
                   pb->cis, pb->cie, pb->cjs, pb->cje, pb->cks, pb->cke);
  }


  if ((NR_RADIATION_ENABLED|| IM_RADIATION_ENABLED)) {
    // copy from pmb block
    AthenaArray<Real> &var_cc = pb->pnrrad->ir;
    AthenaArray<Real> &coarse_cc = pb->pnrrad->coarse_ir_;
    int nu = var_cc.GetDim1() - 1;

    pmr->ProlongateCellCenteredValues(coarse_cc, var_cc, -1, 0, nu,
                   pb->cis, pb->cie, pb->cjs, pb->cje, pb->cks, pb->cke);
  }

  int il = pb->cis, iu = pb->cie+1, jl = pb->cjs, ju = pb->cje + f2,
      kl = pb->cks, ku = pb->cke + f3;
  // Step FFC8. skip the surface fields contacting previously refined MeshBlocks
  if (pmr->flag_ffc_recv_[BoundaryFace::inner_x1]) il++;
  if (pmr->flag_ffc_recv_[BoundaryFace::outer_x1]) iu--;
  if (pmr->flag_ffc_recv_[BoundaryFace::inner_x2]) jl++;
  if (pmr->flag_ffc_recv_[BoundaryFace::outer_x2]) ju--;
  if (pmr->flag_ffc_recv_[BoundaryFace::inner_x3]) kl++;
  if (pmr->flag_ffc_recv_[BoundaryFace::outer_x3]) ku--;

  for (auto fc_pair : pmr->pvars_fc_) {
    FaceField *var_fc = std::get<0>(fc_pair);
    FaceField *coarse_fc = std::get<1>(fc_pair);

    pmr->ProlongateSharedFieldX1((*coarse_fc).x1f, (*var_fc).x1f,
                                 il, iu, pb->cjs, pb->cje, pb->cks, pb->cke);
    pmr->ProlongateSharedFieldX2((*coarse_fc).x2f, (*var_fc).x2f,
                                 pb->cis, pb->cie, jl, ju, pb->cks, pb->cke);
    pmr->ProlongateSharedFieldX3((*coarse_fc).x3f, (*var_fc).x3f,
                                 pb->cis, pb->cie, pb->cjs, pb->cje, kl, ku);
    pmr->ProlongateInternalField(*var_fc, pb->cis, pb->cie,
                                 pb->cjs, pb->cje, pb->cks, pb->cke);
  }
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Mesh::PrepareAndSendFaceFieldCorrection(LogicalLocation *newloc,
//!                              int *ranklist, int *newrank, int *nslist, int nbtold)
//! \brief the first part of the face field correction

void Mesh::PrepareAndSendFaceFieldCorrection(LogicalLocation *newloc,
                                  int *ranklist, int *newrank, int *nslist, int nbtold) {
  const int bnx1 = my_blocks(0)->block_size.nx1;
  const int bnx2 = my_blocks(0)->block_size.nx2;
  const int bnx3 = my_blocks(0)->block_size.nx3;
  const int nf = static_cast<int>(my_blocks(0)->pmr->pvars_fc_.size());
  const int sffc1 = nf*bnx2*bnx3, sffc2 = nf*bnx1*bnx3, sffc3 = nf*bnx2*bnx3;

  // Step FFC1. Construct a hash map of MeshBlocks before refinement
  for (int i = 0; i < nbtold; ++i) {
    LogicalLocation const &loc = loclist[i];
    if (loc.level > root_level) {
      int lev = loc.level - root_level - 1;
      locmap_[lev][loc] = i;
    }
  }

  // Step FFC2. Find pairs that need face field correction
  int s = 0;
  for (int n : refined_) { // loop over refined MeshBlocks
    LogicalLocation const &nloc = newloc[n];
    int xf1 = static_cast<int>(nloc.lx1 & 1LL);
    int xf2 = static_cast<int>(nloc.lx2 & 1LL);
    int xf3 = static_cast<int>(nloc.lx3 & 1LL);
    std::int64_t o1, o2, o3;
    o1 = (xf1 << 1) - 1;
    o2 = (xf2 << 1) - 1;
    o3 = (xf3 << 1) - 1;
    LogicalLocation floc = nloc;
    int lev = floc.level - root_level - 1;
    floc.lx1 += o1;
    bool fskip = false;
    if (floc.lx1 < 0) {
      if (mesh_bcs[BoundaryFace::inner_x1] == BoundaryFlag::periodic)
        floc.lx1 = (nrbx1 << (lev+1)) - 1;
      else
        fskip = true;
    } else if (floc.lx1 >= (nrbx1 << (lev+1))) {
      if (mesh_bcs[BoundaryFace::outer_x1] == BoundaryFlag::periodic)
        floc.lx1 = 0;
      else
        fskip = true;
    }
    if (!fskip) {
      auto itr = locmap_[lev].find(floc);
      if (itr != locmap_[lev].end()) {
        const auto& it = *itr;
        int i = it.second;
        if (ranklist[i] == Globals::my_rank) {
          ffc_send_.emplace_back(i, n, xf1^1, sffc1, -1);
          if (newrank[n] == Globals::my_rank)
            ffc_recv_.emplace_back(i, n, xf1, 0, s);
          s++;
        } else if (newrank[n] == Globals::my_rank) {
          ffc_recv_.emplace_back(i, n, xf1, sffc1, -1);
        }
      }
    }
    if (f2) {
      floc.lx1 = nloc.lx1, floc.lx2 += o2;
      fskip = false;
      if (floc.lx2 < 0) {
        if (mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic)
          floc.lx2 = (nrbx2 << (lev+1)) - 1;
        else
          fskip = true;
      } else if (floc.lx2 >= (nrbx2 << (lev+1))) {
        if (mesh_bcs[BoundaryFace::outer_x2] == BoundaryFlag::periodic)
          floc.lx2 = 0;
        else
          fskip = true;
      }
      if (!fskip) {
        auto itr = locmap_[lev].find(floc);
        if (itr != locmap_[lev].end()) {
          const auto& it = *itr;
          int i = it.second;
          if (ranklist[i] == Globals::my_rank) {
            ffc_send_.emplace_back(i, n, (xf2^1) + 2, sffc2, -1);
            if (newrank[n] == Globals::my_rank)
              ffc_recv_.emplace_back(i, n, xf2 + 2, 0, s);
            s++;
          } else if (newrank[n] == Globals::my_rank) {
            ffc_recv_.emplace_back(i, n, xf2 + 2, sffc2, -1);
          }
        }
      }
    }
    if (f3) {
      floc.lx2 = nloc.lx2, floc.lx3 += o3;
      fskip = false;
      if (floc.lx3 < 0) {
        if (mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic)
          floc.lx3 = (nrbx3 << (lev+1)) - 1;
        else
          fskip = true;
      } else if (floc.lx3 >= (nrbx3 << (lev+1))) {
        if (mesh_bcs[BoundaryFace::outer_x3] == BoundaryFlag::periodic)
          floc.lx3 = 0;
        else
          fskip = true;
      }
      if (!fskip) {
        auto itr = locmap_[lev].find(floc);
        if (itr != locmap_[lev].end()) {
          const auto& it = *itr;
          int i = it.second;
          if (ranklist[i] == Globals::my_rank) {
            ffc_send_.emplace_back(i, n, (xf3^1) + 4, sffc3, -1);
            if (newrank[n] == Globals::my_rank)
              ffc_recv_.emplace_back(i, n, xf3 + 4, 0, s);
            s++;
          } else if (newrank[n] == Globals::my_rank) {
            ffc_recv_.emplace_back(i, n, xf3 + 4, sffc3, -1);
          }
        }
      }
    }
  }

  // Step FFC3. Initiate MPI receive
#ifdef MPI_PARALLEL
  for (FaceFieldCorrection &t : ffc_recv_) {
    if (ranklist[t.from] != Globals::my_rank) {
      int tag = CreateFaceFieldCorrectionMPITag(t.to-nslist[newrank[t.to]], t.face);
      MPI_Irecv(t.buf, t.size, MPI_ATHENA_REAL, ranklist[t.from],
                tag, MPI_COMM_WORLD, &(t.req));
    }
  }
#endif

  // Step FFC4. Pack and send shared face fields
  for (FaceFieldCorrection &f : ffc_send_) {
    MeshBlock *pmb = FindMeshBlock(f.from);
    int p = 0;
    for (FaceField &var_fc : pmb->vars_fc_) {
      switch (f.face) {
        case BoundaryFace::inner_x1:
          BufferUtility::PackData(var_fc.x1f, f.buf, pmb->is, pmb->is,
                                  pmb->js, pmb->je, pmb->ks, pmb->ke, p);
          break;
        case BoundaryFace::outer_x1:
          BufferUtility::PackData(var_fc.x1f, f.buf, pmb->ie+1, pmb->ie+1,
                                  pmb->js, pmb->je, pmb->ks, pmb->ke, p);
          break;
        case BoundaryFace::inner_x2:
          BufferUtility::PackData(var_fc.x2f, f.buf, pmb->is, pmb->ie,
                                  pmb->js, pmb->js, pmb->ks, pmb->ke, p);
          break;
        case BoundaryFace::outer_x2:
          BufferUtility::PackData(var_fc.x2f, f.buf, pmb->is, pmb->ie,
                                  pmb->je+f2, pmb->je+f2, pmb->ks, pmb->ke, p);
          break;
        case BoundaryFace::inner_x3:
          BufferUtility::PackData(var_fc.x3f, f.buf, pmb->is, pmb->ie,
                                  pmb->js, pmb->je, pmb->ks, pmb->ks, p);
          break;
        case BoundaryFace::outer_x3:
          BufferUtility::PackData(var_fc.x3f, f.buf, pmb->is, pmb->ie,
                                  pmb->js, pmb->je, pmb->ke+f3, pmb->ke+f3, p);
          break;
        default:
          break;
      }
    }
#ifdef MPI_PARALLEL
    if (newrank[f.to] != Globals::my_rank) {
      int face = f.face^1;
      int tag = CreateFaceFieldCorrectionMPITag(f.to-nslist[newrank[f.to]], face);
      MPI_Isend(f.buf, f.size, MPI_ATHENA_REAL, newrank[f.to],
                tag, MPI_COMM_WORLD, &(f.req));
    }
#endif
  }
}


//----------------------------------------------------------------------------------------
//! \fn void Mesh::ReceiveAndSetFaceFieldCorrection(int *newrank)
//! \brief the second part of the face field correction

void Mesh::ReceiveAndSetFaceFieldCorrection(int *newrank) {
  // Step FFC5. clear face field correction flags
  for (int n = 0; n < nblocal; ++n) {
    MeshRefinement *pmr = my_blocks(n)->pmr;
    for (int i = 0; i < 6; ++i)
      pmr->flag_ffc_recv_[i] = false;
  }

  // Step FFC6. wait and unpack/set
  for (FaceFieldCorrection &t : ffc_recv_) {
    MeshBlock *pmb = FindMeshBlock(t.to);
    pmb->pmr->flag_ffc_recv_[t.face] = true;
    Real *buf = t.buf;
    if (ranklist[t.from] == Globals::my_rank) // local copy
      buf = ffc_send_[t.src].buf;
#ifdef MPI_PARALLEL
    else // MPI
      MPI_Wait(&(t.req), MPI_STATUS_IGNORE);
#endif
    int p = 0;
    for (FaceField &var_fc : pmb->vars_fc_) {
      switch (t.face) {
        case BoundaryFace::inner_x1:
          BufferUtility::UnpackData(buf, var_fc.x1f, pmb->is, pmb->is,
                                    pmb->js, pmb->je, pmb->ks, pmb->ke, p);
          break;
        case BoundaryFace::outer_x1:
          BufferUtility::UnpackData(buf, var_fc.x1f, pmb->ie+1, pmb->ie+1,
                                    pmb->js, pmb->je, pmb->ks, pmb->ke, p);
          break;
        case BoundaryFace::inner_x2:
          BufferUtility::UnpackData(buf, var_fc.x2f, pmb->is, pmb->ie,
                                    pmb->js, pmb->js, pmb->ks, pmb->ke, p);
          break;
        case BoundaryFace::outer_x2:
          BufferUtility::UnpackData(buf, var_fc.x2f, pmb->is, pmb->ie,
                                    pmb->je+f2, pmb->je+f2, pmb->ks, pmb->ke, p);
          break;
        case BoundaryFace::inner_x3:
          BufferUtility::UnpackData(buf, var_fc.x3f, pmb->is, pmb->ie,
                                    pmb->js, pmb->je, pmb->ks, pmb->ks, p);
          break;
        case BoundaryFace::outer_x3:
          BufferUtility::UnpackData(buf, var_fc.x3f, pmb->is, pmb->ie,
                                    pmb->js, pmb->je, pmb->ke+f3, pmb->ke+f3, p);
          break;
        default:
          break;
      }
    }
  }

  // Step FFC7. Finalize MPI send
#ifdef MPI_PARALLEL
  for (FaceFieldCorrection &f : ffc_send_) {
    if (newrank[f.to] != Globals::my_rank) {
      MPI_Wait(&(f.req), MPI_STATUS_IGNORE);
    }
  }
#endif
}

//----------------------------------------------------------------------------------------
//! \fn int Mesh::CreateAMRMPITag(int lid, int ox1, int ox2, int ox3)
//! \brief calculate an MPI tag for AMR block transfer
//!
//! tag = local id of destination (remaining bits) + ox1(1 bit) + ox2(1 bit) + ox3(1 bit)
//!       + physics(5 bits)
//!
//! See comments on BoundaryBase::CreateBvalsMPITag()

int Mesh::CreateAMRMPITag(int lid, int ox1, int ox2, int ox3) {
  // former "AthenaTagMPI" AthenaTagMPI::amr=8 redefined to 0
  return (lid<<8) | (ox1<<7)| (ox2<<6) | (ox3<<5);
}

//----------------------------------------------------------------------------------------
//! \fn int Mesh::CreateFaceFieldCorrectionMPITag(int lid, int face)
//! \brief calculate an MPI tag for face field correction
//!
//! tag = local id of destination (remaining bits) + ox1(1 bit) + ox2(1 bit) + ox3(1 bit)
//!       + physics(5 bits)
//!

int Mesh::CreateFaceFieldCorrectionMPITag(int lid, int face) {
  // Set the bottom bit = 1
  // It should be OK as this communication does not overlap with te main integrator,
  // and all we need is to distinguish FFC communications from AMR communications.
  return (lid<<8) | (face << 5) | 1;
}

