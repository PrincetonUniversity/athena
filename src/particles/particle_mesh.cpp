//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle_mesh.cpp
//! \brief implements ParticleMesh class used for operations involved in particle-mesh
//!        methods.

// Standard library
#include <algorithm>
#include <cstring>
#include <sstream>

// Athena++ classes headers
#include "../athena.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../utils/buffer_utils.hpp"
#include "particle_mesh.hpp"
#include "particles.hpp"

// Local function prototypes.
static Real _WeightFunction(Real dxi);

//--------------------------------------------------------------------------------------
//! \fn int ParticleMesh::AddMeshAux()
//! \brief adds one auxiliary to the mesh and returns the index.

int ParticleMesh::AddMeshAux() {
  return nmeshaux++;
}

//--------------------------------------------------------------------------------------
//! \fn ParticleMesh::ParticleMesh(Particles *ppar, int nmeshaux)
//! \brief constructs a new ParticleMesh instance.

ParticleMesh::ParticleMesh(Particles *ppar, MeshBlock *pmb) : nmeshaux(0), iweight(-1),
  imom1(-1), imom2(-1), imom3(-1), imass(-1), updated(false),
  is(pmb->is), ie(pmb->ie), js(pmb->js), je(pmb->je), ks(pmb->ks), ke(pmb->ke),
  active1_(ppar->active1_), active2_(ppar->active2_), active3_(ppar->active3_),
  dxi1_(active1_ ? RINF : 0),
  dxi2_(active2_ ? RINF : 0),
  dxi3_(active3_ ? RINF : 0),
  nx1_(pmb->ncells1), nx2_(pmb->ncells2), nx3_(pmb->ncells3),
  ncells_(nx1_ * nx2_ * nx3_),
  npc1_(active1_ ? NPC : 1), npc2_(active2_ ? NPC : 1), npc3_(active3_ ? NPC : 1),
  my_ipar_(ppar->my_ipar_), ppar_(ppar), pmb_(pmb), pmesh_(ppar->pmy_mesh) {
  // Add weight in meshaux.
  iweight = AddMeshAux();
  if (ppar_->imass != -1) imass = AddMeshAux();

  // Add momentum in meshaux
  imom1 = AddMeshAux();
  imom2 = AddMeshAux();
  imom3 = AddMeshAux();

  meshaux.NewAthenaArray(nmeshaux, nx3_, nx2_, nx1_);
  coarse_meshaux_.NewAthenaArray(nmeshaux, pmb->ncc3, pmb->ncc2, pmb->ncc1);

  // Get a shorthand to weights.
  weight.InitWithShallowSlice(meshaux, 4, iweight, 1);
  if (ppar_->imass != -1) density.InitWithShallowSlice(meshaux, 4, imass, 1);

  // Enroll CellCenteredBoundaryVariable object
  AthenaArray<Real> empty_flux[3]=
    {AthenaArray<Real>(),AthenaArray<Real>(), AthenaArray<Real>()};

  pmbvar = new ParticleMeshBoundaryVariable(pmb, &meshaux, &coarse_meshaux_,
                                            empty_flux, this);
  pmbvar->bvar_index = pmb_->pbval->bvars.size();
  pmb_->pbval->bvars.push_back(pmbvar);
  // Add particle mesh boundary variable to the list for main integrator
  // if that particle exert gravity. Otherwise, add it to the list for
  // outputs.

}

//--------------------------------------------------------------------------------------
//! \fn ParticleMesh::~ParticleMesh()
//! \brief destructs a ParticleMesh instance.

ParticleMesh::~ParticleMesh() {
  // Destroy the particle meshblock.
  meshaux.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
//! \fn Real ParticleMesh::FindMaximumWeight()
//! \brief returns the maximum weight in the meshblock.

Real ParticleMesh::FindMaximumWeight() const {
  Real wmax = 0.0;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        wmax = std::max(wmax, weight(k,j,i));
  return wmax;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::InterpolateMeshToParticles(
//!              const AthenaArray<Real>& meshsrc, int ms1,
//!              AthenaArray<Real>& par, int p1, int nprop)
//! \brief interpolates meshsrc from property index ms1 to ms1+nprop-1 onto particle
//!     array par (realprop, auxprop, or work in Particles class) from property index p1
//!     to p1+nprop-1.

void ParticleMesh::InterpolateMeshToParticles(
         const AthenaArray<Real>& meshsrc, int ms1,
         AthenaArray<Real>& par, int p1, int nprop) {
  // Transpose meshsrc.
  int nx1 = meshsrc.GetDim1(), nx2 = meshsrc.GetDim2(), nx3 = meshsrc.GetDim3();
  AthenaArray<Real> u;
  u.NewAthenaArray(nx3,nx2,nx1,nprop);
  for (int n = 0; n < nprop; ++n)
    for (int k = 0; k < nx3; ++k)
      for (int j = 0; j < nx2; ++j)
        for (int i = 0; i < nx1; ++i)
          u(k,j,i,n) = meshsrc(ms1+n,k,j,i);

  // Allocate space for SIMD.
  Real **w1 __attribute__((aligned(CACHELINE_BYTES))) = new Real*[npc1_];
  Real **w2 __attribute__((aligned(CACHELINE_BYTES))) = new Real*[npc2_];
  Real **w3 __attribute__((aligned(CACHELINE_BYTES))) = new Real*[npc3_];
  for (int i = 0; i < npc1_; ++i)
    w1[i] = new Real[SIMD_WIDTH];
  for (int i = 0; i < npc2_; ++i)
    w2[i] = new Real[SIMD_WIDTH];
  for (int i = 0; i < npc3_; ++i)
    w3[i] = new Real[SIMD_WIDTH];
  int imb1v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  int imb2v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  int imb3v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));

  // Loop over each particle.
  int npar = ppar_->npar;
  for (int k = 0; k < npar; k += SIMD_WIDTH) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int kk = 0; kk < std::min(SIMD_WIDTH, npar-k); ++kk) {
      int kkk = k + kk;

      // Find the domain the particle influences.
      Real xi1 = ppar_->xi1(kkk), xi2 = ppar_->xi2(kkk), xi3 = ppar_->xi3(kkk);
      int imb1 = static_cast<int>(xi1 - dxi1_),
          imb2 = static_cast<int>(xi2 - dxi2_),
          imb3 = static_cast<int>(xi3 - dxi3_);
      xi1 = imb1 + 0.5 - xi1;
      xi2 = imb2 + 0.5 - xi2;
      xi3 = imb3 + 0.5 - xi3;

      imb1v[kk] = imb1;
      imb2v[kk] = imb2;
      imb3v[kk] = imb3;

      // Weigh each cell.
#pragma loop count (NPC)
      for (int i = 0; i < npc1_; ++i)
        w1[i][kk] = active1_ ? _WeightFunction(xi1 + i) : 1.0;
#pragma loop count (NPC)
      for (int i = 0; i < npc2_; ++i)
        w2[i][kk] = active2_ ? _WeightFunction(xi2 + i) : 1.0;
#pragma loop count (NPC)
      for (int i = 0; i < npc3_; ++i)
        w3[i][kk] = active3_ ? _WeightFunction(xi3 + i) : 1.0;
    }

#pragma ivdep
    for (int kk = 0; kk < std::min(SIMD_WIDTH, npar-k); ++kk) {
      int kkk = k + kk;

      // Initiate interpolation.
      Real *pd = new Real[nprop];
      for (int i = 0; i < nprop; ++i)
        pd[i] = 0.0;

      int imb1 = imb1v[kk], imb2 = imb2v[kk], imb3 = imb3v[kk];

#pragma loop count (NPC)
      for (int ipc3 = 0; ipc3 < npc3_; ++ipc3) {
#pragma loop count (NPC)
        for (int ipc2 = 0; ipc2 < npc2_; ++ipc2) {
#pragma loop count (NPC)
          for (int ipc1 = 0; ipc1 < npc1_; ++ipc1) {
            Real w = w1[ipc1][kk] * w2[ipc2][kk] * w3[ipc3][kk];

            // Interpolate meshsrc to particles.
            for (int n = 0; n < nprop; ++n)
              pd[n] += w * u(imb3+ipc3,imb2+ipc2,imb1+ipc1,n);
          }
        }
      }

      // Record the final interpolated properties.
      for (int n = 0; n < nprop; ++n)
        par(p1+n,kkk) = pd[n];

      delete [] pd;
    }
  }

  // Release working arrays.
  u.DeleteAthenaArray();
  for (int i = 0; i < npc1_; ++i)
    delete [] w1[i];
  for (int i = 0; i < npc2_; ++i)
    delete [] w2[i];
  for (int i = 0; i < npc3_; ++i)
    delete [] w3[i];
  delete [] w1;
  delete [] w2;
  delete [] w3;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::AssignParticlesToMeshAux(
//       const AthenaArray<Real>& par, int p1, int ma1, int nprop)
//! \brief assigns par (realprop, auxprop, or work in Particles class) from property
//!        index p1 to p1+nprop-1 onto meshaux from property index ma1 and up.

void ParticleMesh::AssignParticlesToMeshAux(
         const AthenaArray<Real>& par, int p1, int ma1, int nprop) {
  // Zero out meshaux.
  std::fill(&weight(0,0,0), &weight(0,0,0) + ncells_, 0.0);
  std::fill(&meshaux(ma1,0,0,0), &meshaux(ma1+nprop,0,0,0), 0.0);

  // Allocate space for SIMD.
  Real **w1 __attribute__((aligned(CACHELINE_BYTES))) = new Real*[npc1_];
  Real **w2 __attribute__((aligned(CACHELINE_BYTES))) = new Real*[npc2_];
  Real **w3 __attribute__((aligned(CACHELINE_BYTES))) = new Real*[npc3_];
  for (int i = 0; i < npc1_; ++i)
    w1[i] = new Real[SIMD_WIDTH];
  for (int i = 0; i < npc2_; ++i)
    w2[i] = new Real[SIMD_WIDTH];
  for (int i = 0; i < npc3_; ++i)
    w3[i] = new Real[SIMD_WIDTH];
  int imb1v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  int imb2v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  int imb3v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));

  // Loop over each particle.
  int npar = ppar_->npar;
  for (int k = 0; k < npar; k += SIMD_WIDTH) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int kk = 0; kk < std::min(SIMD_WIDTH, npar-k); ++kk) {
      int kkk = k + kk;

      // Find the domain the particle influences.
      Real xi1 = ppar_->xi1(kkk), xi2 = ppar_->xi2(kkk), xi3 = ppar_->xi3(kkk);
      int imb1 = static_cast<int>(xi1 - dxi1_),
          imb2 = static_cast<int>(xi2 - dxi2_),
          imb3 = static_cast<int>(xi3 - dxi3_);
      xi1 = imb1 + 0.5 - xi1;
      xi2 = imb2 + 0.5 - xi2;
      xi3 = imb3 + 0.5 - xi3;

      imb1v[kk] = imb1;
      imb2v[kk] = imb2;
      imb3v[kk] = imb3;

      // Weigh each cell.
#pragma loop count (NPC)
      for (int i = 0; i < npc1_; ++i)
        w1[i][kk] = active1_ ? _WeightFunction(xi1 + i) : 1.0;
#pragma loop count (NPC)
      for (int i = 0; i < npc2_; ++i)
        w2[i][kk] = active2_ ? _WeightFunction(xi2 + i) : 1.0;
#pragma loop count (NPC)
      for (int i = 0; i < npc3_; ++i)
        w3[i][kk] = active3_ ? _WeightFunction(xi3 + i) : 1.0;
    }

#pragma ivdep
    for (int kk = 0; kk < std::min(SIMD_WIDTH, npar-k); ++kk) {
      int kkk = k + kk;

      // Fetch particle properties.
      Real *ps = new Real[nprop];
      for (int i = 0; i < nprop; ++i)
        ps[i] = par(p1+i,kkk);

      int imb1 = imb1v[kk], imb2 = imb2v[kk], imb3 = imb3v[kk];

#pragma loop count (NPC)
      for (int ipc3 = 0; ipc3 < npc3_; ++ipc3) {
#pragma loop count (NPC)
        for (int ipc2 = 0; ipc2 < npc2_; ++ipc2) {
#pragma loop count (NPC)
          for (int ipc1 = 0; ipc1 < npc1_; ++ipc1) {
            Real w = w1[ipc1][kk] * w2[ipc2][kk] * w3[ipc3][kk];

            // Record the weights.
            weight(imb3+ipc3,imb2+ipc2,imb1+ipc1) += w;

            // Assign particles to meshaux.
            for (int n = 0; n < nprop; ++n)
              meshaux(ma1+n,imb3+ipc3,imb2+ipc2,imb1+ipc1) += w * ps[n];
          }
        }
      }

      delete [] ps;
    }
  }

  // Release working array.
  for (int i = 0; i < npc1_; ++i)
    delete [] w1[i];
  for (int i = 0; i < npc2_; ++i)
    delete [] w2[i];
  for (int i = 0; i < npc3_; ++i)
    delete [] w3[i];
  delete [] w1;
  delete [] w2;
  delete [] w3;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::InterpolateMeshAndAssignParticles(
//!              const AthenaArray<Real>& meshsrc, int ms1,
//!              AthenaArray<Real>& pardst, int pd1, int ni,
//!              const AthenaArray<Real>& parsrc, int ps1, int ma1, int na)
//! \brief interpolates meshsrc from property index ms1 to ms1 + ni - 1 onto particle
//!     array pardst from index pd1 to pd1 + ni - 1, and assigns parsrc from property
//!     index ps1 to ps1 + na - 1 onto meshaux from ma1 to ma1 + na - 1.  The arrays
//!     parsrc and pardst can be realprop, auxprop, or work in Particles class.

void ParticleMesh::InterpolateMeshAndAssignParticles(
         const AthenaArray<Real>& meshsrc, int ms1,
         AthenaArray<Real>& pardst, int pd1, int ni,
         const AthenaArray<Real>& parsrc, int ps1, int ma1, int na) {
  // Zero out meshaux.
  std::fill(&weight(0,0,0), &weight(0,0,0) + ncells_, 0.0);
  std::fill(&meshaux(ma1,0,0,0), &meshaux(ma1+na,0,0,0), 0.0);

  // Transpose meshsrc.
  int nx1 = meshsrc.GetDim1(), nx2 = meshsrc.GetDim2(), nx3 = meshsrc.GetDim3();
  AthenaArray<Real> u;
  u.NewAthenaArray(nx3,nx2,nx1,ni);
  for (int n = 0; n < ni; ++n)
    for (int k = 0; k < nx3; ++k)
      for (int j = 0; j < nx2; ++j)
        for (int i = 0; i < nx1; ++i)
          u(k,j,i,n) = meshsrc(ms1+n,k,j,i);

  // Allocate space for SIMD.
  Real **w1 __attribute__((aligned(CACHELINE_BYTES))) = new Real*[npc1_];
  Real **w2 __attribute__((aligned(CACHELINE_BYTES))) = new Real*[npc2_];
  Real **w3 __attribute__((aligned(CACHELINE_BYTES))) = new Real*[npc3_];
  for (int i = 0; i < npc1_; ++i)
    w1[i] = new Real[SIMD_WIDTH];
  for (int i = 0; i < npc2_; ++i)
    w2[i] = new Real[SIMD_WIDTH];
  for (int i = 0; i < npc3_; ++i)
    w3[i] = new Real[SIMD_WIDTH];
  int imb1v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  int imb2v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));
  int imb3v[SIMD_WIDTH] __attribute__((aligned(CACHELINE_BYTES)));

  // Loop over each particle.
  int npar = ppar_->npar;
  for (int k = 0; k < npar; k += SIMD_WIDTH) {
#pragma omp simd simdlen(SIMD_WIDTH)
    for (int kk = 0; kk < std::min(SIMD_WIDTH, npar-k); ++kk) {
      int kkk = k + kk;

      // Find the domain the particle influences.
      Real xi1 = ppar_->xi1(kkk), xi2 = ppar_->xi2(kkk), xi3 = ppar_->xi3(kkk);
      int imb1 = static_cast<int>(xi1 - dxi1_),
          imb2 = static_cast<int>(xi2 - dxi2_),
          imb3 = static_cast<int>(xi3 - dxi3_);
      xi1 = imb1 + 0.5 - xi1;
      xi2 = imb2 + 0.5 - xi2;
      xi3 = imb3 + 0.5 - xi3;

      imb1v[kk] = imb1;
      imb2v[kk] = imb2;
      imb3v[kk] = imb3;

      // Weigh each cell.
#pragma loop count (NPC)
      for (int i = 0; i < npc1_; ++i)
        w1[i][kk] = active1_ ? _WeightFunction(xi1 + i) : 1.0;
#pragma loop count (NPC)
      for (int i = 0; i < npc2_; ++i)
        w2[i][kk] = active2_ ? _WeightFunction(xi2 + i) : 1.0;
#pragma loop count (NPC)
      for (int i = 0; i < npc3_; ++i)
        w3[i][kk] = active3_ ? _WeightFunction(xi3 + i) : 1.0;
    }

#pragma ivdep
    for (int kk = 0; kk < std::min(SIMD_WIDTH, npar-k); ++kk) {
      int kkk = k + kk;

      // Initiate interpolation and fetch particle properties.
      Real *pd = new Real[ni];
      Real *ps = new Real[na];
      for (int i = 0; i < ni; ++i)
        pd[i] = 0.0;
      for (int i = 0; i < na; ++i)
        ps[i] = parsrc(ps1+i,kkk);

      int imb1 = imb1v[kk], imb2 = imb2v[kk], imb3 = imb3v[kk];

#pragma loop count (NPC)
      for (int ipc3 = 0; ipc3 < npc3_; ++ipc3) {
#pragma loop count (NPC)
        for (int ipc2 = 0; ipc2 < npc2_; ++ipc2) {
#pragma loop count (NPC)
          for (int ipc1 = 0; ipc1 < npc1_; ++ipc1) {
            Real w = w1[ipc1][kk] * w2[ipc2][kk] * w3[ipc3][kk];

            // Record the weights.
            weight(imb3+ipc3,imb2+ipc2,imb1+ipc1) += w;

            // Interpolate meshsrc to particles.
            for (int n = 0; n < ni; ++n)
              pd[n] += w * u(imb3+ipc3,imb2+ipc2,imb1+ipc1,n);

            // Assign particles to meshaux.
            for (int n = 0; n < na; ++n)
              meshaux(ma1+n,imb3+ipc3,imb2+ipc2,imb1+ipc1) += w * ps[n];
          }
        }
      }

      // Record the final interpolated properties.
      for (int n = 0; n < ni; ++n)
        pardst(pd1+n,kkk) = pd[n];

      delete [] pd;
      delete [] ps;
    }
  }

  // Release working array.
  u.DeleteAthenaArray();
  for (int i = 0; i < npc1_; ++i)
    delete [] w1[i];
  for (int i = 0; i < npc2_; ++i)
    delete [] w2[i];
  for (int i = 0; i < npc3_; ++i)
    delete [] w3[i];
  delete [] w1;
  delete [] w2;
  delete [] w3;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::DepositMeshAux(AthenaArray<Real>& u,
//!                                       int ma1, int mb1, int nprop)
//! \brief deposits data in meshaux from property index ma1 to ma1+nprop-1 to meshblock
//!        data u from property index mb1 and mb1+nprop-1, divided by cell volume.

void ParticleMesh::DepositMeshAux(AthenaArray<Real>& u, int ma1, int mb1, int nprop) {
  Coordinates *pc = pmb_->pcoord;

#pragma ivdep
  for (int n = 0; n < nprop; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i)
          u(mb1+n,k,j,i) += meshaux(ma1+n,k,j,i) / pc->GetCellVolume(k,j,i);
}

//--------------------------------------------------------------------------------------
//! \fn Real _WeightFunction(Real dxi)
//! \brief evaluates the weight function given index distance.

Real _WeightFunction(Real dxi) {
  dxi = std::min(std::abs(dxi), static_cast<Real>(1.5));
  return dxi < 0.5 ? 0.75 - dxi * dxi : 0.5 * ((1.5 - dxi) * (1.5 - dxi));
}
