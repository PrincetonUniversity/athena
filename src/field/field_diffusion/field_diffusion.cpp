//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//  \brief Class to implement diffusion processes in the induction equations

// Athena++ headers
#include "field_diffusion.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../field.hpp"
#include "../../hydro/hydro.hpp"
#include "../../parameter_input.hpp"

// FieldDiffusion constructor

FieldDiffusion::FieldDiffusion(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;
  field_diffusion_defined = false;

  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

  // Check if field diffusion
  coeff_o = pin->GetOrAddReal("problem","coef_o",0.0);
  coeff_h = pin->GetOrAddReal("problem","coef_h",0.0);
  coeff_a = pin->GetOrAddReal("problem","coef_a",0.0);

  if ((coeff_o != 0.0) || (coeff_h != 0.0) || (coeff_a != 0.0)) {
    field_diffusion_defined = true;
    // Allocate memory for scratch vectors
    etaB.NewAthenaArray(3,ncells3,ncells2,ncells1);
    e_oa.x1e.NewAthenaArray(ncells3+1,ncells2+1,ncells1);
    e_oa.x2e.NewAthenaArray(ncells3+1,ncells2,ncells1+1);
    e_oa.x3e.NewAthenaArray(ncells3,ncells2+1,ncells1+1);
    e_h.x1e.NewAthenaArray(ncells3+1,ncells2+1,ncells1);
    e_h.x2e.NewAthenaArray(ncells3+1,ncells2,ncells1+1);
    e_h.x3e.NewAthenaArray(ncells3,ncells2+1,ncells1+1);
    pflux.x1f.NewAthenaArray(ncells3,ncells2,ncells1+1);
    pflux.x2f.NewAthenaArray(ncells3,ncells2+1,ncells1);
    pflux.x3f.NewAthenaArray(ncells3+1,ncells2,ncells1);

    jfx.NewAthenaArray(3,ncells3+1,ncells2+1,ncells1+1);
    jfy.NewAthenaArray(3,ncells3+1,ncells2+1,ncells1+1);
    jfz.NewAthenaArray(3,ncells3+1,ncells2+1,ncells1+1);
    jcc.NewAthenaArray(3,ncells3+1,ncells2+1,ncells1+1);

    eta_tot_.NewAthenaArray(nthreads,ncells1);
    bmag_.NewAthenaArray(ncells3,ncells2,ncells1);

    jedge_.x1e.NewAthenaArray(ncells3+1,ncells2+1,ncells1);
    jedge_.x2e.NewAthenaArray(ncells3+1,ncells2,ncells1+1);
    jedge_.x3e.NewAthenaArray(ncells3,ncells2+1,ncells1+1);

    face_area_.NewAthenaArray(nthreads,ncells1);
    face_area_p1_.NewAthenaArray(nthreads,ncells1);
    edge_length_.NewAthenaArray(nthreads,ncells1);
    edge_length_m1_.NewAthenaArray(nthreads,ncells1);
    cell_volume_.NewAthenaArray(nthreads,ncells1);
    dx1_.NewAthenaArray(nthreads,ncells1);
    dx2_.NewAthenaArray(nthreads,ncells1);
    dx3_.NewAthenaArray(nthreads,ncells1);

    if(pmb->pmy_mesh->FieldDiffusivity_==NULL)
      CalcMagDiffCoeff_ = ConstDiffusivity;
    else
      CalcMagDiffCoeff_ = pmb->pmy_mesh->FieldDiffusivity_;
  }
}

// destructor

FieldDiffusion::~FieldDiffusion()
{
  if ((coeff_o != 0.0) || (coeff_h != 0.0) || (coeff_a != 0.0)) {
    etaB.DeleteAthenaArray();
    e_oa.x1e.DeleteAthenaArray();
    e_oa.x2e.DeleteAthenaArray();
    e_oa.x3e.DeleteAthenaArray();
    e_h.x1e.DeleteAthenaArray();
    e_h.x2e.DeleteAthenaArray();
    e_h.x3e.DeleteAthenaArray();
    pflux.x1f.DeleteAthenaArray();
    pflux.x2f.DeleteAthenaArray();
    pflux.x3f.DeleteAthenaArray();

    jfx.DeleteAthenaArray();
    jfy.DeleteAthenaArray();
    jfz.DeleteAthenaArray();
    jcc.DeleteAthenaArray();

    eta_tot_.DeleteAthenaArray();
    bmag_.DeleteAthenaArray();

    jedge_.x1e.DeleteAthenaArray();
    jedge_.x2e.DeleteAthenaArray();
    jedge_.x3e.DeleteAthenaArray();

    face_area_.DeleteAthenaArray();
    face_area_p1_.DeleteAthenaArray();
    edge_length_.DeleteAthenaArray();
    edge_length_m1_.DeleteAthenaArray();
    cell_volume_.DeleteAthenaArray();
    dx1_.DeleteAthenaArray();
    dx2_.DeleteAthenaArray();
    dx3_.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::CalcFieldDiffusionEMF
//  \brief Calculate diffusion EMF(Ohmic & Ambipolar for now)

void FieldDiffusion::CalcFieldDiffusionEMF(FaceField &bi,
     const AthenaArray<Real> &bc, EdgeField &e)
{
  Field *pf = pmy_block->pfield;
  Hydro *ph = pmy_block->phydro;
  Mesh  *pm = pmy_block->pmy_mesh;

  if((coeff_o==0.0) && (coeff_a==0.0)) return;

  SetFieldDiffusivity(ph->w,pf->bcc);

  //CalcCurrent(pfield->b);
  CalcCurrent(bi);
  ClearEMF(e_oa);
  if (coeff_o != 0.0)
    //OhmicEMF(pfield->b, pfield->bcc, e_oa);
    OhmicEMF(bi, bc, e_oa);
  if (coeff_a != 0.0)
    //AmbipolarEMF(pfield->b, pfield->bcc, e_oa);
    AmbipolarEMF(bi, bc, e_oa);

  // calculate the Poynting flux pflux and add to energy flux in Hydro class
  if (NON_BAROTROPIC_EOS) {
    //PoyntingFlux(e_oa, pfield->bcc);
    PoyntingFlux(e_oa, bc);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::AddEMF(EdgeField &e_src, EdgeField &e_des)
//  \brief Add source EMF to destination EMF

void FieldDiffusion::AddEMF(const EdgeField &e_src, EdgeField &e_des)
{
  int size1 = e_src.x1e.GetSize();
  int size2 = e_src.x2e.GetSize();
  int size3 = e_src.x3e.GetSize();

#pragma simd
  for (int i=0; i<size1; ++i)
    e_des.x1e(i) += e_src.x1e(i);

#pragma simd
  for (int i=0; i<size2; ++i)
    e_des.x2e(i) += e_src.x2e(i);

#pragma simd
  for (int i=0; i<size3; ++i)
    e_des.x3e(i) += e_src.x3e(i);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FieldDiffusion::ClearFieldDiffusionEMF(EdgeField &e)
//  \brief Clear EMF

void FieldDiffusion::ClearEMF(EdgeField &e)
{
  int size1 = e.x1e.GetSize();
  int size2 = e.x2e.GetSize();
  int size3 = e.x3e.GetSize();

#pragma simd
  for (int i=0; i<size1; ++i)
    e.x1e(i) = 0.0;

#pragma simd
  for (int i=0; i<size2; ++i)
    e.x2e(i) = 0.0;

#pragma simd
  for (int i=0; i<size3; ++i)
    e.x3e(i) = 0.0;

  return;
}


//--------------------------------------------------------------------------------------
// Set magnetic diffusion coefficients

void FieldDiffusion::SetFieldDiffusivity(const AthenaArray<Real> &w, const AthenaArray<Real> &bc)
{
  MeshBlock *pmb = pmy_block;
  int il = pmb->is-NGHOST; int jl = pmb->js; int kl = pmb->ks;
  int iu = pmb->ie+NGHOST; int ju = pmb->je; int ku = pmb->ke;
  if (pmb->block_size.nx2 > 1) {
    jl -= NGHOST; ju += NGHOST;
  }
  if (pmb->block_size.nx3 > 1) {
    kl -= NGHOST; ku += NGHOST;
  }

  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{

  for (int k=kl; k<=ku; ++k){
#pragma omp for schedule(static)
    for (int j=jl; j<=ju; ++j){
#pragma simd
      for (int i=il; i<=iu; ++i){
        Real Bsq = SQR(bc(IB1,k,j,i)) + SQR(bc(IB2,k,j,i)) + SQR(bc(IB3,k,j,i));
        bmag_(k,j,i) = sqrt(Bsq);
      }
    }
  }
  // set diffusivities
  CalcMagDiffCoeff_(this, w, bmag_, il, iu, jl, ju, kl, ku);

} // end of omp parallel region

  return;
}
