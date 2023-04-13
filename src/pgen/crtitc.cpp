//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min

// Athena++ headers
#include "../globals.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../cr/cr.hpp"
#include "../cr/integrators/cr_integrators.hpp"
#include "../thermal_conduction/tc.hpp"
#include <stdio.h>  // fopen and fwrite

// data for the cooling function


static int ninputline=100001;
static int list_len=1048577;



//temperature unit at 2e5K, cool_0=0.1414
//kappa_0=0.0414
//
static Real tunit;
static Real rhounit;
static Real time_unit;
static Real gamma_idx;
static Real totP=1.06;

static AthenaArray<Real> cool_t;
static AthenaArray<Real> cool_coef;
static AthenaArray<Real> cool_index;
static AthenaArray<Real> cool_tef;
static AthenaArray<Real> cell_vol;
static int nline=4;

static Real const_pb=0.01;

static AthenaArray<Real> input_random_pot;


static Real cool_coef1=0.1414;
static Real heat_coef=0.1;
static Real cool_alpha1=-1.0;

static int n_user_var=7;

//======================================================================================
/*! \file beam.cpp
 *  \brief Beam test for the radiative transfer module
 *
 *====================================================================================*/


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================

static Real sigma=1.e8;

static Real pot_B=5.0;
//static Real pot_B=0.0;

void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc);

void TCKappa(MeshBlock *pmb, 
	         AthenaArray<Real> &prim, AthenaArray<Real> &bcc);

void ExactCooling(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

void ExactCooling2(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

void FixCRsourceLeft(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke);
void FixCRsourceRight(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, AthenaArray<Real> &u_cr, Real time, Real dt, 
     int is, int ie, int js, int je, int ks, int ke);



void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // isothermal sound speed for 1.e4 K is 1.16061e6 cm/s
  // length unit is pc=3.086*10^18 cm
  rhounit = pin->GetOrAddReal("problem", "rhounit", 1.25e-25);
  tunit = pin->GetOrAddReal("problem", "Tunit", 2.0e5);
 // time_unit = pin->GetOrAddReal("problem", "Timeunit", 3.3497e11);

  gamma_idx = pin->GetOrAddReal("hydro", "gamma", 5.0/3.0);

  turb_flag = pin->GetInteger("problem","turb_flag");    
 
  EnrollUserExplicitSourceFunction(ExactCooling2);
//  EnrollUserBoundaryFunction(INNER_X1, FixCRsourceLeft);
//  EnrollUserBoundaryFunction(OUTER_X1, FixCRsourceRight);


  // the picesewise power law cooling

  cool_t.NewAthenaArray(nline);

  cool_coef.NewAthenaArray(nline);
  cool_index.NewAthenaArray(nline);
  cool_tef.NewAthenaArray(nline);

  cool_t(0) = 1.e-10;
  cool_t(1) = 0.02;
  cool_t(2) = 20.0;
  cool_t(3) = 1.e8;

  cool_coef(0) = 0.3/pow(0.02,7.0);
  cool_coef(1) = 0.3;
  cool_coef(2) = 0.3/pow(20.0,1.5);

  cool_index(0) = 6.0;
  cool_index(1) = -1.0;
  cool_index(2) = 0.5;

  cool_tef(nline-1) = 0.0;

  Real lambda_t_n=cool_coef(nline-2)*pow(cool_t(nline-1),cool_index(nline-2))
                                                            /cool_t(nline-1);

  for(int i=nline-2; i>=0; i--){
    Real slope = cool_index(i);
    Real coef = cool_coef(i);
    if(fabs(slope-1.0) < TINY_NUMBER){
      cool_tef(i) = cool_tef(i+1) + lambda_t_n*log(cool_t(i+1)/cool_t(i))/coef;
    }else{
      cool_tef(i) = cool_tef(i+1) + lambda_t_n*(pow(cool_t(i+1),1.0-slope) - 
                                     pow(cool_t(i),1.0-slope))/(coef*(1.0-slope));        
    }

  }


  // read in the random vector potential list
/*
  input_random_pot.NewAthenaArray(list_len);
  FILE *finput;
  if ( (finput=fopen("./ini_vect.txt","r"))==NULL )
  {
      printf("Open input file error");
      return;
  }

  for(int i=0; i<list_len; ++i)
  	fscanf(finput,"%lf",&(input_random_pot(i)));

  fclose(finput);
*/


}



void MeshBlock::UserWorkInLoop(void)
{

  
/*
  if(MAGNETIC_FIELDS_ENABLED){

    // initialie vector potential (z component) 
    // Then Bx=partial A/\partial y
    // By=-partial A/\partial x
 // logical block location
 
    int block_i = loc.lx1;
    int block_j = loc.lx2;
    int block_k = loc.lx3;
 

    int totx = pmy_mesh->mesh_size.nx1;
    int toty = pmy_mesh->mesh_size.nx2;
    int totz = pmy_mesh->mesh_size.nx3;


 

    AthenaArray<Real> area, len, len_p1, vec_pot;
    int nz1 = block_size.nx1 + 2*(NGHOST);
    int nz2 = block_size.nx2 + 2*(NGHOST);
    int nz3 = block_size.nx3;
    if(nz3 > 1) nz3 += 2*(NGHOST);

    area.NewAthenaArray(nz1);
    len.NewAthenaArray(nz1);
    len_p1.NewAthenaArray(nz1);

    vec_pot.NewAthenaArray(nz3,nz2,nz1);


    for(int k=0; k<nz3; ++k)
    for(int j=0; j<nz2; ++j)
    for(int i=0; i<nz1; ++i){
        // need to know the relative size of the whole mesh
        int globalk=0;
        if(nz3 > 1){
          globalk=block_k*block_size.nx3+k-ks;
          if(globalk <0)  globalk += totz;
          else if(globalk >  totz-1) globalk -= totz;
        }

        int globalj=block_j*block_size.nx2+j-js;
        if (globalj < 0) globalj += toty;
        else if(globalj >  toty-1) globalj -= toty;

        int globali=block_i*block_size.nx1+i-is;
        if (globali < 0) globali += totx;
        else if(globali > totx - 1) globali -= totx;

        int index = globalk * totx*toty + globalj*totx + globali;

//      vec_pot(k,j,i) = input_random_pot(index);
        vec_pot(k,j,i) = pot_B * phydro->u(IDN,k,j,i); 
    }    


    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        pcoord->Face1Area(k,j,is,ie+1,area);
        pcoord->Edge3Length(k,j  ,is,ie+1,len);
        pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);

        for (int i=is; i<=ie+1; ++i) {
          pfield->b.x1f(k,j,i) += (len_p1(i)*vec_pot(k,j+1,i) -
                                    len(i)*vec_pot(k,j,i))/area(i);
        }
      }
    }

    if(block_size.nx2 > 1){

      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          pcoord->Face2Area(k,j,is,ie,area);
          pcoord->Edge3Length(k,j,is,ie+1,len);         
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) += -(len(i+1)*vec_pot(k,j,i+1)
                               -len(i)*vec_pot(k,j,i))/area(i);
          }
        }
      }

    }

    if(block_size.nx3 > 1){

      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
    }// end nx3

    // set cell centerd magnetic field
    // Add magnetic energy density to the total energy
    pfield->CalculateCellCenteredField(pfield->b,pfield->bcc,pcoord,is,ie,js,je,ks,ke);

    Real local_pb=0.0;
    Real global_pb=0.0;
    int local_vol=0.0;
    int global_vol=0.0;

    //Add the magnetic energy back
    for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
        for(int i=is; i<=ie; ++i){
          phydro->u(IEN,k,j,i) +=
            0.5*(SQR((pfield->bcc(IB1,k,j,i)))
               + SQR((pfield->bcc(IB2,k,j,i)))
               + SQR((pfield->bcc(IB3,k,j,i))));

            local_pb += 0.5*(SQR((pfield->bcc(IB1,k,j,i)))
               + SQR((pfield->bcc(IB2,k,j,i)))
               + SQR((pfield->bcc(IB3,k,j,i))));
            local_vol++;
      
        }
      }
    }

    MPI_Allreduce(&local_pb, &global_pb, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_vol, &global_vol, 1, MPI_INT,MPI_SUM, MPI_COMM_WORLD);

    printf("Average PB: %e\n",global_pb/global_vol);

    area.DeleteAthenaArray();
    len.DeleteAthenaArray();
    len_p1.DeleteAthenaArray();
    vec_pot.DeleteAthenaArray();


  }// end MHD
*/

  return;

}

// Apply over energy balance due to heating and cooling
void Mesh::UserWorkInLoop(void)
{
  // calculate the overall heating and cooling rate
  MeshBlock *pmb = pblock;
  Real tot_heating = 0.0;
  Real tot_cooling = 0.0;
  Real tot_vol = 0.0;
  Real tot_mass = 0.0;
  while(pmb != nullptr){
    int kl = pmb->ks, ku = pmb->ke;
    int jl = pmb->js, ju = pmb->je;
    int il = pmb->is, iu = pmb->ie;

    int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
    cell_vol.NewAthenaArray(ncells1);
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
        pmb->pcoord->CellVolume(k,j,il,iu,cell_vol);
        for(int i=il; i<=iu; ++i){
          tot_cooling += pmb->user_out_var(0,k,j,i) * cell_vol(i);
          tot_heating += pmb->user_out_var(1,k,j,i) * cell_vol(i);
          tot_vol += cell_vol(i);
          tot_mass += pmb->phydro->u(IDN,k,j,i) * cell_vol(i);
        }
      }
    }

    pmb = pmb->next;
  }

  // Now do MPI sum over different cores
#ifdef MPI_PARALLEL

  Real global_cooling=0.0;
  Real global_heating=0.0;
  Real global_vol=0.0;
  Real global_mass = 0.0;
  MPI_Allreduce(&tot_heating, &global_heating, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
  
  MPI_Allreduce(&tot_cooling, &global_cooling, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  
  MPI_Allreduce(&tot_vol, &global_vol, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

  MPI_Allreduce(&tot_mass, &global_mass, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);

  tot_cooling = global_cooling;
  tot_heating = global_heating;
  tot_vol = global_vol;
  tot_mass = global_mass;

#endif  
  //this is the amount of energy we need to add back per unit volume, the same for all cells
  Real diff = (tot_heating + tot_cooling)/tot_mass;
  pmb = pblock;
  while(pmb != nullptr){
    int kl = pmb->ks, ku = pmb->ke;
    int jl = pmb->js, ju = pmb->je;
    int il = pmb->is, iu = pmb->ie;

    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
        for(int i=il; i<=iu; ++i){
          pmb->user_out_var(5,k,j,i) = diff * pmb->phydro->u(IDN,k,j,i);
        }
      }
    }

    pmb = pmb->next;
  }

  // apply this to all the region, including the ghost zones
  pmb=pblock;
  while(pmb != nullptr){
    Real gamma_1 = pmb->peos->GetGamma() - 1.0;
    int kl=0, ku=pmb->ncells3-1;
    int jl=0, ju=pmb->ncells2-1;
    int il=0, iu=pmb->ncells1-1;
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
        for(int i=il; i<=iu; ++i){
          pmb->phydro->u(IEN,k,j,i) -= diff * pmb->phydro->u(IDN,k,j,i);
          pmb->phydro->w(IPR,k,j,i) -= diff* pmb->phydro->u(IDN,k,j,i) * gamma_1;
        }
      }
    }

   // transform conservative to privimative variables
//    Hydro *ph=pmb->phydro;
//    Field *pf=pmb->pfield;
//    pmb->peos->ConservedToPrimitiveCellAverage(ph->u, ph->w, pf->b,
//                                                 ph->w, pf->bcc, pmb->pcoord,
//                                                 il, iu, jl, ju, kl, ku);

    pmb = pmb->next;
  }


 // convert conservitive to primitive variables

}

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{ 

  cool_t.DeleteAthenaArray();
  cool_coef.DeleteAthenaArray();
  cool_index.DeleteAthenaArray();
  cool_tef.DeleteAthenaArray();  

//  input_random_pot.DeleteAthenaArray();

}


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{



  AllocateUserOutputVariables(n_user_var);

  if(CR_ENABLED)
    pcr->EnrollOpacityFunction(Diffusion);

  if(TC_ENABLED){
      ptc->EnrollOpacityFunction(TCKappa);
  }

  // the velocity components will be stored 
  // user_out_var[6],[7],[8]

}



void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  std::srand(gid);
  Real amp = 0.1;

  int kl=ks, ku=ke;
  if(ku > kl){
    ku += NGHOST;
    kl -= NGHOST;
  }

  int jl=js, ju=je;
  if(ju > jl){
    ju += NGHOST;
    jl -= NGHOST;
  }
  int il = is-NGHOST, iu=ie+NGHOST;

  // logical block location
//  int block_i = loc.lx1;
//  int block_j = loc.lx2;
//  int block_k = loc.lx3;

  int totx = pin->GetInteger("mesh", "nx1");
  int toty = pin->GetInteger("mesh", "nx2");
  int totz = pin->GetInteger("mesh", "nx3");

 // Initialize hydro variable

  Real rho = 0.57735;

//  Real tem = 0.3276;
  Real tem = 1.73205;
  Real pc = rho*tem*0.1;

  Real va = sqrt(2.0*const_pb/rho);
  Real vm = 100.0;
  if(CR_ENABLED) vm = pcr->vmax;
  Real fc = va * (pc*3.0 + pc)/vm;


  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for(int i=is; i<=ie; ++i) {

        phydro->u(IDN,k,j,i) = rho * (1.0 + amp 
        	                       * ((double)rand()/(double)RAND_MAX-0.5));
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = tem * rho/(gamma_idx-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
        
        if(CR_ENABLED){
            pcr->u_cr(CRE,k,j,i) = pc*3.0 * (1.0 + amp 
        	                       * ((double)rand()/(double)RAND_MAX-0.5));
            pcr->u_cr(CRF1,k,j,i) = 0.0;
            pcr->u_cr(CRF2,k,j,i) = 0.0;
            pcr->u_cr(CRF3,k,j,i) = 0.0;
        }
      }// end i
    }
  }
  //Need to set opactiy sigma in the ghost zones
  if(CR_ENABLED){

  // Default values are 1/3
    int nz1 = block_size.nx1 + 2*(NGHOST);
    int nz2 = block_size.nx2;
    if(nz2 > 1) nz2 += 2*(NGHOST);
    int nz3 = block_size.nx3;
    if(nz3 > 1) nz3 += 2*(NGHOST);
    for(int k=0; k<nz3; ++k){
      for(int j=0; j<nz2; ++j){
        for(int i=0; i<nz1; ++i){
          pcr->sigma_diff(0,k,j,i) = sigma;
          pcr->sigma_diff(1,k,j,i) = sigma;
          pcr->sigma_diff(2,k,j,i) = sigma;
        }
      }
    }// end k,j,i

  }// End CR



  if(MAGNETIC_FIELDS_ENABLED){

    // initialie vector potential (z component) 
    // Then Bx=partial A/\partial y
    // By=-partial A/\partial x


    AthenaArray<Real> area, len, len_p1, vec_pot;
    int nz1 = block_size.nx1 + 2*(NGHOST);
    int nz2 = block_size.nx2 + 2*(NGHOST);
    int nz3 = block_size.nx3;
    if(nz3 > 1) nz3 += 2*(NGHOST);

    area.NewAthenaArray(nz1);
    len.NewAthenaArray(nz1);
    len_p1.NewAthenaArray(nz1);

    vec_pot.NewAthenaArray(nz3,nz2,nz1);


    for(int k=0; k<nz3; ++k)
    for(int j=0; j<nz2; ++j)
    for(int i=0; i<nz1; ++i){
        // need to know the relative size of the whole mesh
//        int globalk=0;
//        if(nz3 > 1){
//        	globalk=block_k*block_size.nx3+k-ks;
//        	if(globalk <0)  globalk += totz;
//        	else if(globalk >  totz-1) globalk -= totz;
//        }

//        int globalj=block_j*block_size.nx2+j-js;
//        if (globalj < 0) globalj += toty;
//        else if(globalj >  toty-1) globalj -= toty;

//       int globali=block_i*block_size.nx1+i-is;
//        if (globali < 0) globali += totx;
//        else if(globali > totx - 1) globali -= totx;

//        int index = globalk * totx*toty + globalj*totx + globali;

//    	vec_pot(k,j,i) = input_random_pot(index);
       vec_pot(k,j,i) = 0.0;
    }    


    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        pcoord->Face1Area(k,j,is,ie+1,area);
        pcoord->Edge3Length(k,j  ,is,ie+1,len);
        pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);

        for (int i=is; i<=ie+1; ++i) {
//          pfield->b.x1f(k,j,i) = pot_B*(len_p1(i)*vec_pot(k,j+1,i) -
//                                    len(i)*vec_pot(k,j,i))/area(i);
          pfield->b.x1f(k,j,i) = sqrt(2.0*const_pb);
        }
      }
    }

    if(block_size.nx2 > 1){

      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je+1; ++j) {
          pcoord->Face2Area(k,j,is,ie,area);
          pcoord->Edge3Length(k,j,is,ie+1,len);         
          for (int i=is; i<=ie; ++i) {
            pfield->b.x2f(k,j,i) = -0.0 * pot_B*(len(i+1)*vec_pot(k,j,i+1)
                               -0.0 * len(i)*vec_pot(k,j,i))/area(i);
          }
        }
      }

    }

    if(block_size.nx3 > 1){

      for (int k=ks; k<=ke+1; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
            pfield->b.x3f(k,j,i) = 0.0;
          }
        }
      }
    }// end nx3

    // set cell centerd magnetic field
    // Add magnetic energy density to the total energy
    pfield->CalculateCellCenteredField(pfield->b,pfield->bcc,pcoord,is,ie,js,je,ks,ke);

    Real local_pb=0.0;
    Real global_pb=0.0;
    int local_vol=0.0;
    int global_vol=0.0;

    for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
        for(int i=is; i<=ie; ++i){
          phydro->u(IEN,k,j,i) +=
            0.5*(SQR((pfield->bcc(IB1,k,j,i)))
               + SQR((pfield->bcc(IB2,k,j,i)))
               + SQR((pfield->bcc(IB3,k,j,i))));

//            local_pb += 0.5*(SQR((pfield->bcc(IB1,k,j,i)))
//               + SQR((pfield->bcc(IB2,k,j,i)))
//               + SQR((pfield->bcc(IB3,k,j,i))));
//            local_vol++;
      
        }
      }
    }

//    MPI_Allreduce(&local_pb, &global_pb, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
//    MPI_Allreduce(&local_vol, &global_vol, 1, MPI_INT,MPI_SUM, MPI_COMM_WORLD);

//    printf("Average PB: %e\n",global_pb/global_vol);

    area.DeleteAthenaArray();
    len.DeleteAthenaArray();
    len_p1.DeleteAthenaArray();
    vec_pot.DeleteAthenaArray();

  }// end MHD
  



  
  return;
}



void TCKappa(MeshBlock *pmb, AthenaArray<Real> &prim, AthenaArray<Real> &bcc)
{ 
  // set the default opacity to be a large value in the default hydro case
  ThermalConduction *ptc=pmb->ptc;
  int kl=pmb->ks, ku=pmb->ke;
  int jl=pmb->js, ju=pmb->je;
  int il=pmb->is-1, iu=pmb->ie+1;
  if(pmb->block_size.nx2 > 1){
    jl -= 1;
    ju += 1;
  }
  if(pmb->block_size.nx3 > 1){
    kl -= 1;
    ku += 1;
  }

  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
#pragma omp simd
      for(int i=il; i<=iu; ++i){
        Real tgas = prim(IEN,k,j,i)/prim(IDN,k,j,i);
        Real kappa_parallel = 0.01;
        Real kappa_pernd = 1.e-8;
//        kappa_parallel = std::max(kappa_parallel, 10.0*kappa_pernd);

        ptc->kappa(0,k,j,i) = kappa_parallel;
        ptc->kappa(1,k,j,i) = kappa_pernd;
        ptc->kappa(2,k,j,i) = kappa_pernd;

      }
    }
  }

  if(MAGNETIC_FIELDS_ENABLED){
    //First, calculate B_dot_grad_Pc
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
      // diffusion coefficient is calculated with respect to B direction
      // Use a simple estimate of Grad Pc
        for(int i=il; i<=iu; ++i){
          // Now calculate the angles of B
          Real bxby = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i));
          Real btot = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i) +
                           bcc(IB3,k,j,i)*bcc(IB3,k,j,i));

          if(btot > TINY_NUMBER){
            ptc->b_angle(0,k,j,i) = bxby/btot;
            ptc->b_angle(1,k,j,i) = bcc(IB3,k,j,i)/btot;
          }else{
            ptc->b_angle(0,k,j,i) = 1.0;
            ptc->b_angle(1,k,j,i) = 0.0;
          }
          if(bxby > TINY_NUMBER){
            ptc->b_angle(2,k,j,i) = bcc(IB2,k,j,i)/bxby;
            ptc->b_angle(3,k,j,i) = bcc(IB1,k,j,i)/bxby;
          }else{
            ptc->b_angle(2,k,j,i) = 0.0;
            ptc->b_angle(3,k,j,i) = 1.0;
          }


        }//end i        

      }// end j
    }// end k
  }// end MHD



}


void Diffusion(MeshBlock *pmb, AthenaArray<Real> &u_cr, 
        AthenaArray<Real> &prim, AthenaArray<Real> &bcc)
{ 


  // set the default opacity to be a large value in the default hydro case
  CosmicRay *pcr=pmb->pcr;
  int kl=pmb->ks, ku=pmb->ke;
  int jl=pmb->js, ju=pmb->je;
  int il=pmb->is-1, iu=pmb->ie+1;
  if(pmb->block_size.nx2 > 1){
    jl -= 1;
    ju += 1;
  }
  if(pmb->block_size.nx3 > 1){
    kl -= 1;
    ku += 1;
  }

  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
#pragma omp simd
      for(int i=il; i<=iu; ++i){

        pcr->sigma_diff(0,k,j,i) = sigma;
        pcr->sigma_diff(1,k,j,i) = sigma;
        pcr->sigma_diff(2,k,j,i) = sigma;  

      }
    }
  }

  Real invlim=1.0/pcr->vmax;

  // The information stored in the array
  // b_angle is
  // b_angle[0]=sin_theta_b
  // b_angle[1]=cos_theta_b
  // b_angle[2]=sin_phi_b
  // b_angle[3]=cos_phi_b




  if(MAGNETIC_FIELDS_ENABLED){
    //First, calculate B_dot_grad_Pc
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
    // x component
        pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
        for(int i=il; i<=iu; ++i){
          Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                         + pcr->cwidth(i);
          Real dprdx=(u_cr(CRE,k,j,i+1) -  u_cr(CRE,k,j,i-1))/3.0;
          dprdx /= distance;
          pcr->b_grad_pc(k,j,i) = bcc(IB1,k,j,i) * dprdx;
        }
    //y component
        if(pmb->block_size.nx2 > 1){
          pmb->pcoord->CenterWidth2(k,j-1,il,iu,pcr->cwidth1);       
          pmb->pcoord->CenterWidth2(k,j,il,iu,pcr->cwidth);
          pmb->pcoord->CenterWidth2(k,j+1,il,iu,pcr->cwidth2);

          for(int i=il; i<=iu; ++i){
            Real distance = 0.5*(pcr->cwidth1(i) + pcr->cwidth2(i))
                           + pcr->cwidth(i);
            Real dprdy=(u_cr(CRE,k,j+1,i) - u_cr(CRE,k,j-1,i))/3.0;
            dprdy /= distance;
            pcr->b_grad_pc(k,j,i) += bcc(IB2,k,j,i) * dprdy;

          }
        }
    // z component
        if(pmb->block_size.nx3 > 1){
          pmb->pcoord->CenterWidth3(k-1,j,il,iu,pcr->cwidth1);       
          pmb->pcoord->CenterWidth3(k,j,il,iu,pcr->cwidth);
          pmb->pcoord->CenterWidth3(k+1,j,il,iu,pcr->cwidth2);

          for(int i=il; i<=iu; ++i){
            Real distance = 0.5*(pcr->cwidth1(i) + pcr->cwidth2(i))
                            + pcr->cwidth(i);
            Real dprdz=(u_cr(CRE,k+1,j,i) -  u_cr(CRE,k-1,j,i))/3.0;
            dprdz /= distance;
            pcr->b_grad_pc(k,j,i) += bcc(IB3,k,j,i) * dprdz;

            // now only get the sign
  //          if(pcr->b_grad_pc(k,j,i) > TINY_NUMBER) pcr->b_grad_pc(k,j,i) = 1.0;
  //          else if(-pcr->b_grad_pc(k,j,i) > TINY_NUMBER) pcr->b_grad_pc(k,j,i) 
  //            = -1.0;
  //          else pcr->b_grad_pc(k,j,i) = 0.0;
          }
        }
      // now calculate the streaming velocity
      // streaming velocity is calculated with respect to the current coordinate 
      //  system
      // diffusion coefficient is calculated with respect to B direction
        for(int i=il; i<=iu; ++i){
          Real pb= bcc(IB1,k,j,i)*bcc(IB1,k,j,i)
                  +bcc(IB2,k,j,i)*bcc(IB2,k,j,i)
                  +bcc(IB3,k,j,i)*bcc(IB3,k,j,i);
          Real inv_sqrt_rho = 1.0/sqrt(prim(IDN,k,j,i));
          Real va1 = bcc(IB1,k,j,i)*inv_sqrt_rho;
          Real va2 = bcc(IB2,k,j,i)*inv_sqrt_rho;
          Real va3 = bcc(IB3,k,j,i)*inv_sqrt_rho;

          Real va = sqrt(pb/prim(IDN,k,j,i));

          Real dpc_sign = 0.0;
          if(pcr->b_grad_pc(k,j,i) > TINY_NUMBER) dpc_sign = 1.0;
          else if(-pcr->b_grad_pc(k,j,i) > TINY_NUMBER) dpc_sign = -1.0;
          
          pcr->v_adv(0,k,j,i) = -va1 * dpc_sign;
          pcr->v_adv(1,k,j,i) = -va2 * dpc_sign;
          pcr->v_adv(2,k,j,i) = -va3 * dpc_sign;

          // now the diffusion coefficient

          if(va < TINY_NUMBER){
            pcr->sigma_adv(0,k,j,i) = pcr->max_opacity;
          }else{
            pcr->sigma_adv(0,k,j,i) = fabs(pcr->b_grad_pc(k,j,i))
                                   /(sqrt(pb) * va * (4.0/3.0) 
                                    * invlim * u_cr(CRE,k,j,i)); 
          }

          pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
          pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;  

          // Now calculate the angles of B
          Real bxby = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) +
                           bcc(IB2,k,j,i)*bcc(IB2,k,j,i));
          Real btot = sqrt(pb);
          if(btot > TINY_NUMBER){
            pcr->b_angle(0,k,j,i) = bxby/btot;
            pcr->b_angle(1,k,j,i) = bcc(IB3,k,j,i)/btot;
          }else{
            pcr->b_angle(0,k,j,i) = 1.0;
            pcr->b_angle(1,k,j,i) = 0.0;
          }
          if(bxby > TINY_NUMBER){
            pcr->b_angle(2,k,j,i) = bcc(IB2,k,j,i)/bxby;
            pcr->b_angle(3,k,j,i) = bcc(IB1,k,j,i)/bxby;
          }else{
            pcr->b_angle(2,k,j,i) = 0.0;
            pcr->b_angle(3,k,j,i) = 1.0;            
          }

        }//        

      }// end j
    }// end k

  }// End MHD  
  else{



    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
    // x component
        pmb->pcoord->CenterWidth1(k,j,il-1,iu+1,pcr->cwidth);
        for(int i=il; i<=iu; ++i){
           Real distance = 0.5*(pcr->cwidth(i-1) + pcr->cwidth(i+1))
                          + pcr->cwidth(i);
           Real grad_pr=(u_cr(CRE,k,j,i+1) -  u_cr(CRE,k,j,i-1))/3.0;
           grad_pr /= distance;

           Real va = 0.0;

           if(va < TINY_NUMBER){
             pcr->sigma_adv(0,k,j,i) = pcr->max_opacity;
             pcr->v_adv(0,k,j,i) = 0.0;
           }else{
             Real sigma2 = fabs(grad_pr)/(va * (4.0/3.0) 
                               * invlim * u_cr(CRE,k,j,i)); 
             if(fabs(grad_pr) < TINY_NUMBER){
               pcr->sigma_adv(0,k,j,i) = 0.0;
               pcr->v_adv(0,k,j,i) = 0.0;
             }else{
               pcr->sigma_adv(0,k,j,i) = sigma2;
               pcr->v_adv(0,k,j,i) = -va * grad_pr/fabs(grad_pr);     
             }
          }

          pcr->sigma_adv(1,k,j,i) = pcr->max_opacity;
          pcr->sigma_adv(2,k,j,i) = pcr->max_opacity;
         

          pcr->v_adv(1,k,j,i) = 0.0;
          pcr->v_adv(2,k,j,i) = 0.0;




        }

      }
    }

  }
}



void ExactCooling(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

    // scale the temperature unit and cooling time
    // so the cooling will be just 
    // dT/dt =  -cool_coef *rho * T^alpha
  // for the current unit, cool_coef=1.4142


  int kl=pmb->ks, ku=pmb->ke;
  int jl=pmb->js, ju=pmb->je;
  int il=pmb->is, iu=pmb->ie;



  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
#pragma omp simd
      for(int i=il; i<=iu; ++i){
        // Need to take density and temperature at time step n from 
        // cons, not from prim
        // Because we do not need intermediate step to calculate the cooling 
        // function
        Real rho = cons(IDN,k,j,i);
        Real eint = cons(IEN,k,j,i) 
                    - 0.5 *(cons(IM1,k,j,i)*cons(IM1,k,j,i)
                      + cons(IM2,k,j,i)*cons(IM2,k,j,i)
                      + cons(IM3,k,j,i)*cons(IM3,k,j,i))/rho;
        if(MAGNETIC_FIELDS_ENABLED){
             eint -= 0.5 * (bcc(IB1,k,j,i) * bcc(IB1,k,j,i)
                      + bcc(IB2,k,j,i) * bcc(IB2,k,j,i)
                      + bcc(IB3,k,j,i) * bcc(IB3,k,j,i));
        }
        Real t_i = eint *(gamma_idx - 1.0)/rho;
          Real newt= pow(t_i,1.0-cool_alpha1)-(1.0-cool_alpha1)*rho*cool_coef1*dt;
          newt = pow(newt,1.0/(1.0-cool_alpha1));

//        Real newt2 = t_i - cool_coef * rho * pow(t_i,cool_alpha) * dt;
//        newt = std::max(newt2,newt);

          cons(IEN,k,j,i) += (newt - t_i) * rho/(gamma_idx - 1.0); 

          pmb->user_out_var(0,k,j,i) = (newt - t_i) * rho/(gamma_idx - 1.0);
 
 //       pmb->user_out_var(0,k,j,i) = ((newt - t_i) * rho/(gamma_idx - 1.0))/dt;
 //       pmb->user_out_var(1,k,j,i) = -cool_coef1 * rho * pow(t_i, cool_alpha1) * rho/(gamma_idx - 1.0);  
        // add a constant heating rate
          cons(IEN,k,j,i) += dt * heat_coef/(gamma_idx-1.0);     
          pmb->user_out_var(1,k,j,i) = dt * heat_coef/(gamma_idx-1.0);
      }     
    }
  }

}




void ExactCooling2(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

    // scale the temperature unit and cooling time
    // so the cooling will be just 
    // dT/dt = coef() Lambda(T)


  int kl=pmb->ks, ku=pmb->ke;
  int jl=pmb->js, ju=pmb->je;
  int il=pmb->is, iu=pmb->ie;


  Real lambda_t_n=cool_coef(nline-2)*pow(cool_t(nline-1),cool_index(nline-2))
                                                            /cool_t(nline-1);


  Real ncells1=pmb->block_size.nx1 + 2*(NGHOST);
    
//  cell_vol.NewAthenaArray(ncells1);

    
//  Real tot_heating = 0.0;
//  Real tot_cooling = 0.0;
//  Real tot_vol=0.0;

    
  for(int k=kl; k<=ku; ++k){
    for(int j=jl; j<=ju; ++j){
//      pmb->pcoord->CellVolume(k,j,il,iu,cell_vol);
      for(int i=il; i<=iu; ++i){
        // Need to take density and temperature at time step n from 
        // cons, not from prim
        // Because we do not need intermediate step to calculate the cooling 
        // function
        Real rho = cons(IDN,k,j,i);
        Real eint = cons(IEN,k,j,i) 
                    - 0.5 *(cons(IM1,k,j,i)*cons(IM1,k,j,i)
                      + cons(IM2,k,j,i)*cons(IM2,k,j,i)
                      + cons(IM3,k,j,i)*cons(IM3,k,j,i))/rho;
        if(MAGNETIC_FIELDS_ENABLED){
             eint -= 0.5 * (bcc(IB1,k,j,i) * bcc(IB1,k,j,i)
                      + bcc(IB2,k,j,i) * bcc(IB2,k,j,i)
                      + bcc(IB3,k,j,i) * bcc(IB3,k,j,i));
        }
        Real t_i = eint *(gamma_idx - 1.0)/rho;


        if(t_i > cool_t(0)){

          if(t_i > cool_t(nline-1)) t_i = cool_t(nline-1);

          int t_loc=0;
          while((t_loc < nline-2) && (cool_t(t_loc+1) < t_i) ){
            ++t_loc;
          }

          Real slope = cool_index(t_loc);
          Real coef = cool_coef(t_loc);

          Real tef = cool_tef(t_loc+1);
          if(fabs(slope-1.0) < TINY_NUMBER){
            tef += lambda_t_n*log(cool_t(t_loc+1)/t_i)/coef;
          }else{
            tef += lambda_t_n*(pow(cool_t(t_loc+1),1.0-slope) - 
                            pow(t_i,1.0-slope))/(coef*(1.0-slope));        
          }

          Real new_tef = tef + rho * dt * lambda_t_n;
          // Now invert TEF to get the current temperature
          // new_tef > tef
          int tef_loc=t_loc+1;
          while((tef_loc > 0) && (new_tef > cool_tef(tef_loc))){
            --tef_loc;
          }

          Real diff_tef = (new_tef - cool_tef(tef_loc+1))/lambda_t_n;
          slope = cool_index(tef_loc);
          coef = cool_coef(tef_loc);

          Real tnew = t_i;
          if(fabs(slope-1.0) < TINY_NUMBER){
            tnew = exp(log(cool_t(tef_loc+1))-(coef*diff_tef));
          }else{
            tnew = pow(cool_t(tef_loc+1),1.0-slope) 
                               - (1.0-slope) * coef * diff_tef;
            tnew = pow(tnew,1.0/(1.0-slope));
          }



          cons(IEN,k,j,i) += (tnew - t_i) * rho/(gamma_idx - 1.0);  

          pmb->user_out_var(0,k,j,i) = (tnew - t_i) * rho/(gamma_idx - 1.0);
            
//          tot_cooling += pmb->user_out_var(0,k,j,i) * cell_vol(i);

          //add a constant heating rate
          cons(IEN,k,j,i) += dt * rho* heat_coef/(gamma_idx-1.0);     

          pmb->user_out_var(1,k,j,i) = dt * rho * heat_coef/(gamma_idx-1.0);
//            tot_heating += cell_vol(i) * pmb->user_out_var(1,k,j,i);
              
        }else{
            pmb->user_out_var(0,k,j,i) = 0.0;
            pmb->user_out_var(1,k,j,i) = 0.0;          	
        }
            
//          tot_vol += cell_vol(i);

        
      }
    }
  }
 /*   
    Real global_cooling=0.0;
    Real global_heating=0.0;
    Real global_vol=0.0;
    MPI_Allreduce(&tot_heating, &global_heating, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
    
    MPI_Allreduce(&tot_cooling, &global_cooling, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    MPI_Allreduce(&tot_vol, &global_vol, 1, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
    
    //this is the amount of energy we need to add back per unit volume, the same for all cells
    Real diff = -(global_heating + global_cooling)/global_vol;
    
    
    for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
        for(int i=il; i<=iu; ++i){
            cons(IEN,k,j,i) += diff;
        }
      }
    }

*/
//    cell_vol.DeleteAthenaArray();

}




