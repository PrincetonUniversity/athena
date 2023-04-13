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
//! \file rad_source.cpp
//  \brief Add radiation source terms to both radiation and gas
//======================================================================================




#include <algorithm>
#include <string>
#include <vector>

//Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../radiation.hpp"
#include "../implicit/radiation_implicit.hpp"
#include "../../coordinates/coordinates.hpp" //
#include "../../hydro/hydro.hpp"
#include "../../field/field.hpp"
#include "../../eos/eos.hpp"

// class header
#include "rad_integrators.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


void RadIntegrator::CalSourceTerms(MeshBlock *pmb, const Real dt, 
        AthenaArray<Real> &u, 
        AthenaArray<Real> &ir_ini, AthenaArray<Real> &ir)
{
  Radiation *prad=pmb->prad;
  Coordinates *pco = pmb->pcoord;
  
  Real& prat = prad->prat;
  Real invcrat = 1.0/prad->crat;

  
  Real *sigma_at, *sigma_aer, *sigma_s, *sigma_p;
  Real *lab_ir;
  
  int &nang =prad->nang;
  int &nfreq=prad->nfreq;
  
  
  // Get the temporary array
  AthenaArray<Real> &wmu_cm = wmu_cm_;
  AthenaArray<Real> &tran_coef = tran_coef_;
  AthenaArray<Real> &ir_cm = ir_cm_;
  AthenaArray<Real> &cm_to_lab = cm_to_lab_;


  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
 
  for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=is; i<=ie; ++i){

        // for implicit update, using the quantities from the partially 
        // updated u, not from w 

        Real rho = u(IDN,k,j,i);
        Real vx = vel_source_(k,j,i,0);
        Real vy = vel_source_(k,j,i,1);
        Real vz = vel_source_(k,j,i,2);
        Real vel = vx*vx + vy*vy + vz*vz;
        
        
        Real lorzsq = 1.0/(1.0 - vel  * invcrat * invcrat);
        Real lorz = sqrt(lorzsq);
        
        
        sigma_at=&(prad->sigma_a(k,j,i,0));
        sigma_aer=&(prad->sigma_ae(k,j,i,0));
        sigma_s=&(prad->sigma_s(k,j,i,0));
        sigma_p=&(prad->sigma_planck(k,j,i,0));

         // Prepare the transformation coefficients
        Real numsum = 0.0;

#pragma omp simd reduction(+:numsum)
        for(int n=0; n<nang; ++n){
           Real vdotn = vx * prad->mu(0,k,j,i,n) + vy * prad->mu(1,k,j,i,n)
                        + vz * prad->mu(2,k,j,i,n);
           Real vnc = 1.0 - vdotn * invcrat;
           tran_coef(n) = lorz * vnc;
           wmu_cm(n) = prad->wmu(n)/(tran_coef(n) * tran_coef(n));
           numsum += wmu_cm(n);
           cm_to_lab(n) = tran_coef(n)*tran_coef(n)*tran_coef(n)*tran_coef(n);
           
        }
           // Normalize weight in co-moving frame to make sure the sum is one
        numsum = 1.0/numsum;
#pragma omp simd
        for(int n=0; n<nang; ++n){
           wmu_cm(n) *= numsum;
        }
                
        for(int ifr=0; ifr<nfreq; ++ifr){
          lab_ir=&(ir_ini(k,j,i,ifr*nang));
          for(int n=0; n<nang; ++n)
            ir_cm(n+ifr*nang) = lab_ir[n];
          if(IM_RADIATION_ENABLED){
            // add the explicit flux divergence term
            Real *divflux = &(divflx_(k,j,i,ifr*nang));
            for(int n=0; n<nang; ++n)
              ir_cm(n+ifr*nang) += divflux[n];
            // for Gauss-Seidel iteration, add flux from left hand side
            if(pmb->pmy_mesh->pimrad->ite_scheme == 1){
              for(int n=0; n<nang; ++n)
                ir_cm(n+ifr*nang) -= left_coef1_(k,j,i,n+ifr*nang) 
                                   * ir(k,j,i-1,n+ifr*nang);
              if(je > js){
                for(int n=0; n<nang; ++n)
                  ir_cm(n+ifr*nang) -= left_coef2_(k,j,i,n+ifr*nang) 
                                     * ir(k,j-1,i,n+ifr*nang);                
              }
              if(ke > ks){
                for(int n=0; n<nang; ++n)
                  ir_cm(n+ifr*nang) -= left_coef3_(k,j,i,n+ifr*nang) 
                                     * ir(k-1,j,i,n+ifr*nang);                 
              }

            }
            // add angular flux
            if(prad->angle_flag == 1){
              Real *p_angflx  = &(ang_flx_(k,j,i,ifr*nang));
              for(int n=0; n<nang; ++n)
                ir_cm(n+ifr*nang) += p_angflx[n];
            }// end angle_flag == 1

          }
#pragma omp simd
          for(int n=0; n<nang; ++n){
            ir_cm(n+ifr*nang) *= cm_to_lab(n);
          }
          // apply floor to ir_cm
#pragma omp simd
          for(int n=0; n<nang; ++n){
            ir_cm(n+ifr*nang) = std::max(ir_cm(n+ifr*nang),TINY_NUMBER);
          }
    
          for(int n=0; n<nang; n++){
            implicit_coef_(n+ifr*nang) = 1.0;
            if(IM_RADIATION_ENABLED){
              implicit_coef_(n+ifr*nang) += const_coef_(k,j,i,n+ifr*nang);
              if(prad->angle_flag == 1){
                implicit_coef_(n+ifr*nang) += imp_ang_coef_(k,j,i,n);
              }
            }
          }// end ang

        }// End frequency
        
        if(nfreq == 1){

         // Add absorption and scattering opacity source
          tgas_new_(k,j,i) = AbsorptionScattering(wmu_cm,tran_coef, sigma_at, sigma_p, 
                             sigma_aer, sigma_s, dt, lorz, rho, tgas_(k,j,i), 
                             implicit_coef_,ir_cm);
        
         // Add compton scattering
          if(compton_flag_ > 0)
            Compton(wmu_cm,tran_coef, sigma_s, dt, lorz, rho, tgas_new_(k,j,i), ir_cm);
        }else{
        
          // get monochromatic specific intensity 

          GetCmMCIntensity(ir_cm, tran_coef, ir_cen_, ir_slope_);
          // calculate the shift ratio
          ForwardSplitting(tran_coef, ir_cm, ir_slope_, split_ratio_,
                                         map_bin_start_,map_bin_end_);
          MapIrcmFrequency(ir_cm,ir_shift_);
          
          DetermineShiftRatio(ir_cm,ir_shift_,delta_ratio_);

          // calculate the source term 
          tgas_new_(k,j,i) = MultiGroupAbsScat(wmu_cm,tran_coef, sigma_at, sigma_p, sigma_aer,
                              sigma_s, dt, lorz, rho, tgas_(k,j,i), implicit_coef_,ir_shift_);

          // Add compton scattering 
          // Compton scattering for implicit scheme is added separately
          if(compton_flag_ > 0)
            MultiGroupCompton(wmu_cm,tran_coef,dt,lorz,rho,tgas_new_(k,j,i),ir_shift_);

          // inverseshift
          InverseMapFrequency(ir_shift_,ir_cm);
       
        }
       
         //update specific intensity in the lab frame
         // do not modify ir_ini
         for(int ifr=0; ifr<nfreq; ++ifr){
             lab_ir = &(ir(k,j,i,nang*ifr));
           for(int n=0; n<nang; ++n){
             lab_ir[n] = std::max(ir_cm(n+ifr*nang)/cm_to_lab(n), TINY_NUMBER);
           }
         }

      }// end i
    }// end j
  }// end k

}



void RadIntegrator::AddMultiGroupCompt(MeshBlock *pmb, const Real dt, 
        AthenaArray<Real> &u, AthenaArray<Real> &ir)
{

  // need to transform lab frame ir to co-moving frame
  Radiation *prad=pmb->prad;
  Coordinates *pco = pmb->pcoord;
  
  Real& prat = prad->prat;
  Real invcrat = 1.0/prad->crat;


  Real *lab_ir;
  
  int &nang =prad->nang;
  int &nfreq=prad->nfreq;

  // only apply for multi-grou case
  if((nfreq > 1) && (compton_flag_ > 0)){
  
  
  // Get the temporary array
    AthenaArray<Real> &wmu_cm = wmu_cm_;
    AthenaArray<Real> &tran_coef = tran_coef_;
    AthenaArray<Real> &ir_cm = ir_cm_;
    AthenaArray<Real> &cm_to_lab = cm_to_lab_;


    int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
    int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
 
    for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
        for(int i=is; i<=ie; ++i){

          Real rho = u(IDN,k,j,i);
          Real vx = vel_source_(k,j,i,0);
          Real vy = vel_source_(k,j,i,1);
          Real vz = vel_source_(k,j,i,2);
          Real vel = vx*vx + vy*vy + vz*vz;
        
        
          Real lorzsq = 1.0/(1.0 - vel  * invcrat * invcrat);
          Real lorz = sqrt(lorzsq);
        


         // Prepare the transformation coefficients
          Real numsum = 0.0;

          for(int n=0; n<nang; ++n){
             Real vdotn = vx * prad->mu(0,k,j,i,n) + vy * prad->mu(1,k,j,i,n)
                        + vz * prad->mu(2,k,j,i,n);
             Real vnc = 1.0 - vdotn * invcrat;
             tran_coef(n) = lorz * vnc;
             wmu_cm(n) = prad->wmu(n)/(tran_coef(n) * tran_coef(n));
             numsum += wmu_cm(n);
             cm_to_lab(n) = tran_coef(n)*tran_coef(n)*tran_coef(n)*tran_coef(n);
           
          }
           // Normalize weight in co-moving frame to make sure the sum is one
          numsum = 1.0/numsum;
#pragma omp simd
          for(int n=0; n<nang; ++n){
            wmu_cm(n) *= numsum;
          }
                
          for(int ifr=0; ifr<nfreq; ++ifr){
            lab_ir=&(ir(k,j,i,ifr*nang));
            for(int n=0; n<nang; ++n)
              ir_cm(n+ifr*nang) = std::max(lab_ir[n] * cm_to_lab(n), TINY_NUMBER);
          }// End frequency


          GetCmMCIntensity(ir_cm, tran_coef, ir_cen_, ir_slope_);
          // calculate the shift ratio
          ForwardSplitting(tran_coef, ir_cm, ir_slope_, split_ratio_,
                                         map_bin_start_,map_bin_end_);
          MapIrcmFrequency(ir_cm,ir_shift_);
          
          DetermineShiftRatio(ir_cm,ir_shift_,delta_ratio_);


          // Add compton scattering 
          MultiGroupCompton(wmu_cm,tran_coef,dt,lorz,rho,tgas_new_(k,j,i),ir_shift_);

          // inverseshift
          InverseMapFrequency(ir_shift_,ir_cm);

          for(int ifr=0; ifr<nfreq; ++ifr){
            lab_ir = &(ir(k,j,i,nang*ifr));
            for(int n=0; n<nang; ++n){
              lab_ir[n] = std::max(ir_cm(n+ifr*nang)/cm_to_lab(n), TINY_NUMBER);
            }
          }// end ifr


        }// end i
      }// end j
    }// end k
  }// end nfreq > 1
}

// ir_ini and ir only differ by the source term for explicit scheme
// for implicit scheme, ir also includes the flux divergence term

void RadIntegrator::AddSourceTerms(MeshBlock *pmb, AthenaArray<Real> &u, 
                       AthenaArray<Real> &ir_ini, AthenaArray<Real> &ir)
{

  Radiation *prad=pmb->prad;
  Field *pfield = pmb->pfield;

  Real& prat = prad->prat;
  Real invcrat = 1.0/prad->crat;
  Real invredc = 1.0/prad->reduced_c;
  Real invredfactor = invredc/invcrat;

  Real gm1 = pmb->peos->GetGamma() - 1.0;

  Real rho_floor = pmb->peos->GetDensityFloor();
  
  
  Real *sigma_at, *sigma_aer, *sigma_s, *sigma_p;

  
  int &nang =prad->nang;
  int &nfreq=prad->nfreq;


  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
 
  for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=is; i<=ie; ++i){
      
        // first, calculate Er and Fr in lab frame before the step
        Real er0 = 0.0;
        Real frx0 = 0.0;
        Real fry0 = 0.0;
        Real frz0 = 0.0;
        
        for(int ifr=0; ifr<nfreq; ++ifr){
          Real *p_ir0 =  &(ir_ini(k,j,i,ifr*nang));
          Real er_fr = 0.0;
          Real frx_fr = 0.0;
          Real fry_fr = 0.0;
          Real frz_fr = 0.0;

          if(IM_RADIATION_ENABLED){

            if(pmb->pmy_mesh->pimrad->ite_scheme == 0 || 
               pmb->pmy_mesh->pimrad->ite_scheme == 2){ 
              for(int n=0; n<nang; ++n){
                Real ir_weight = p_ir0[n];
                ir_weight += divflx_(k,j,i,ifr*nang+n);
                if(prad->angle_flag == 1){
                  ir_weight += ang_flx_(k,j,i,ifr*nang+n);
                  ir_weight -= ((imp_ang_coef_(k,j,i,n)) 
                           * ir(k,j,i,ifr*nang+n));
                }
                ir_weight -= (const_coef_(k,j,i,ifr*nang+n) 
                           * ir(k,j,i,ifr*nang+n));
                ir_weight *= prad->wmu(n);
                er_fr  += ir_weight;
                frx_fr += ir_weight * prad->mu(0,k,j,i,n);
                fry_fr += ir_weight * prad->mu(1,k,j,i,n);
                frz_fr += ir_weight * prad->mu(2,k,j,i,n);
              }// end Jacobi iteration
            }else if(pmb->pmy_mesh->pimrad->ite_scheme == 1){
              for(int n=0; n<nang; ++n){
                Real ir_weight = p_ir0[n];
                ir_weight += divflx_(k,j,i,ifr*nang+n);
                if(prad->angle_flag == 1){
                  ir_weight += ang_flx_(k,j,i,ifr*nang+n);
                  ir_weight -= ((imp_ang_coef_(k,j,i,n)) 
                           * ir(k,j,i,ifr*nang+n));                  
                }
                ir_weight -= (const_coef_(k,j,i,ifr*nang+n) 
                           * ir(k,j,i,ifr*nang+n));
                ir_weight -= (left_coef1_(k,j,i,ifr*nang+n)
                           * ir(k,j,i-1,ifr*nang+n));
                if(je > js){
                  ir_weight -= (left_coef2_(k,j,i,ifr*nang+n)
                           * ir(k,j-1,i,ifr*nang+n));
                }
                if(ke > ks){
                  ir_weight -= (left_coef3_(k,j,i,ifr*nang+n)
                           * ir(k-1,j,i,ifr*nang+n));
                }
                ir_weight *= prad->wmu(n);
                er_fr  += ir_weight;
                frx_fr += ir_weight * prad->mu(0,k,j,i,n);
                fry_fr += ir_weight * prad->mu(1,k,j,i,n);
                frz_fr += ir_weight * prad->mu(2,k,j,i,n);
              }
            }// end Gauss-Seidel iteration
          }else{
#pragma omp simd reduction (+:er_fr,frx_fr,fry_fr,frz_fr)
            for(int n=0; n<nang; ++n){
              Real ir_weight = p_ir0[n];
              ir_weight *= prad->wmu(n);
              er_fr  += ir_weight;
              frx_fr += ir_weight * prad->mu(0,k,j,i,n);
              fry_fr += ir_weight * prad->mu(1,k,j,i,n);
              frz_fr += ir_weight * prad->mu(2,k,j,i,n);
            }
          }
          er0 += er_fr;
          frx0 += frx_fr;
          fry0 += fry_fr;
          frz0 += frz_fr;
        }

        Real er = 0.0;
        Real frx = 0.0;
        Real fry = 0.0;
        Real frz = 0.0;
        
        for(int ifr=0; ifr<nfreq; ++ifr){
          Real *p_ir =  &(ir(k,j,i,ifr*nang));
          Real er_fr = 0.0;
          Real frx_fr = 0.0;
          Real fry_fr = 0.0;
          Real frz_fr = 0.0;
#pragma omp simd reduction (+:er_fr,frx_fr,fry_fr,frz_fr)
          for(int n=0; n<nang; ++n){
            Real ir_weight = p_ir[n] * prad->wmu(n);
            er_fr  += ir_weight;
            frx_fr += ir_weight * prad->mu(0,k,j,i,n);
            fry_fr += ir_weight * prad->mu(1,k,j,i,n);
            frz_fr += ir_weight * prad->mu(2,k,j,i,n);
          }
          er += er_fr;
          frx += frx_fr;
          fry += fry_fr;
          frz += frz_fr;
        }// end ifr

        
        // Now apply the radiation source terms to gas with energy and
        // momentum conservation
        u(IM1,k,j,i) += (-prat*(frx- frx0) * invredc);
        u(IM2,k,j,i) += (-prat*(fry- fry0) * invredc);
        u(IM3,k,j,i) += (-prat*(frz- frz0) * invredc);

        //limit the velocity by speed of light
        Real vx = u(IM1,k,j,i)/u(IDN,k,j,i);
        Real vy = u(IM2,k,j,i)/u(IDN,k,j,i);
        Real vz = u(IM3,k,j,i)/u(IDN,k,j,i);
        Real vel = vx*vx+vy*vy+vz*vz;
        vel = sqrt(vel);
        Real ratio = vel * invcrat;
        if(ratio > prad->vmax){
          Real factor = prad->vmax/ratio;
          u(IM1,k,j,i) *= factor;
          u(IM2,k,j,i) *= factor;
          u(IM3,k,j,i) *= factor;
        }

        Real ekin = 0.5 *(SQR(u(IM1,k,j,i))+SQR(u(IM2,k,j,i))
                  +SQR(u(IM3,k,j,i)))/u(IDN,k,j,i);
        Real pb = 0.0;
        if(MAGNETIC_FIELDS_ENABLED){
          pb = 0.5*(SQR(pfield->bcc(IB1,k,j,i))+SQR(pfield->bcc(IB2,k,j,i))
               +SQR(pfield->bcc(IB3,k,j,i)));
        }       

        if(prad->set_source_flag == 2){
          Real eint = tgas_new_(k,j,i) * u(IDN,k,j,i)/gm1;
          u(IEN,k,j,i) = eint + pb + ekin;         
        }else{
          Real e_source = (-prat*(er - er0 ) * invredfactor);  
          // first check that gas internal energy will not become negative
          Real eint = u(IEN,k,j,i) + e_source - ekin - pb;
          Real tgas = eint * gm1/u(IDN,k,j,i);
          if(eint < 0.0){
            eint = tgas_new_(k,j,i) * u(IDN,k,j,i)/gm1;
            u(IEN,k,j,i) = eint + pb + ekin;
          }else if(tgas > prad->t_ceiling_(k,j,i)){
            eint = prad->t_ceiling_(k,j,i) * u(IDN,k,j,i)/gm1;
            u(IEN,k,j,i) = ekin + pb + eint;
          }else{
            u(IEN,k,j,i) += e_source;   
          }       
        }// end source_flag

      }// end i
    }// end j
  }// end k
        
        // check internal energy
//        Real ekin=0.5 * (u(IM1,k,j,i) * u(IM1,k,j,i)
//                       + u(IM2,k,j,i) * u(IM2,k,j,i)
//                       + u(IM3,k,j,i) * u(IM3,k,j,i))/u(IDN,k,j,i);
//        Real pb=0.0;
//        if(MAGNETIC_FIELDS_ENABLED){
//          if(step==1){
//            const Real& bcc1 = pmb->pfield->bcc1(IB1,k,j,i);
//            const Real& bcc2 = pmb->pfield->bcc1(IB2,k,j,i);
//            const Real& bcc3 = pmb->pfield->bcc1(IB3,k,j,i);
//            pb=0.5*(SQR(bcc1) + SQR(bcc2) + SQR(bcc3));
//          }else{
//            const Real& bcc1 = pmb->pfield->bcc(IB1,k,j,i);
//            const Real& bcc2 = pmb->pfield->bcc(IB2,k,j,i);
//            const Real& bcc3 = pmb->pfield->bcc(IB3,k,j,i);
//            pb=0.5*(SQR(bcc1) + SQR(bcc2) + SQR(bcc3));
//          }
//        }
//        Real eint=u(IEN,k,j,i) - pb - ekin;
//        if(eint < 0.0){
//          Real gm1 = pmb->phydro->peos->GetGamma() - 1.0;
//          eint = u(IDN,k,j,i) * tgas /gm1;
//          u(IEN,k,j,i) = eint  + pb + ekin;
//        }
        

}





