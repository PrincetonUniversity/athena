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

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>      // sqrt
#include <fstream>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../parameter_input.hpp"

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real tgas, er1, er2, er3;
  Real sigma1,sigma2,sigma3;
  Real kappa_es, bd_flag;

  er1 = pin->GetOrAddReal("problem","er_1",10.0);
  er2 = pin->GetOrAddReal("problem","er_2",20.0);
  er3 = pin->GetOrAddReal("problem","er_3",30.0);
  tgas = pin->GetOrAddReal("problem","tgas",1.0);
  sigma1 = pin->GetOrAddReal("problem","sigma_1",100.0);
  sigma2 = pin->GetOrAddReal("problem","sigma_2",200.0);
  sigma3 = pin->GetOrAddReal("problem","sigma_3",300.0);
  kappa_es = pin->GetOrAddReal("problem","kappa_es",0.0);
  bd_flag = pin->GetOrAddInteger("problem","black_body",0.0);

  Real gamma = peos->GetGamma();
  Real tr_ini = std::pow(er1,0.25);

  // Initialize hydro variable
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = tgas/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  // Now initialize opacity and specific intensity
  if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
    int nfreq = pnrrad->nfreq;
    int nang = pnrrad->nang;
    AthenaArray<Real> ir_cm;
    ir_cm.NewAthenaArray(pnrrad->n_fre_ang);

    pnrrad->kappa_es=kappa_es;

    Real *ir_lab;

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          for (int ifr=0; ifr<pnrrad->nfreq; ++ifr) {
            ir_lab = &(pnrrad->ir(k,j,i,ifr*nang));
            if (bd_flag == 0) {
              Real er_ini=1.0;
              if (ifr == 0) {
                er_ini = er1;
              } else if (ifr == 1) {
                er_ini = er2;
              } else {
                er_ini = er3;
              }
              for (int n=0; n<pnrrad->nang; n++) {
                ir_lab[n] = er_ini;
              }
            } else {
              Real emission = er1;
              // Initialize with blackbody spectrum
              if (ifr == nfreq-1) {
                emission *= (1.0-pnrrad->FitIntPlanckFunc(pnrrad->nu_grid(ifr)/tr_ini));
              } else {
                emission *= pnrrad->IntPlanckFunc(pnrrad->nu_grid(ifr)/tr_ini,
                                                  pnrrad->nu_grid(ifr+1)/tr_ini);
              }
              for (int n=0; n<pnrrad->nang; n++) {
                ir_lab[n] = emission;
              }
            }
          }
        }
      }
    }

    for (int k=0; k<ncells3; ++k) {
      for (int j=0; j<ncells2; ++j) {
        for (int i=0; i<ncells1; ++i) {
          for (int ifr=0; ifr < nfreq; ++ifr) {
            if (ifr == 0) {
              pnrrad->sigma_s(k,j,i,ifr) = kappa_es;
              pnrrad->sigma_a(k,j,i,ifr) = sigma1;
              pnrrad->sigma_pe(k,j,i,ifr) = sigma1;
              pnrrad->sigma_p(k,j,i,ifr) = sigma1;
            } else if (ifr == 1) {
              pnrrad->sigma_s(k,j,i,ifr) = kappa_es;
              pnrrad->sigma_a(k,j,i,ifr) = sigma2;
              pnrrad->sigma_pe(k,j,i,ifr) = sigma2;
              pnrrad->sigma_p(k,j,i,ifr) = sigma2;
            } else {
              pnrrad->sigma_s(k,j,i,ifr) = kappa_es;
              pnrrad->sigma_a(k,j,i,ifr) = sigma3;
              pnrrad->sigma_pe(k,j,i,ifr) = sigma3;
              pnrrad->sigma_p(k,j,i,ifr) = sigma3;
            }
          }
        }
      }
    }
  }
  return;
}


// calculate the sum of radiation energy density across the mesh for each frequency groups
void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  MeshBlock *pmb = my_blocks(0);
  int  totnum = pmb->pnrrad->nfreq+1;
  AthenaArray<Real> sum_var;
  sum_var.NewAthenaArray(totnum);

  for (int nb=0; nb<nblocal; ++nb) {
    pmb = my_blocks(nb);
    if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED)
      pmb->pnrrad->CalculateMoment(pmb->pnrrad->ir);
    //  Compute the sum
    for (int k=pmb->ks; k<=pmb->ke; k++) {
      for (int j=pmb->js; j<=pmb->je; j++) {
        for (int i=pmb->is; i<=pmb->ie; i++) {
          // get the initial solution
          Real dx = pmb->pcoord->x1f(i+1) - pmb->pcoord->x1f(i);

          sum_var(0) += dx*pmb->phydro->w(IPR,k,j,i)/pmb->phydro->w(IDN,k,j,i);

          for (int ifr=0; ifr<pmb->pnrrad->nfreq; ++ifr) {
            sum_var(ifr+1) += pmb->pnrrad->rad_mom_nu(ifr*13+IER,k,j,i) * dx;
          }
        }
      }
    }
  }

#ifdef MPI_PARALLEL
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE,&sum_var(0),totnum,MPI_ATHENA_REAL,MPI_SUM,0,
               MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&sum_var(0),&sum_var(0),totnum,MPI_ATHENA_REAL,MPI_SUM,0,
               MPI_COMM_WORLD);
  }
#endif

  // only the root process outputs the data
  if (Globals::my_rank == 0) {
    // divide by the mesh size
    for (int i=0; i<totnum; ++i) {
      sum_var(i) /= (mesh_size.x1max - mesh_size.x1min);
    }
    // open output file and write out errors
    std::string fname;
    fname.assign("Averaged_quantity.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if ((pfile = fopen(fname.c_str(),"r")) != NULL) {
      if ((pfile = freopen(fname.c_str(),"a",pfile)) == NULL) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Error output file could not be opened" <<std::endl;
        throw std::runtime_error(msg.str().c_str());
      }

      // The file does not exist -- open the file in write mode and add headers
    } else {
      if ((pfile = fopen(fname.c_str(),"w")) == NULL) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Error output file could not be opened" <<std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      fprintf(pfile,"# T  Er  ");
      fprintf(pfile,"\n");
    }

    // write errors
    for (int n=0; n<totnum; ++n)
      fprintf(pfile,"%e  ",sum_var(n));
    fprintf(pfile,"\n");
    fclose(pfile);
  }
}
