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
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../nr_radiation/integrators/rad_integrators.hpp"
#include "../nr_radiation/radiation.hpp"
#include "../parameter_input.hpp"

static Real amp = 1.e-6;


// The real and imaginary parts of the linear waves
// used as initial condition as well as test solution to compare with
static Real v_r, v_i, p_r, p_i, er_r, er_i, fr_r, fr_i;
static Real omegareal, omegaimg;

// AMR refinement condition
int RefinementCondition(MeshBlock *pmb);

//======================================================================================
/*! \file rad_linearwave.cpp
 *  \brief Linear wave test for the radiative transfer module
 *
 *====================================================================================*/


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  if (!pin->GetOrAddBoolean("problem","compute_error",false)) return;

  int totnum=5;
  // Initialize errors to zero
  Real * l1_err = new Real[totnum];
  Real * max_err = new Real[totnum];
  for (int i=0; i<totnum; ++i) {
    l1_err[i]=0.0;
    max_err[i]=0.0;
  }
  Real knum = 2.0 * PI;

  for (int nb=0; nb<nblocal; ++nb) {
    MeshBlock *pmb = my_blocks(nb);
    if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
      pmb->pnrrad->CalculateMoment(pmb->pnrrad->ir);
    }
    //  Compute errors
    for (int k=pmb->ks; k<=pmb->ke; k++) {
      for (int j=pmb->js; j<=pmb->je; j++) {
        for (int i=pmb->is; i<=pmb->ie; i++) {
        // get the initial solution
        Real const &x1 = pmb->pcoord->x1v(i);
        Real theta = knum * x1;
        Real delrho = amp * cos(theta);
        Real delv = amp* (v_r * cos(theta) + v_i * sin(theta));
        Real delp = amp * (p_r * cos(theta) + p_i * sin(theta));
        Real der = amp * (er_r * cos(theta) + er_i * sin(theta));
        Real dfr = amp * (fr_r * cos(theta) + fr_i * sin(theta));

        Real diffrho = (pmb->phydro->u(IDN,k,j,i) - 1.0)
                       - delrho * std::exp(-omegaimg*time);
        Real diffvel= pmb->phydro->w(IVX,k,j,i)
                      - delv*std::exp(-omegaimg*time);
        Real diffpre = (pmb->phydro->w(IEN,k,j,i) - 1.0)
                       - delp *std::exp(-omegaimg*time);
        Real differ= (pmb->pnrrad->rad_mom(IER,k,j,i) - 1.0)
                     - der*std::exp(-omegaimg*time);
        Real difffr= pmb->pnrrad->rad_mom(IFR1,k,j,i)
                     -  dfr *std::exp(-omegaimg*time);

        l1_err[0] += std::abs(diffrho);
        max_err[0] = std::max(std::abs(diffrho), max_err[0]);

        l1_err[1] += std::abs(diffvel);
        max_err[1] = std::max(std::abs(diffvel), max_err[1]);

        l1_err[2] += std::abs(diffpre);
        max_err[2] = std::max(std::abs(diffpre), max_err[2]);

        l1_err[3] += std::abs(differ);
        max_err[3] = std::max(std::abs(differ), max_err[3]);

        l1_err[4] += std::abs(difffr);
        max_err[4] = std::max(std::abs(difffr), max_err[4]);
        }
      }
    }
  }
  // normalize errors by number of cells
  for (int i=0; i<totnum; ++i) {
    l1_err[i] = l1_err[i]/ static_cast<Real>(GetTotalCells());
  }
  Real rms_err = 0.0, max_max_over_l1=0.0;

#ifdef MPI_PARALLEL
  if (Globals::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE,&l1_err,totnum,MPI_ATHENA_REAL,MPI_SUM,0,
               MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE,&max_err,totnum,MPI_ATHENA_REAL,MPI_MAX,0,
               MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&l1_err,&l1_err,totnum,MPI_ATHENA_REAL,MPI_SUM,0,
               MPI_COMM_WORLD);
    MPI_Reduce(&max_err,&max_err,totnum,MPI_ATHENA_REAL,MPI_MAX,0,
               MPI_COMM_WORLD);
  }
#endif

  // only the root process outputs the data
  if (Globals::my_rank == 0) {
    // compute rms error
    for (int i=0; i<totnum; ++i) {
       rms_err += SQR(l1_err[i]);
       max_max_over_l1 = std::max(max_max_over_l1, (max_err[i]/l1_err[i]));
    }
    rms_err = std::sqrt(rms_err);

    // open output file and write out errors
    std::string fname;
    fname.assign("linearwave-errors.dat");
    std::stringstream msg;
    FILE *pfile;

    // The file exists -- reopen the file in append mode
    if((pfile = fopen(fname.c_str(),"r")) != NULL) {
      if((pfile = freopen(fname.c_str(),"a",pfile)) == NULL) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Error output file could not be opened" <<std::endl;
        throw std::runtime_error(msg.str().c_str());
      }

    // The file does not exist -- open the file in write mode and add headers
    } else {
      if((pfile = fopen(fname.c_str(),"w")) == NULL) {
        msg << "### FATAL ERROR in function [Mesh::UserWorkAfterLoop]"
            << std::endl << "Error output file could not be opened" <<std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      fprintf(pfile,"# Nx1  Nx2  Nx3  Ncycle  ");
      fprintf(pfile,"RMS-L1-Error  d_L1  M1_L1  P_L1  Er_L1  Fr_L1");
      fprintf(pfile,"  Largest-Max/L1  d_max  M1_max  P_max  Er_max  Fr_max");
      fprintf(pfile,"\n");
    }

    // write errors
    fprintf(pfile,"%d  %d",mesh_size.nx1,mesh_size.nx2);
    fprintf(pfile,"  %d  %d",mesh_size.nx3,ncycle);
    fprintf(pfile,"  %e  %e  %e  %e",rms_err,l1_err[0],l1_err[1],l1_err[2]);
    fprintf(pfile,"  %e  %e",l1_err[3],l1_err[4]);
    fprintf(pfile,"  %e  %e  ",max_max_over_l1,max_err[0]);
    fprintf(pfile,"%e  %e  %e  %e",max_err[1],max_err[2],max_err[3],max_err[4]);
    fprintf(pfile,"\n");
    fclose(pfile);
  }
  delete[] l1_err;
  delete[] max_err;
}



void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  int regime = pin->GetOrAddInteger("problem","regime",8);
  Real sigma0;
  switch(regime) {
    case 1:
      v_r= 1.2909945401100994;
      v_i= 0.00008594287531051318;
      p_r= 1.666666822510812;
      p_i= 0.00021916408244003488;
      er_r= 0.000023938771551742892;
      er_i= 0.002595343282507656;
      fr_r= 5.803552452185792e-6;
      fr_i= 0.004579152314659846;
      omegareal= 8.111557926068842;
      omegaimg=  0.0005399950114077836;
      sigma0 = 0.01;
      break;
    case 2:
      v_r= 1.2910035706152658;
      v_i= 0.0008586435251690642;
      p_r= 1.666682219042423;
      p_i= 0.002189719996892813;
      er_r= 0.0023915640535869353;
      er_i= 0.025910075592072698;
      fr_r= 0.0005794967774284472;
      fr_i= 0.045749450689271914;
      omegareal= 8.111614666406222;
      omegaimg=  0.00539501638144715;
      sigma0 = 0.1;
      break;
    case 3:
      v_r= 1.2917695143401924;
      v_i= 0.00787279209054369;
      p_r= 1.6679409990335867;
      p_i= 0.0201456313053114;
      er_r= 0.2181601920018504;
      er_i= 0.22077999320680775;
      fr_r= 0.050320573067502435;
      fr_i= 0.4191946132602506;
      omegareal= 8.116427232764806;
      omegaimg=  0.04946621158978377;
      sigma0 = 1.0;
      break;
    case 4:
      v_r= 1.2920340482874326;
      v_i= 0.009218020632790785;
      p_r= 1.6616893859774136;
      p_i= 0.02432338104787091;
      er_r= 2.35100088818539;
      er_i= 0.03457134839094067;
      fr_r= 0.20390000836732208;
      fr_i= 0.47734523049957645;
      omegareal= 8.11808934857535;
      omegaimg=  0.057918531801229335;
      sigma0 = 10.0;
      break;
    case 5:
      v_r= 1.2909668031826855;
      v_i= 0.001045499271992685;
      p_r= 1.6580540395720993;
      p_i= 0.002767087292105325;
      er_r= 2.6288579088325177;
      er_i= 0.0005865742021802234;
      fr_r= 0.1725541379069982;
      fr_i= 0.05379877178913972;
      omegareal= 8.11138364981405;
      omegaimg=  0.006569065664451391;
      sigma0 = 100.0;
      break;
    case 6:
      v_r= 1.290810912416655;
      v_i= 0.008591414402392937;
      p_r= 1.6661135402797116;
      p_i= 0.0219056348896828;
      er_r= -0.000051245227843874086;
      er_i= 0.0025938927144686375;
      fr_r= -0.0001441706051503446;
      fr_i= 0.004575470538400113;
      omegareal= 8.110404159243402;
      omegaimg=  0.05398144874100638;
      sigma0 = 0.01;
      break;
    case 7:
      v_r= 1.2728792327581575;
      v_i= 0.08294778800855061;
      p_r= 1.612846028223176;
      p_i= 0.2082461771006987;
      er_r= -0.0047877737151852395;
      er_i= 0.024520349202385106;
      fr_r= -0.013679926554110659;
      fr_i= 0.04217262048775258;
      omegareal= 7.997736093080079;
      omegaimg=  0.5211763228783722;
      sigma0 = 0.1;
      break;
    case 8:
      v_r= 1.0183496112257058;
      v_i= 0.1121215711780068;
      p_r= 1.0220380692314723;
      p_i= 0.18993018794365163;
      er_r= -0.026018127896336885;
      er_i= 0.12095401964915764;
      fr_r= -0.10566859341556321;
      fr_i= 0.030196412832965945;
      omegareal= 6.398479314825398;
      omegaimg=  0.7044806086435424;
      sigma0 = 1.0;
      break;
    case 9:
      v_r= 1.2092908080104792;
      v_i= 0.25516371383473957;
      p_r= 1.1686533807826482;
      p_i= 0.28223448062039036;
      er_r= 0.6544229491473885;
      er_i= 1.0653905194884314;
      fr_r= -0.04918496753828805;
      fr_i= 0.17766949176109614;
      omegareal= 7.598198236998773;
      omegaimg=  1.6032408976918122;
      sigma0 = 10.0;
      break;
    case 10:
      v_r= 1.3993585014150265;
      v_i= 0.03450304866280985;
      p_r= 1.4211576732442506;
      p_i= 0.04441546638644912;
      er_r= 1.6841247487314883;
      er_i= 0.17440949185811563;
      fr_r= 0.18330452762195915;
      fr_i= 0.038269237092134034;
      omegareal= 8.792428775567739;
      omegaimg=  0.21678904841106914;
      sigma0 = 100.0;
      break;
    case 11:
      v_r= 1.271935073054169;
      v_i= 0.08266820679991217;
      p_r= 1.6110935636801302;
      p_i= 0.20757539251245596;
      er_r= -0.0006798060846963105;
      er_i= 0.0023794592885675314;
      fr_r= -0.0014238164476813107;
      fr_i= 0.004188456336444625;
      omegareal= 7.991803762700347;
      omegaimg=  0.5194196623360917;
      sigma0 = 0.01;
      break;
    case 12:
      v_r= 1.014892795882639;
      v_i= 0.08976614261743401;
      p_r= 1.0234931616949725;
      p_i= 0.15900654619222715;
      er_r= -0.003705180537531039;
      er_i= 0.007369404109744249;
      fr_r= -0.010447571816544906;
      fr_i= 0.0022692464424416106;
      omegareal= 6.376759503452208;
      omegaimg=  0.5640173083760487;
      sigma0 = 0.1;
      break;
    case 13:
      v_r= 0.9977393625226653;
      v_i= 0.1324385410980396;
      p_r= 0.9963589910920785;
      p_i= 0.03350714894270625;
      er_r= -0.009348064816826119;
      er_i= 0.07057816396845434;
      fr_r= -0.011965869992714585;
      fr_i=  0.006087909287739512;
      omegareal= 6.268981302997137;
      omegaimg= 0.8321358955315021;
      sigma0 = 1.0;
      break;
    case 14:
      v_r= 2.1806327929670344;
      v_i= 0.6503652227748439;
      p_r= 1.2398074098916425;
      p_i= 0.24885283255343593;
      er_r= 0.9567316753829218;
      er_i= 0.9851131720168348;
      fr_r= 0.12816968396550954;
      fr_i= 0.2810151431629262;
      omegareal= 13.701319925124455;
      omegaimg=  4.086365212039477;
      sigma0 = 10.0;
      break;
    case 15:
      v_r= 2.2904573528686267;
      v_i= 0.0604384248330597;
      p_r= 1.344412710622572;
      p_i= 0.025960902814908177;
      er_r= 1.3776131567434966;
      er_i= 0.10314648118800326;
      fr_r= 0.3038178476543796;
      fr_i= 0.032551125836846495;
      omegareal= 14.391367986265605;
      omegaimg=  0.37974582290015857;
      sigma0 = 100.0;
      break;
    default:
    std::stringstream msg;
    msg << "### FATAL ERROR in ProblemGenerator" << std::endl
        << "Invalid Linear wave mode "<< regime << " is specified" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  Real gamma = peos->GetGamma();
  Real knum = 2.0 * PI;
  Real rho0 = 1.0, p0 = 1.0;
  //Real e0 = p0/(gamma-1.0);

  // Initialize hydro variable
  for (int k=0; k<ncells3; ++k) {
    for (int j=0; j<ncells2; ++j) {
      for (int i=0; i<ncells1; ++i) {
        Real const &x1 = pcoord->x1v(i);
        // Real const &x2 = pcoord->x2v(j);
        // Real const &x3 = pcoord->x3v(k);
        Real theta = knum * x1;
        Real delv = amp* (v_r * cos(theta) + v_i * sin(theta));
        Real delp = amp * (p_r * cos(theta) + p_i * sin(theta));

        phydro->u(IDN,k,j,i) = rho0 + amp * cos(theta);
        phydro->u(IM1,k,j,i) = delv;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = (p0 + delp)/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
        if (NR_RADIATION_ENABLED || IM_RADIATION_ENABLED) {
          Real der = amp * (er_r * cos(theta) + er_i * sin(theta));
          Real dfr = amp * (fr_r * cos(theta) + fr_i * sin(theta));

          Real jr = (1.0+der);
          Real hr = dfr;
          for (int ifr=0; ifr<pnrrad->nfreq; ++ifr) {
            for (int n=0; n<pnrrad->nang; ++n) {
               Real const &weight = pnrrad->wmu(n);
               Real const &miux = pnrrad->mu(0,k,j,i,n);
               pnrrad->ir(k,j,i,ifr*pnrrad->nang+n) = (jr/(4.0*weight)
                                                       + hr/(4.0*weight*miux));
            }
            pnrrad->sigma_s(k,j,i,ifr) = 0.0;
            pnrrad->sigma_a(k,j,i,ifr) = sigma0;
            pnrrad->sigma_pe(k,j,i,ifr) = sigma0;
            pnrrad->sigma_p(k,j,i,ifr) = sigma0;
          }
        }
      }
    }
  }
  return;
}

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
  }
}


// refinement condition: density curvature
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> const &w = pmb->phydro->w;
  Real dmax = 0.0, dmin = 2.0;  // max and min densities
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
        if (w(IDN,k,j,i) > dmax) dmax = w(IDN,k,j,i);
        if (w(IDN,k,j,i) < dmin) dmin = w(IDN,k,j,i);
      }
    }
  }
  // refine : delta rho > 0.9*amp
  if (dmax-1.0 > 0.9*amp) return 1;
  return -1;
}
