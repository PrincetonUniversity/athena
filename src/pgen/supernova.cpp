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
//! \file supernova.cpp
//  \brief Problem generator for shock tube problems.  
//
// Problem generator for a supernova explosion
//======================================================================================

// C++ headers
#include <fstream>    // ifstream
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>
#include <stdio.h>
#include <math.h>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../field/field.hpp"
#include "../eos/eos.hpp"
#include "../coordinates/coordinates.hpp"
#include "../outputs/outputs.hpp"


Real efin;  //energy of the thermal bomb
Real tbomb;  //time duration of the thermal bomb
Real mdot;   //mass loss in 1e-5 Msol/yr
int nbomb;   //number of grid point where the energy is injected
Real global_M_enc; // set the mass to mass interior to inner boundary
Real last_dt;

//constants
Real m_amu = 1.67377e-24;
Real m_p = 1.6726e-24;
Real kboltz = 1.381e-16;

//composition dependent parameters
Real mu = 0.60*m_amu;
Real muA = 1.24*m_amu;
Real qA = muA/mu - 1.0;



void Outflow_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);




void thermal_bomb(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

    Real volume;
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
            for (int i=pmb->is; i<=pmb->is+nbomb; ++i) {
                if (time <= tbomb)
                {
                volume = pmb->pcoord->GetCellVolume(k, j, i);
                std::cout << cons(IEN,k,j,i) << "," << dt*efin/(volume*tbomb*nbomb) << '\n';
                cons(IEN,k,j,i) += dt*efin/(volume*tbomb*nbomb);
                }
            }
        }
    }
}


void Newtonian( MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim,
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(2);
  return;
}


void Mesh::InitUserMeshData(ParameterInput *pin)
{
    efin = pin->GetReal("problem","efin");
    tbomb = pin->GetReal("problem","tbomb");
    nbomb = pin->GetInteger("problem","nbomb");
    global_M_enc = pin->GetReal("problem","M_enc");


    // AllocateRealUserMeshDataField(1);
    
    // EnrollUserExplicitSourceFunction(thermal_bomb);


   EnrollUserBoundaryFunction(OUTER_X1, Outflow_X2);

    EnrollUserExplicitSourceFunction(Newtonian);
    
}


void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
    
    mdot = pin->GetReal("problem","mdot");
    
    int number_of_lines = 0;
    std::string line;
    std::string fname_p;
    fname_p = pin->GetString("problem", "profile");
    
    std::ifstream file(fname_p.c_str());
    if (not file.is_open())
    {
        std::stringstream msg;
        msg << "### FATAL ERROR in Problem Generator" << std::endl
        << "file " << fname_p << " cannot be opened" << std::endl;
        throw std::runtime_error(msg.str().c_str());
    }
    while (std::getline(file, line))
        ++number_of_lines;
    std::cout << "Number of lines in text file: " << number_of_lines;
    
    int imax = block_size.nx1;

    // reading profile from file
    AthenaArray<Real> rad_from_file;
    rad_from_file.NewAthenaArray(number_of_lines);
    
    AthenaArray<Real> rho_from_file;
    rho_from_file.NewAthenaArray(number_of_lines);
    
    AthenaArray<Real> rho_mapped;
    rho_mapped.NewAthenaArray(imax);
    
    AthenaArray<Real> vel_from_file;
    vel_from_file.NewAthenaArray(number_of_lines);
    
    AthenaArray<Real> vel_mapped;
    vel_mapped.NewAthenaArray(imax);
    
    AthenaArray<Real> press_from_file;
    press_from_file.NewAthenaArray(number_of_lines);
    
    AthenaArray<Real> press_mapped;
    press_mapped.NewAthenaArray(imax);
    
    std::ifstream file1(fname_p.c_str());
    for (int n = 0; n < number_of_lines; ++n)
    {
        Real massx, radx, velx, rhox, tempx, pressx;
        file1 >> massx >> radx >> velx >> rhox >> tempx >> pressx;
        rad_from_file(n) = radx;
        rho_from_file(n) = rhox;
        vel_from_file(n) = velx;
        press_from_file(n) = pressx;
    }
    
    int lower_index;
    int upper_index;
    int middle_index;
    
    Real rho_current;
    Real vel_current;
    Real press_current;
    
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
            if (pcoord->x1v(i) < rad_from_file(number_of_lines-1)) {
                lower_index = 0;
                upper_index = number_of_lines-1;
                while ((upper_index - lower_index)>1){
                    middle_index = (lower_index + upper_index)/2;
                    if (pcoord->x1v(i) < rad_from_file(middle_index)){
                        upper_index = middle_index;
                    }
                    else {
                        lower_index = middle_index;
                    }
                }
                rho_current = (pcoord->x1v(i) - rad_from_file(lower_index))*(rho_from_file(upper_index)-rho_from_file(lower_index))/(rad_from_file(upper_index)-rad_from_file(lower_index)) + rho_from_file(lower_index);
                vel_current = (pcoord->x1v(i) - rad_from_file(lower_index))*(vel_from_file(upper_index)-vel_from_file(lower_index))/(rad_from_file(upper_index)-rad_from_file(lower_index)) + vel_from_file(lower_index);
                press_current = (pcoord->x1v(i) - rad_from_file(lower_index))*(press_from_file(upper_index)-press_from_file(lower_index))/(rad_from_file(upper_index)-rad_from_file(lower_index)) + press_from_file(lower_index);
                
                phydro->u(IDN,k,j,i) = rho_current;
                phydro->u(IM1,k,j,i) = vel_current*rho_current;
                phydro->u(IM2,k,j,i) = 0.0;
                phydro->u(IM3,k,j,i) = 0.0;
                phydro->u(IEN,k,j,i) = press_current/(peos->GetGamma() - 1.0) + 0.5*rho_current*(vel_current*vel_current);
            } else {
                rho_current = mdot*5.016e13/(pcoord->x1v(i)*pcoord->x1v(i));
                vel_current = 1.0e6;
                press_current = rho_current*kboltz*5000.0/mu;
                phydro->u(IDN,k,j,i) = rho_current;
                phydro->u(IM1,k,j,i) = vel_current*rho_current;
                phydro->u(IM2,k,j,i) = 0.0;
                phydro->u(IM3,k,j,i) = 0.0;
                phydro->u(IEN,k,j,i) = press_current/(peos->GetGamma() - 1.0) + 0.5*rho_current*(vel_current*vel_current);
            }
            
            if (i == is) {
                phydro->u(IDN,k,j,i-1) = phydro->u(IDN,k,j,i);
                phydro->u(IM1,k,j,i-1) = phydro->u(IM1,k,j,i);
                phydro->u(IEN,k,j,i-1) = phydro->u(IEN,k,j,i);
                phydro->u(IDN,k,j,i-2) = phydro->u(IDN,k,j,i);
                phydro->u(IM1,k,j,i-2) = phydro->u(IM1,k,j,i);
                phydro->u(IEN,k,j,i-2) = phydro->u(IEN,k,j,i);
            }
            
        }
    }}
    
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{

  Real G = 6.67428e-8;

  Real M_enc = global_M_enc;

  for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      M_enc = global_M_enc; //reset mass every radial coordinate
      for(int i=is; i<=ie; ++i){
        Real rho = phydro->u(IDN,k,j,i);
        Real rr = pcoord->x1f(i+1); //right face radial coordinate
        Real rl = pcoord->x1f(i); //left face radial coordinate
 
        Real DeltaVol = pcoord->GetCellVolume(k, j, i);
        Real DeltaM = rho * DeltaVol;
        Real phil = -G * M_enc/rl;
        Real phir = -G * (M_enc+DeltaM)/rr;
        M_enc += DeltaM;
         
        //                                dPhi / dr
        Real src = - last_dt * rho * G * M_enc / (pcoord->x1v(i)*pcoord->x1v(i));
        phydro->u(IM1,k,j,i) += src;
        phydro->u(IEN,k,j,i) += src * phydro->u(IVX,k,j,i);

        user_out_var(0,k,j,i) = G * M_enc / (pcoord->x1v(i)*pcoord->x1v(i));
        user_out_var(1,k,j,i) = M_enc;
        
      }
    }
  }


}


void Newtonian(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, 
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real G = 6.67428e-8;
  last_dt = dt;

  Real M_enc = global_M_enc;

  for(int k=pmb->ks; k<=pmb->ke; ++k){
    for(int j=pmb->js; j<=pmb->je; ++j){
      M_enc = global_M_enc; //reset mass every radial coordinate
      for(int i=pmb->is; i<=pmb->ie; ++i){
        Real rho = prim(IDN,k,j,i);
        Real rr = pmb->pcoord->x1f(i+1); //right face radial coordinate
        Real rl = pmb->pcoord->x1f(i); //left face radial coordinate
 
        Real DeltaVol = pmb->pcoord->GetCellVolume(k, j, i);
        Real DeltaM = rho * DeltaVol;
        // Real phil = -G * M_enc/rl;
        // Real phir = -G * (M_enc+DeltaM)/rr;
        M_enc += DeltaM;
         
        //                                dPhi / dr
        Real src = - dt * rho * G * M_enc / (pmb->pcoord->x1v(i)*pmb->pcoord->x1v(i));
        cons(IM1,k,j,i) += src;
        cons(IEN,k,j,i) += src * prim(IVX,k,j,i);

        // custom0(i) = src;
        // custom1(i) = M_enc;
        
      }
    }
  }

}




void Outflow_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=NGHOST; ++i) {
          Real &x1g = pco->x1v(ie+i);
          Real &x1 = pco->x1v(ie+i-1);
          if(a(IVX,k,j,ie) < 0.0){
            a(IDN,k,j,ie+i) = a(IDN,k,j,ie);
            a(IVX,k,j,ie+i) = 0.0;
          }else{
            a(IDN,k,j,ie+i) = a(IDN,k,j,ie+i-1) * x1*x1/(x1g*x1g);
            a(IVX,k,j,ie+i) = a(IVX,k,j,ie);
          }
          a(IVY,k,j,ie+i) = a(IVY,k,j,ie);
          a(IVZ,k,j,ie+i) = a(IVZ,k,j,ie);
          a(IEN,k,j,ie+i) = a(IEN,k,j,ie);
        
      }
    }
  }
  return;
}


