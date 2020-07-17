//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file geodesic_grid.cpp
//  \brief implements GeodesicGrid, a geodesic grid on the unit sphere

#include <cassert>
#include <algorithm>  // max
#include <cmath>      // acos, cos, sin, sqrt
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str, string
#include <iostream>

#include "geodesic_grid.hpp"
#include "../athena_arrays.hpp"            // AthenaArray
#include "../athena.hpp"

GeodesicGrid::GeodesicGrid(int nlev) {

  // Some basic vertices of the icosahedron that are needed to build the grid
  Real SinAng = 2.0/std::sqrt(5.0);
  Real CosAng = 1.0/std::sqrt(5.0);

  Real P1[3] = {0.0, 0.0, 1.0};
  Real P2[3] = {SinAng, 0.0, CosAng};
  Real P3[3] = {SinAng*std::cos(0.2*PI), SinAng*std::sin(0.2*PI), -CosAng};
  Real P4[3] = {SinAng*std::cos(-0.4*PI), SinAng*std::sin(-0.4*PI), CosAng};
  Real P5[3] = {SinAng*std::cos(-0.2*PI), SinAng*std::sin(-0.2*PI), -CosAng};
  Real P6[3] = {0.0, 0.0, -1.0};

  geo_grid_xcomp.NewAthenaArray(5,nlev+2*NRGHOST,2*nlev+2*NRGHOST);
  geo_grid_ycomp.NewAthenaArray(5,nlev+2*NRGHOST,2*nlev+2*NRGHOST);
  geo_grid_zcomp.NewAthenaArray(5,nlev+2*NRGHOST,2*nlev+2*NRGHOST);

  geo_grid_pol_xcomp.NewAthenaArray(2);
  geo_grid_pol_ycomp.NewAthenaArray(2);
  geo_grid_pol_zcomp.NewAthenaArray(2);

  /// Here we find the unit normal components at the cell centers for one block (out of five)
  int row_index = NRGHOST;
  for (int l = 0; l < nlev; ++l){
    int col_index = NRGHOST;

    for (int m = l; m < nlev; ++m){
      Real x_point = ((m-l+1)*P2[0] + (nlev-m-1)*P1[0] + l*P4[0])/(Real)(nlev);
      Real y_point = ((m-l+1)*P2[1] + (nlev-m-1)*P1[1] + l*P4[1])/(Real)(nlev);
      Real z_point = ((m-l+1)*P2[2] + (nlev-m-1)*P1[2] + l*P4[2])/(Real)(nlev);

      Real norm = std::sqrt(x_point*x_point+y_point*y_point+z_point*z_point);

      geo_grid_xcomp(0,row_index,col_index) = x_point/norm;
      geo_grid_ycomp(0,row_index,col_index) = y_point/norm;
      geo_grid_zcomp(0,row_index,col_index) = z_point/norm;

      col_index += 1;
    }

    for (int m = nlev-l; m < nlev; ++m){
      Real x_point = ((nlev-l)*P2[0] + (m-nlev+l+1)*P5[0] + (nlev-m-1)*P4[0])/(Real)(nlev);
      Real y_point = ((nlev-l)*P2[1] + (m-nlev+l+1)*P5[1] + (nlev-m-1)*P4[1])/(Real)(nlev);
      Real z_point = ((nlev-l)*P2[2] + (m-nlev+l+1)*P5[2] + (nlev-m-1)*P4[2])/(Real)(nlev);

      Real norm = std::sqrt(x_point*x_point+y_point*y_point+z_point*z_point);

      geo_grid_xcomp(0,row_index,col_index) = x_point/norm;
      geo_grid_ycomp(0,row_index,col_index) = y_point/norm;
      geo_grid_zcomp(0,row_index,col_index) = z_point/norm;

      col_index += 1;
    }

    for (int m = l; m < nlev; ++m){
      Real x_point = ((m-l+1)*P3[0] + (nlev-m-1)*P2[0] + l*P5[0])/(Real)(nlev);
      Real y_point = ((m-l+1)*P3[1] + (nlev-m-1)*P2[1] + l*P5[1])/(Real)(nlev);
      Real z_point = ((m-l+1)*P3[2] + (nlev-m-1)*P2[2] + l*P5[2])/(Real)(nlev);

      Real norm = std::sqrt(x_point*x_point+y_point*y_point+z_point*z_point);

      geo_grid_xcomp(0,row_index,col_index) = x_point/norm;
      geo_grid_ycomp(0,row_index,col_index) = y_point/norm;
      geo_grid_zcomp(0,row_index,col_index) = z_point/norm;

      col_index += 1;
    }

    for (int m = nlev-l; m < nlev; ++m){
      Real x_point = ((nlev-l)*P3[0] + (m-nlev+l+1)*P6[0] + (nlev-m-1)*P5[0])/(Real)(nlev);
      Real y_point = ((nlev-l)*P3[1] + (m-nlev+l+1)*P6[1] + (nlev-m-1)*P5[1])/(Real)(nlev);
      Real z_point = ((nlev-l)*P3[2] + (m-nlev+l+1)*P6[2] + (nlev-m-1)*P5[2])/(Real)(nlev);

      Real norm = std::sqrt(x_point*x_point+y_point*y_point+z_point*z_point);

      geo_grid_xcomp(0,row_index,col_index) = x_point/norm;
      geo_grid_ycomp(0,row_index,col_index) = y_point/norm;
      geo_grid_zcomp(0,row_index,col_index) = z_point/norm;

      col_index += 1;
    }

    row_index += 1;

  }

  // Rotating the block 4 times by 2*PI/5, we cover the entire grid (apart from poles!)
  for (int bl = 1; bl < 5; ++bl){
    for (int l = NRGHOST; l < nlev+NRGHOST; ++l){
      for (int m = NRGHOST; m < 2*nlev+NRGHOST; ++m){
        Real x_point = geo_grid_xcomp(0,l,m)*std::cos(bl*0.4*PI)+geo_grid_ycomp(0,l,m)*std::sin(bl*0.4*PI);
        Real y_point = geo_grid_ycomp(0,l,m)*std::cos(bl*0.4*PI)-geo_grid_xcomp(0,l,m)*std::sin(bl*0.4*PI);
        Real z_point = geo_grid_zcomp(0,l,m);

        geo_grid_xcomp(bl,l,m) = x_point;
        geo_grid_ycomp(bl,l,m) = y_point;
        geo_grid_zcomp(bl,l,m) = z_point;
      }
    }
  }

  // Coordinates of the poles
  geo_grid_pol_xcomp(0) = 0.0; //north pole
  geo_grid_pol_ycomp(0) = 0.0;
  geo_grid_pol_zcomp(0) = 1.0;

  geo_grid_pol_xcomp(1) = 0.0; //south pole
  geo_grid_pol_ycomp(1) = 0.0;
  geo_grid_pol_zcomp(1) = -1.0;

  // Fill in the ghost regions of the blocks
  FillInGhostRegions(nlev, geo_grid_xcomp, geo_grid_pol_xcomp);
  FillInGhostRegions(nlev, geo_grid_ycomp, geo_grid_pol_ycomp);
  FillInGhostRegions(nlev, geo_grid_zcomp, geo_grid_pol_zcomp);


  // Compute total number of points in the grid
  this->nlev_ = nlev;
  this->numpoints_ = 10*nlev_*nlev_ + 2;


  // Arrange the indices of all points in a 5-block structure (plus 2 poles)
  // This will be used in order to find the indices of the neighbors
  geo_grid_index.NewAthenaArray(5,nlev+2*NRGHOST,2*nlev+2*NRGHOST);
  geo_grid_pol_index.NewAthenaArray(2);

  for (int bl = 0; bl < 5; ++bl){
    for (int l = 0; l < nlev; ++l){
      for (int m = 0; m < 2*nlev; ++m){
        int lm = bl*2*nlev*nlev + l*2*nlev + m;
        geo_grid_index(bl,l+NRGHOST,m+NRGHOST) = lm;
      }
    }
  }
  for (int pl = 0; pl < 2; ++pl){
    int lm = 5*2*nlev*nlev + pl;
    geo_grid_pol_index(pl) = lm;
  }
  FillInGhostRegions(nlev, geo_grid_index, geo_grid_pol_index);

}

//_______________________________________________________________________
///////////////////////// Number of vertices ////////////////////////////

int GeodesicGrid::NumVertices() const {
  return numpoints_;
};

//_______________________________________________________________________
//////////////////////// Indices of neighbors ///////////////////////////

int GeodesicGrid::NumNeighbors(int ic, int neighbors[6]) const {
  assert(ic >= 0 && ic < NumVertices());

  int num_neigh;

  if(ic == 10*nlev_*nlev_){ //north pole

    for(int bl = 0; bl < 5; ++bl){
      neighbors[bl] = geo_grid_index(bl,1,1);
    }
    neighbors[5] = EMPTY;
    num_neigh = 5;

  } else if (ic == 10*nlev_*nlev_+1){ //south pole

    for(int bl = 0; bl < 5; ++bl){
      neighbors[bl] = geo_grid_index(bl,nlev_,2*nlev_);
    }
    neighbors[5] = EMPTY;
    num_neigh = 5;

  } else {

    int ibl0 =  ic / (2*nlev_*nlev_);
    int ibl1 = (ic % (2*nlev_*nlev_)) / (2*nlev_);
    int ibl2 = (ic % (2*nlev_*nlev_)) % (2*nlev_);

    neighbors[0] = geo_grid_index(ibl0,ibl1+NRGHOST,  ibl2+NRGHOST+1);
    neighbors[1] = geo_grid_index(ibl0,ibl1+NRGHOST+1,ibl2+NRGHOST);
    neighbors[2] = geo_grid_index(ibl0,ibl1+NRGHOST+1,ibl2+NRGHOST-1);
    neighbors[3] = geo_grid_index(ibl0,ibl1+NRGHOST,  ibl2+NRGHOST-1);
    neighbors[4] = geo_grid_index(ibl0,ibl1+NRGHOST-1,ibl2+NRGHOST);

    if(ic % (2*nlev_*nlev_) == nlev_-1 || ic % (2*nlev_*nlev_) == 2*nlev_-1){
      neighbors[5] = EMPTY; //vertices of the orig. icosahedron apart from poles
      num_neigh = 5;
    } else {                //all points that have 6 neighbors
      neighbors[5] = geo_grid_index(ibl0,ibl1+NRGHOST-1,ibl2+NRGHOST+1);
      num_neigh = 6;
    }

  }

  return num_neigh;
};

//___________________________________________________________________________
/////////////////// Polar coordinates of the vertices ///////////////////////

void GeodesicGrid::PositionPolar(int ic, Real * theta, Real * phi) const {
  assert(ic >= 0 && ic < NumVertices());

  Real x_, y_, z_;
  Position(ic,&x_,&y_,&z_);

  *theta = std::acos(z_);
  *phi   = std::atan2(y_,x_);

};

//___________________________________________________________________________
/////////////////// Cartesian coordinates of the vertices ///////////////////

void GeodesicGrid::Position(int ic, Real * x, Real * y, Real * z) const {
  assert(ic >= 0 && ic < NumVertices());

  int ibl0 =  ic / (2*nlev_*nlev_);
  int ibl1 = (ic % (2*nlev_*nlev_)) / (2*nlev_);
  int ibl2 = (ic % (2*nlev_*nlev_)) % (2*nlev_);

  if(ibl0 == 5){ //north or south poles
    *x=geo_grid_pol_xcomp(ibl2);
    *y=geo_grid_pol_ycomp(ibl2);
    *z=geo_grid_pol_zcomp(ibl2);
  } else {       //all other points apart from the poles
    *x=geo_grid_xcomp(ibl0,ibl1+NRGHOST,ibl2+NRGHOST);
    *y=geo_grid_ycomp(ibl0,ibl1+NRGHOST,ibl2+NRGHOST);
    *z=geo_grid_zcomp(ibl0,ibl1+NRGHOST,ibl2+NRGHOST);
  }
};

//___________________________________________________________________________
////////////////// Cartesian coordinates of the midpoints ///////////////////

void GeodesicGrid::PositionMid(int ic1, int ic2, Real * x, Real * y, Real * z) const {
  assert(ic1 >= 0 && ic1 < NumVertices());
  assert(ic2 >= 0 && ic2 < NumVertices());
  assert(AreNeighbors(ic1, ic2));

  Real x1, y1, z1;
  Real x2, y2, z2;

  Position(ic1,&x1,&y1,&z1);
  Position(ic2,&x2,&y2,&z2);

  Real xm = 0.5*(x1+x2);
  Real ym = 0.5*(y1+y2);
  Real zm = 0.5*(z1+z2);

  Real norm = std::sqrt(xm*xm+ym*ym+zm*zm);

  *x = xm/norm;
  *y = ym/norm;
  *z = zm/norm;
};

//___________________________________________________________________________
//////////////////////////// Unit flux directions ///////////////////////////

void GeodesicGrid::UnitFluxDir(int ic1, int ic2, Real * dtheta, Real * dphi) const {
  assert(ic1 >= 0 && ic1 < NumVertices());
  assert(ic2 >= 0 && ic2 < NumVertices());
  assert(AreNeighbors(ic1, ic2));

  int ibl0 =  ic1 / (2*nlev_*nlev_);

  Real zeta1, psi1;

  PositionPolar(ic1,&zeta1,&psi1);

  Real xm, ym, zm;
  PositionMid(ic1,ic2,&xm,&ym,&zm);

  Real zetam = std::acos(zm);
  Real psim  = std::atan2(ym,xm);

  Real a_par, p_par;

  if (std::abs(psim-psi1) < 1.0e-10 || ibl0 == 5){
    *dtheta = std::copysign(1.0,zetam-zeta1);
    *dphi = 0.0;
  }
  else{
    GreatCircleParam(zeta1,zetam,psi1,psim,&a_par,&p_par);
    Real zeta_deriv = a_par*std::sin(psim-p_par)/(1.0+a_par*a_par*std::cos(psim-p_par)*std::cos(psim-p_par));
    Real denom = 1.0/std::sqrt(zeta_deriv*zeta_deriv+std::sin(zetam)*std::sin(zetam));
    Real signfactor = std::copysign(1.0,psim-psi1)*std::copysign(1.0,0.5*M_PI-std::abs(psim-psi1));
    *dtheta = signfactor * zeta_deriv * denom;
    *dphi   = signfactor * denom;
  }

};

//___________________________________________________________________________
//////////////// Compute weight and lengths of dual edges ///////////////////

Real GeodesicGrid::ComputeWeightAndDualEdges(int ic, Real length[6]) const {
  assert(ic >= 0 && ic < NumVertices());

  int nvec[6];
  int nnum = NumNeighbors(ic, nvec);

  Real x0, y0, z0;
  Position(ic,&x0,&y0,&z0);

  Real weight_current = 0.0;

  for (int nb = 0; nb < nnum; ++nb){

    Real xn1, yn1, zn1;
    Real xn2, yn2, zn2;
    Real xn3, yn3, zn3;

    Position(nvec[(nb + nnum - 1)%nnum],&xn1,&yn1,&zn1);
    Position(nvec[nb],                  &xn2,&yn2,&zn2);
    Position(nvec[(nb + 1)%nnum],       &xn3,&yn3,&zn3);

    Real xc1, yc1, zc1;
    Real xc2, yc2, zc2;

    CircumcenterNormalized(x0,xn1,xn2,y0,yn1,yn2,z0,zn1,zn2,&xc1,&yc1,&zc1);
    CircumcenterNormalized(x0,xn2,xn3,y0,yn2,yn3,z0,zn2,zn3,&xc2,&yc2,&zc2);

    Real scalprod_c1 = x0*xc1 + y0*yc1 + z0*zc1;
    Real scalprod_c2 = x0*xc2 + y0*yc2 + z0*zc2;
    Real scalprod_12 = xc1*xc2 + yc1*yc2 + zc1*zc2;

    Real numerator = std::abs(x0*(yc1*zc2-yc2*zc1)+y0*(xc2*zc1-xc1*zc2)+z0*(xc1*yc2-yc1*xc2));
    Real denominator = 1.0+scalprod_c1+scalprod_c2+scalprod_12;

    weight_current += 2.0*std::atan(numerator/denominator);

    length[nb] = std::acos(scalprod_12);

  }

  if (nnum == 5){
    length[5] = std::numeric_limits<Real>::quiet_NaN();
  }

  return weight_current;
};
//___________________________________________________________________________
/////////////////////////// Compute weight only /////////////////////////////

Real GeodesicGrid::ComputeWeight(int ic) const {
  assert(ic >= 0 && ic < NumVertices());

  Real dum_array[6];

  Real weight = ComputeWeightAndDualEdges(ic, dum_array);

  return weight;
};

//___________________________________________________________________________
/////////////////// Check if the points are neighbors ///////////////////////

bool GeodesicGrid::AreNeighbors(int ic1, int ic2) const {
  assert(ic1 >= 0 && ic1 <= NumVertices());
  assert(ic2 >= 0 && ic2 <= NumVertices());

  int nvec[6];
  int nnum = NumNeighbors(ic1, nvec);

  bool out = false;
  for (int i = 0; i < nnum; ++i) {
    out |= ic2 == nvec[i];
  }
  return out;
}

//_______________________________________________________________________
////////////// Compute xi and eta coordinates for WENO scheme ///////////

void GeodesicGrid::ComputeXiEta(int ic, Real xi[6], Real eta[6]) const {
  assert(ic >= 0 && ic < NumVertices());

  Real x0, y0, z0;
  Position(ic,&x0,&y0,&z0);

  int nvec[6];
  int nnum = NumNeighbors(ic, nvec);

  Real A_angle = 0;

  for (int nb = 0; nb < nnum; ++nb){

    Real xn1, yn1, zn1;
    Real xn2, yn2, zn2;

    Position(nvec[nb],           &xn1,&yn1,&zn1);
    Position(nvec[(nb + 1)%nnum],&xn2,&yn2,&zn2);

    Real n1_x = y0*zn1 - yn1*z0;
    Real n1_y = z0*xn1 - zn1*x0;
    Real n1_z = x0*yn1 - xn1*y0;

    Real n2_x = y0*zn2 - yn2*z0;
    Real n2_y = z0*xn2 - zn2*x0;
    Real n2_z = x0*yn2 - xn2*y0;

    Real norm1 = std::sqrt(n1_x*n1_x+n1_y*n1_y+n1_z*n1_z);
    Real norm2 = std::sqrt(n2_x*n2_x+n2_y*n2_y+n2_z*n2_z);

    Real cosA = std::min((n1_x*n2_x+n1_y*n2_y+n1_z*n2_z)/(norm1*norm2),1.0);

    Real scalprod_c1 = x0*xn1 + y0*yn1 + z0*zn1;

    Real c_len = std::acos(scalprod_c1);

    xi[nb] = c_len*std::cos(A_angle);
    eta[nb] = c_len*std::sin(A_angle);

    A_angle += std::acos(cosA);

  }

  if (nnum == 5){
    xi[5] = std::numeric_limits<Real>::quiet_NaN();
    eta[5] = std::numeric_limits<Real>::quiet_NaN();
  }


};

//_______________________________________________________________________
/////////////// Compute arc length between two grid points //////////////

Real GeodesicGrid::ArcLength(int ic1, int ic2) const {
  assert(ic1 >= 0 && ic1 <= NumVertices());
  assert(ic2 >= 0 && ic2 <= NumVertices());
  assert(AreNeighbors(ic1, ic2));

  Real x1, y1, z1;
  Real x2, y2, z2;

  Position(ic1,&x1,&y1,&z1);
  Position(ic2,&x2,&y2,&z2);

  Real arclen = std::acos(x1*x2+y1*y2+z1*z2);

  return arclen;
}



////////////////////////// additional functions ///////////////////////////////

void GeodesicGrid::CircumcenterNormalized(Real x1, Real x2, Real x3, Real y1, Real y2, Real y3, Real z1, Real z2, Real z3, Real * x_cc, Real * y_cc, Real * z_cc) const {
  Real a = std::sqrt((x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2));
  Real b = std::sqrt((x1-x3)*(x1-x3)+(y1-y3)*(y1-y3)+(z1-z3)*(z1-z3));
  Real c = std::sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
  Real denom = 1.0/((a+c+b)*(a+c-b)*(a+b-c)*(b+c-a));
  Real x_c = (x1*a*a*(b*b+c*c-a*a)+x2*b*b*(c*c+a*a-b*b)+x3*c*c*(a*a+b*b-c*c))*denom;
  Real y_c = (y1*a*a*(b*b+c*c-a*a)+y2*b*b*(c*c+a*a-b*b)+y3*c*c*(a*a+b*b-c*c))*denom;
  Real z_c = (z1*a*a*(b*b+c*c-a*a)+z2*b*b*(c*c+a*a-b*b)+z3*c*c*(a*a+b*b-c*c))*denom;
  Real norm_c = std::sqrt(x_c*x_c+y_c*y_c+z_c*z_c);
  *x_cc = x_c/norm_c;
  *y_cc = y_c/norm_c;
  *z_cc = z_c/norm_c;
}

void GeodesicGrid::GreatCircleParam(Real zeta1, Real zeta2, Real psi1, Real psi2, Real * apar, Real * psi0)const {
  Real atilde = (std::sin(psi2)/std::tan(zeta1)-std::sin(psi1)/std::tan(zeta2))/std::sin(psi2-psi1);
  Real btilde = (std::cos(psi2)/std::tan(zeta1)-std::cos(psi1)/std::tan(zeta2))/std::sin(psi1-psi2);
  *psi0 = std::atan2(btilde, atilde);
  *apar = std::sqrt(atilde*atilde+btilde*btilde);
}




