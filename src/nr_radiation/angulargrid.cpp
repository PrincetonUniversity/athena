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
//! \file angulargrid.cpp
//  \brief implementation of angular grid in class Radiation
//======================================================================================

// C headers

// C++ headers
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "./radiation.hpp"

//--------------------------------------------------------------------------------------
// \!fn void Radiation::AngularGrid(int angle_flag, int nmu)
// \brief function to create the angular grid

void NRRadiation::AngularGrid(int angle_flag, int nmu) {
  std::stringstream msg;
  MeshBlock *pmb=pmy_block;
  // allocate some temporaray arrays
  AthenaArray<Real> mu2tmp, mutmp, wtmp2, wtmp;
  AthenaArray<Real> pmat, pinv, wpf;
  AthenaArray<int> plab, pl;

  int n_ang = nang/noct;

  mu2tmp.NewAthenaArray(nmu);
  mutmp.NewAthenaArray(n_ang,3);
  wtmp2.NewAthenaArray(nmu-1);
  wtmp.NewAthenaArray(nmu);

  pmat.NewAthenaArray(nmu,nmu);
  pinv.NewAthenaArray(nmu-1,nmu-1);
  plab.NewAthenaArray(n_ang);
  pl.NewAthenaArray(nmu,3);
  wpf.NewAthenaArray(nmu-1);


  // initialize coordinate direction
  int axisx=0, axisy=1, axisz=2;
  pmy_block->pcoord->AxisDirection(&axisx, &axisy, &axisz);

  // check the dimension of the problem

  int n1z = pmb->ncells1, n2z = pmb->ncells2, n3z = pmb->ncells3;

  int ndim = 1;
  if (n2z > 1) ndim = 2;
  if (n3z > 1) ndim = 3;

  if (angle_flag == 0) {
    if (ndim > 1) {
      Real deltamu = 2.0 / (2 * nmu - 1);
      mu2tmp(0) = 1.0 / (3.0 * (2 * nmu - 1));
      for (int i=1; i<nmu; i++) {
        mu2tmp(i) = mu2tmp(i-1) + deltamu;
      }

      Real w2 = 4.0 * mu2tmp(0);
      Real wsum2 = std::sqrt(w2);
      wtmp2(0) = wsum2;

      for (int i=1; i<nmu-2; i++) {
        w2 += deltamu;
        wtmp2(i) = std::sqrt(w2);
        wsum2 += wtmp2(i);
      }

      if (nmu > 2)
        wtmp2(nmu-2) = 2.0*(nmu-1)/3.0 - wsum2;

      wtmp(0) = wtmp2(0);
      Real wsum = wtmp(0);
      for (int i=1; i<nmu-1; ++i) {
        wtmp(i) = wtmp2(i) - wtmp2(i-1);
        wsum += wtmp(i);
      }
      wtmp(nmu-1) = 1.0 - wsum;

      int np = 0;
      int iang = 0;

      // initialize to be zero
      for (int i=0; i<nmu; ++i)
        for (int j=0; j<nmu; ++j)
          pmat(i,j) = 0.0;

      for (int i=0; i<nmu; i++) {
        for (int j=0; j<nmu; j++) {
          for (int k=0; k<nmu; k++) {
            if (i + j + k == nmu - 1) {
              // assign cosines to temporary array grid
              mutmp(iang,0) = std::sqrt(mu2tmp(j));
              mutmp(iang,1) = std::sqrt(mu2tmp(k));
              mutmp(iang,2) = std::sqrt(mu2tmp(i));

              int ip=Permutation(i,j,k,np,pl);
              if (ip == -1) {
                pl(np,0) = i;
                pl(np,1) = j;
                pl(np,2) = k;
                pmat(i,np) += 1.0;
                plab(iang) = np;
                np++;
              } else {
                pmat(i,ip) += 1.0;
                plab(iang) = ip;
              }

              iang++;
            }
          }
        }
      }
      if (nmu > 1) {
        //  Invert matrix of Permutations families */
        InverseMatrix(nmu-1, pmat,pinv);
        // Solve for and assign weights for each Permutation family
        MatrixMult(nmu-1,nmu-1,pinv,wtmp,wpf);

        for (int l=0; l<noct; ++l) {
          for (int i=0; i<n_ang; ++i) {
              wmu(l*n_ang+i) = wpf(plab(i));
          } // end nang
        } // end noct
      } else {
        for (int l=0; l<noct; ++l) {
            wmu(l) = 1.0;
        }
      }
    }

    if (ndim == 1) {
      AthenaArray<Real> mutmp1d, wtmp1d;
      mutmp1d.NewAthenaArray(2*nmu);
      wtmp1d.NewAthenaArray(2*nmu);

      Gauleg(2*nmu, -1.0, 1.0, mutmp1d, wtmp1d);
      for (int i=0; i<2*nmu; ++i) {
        wmu(i) = 0.5 * wtmp1d(i);
      }
      for (int n1=0; n1<n1z; ++n1) {
        for (int i=0; i<2*nmu; ++i) {
          mu(0,0,0,n1,i)=mutmp1d(i);
        }
      }
    } else if (ndim == 2) {
      // for spherical coordinate system, it should be r-theta-phi
      for (int n2=0; n2<n2z; ++n2) {
        for (int n1=0; n1<n1z; ++n1) {
          for (int j=0; j<2; ++j) {
            for (int k=0; k<2; ++k) {
              int l=2*j+k;
              for (int i=0; i<n_ang; ++i) {
                int mi = l*n_ang + i;
                if (k == 0) {
                  mu(axisx,0,n2,n1,mi) =  mutmp(i,0);
                } else {
                  mu(axisx,0,n2,n1,mi) = -mutmp(i,0);
                }
                if (j == 0) {
                  mu(axisy,0,n2,n1,mi) =  mutmp(i,1);
                } else {
                  mu(axisy,0,n2,n1,mi) = -mutmp(i,1);
                }
              }
            }
          }
        }
      }
      // for the angular weight
      for (int i=0; i<noct * n_ang; ++i) {
          wmu(i) *= 0.25;
      }
    } else if (ndim == 3) {
      for (int n3=0; n3<n3z; ++n3) {
        for (int n2=0; n2<n2z; ++n2) {
          for (int n1=0; n1<n1z; ++n1) {
            for (int j=0; j<2; ++j) {
              for (int k=0; k<2; ++k) {
                for (int l=0; l<2; ++l) {
                  int m=4*j+2*k+l;

                  for (int i=0; i<n_ang; ++i) {
                    int mi = m*n_ang + i;
                    if (l == 0) {
                      mu(axisx,n3,n2,n1,mi) =  mutmp(i,0);
                    } else {
                      mu(axisx,n3,n2,n1,mi) = -mutmp(i,0);
                    }
                    if (k == 0) {
                      mu(axisy,n3,n2,n1,mi) =  mutmp(i,1);
                    } else {
                      mu(axisy,n3,n2,n1,mi) = -mutmp(i,1);
                    }
                    if (j == 0) {
                      mu(axisz,n3,n2,n1,mi) =  mutmp(i,2);
                    } else {
                      mu(axisz,n3,n2,n1,mi) = -mutmp(i,2);
                    }
                  } // end i
                } // end l
              } // end k
            } // end j
          } // end n1
        } // end n2
      } // end n3

      for (int i=0; i<noct * n_ang; ++i) {
        wmu(i) *= 0.125;
      }
    } // end nDIM = 3
  } else if (angle_flag == 10) {
    if ((ndim == 1) || (ndim == 3)) {
      msg <<"[Warning]: ang_quad = 10 should be used" <<
      "only for 2D problems. \n" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    Real delmu=1.0/std::sqrt((Real)3.0);
    Real sintheta = std::sqrt((Real)2.0)/std::sqrt((Real)3.0);
    Real phi = 0.5 * PI / (Real) (2*nmu);

    for (int i=0; i<nmu; ++i) {
      mutmp(i,0) = sintheta * cos(phi*(Real)(2*i+1));
      mutmp(i,1) = sintheta * sin(phi*(Real)(2*i+1));
      mutmp(i,2) = delmu;
    }
    for (int n2=0; n2<n2z; ++n2) {
      for (int n1=0; n1<n1z; ++n1) {
        for (int j=0; j<2; ++j) {
          for (int k=0; k<2; ++k) {
            int l=2*j+k;

            for (int i=0; i<n_ang; ++i) {
              int mi = l*n_ang + i;
              if (k == 0) {
                mu(axisx,0,n2,n1,mi) =  mutmp(i,0);
              } else {
                mu(axisx,0,n2,n1,mi) = -mutmp(i,0);
              }

              if (j == 0) {
                mu(axisy,0,n2,n1,mi) =  mutmp(i,1);
              } else {
                mu(axisy,0,n2,n1,mi) = -mutmp(i,1);
              }
            }
          }
        }
      }
    }
    for (int i=0; i<noct*n_ang; ++i) {
      wmu(i) = 0.25/((Real)nmu);
    }
  } else {
    msg << "### FATAL ERROR in function [InitialAngle]" << std::endl
    << "Type of angular discretization unknow: "<< angle_flag << "\n ";
    throw std::runtime_error(msg.str().c_str());
  }

  // Now change the angle cosines to different coordinate systems
  pmy_block->pcoord->ConvertAngle(pmy_block,nang, mu);
  return;
}

// overwrite function to create angles that are defined with respect
// to local coordinate system

// This is the so-called spherical polar tetrad,
// which is ideal for spherical polar coordinate

// nzeta is number of zeta angles in one octant
// npsi is the number of psi angles between 0 and pi
void NRRadiation::AngularGrid(int angle_flag, int nzeta, int npsi) {
  MeshBlock *pmb=pmy_block;
  if (angle_flag == 1) {
    // initialize coordinate direction
    int axisx=0, axisy=1, axisz=2;
    pmy_block->pcoord->AxisDirection(&axisx, &axisy, &axisz);

    // in spherical polar coordinate
    //  *axisx = 1;
    //  *axisy = 2;
    //  *axisz = 0;
    // check the dimension of the problem

    int n1z = pmb->ncells1, n2z = pmb->ncells2, n3z = pmb->ncells3;

    int ndim = 1;
    if (n2z > 1) ndim = 2;
    if (n3z > 1) ndim = 3;

    // separate ghost zones and active zones
    // so that they can be compatible with different angular scheme
    if (nzeta > 0) {
      coszeta_v.NewAthenaArray(2*nzeta);
      zeta_v_full.NewAthenaArray(2*nzeta+2*NGHOST);
      zeta_f_full.NewAthenaArray(2*nzeta+1+2*NGHOST);
      dzeta_v.NewAthenaArray(2*nzeta+2*NGHOST);
      dzeta_f.NewAthenaArray(2*nzeta+1+2*NGHOST);
      coszeta_f.NewAthenaArray(2*nzeta+1);
      len_zeta.NewAthenaArray(2*nzeta); // This id Delta (cos\theta)

      int zs = 0; // ze = 2*nzeta - 1;

      Real dcoszeta = 1.0/nzeta;
      coszeta_f(zs) = 1.0;
      for (int i=1; i<nzeta; ++i) {
        coszeta_f(i+zs) = 1.0 - i *dcoszeta;
      }
      coszeta_f(nzeta+zs) = 0.0;
      for (int i=nzeta+1; i<2*nzeta+1; ++i)
        coszeta_f(i+zs) = -coszeta_f(2*nzeta-i+zs);

      for (int i=0; i<nzeta; ++i) {
        coszeta_v(i+zs) = 0.5*(coszeta_f(i+zs)+coszeta_f(i+zs+1));
      }
    // re-normalize
      Real normalization = 2*nzeta/std::sqrt(4*nzeta*nzeta-1);

      for (int i=0; i<nzeta; ++i) {
        coszeta_v(i+zs) *= normalization;
      }
      for (int i=nzeta; i<2*nzeta; ++i) {
        coszeta_v(i+zs) = -coszeta_v(2*nzeta-i-1+zs);
      }

      for (int i=0; i<nzeta*2; ++i) {
        len_zeta(i) = coszeta_f(i) - coszeta_f(i+1);
      }


      // the following arrays include ghost zones
      // this is used for reconstruction
      for (int i=0; i<2*nzeta; ++i) {
        zeta_v_full(i+NGHOST) = acos(coszeta_v(i));
      }
      for (int i=0; i<2*nzeta+1; ++i) {
        zeta_f_full(i+NGHOST) = acos(coszeta_f(i));
      }
      for (int i=1; i<=NGHOST; ++i) {
        Real delta = zeta_v_full(i+NGHOST) - zeta_v_full(NGHOST+i-1);
        zeta_v_full(NGHOST-i) = zeta_v_full(NGHOST-i+1) - delta;
      }
      for (int i=1; i<=NGHOST; ++i) {
        Real delta = zeta_v_full(2*nzeta+NGHOST-i) - zeta_v_full(2*nzeta+NGHOST-i-1);
        zeta_v_full(2*nzeta+NGHOST+i-1) = zeta_v_full(2*nzeta+NGHOST+i-2) + delta;
      }

      for (int i=0; i<2*nzeta+2*NGHOST-1; ++i) {
        dzeta_v(i) = zeta_v_full(i+1) - zeta_v_full(i);
      }
      for (int i=0; i<2*nzeta+1; ++i) {
        zeta_f_full(i+NGHOST) = acos(coszeta_f(i));
      }

      for (int i=1; i<=NGHOST; ++i) {
        Real delta = zeta_f_full(i+NGHOST) - zeta_f_full(NGHOST+i-1);
        zeta_f_full(NGHOST-i) = zeta_f_full(NGHOST-i+1) - delta;
      }

      for (int i=1; i<=NGHOST; ++i) {
        Real delta = zeta_f_full(2*nzeta+NGHOST-i+1) - zeta_f_full(2*nzeta+NGHOST-i);
        zeta_f_full(2*nzeta+NGHOST+i) = zeta_f_full(2*nzeta+NGHOST+i-1) + delta;
      }

      for (int i=0; i<2*nzeta+2*NGHOST; ++i) {
        dzeta_f(i) = zeta_f_full(i+1)-zeta_f_full(i);
      }

    // end if nzeta > 0
    } else {
      coszeta_v.NewAthenaArray(1);
      zeta_v_full.NewAthenaArray(1);
      zeta_f_full.NewAthenaArray(1);
      dzeta_v.NewAthenaArray(1);
      dzeta_f.NewAthenaArray(1);
      coszeta_f.NewAthenaArray(1);
      len_zeta.NewAthenaArray(1); // This id Delta (cos\theta)

      coszeta_f(0) = 1.0;
      coszeta_v(0) = 1.0;

      len_zeta(0) = 1.0;

      zeta_v_full(0) = 0.0;
      zeta_f_full(0) = 0.0;

      dzeta_v(0) = 1.0;
      dzeta_f(0) = 1.0;
    }


    // construct psi angles
    // for 2D problem in spherical polar, npsi = 1
    if (npsi > 0) {
      psi_v.NewAthenaArray(2*npsi);
      psi_f.NewAthenaArray(2*npsi+1);
      sin_psi_f.NewAthenaArray(2*npsi+1);
      len_psi.NewAthenaArray(2*npsi);
      psi_v_full.NewAthenaArray(2*npsi+2*NGHOST);
      psi_f_full.NewAthenaArray(2*npsi+1+2*NGHOST);
      dpsi_v.NewAthenaArray(2*npsi+2*NGHOST);
      dpsi_f.NewAthenaArray(2*npsi+2*NGHOST+1);

      Real dpsi=PI/npsi;

      for (int i=0; i<npsi; ++i) {
        psi_v(i) = (i+0.5)*dpsi;
      }
      for (int i=npsi; i<2*npsi; ++i)
        psi_v(i)=2.0*PI-psi_v(2*npsi-i-1);

      psi_f(0) = 0.0;
      for (int i=1; i<2*npsi; ++i)
        psi_f(i)=2*psi_v(i-1) - psi_f(i-1);
      psi_f(2*npsi) = 2*PI;

      for (int i=0; i<2*npsi+1; ++i)
        sin_psi_f(i) = sin(psi_f(i));

      for (int i=0; i<2*npsi; ++i)
        len_psi(i) = psi_f(i+1) - psi_f(i);

      //  now the arrays with ghost zones

      for (int i=0; i<2*npsi; ++i) {
        psi_v_full(i+NGHOST) = psi_v(i);
      }
      for (int i=0; i<2*npsi+1; ++i) {
        psi_f_full(i+NGHOST) = psi_f(i);
      }
      for (int i=1; i<=NGHOST; ++i) {
        psi_v_full(NGHOST-i) = psi_v(2*npsi-i)-2*PI;
      }
      for (int i=1; i<=NGHOST; ++i) {
        psi_v_full(2*npsi+NGHOST+i-1) = psi_v(i-1) + 2*PI;
      }

      for (int i=0; i<2*npsi+2*NGHOST-1; ++i) {
        dpsi_v(i) = psi_v_full(i+1) - psi_v_full(i);
      }

      for (int i=0; i<2*npsi+1; ++i) {
        psi_f_full(i+NGHOST) = psi_f(i);
      }

      for (int i=1; i<=NGHOST; ++i) {
        psi_f_full(NGHOST-i) = psi_f(2*npsi-i)-2*PI;
      }

      for (int i=1; i<=NGHOST; ++i) {
        psi_f_full(2*npsi+NGHOST+i) = psi_f(i) + 2*PI;
      }

      for (int i=0; i<2*npsi+2*NGHOST; ++i) {
        dpsi_f(i) = psi_f_full(i+1)-psi_f_full(i);
      }
    } // end if npsi > 0

    // for 1D problem
    if (ndim == 1) {
      for (int i=0; i<n1z; ++i) {
         for (int n=0; n<2*nzeta; ++n) {
            //x,k,j,i,n
            mu(0,0,0,i,n) = coszeta_v(n);
         }
      }
    } else if (ndim == 2) {
      for (int j=0; j<n2z; ++j) {
        for (int i=0; i<n1z; ++i) {
          if (npsi == 1) {
            for (int n=0; n<2*nzeta; ++n) {
              for (int m=0; m<2*npsi; ++m) {
                int ang_num = n*(2*npsi)+m;
                Real sinzeta_v = std::sqrt(1.0 - coszeta_v(n)
                                      * coszeta_v(n));
                mu(axisz,0,j,i,ang_num) = coszeta_v(n);
                if (m==0)
                  mu(axisx,0,j,i,ang_num) = sinzeta_v;
                else
                  mu(axisx,0,j,i,ang_num) = -sinzeta_v;
              }
            }
          } else {// the case in x -y plane
            if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
              // in spherical polar, 2D, we still need 3D angular grid
              for (int n=0; n<2*nzeta; ++n) {
                for (int m=0; m<2*npsi; ++m) {
                  int ang_num = n*(2*npsi)+m;
                  Real sinzeta_v = std::sqrt(1.0 - coszeta_v(n)
                                      * coszeta_v(n));
                  mu(axisx,0,j,i,ang_num) = sinzeta_v * cos(psi_v(m));
                  mu(axisy,0,j,i,ang_num) = sinzeta_v * sin(psi_v(m));
                  mu(axisz,0,j,i,ang_num) = coszeta_v(n);
                }
              }
            } else {
              for (int m=0; m<2*npsi; ++m) {
                mu(0,0,j,i,m) = cos(psi_v(m));
                mu(1,0,j,i,m) = sin(psi_v(m));
              }
            }
          }
        } // end i
      } // end j

    } else {
      // 3D case
      for (int k=0; k<n3z; ++k) {
        for (int j=0; j<n2z; ++j) {
          for (int i=0; i<n1z; ++i) {
            for (int n=0; n<2*nzeta; ++n) {
              for (int m=0; m<2*npsi; ++m) {
                int ang_num = n*(2*npsi)+m;
                Real sinzeta_v = std::sqrt(1.0 - coszeta_v(n)
                                      * coszeta_v(n));
                mu(axisx,k,j,i,ang_num) = sinzeta_v * cos(psi_v(m));
                mu(axisy,k,j,i,ang_num) = sinzeta_v * sin(psi_v(m));
                mu(axisz,k,j,i,ang_num) = coszeta_v(n);
              }
            }
          }
        }
      }
    }
    // equal weight for all cases
    for (int n=0; n<nang; ++n) {
      wmu(n) = 1.0/nang;
    }
  } // end ang_flag==1
}
