//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file wave_extract.cpp
//  \brief implementation of functions in the WaveExtract classes

#include <cstdio>
#include <stdexcept>
#include <sstream>

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#include "wave_extract.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/spherical_grid.hpp"
#include "z4c.hpp"

WaveExtract::WaveExtract(Mesh * pmesh, ParameterInput * pin, int n, int res_flag):
    pmesh(pmesh), pofile(NULL) {
  int nlev = pin->GetOrAddInteger("z4c", "extraction_nlev", 3);
  Real rad;
  std::string rad_parname;
  rad_parname = "extraction_radius_";
  std::string n_str = std::to_string(n);
  rad_parname += n_str;
  rad = pin->GetOrAddReal("z4c", rad_parname, 10.0);
  rad_id = n;
  ofname = pin->GetOrAddString("z4c", "extract_filename", "wave");
  root = pin->GetOrAddInteger("z4c", "mpi_root", 0);
  lmax = pin->GetOrAddInteger("z4c", "lmax", 2);
  psphere = new SphericalGrid(nlev, rad);
  ofname += n_str;
  ofname += ".txt";
 int np = psphere->NumVertices();
//  Real theta, phi, x, y, z;
  //printf("np = %d\n",np);
//  for(int ip=0;ip<np;++ip){
//    psphere->SphericalGrid::PositionPolar(ip,&theta,&phi);
//    psphere->SphericalGrid::Position(ip,&x,&y,&z);
//    printf("theta = %.16f, phi = %.16f, ip = %d, x = %.16f, y = %.16f, z = %.16f\n",theta,phi,ip,x,y,z);
//}


#ifdef MPI_PARALLEL
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ioproc = (root == rank);
#else
  ioproc = true;
#endif
  if(res_flag == 0){
    if (ioproc) {
      pofile = fopen(ofname.c_str(), "w");
      if (NULL == pofile) {
        std::stringstream msg;
        msg << "### FATAL ERROR in WaveExtract constructor" << std::endl;
        msg << "Could not open file '" << ofname << "' for writing!";
        throw std::runtime_error(msg.str().c_str());
      }
    fprintf(pofile, "# 1:iter 2:time 3:l=2 m=-2 R 4: l=2 m=-2 I 5: l=2 m=-1 R 6: l=2 m=-1 I 7: l=2 m=0 R 8: l=2 m=0 I 9: l=2 m=1 R 10: l=2 m=1 I 11: l=2 m=2 R 12: l=2 m=2 I\n");
    }
   } else if(res_flag == 1){
     if (ioproc) {
       pofile = fopen(ofname.c_str(), "a");
       if (NULL == pofile) {
         std::stringstream msg;
         msg << "### FATAL ERROR in WaveExtract constructor" << std::endl;
         msg << "Could not open file '" << ofname << "' for writing!";
         throw std::runtime_error(msg.str().c_str());
       }
     }
   }
}

WaveExtract::~WaveExtract() {
  delete psphere;
  if (ioproc) {
    fclose(pofile);
  }
}

void WaveExtract::ReduceMultipole() {
  psi.NewAthenaArray(lmax-1,2*(lmax)+1,2);
  psi.ZeroClear();
  MeshBlock const * pmb = pmesh->pblock;
  while (pmb != NULL) {
    for(int l=2;l<lmax+1;++l){
      for(int m=-l;m<l+1;++m){
        psi(l-2,m+l,0) +=pmb->pwave_extr_loc[rad_id]->psi(l-2,m+l,0);
        psi(l-2,m+l,1) +=pmb->pwave_extr_loc[rad_id]->psi(l-2,m+l,1);
      }
    }
    pmb = pmb->next;
  }
#ifdef MPI_PARALLEL
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (root == rank) {
    for(int l=2;l<lmax+1;++l){
      for(int m=-l;m<l+1;++m){
      MPI_Reduce(MPI_IN_PLACE, &psi(l-2,m+l,0), 1, MPI_ATHENA_REAL, MPI_SUM, root, MPI_COMM_WORLD);
      MPI_Reduce(MPI_IN_PLACE, &psi(l-2,m+l,1), 1, MPI_ATHENA_REAL, MPI_SUM, root, MPI_COMM_WORLD);
      }
    }
  }
  else {
    for(int l=2;l<lmax+1;++l){
      for(int m=-l;m<l+1;++m){
      MPI_Reduce(&psi(l-2,m+l,0), &psi(l-2,m+l,0), 1, MPI_ATHENA_REAL, MPI_SUM, root, MPI_COMM_WORLD);
      MPI_Reduce(&psi(l-2,m+l,1), &psi(l-2,m+l,1), 1, MPI_ATHENA_REAL, MPI_SUM, root, MPI_COMM_WORLD);
      }
    }
  }

#endif
}

void WaveExtract::Write(int iter, Real time) const {
  if (ioproc) {
      fprintf(pofile, "%d %.*g ", iter, FPRINTF_PREC, time);
      for(int l=2;l<lmax+1;++l){
        for(int m=-l;m<l+1;++m){
          fprintf(pofile, "%.*g %.*g ",
                  FPRINTF_PREC, psi(l-2,m+l,0),
                  FPRINTF_PREC, psi(l-2,m+l,1));
        }
      }
      fprintf(pofile, "\n");
      fflush(pofile);
  }
}

WaveExtractLocal::WaveExtractLocal(SphericalGrid * psphere, MeshBlock * pmb, ParameterInput * pin, int n) {
  std::string rad_parname;
  rad_parname = "extraction_radius_";
  std::string n_str = std::to_string(n);
  rad_parname += n_str;
  rad = pin->GetOrAddReal("z4c", rad_parname.c_str(), 10.0);
  lmax = pin->GetOrAddInteger("z4c", "lmax",2);
  ppatch = new SphericalPatch(psphere, pmb, SphericalPatch::vertex);
  datareal.NewAthenaArray(ppatch->NumPoints());
  dataim.NewAthenaArray(ppatch->NumPoints());
  weight.NewAthenaArray(ppatch->NumPoints());
  for (int ip = 0; ip < ppatch->NumPoints(); ++ip) {
    weight(ip) = ppatch->psphere->ComputeWeight(ppatch->idxMap(ip));
    weight(ip) /= rad*rad;
  }
}

WaveExtractLocal::~WaveExtractLocal() {
  delete ppatch;
}

void WaveExtractLocal::Decompose_multipole(AthenaArray<Real> const & u_R, AthenaArray<Real> const & u_I) {
    ppatch->InterpToSpherical(u_R, &datareal);
    ppatch->InterpToSpherical(u_I, &dataim);
    Real theta, phi, ylmR, ylmI,x,y,z;
    psi.NewAthenaArray(lmax-1,2*(lmax)+1,2);
    psi.ZeroClear();
//        for (int ip = 0; ip < ppatch->NumPoints(); ++ip) {
//     ppatch->psphere->GeodesicGrid::PositionPolar(ppatch->idxMap(ip),&theta,&phi);
//    ppatch->psphere->SphericalGrid::Position(ppatch->idxMap(ip),&x,&y,&z);
//   printf("theta = %.16f, phi = %.16f, ip_patch = %d, ip_global = %d, x = %.16f, y = %.16f, z = %.16f\n",theta,phi,ip,ppatch->idxMap(ip),x,y,z);
//}
    for (int l = 2; l < lmax+1; ++l){
      for (int m = -l; m < l+1 ; ++m){
        psilmR=0.0;
        psilmI=0.0;
          for (int ip = 0; ip < ppatch->NumPoints(); ++ip) {
            ppatch->psphere->GeodesicGrid::PositionPolar(ppatch->idxMap(ip),&theta,&phi);
            swsh(&ylmR,&ylmI,l,m,theta,phi);
            psilmR += datareal(ip)*weight(ip)*ylmR + dataim(ip)*weight(ip)*ylmI;
            psilmI += dataim(ip)*weight(ip)*ylmR -datareal(ip)*weight(ip)*ylmI;
          }
        psi(l-2,m+l,0) = psilmR;
        psi(l-2,m+l,1) = psilmI;
      }
   }
}

//Factorial
Real WaveExtractLocal::fac(Real n){
 if(n==0 || n==1){
   return 1.0;
 }
 else{
   n=n*fac(n-1);
   return n;
 }
}

//Calculate spin weighted spherical harmonics sw=-2 using Wigner-d matrix notation see e.g. Eq II.7, II.8 in 0709.0093
void WaveExtractLocal::swsh(Real * ylmR, Real * ylmI, int l, int m, Real theta, Real phi){
  Real wignerd = 0;
  int k1,k2,k;
  k1 = std::max(0, m-2);
  k2 = std::min(l+m,l-2);
  for (k = k1; k<k2+1; ++k){
    wignerd += pow((-1),k)*sqrt(fac(l+m)*fac(l-m)*fac(l+2)*fac(l-2))*pow(std::cos(theta/2.0),2*l+m-2-2*k)*pow(std::sin(theta/2.0),2*k+2-m)/(fac(l+m-k)*fac(l-2-k)*fac(k)*fac(k+2-m));
  }
  *ylmR = sqrt((2*l+1)/(4*M_PI))*wignerd*std::cos(m*phi);
  *ylmI = sqrt((2*l+1)/(4*M_PI))*wignerd*std::sin(m*phi);
}

