//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file gravity_fluxes.cpp
//  \brief Calculate gravity fluxes

// Athena++ headers
#include "hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../gravity/gravity.hpp"


//----------------------------------------------------------------------------------------
//! \fn void Hydro::AddGravityFlux(void)
//  \brief Add gravity flux using new potential directly

void Hydro::AddGravityFlux(void) {
  MeshBlock *pmb=pmy_block;
  Coordinates *pco=pmb->pcoord;
  Real four_pi_G=pmb->pgrav->four_pi_G, grav_mean_rho=pmb->pgrav->grav_mean_rho;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  AthenaArray<Real> &x1flux=flux[X1DIR];
  AthenaArray<Real> &x2flux=flux[X2DIR];
  AthenaArray<Real> &x3flux=flux[X3DIR];
  AthenaArray<Real> &phi=pmb->pgrav->phi;
  Real phil;
  Real gxl,gyl,gzl;

  // i-direction
  for (int k=ks; k<=ke; ++k) {
    Real dx3 = pco->dx3v(k);
    for (int j=js; j<=je; ++j) {
      Real dx2 = pco->dx2v(j);
      for (int i=is; i<=ie+1; ++i) {
        Real dx1 = pco->dx1v(i);
        phil = 0.5*(phi(k,j,i-1)+phi(k,j,i));
        // gx, gy, and gz centered at L and R x1-faces
        gxl =       (phi(k,j  ,i-1) - phi(k,j  ,i  ))/dx1;
        if (pmb->block_size.nx2 > 1) { // 2D or 3D
          gyl = 0.25*((phi(k,j-1,i-1) - phi(k,j+1,i-1)) +
                      (phi(k,j-1,i  ) - phi(k,j+1,i  )))/dx2;
          if (pmb->block_size.nx3 > 1) { // 3D
            gzl = 0.25*((phi(k-1,j,i-1) - phi(k+1,j,i-1)) +
                        (phi(k-1,j,i  ) - phi(k+1,j,i  )))/dx3;
          }
        }
        // momentum fluxes in x1-dir.
        // 2nd term is needed only if Jean's swindle used
        x1flux(IM1,k,j,i) += 0.5*(gxl*gxl-gyl*gyl-gzl*gzl)/four_pi_G
                          +  grav_mean_rho*phil;
        x1flux(IM2,k,j,i) += gxl*gyl/four_pi_G;
        x1flux(IM3,k,j,i) += gxl*gzl/four_pi_G;
        // energy source term is included as a source term separately
        // see hydro/srctrm/gravity.cpp
      }
    }
  }
  // j-direction
  if (pmb->block_size.nx2 > 1) {
    for (int k=ks; k<=ke; ++k) {
      Real dx3 = pco->dx3v(k);
      for (int j=js; j<=je+1; ++j) {
        Real dx2 = pco->dx2v(j);
        for (int i=is; i<=ie; i++) {
          Real dx1 = pco->dx1v(i);
          phil = 0.5*(phi(k,j-1,i)+phi(k,j,i));
          // gx, gy, and gz centered at L and R x2-faces
          gxl = 0.25*((phi(k,j-1,i-1) - phi(k,j-1,i+1)) +
                      (phi(k,j  ,i-1) - phi(k,j  ,i+1)))/dx1;
          gyl =       (phi(k,j-1,i  ) - phi(k,j  ,i  ))/dx2;
          if (pmb->block_size.nx3 > 1) { // 3D
            gzl = 0.25*((phi(k-1,j-1,i) - phi(k+1,j-1,i)) +
                        (phi(k-1,j  ,i) - phi(k+1,j  ,i)))/dx3;
          }
          // momentum fluxes in x2-dir.
          // 2nd term is needed only if Jean's swindle used
          x2flux(IM1,k,j,i) += gyl*gxl/four_pi_G;
          x2flux(IM2,k,j,i) += 0.5*(gyl*gyl-gxl*gxl-gzl*gzl)/four_pi_G
                            +  grav_mean_rho*phil;
          x2flux(IM3,k,j,i) += gyl*gzl/four_pi_G;
          // energy source term is included as a source term separately
          // see hydro/srctrm/gravity.cpp
        }
      }
    }
  }
  // k-direction
  if (pmb->block_size.nx3 > 1) {
    for (int k=ks; k<=ke+1; ++k) {
      Real dx3 = pco->dx3v(k);
      for (int j=js; j<=je; ++j) {
        Real dx2 = pco->dx2v(j);
        for (int i=is; i<=ie; i++) {
          Real dx1 = pco->dx1v(i);
          phil = 0.5*(phi(k-1,j,i)+phi(k,j,i));
          // gx, gy, and gz centered at L and R x3-faces
          gxl = 0.25*((phi(k-1,j,i-1) - phi(k-1,j,i+1)) +
                      (phi(k  ,j,i-1) - phi(k  ,j,i+1)))/dx1;
          gyl = 0.25*((phi(k-1,j-1,i) - phi(k-1,j+1,i)) +
                      (phi(k  ,j-1,i) - phi(k  ,j+1,i)))/dx2;
          gzl =       (phi(k-1,j  ,i) - phi(k  ,j  ,i))/dx3;
          // momentum fluxes in x3-dir.
          // 2nd term is needed only if Jean's swindle used
          x3flux(IM1,k,j,i) += gzl*gxl/four_pi_G;
          x3flux(IM2,k,j,i) += gzl*gyl/four_pi_G;
          x3flux(IM3,k,j,i) += 0.5*(gzl*gzl-gxl*gxl-gyl*gyl)/four_pi_G
                            +  grav_mean_rho*phil;
          // energy source term is included as a source term separately
          // see hydro/srctrm/gravity.cpp
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::AddGravityFluxWithGflx(void)
//  \brief Adds gravity flux to hydro flux using pre-calculated gflx.
//         Implemented for one poisson solver call. (old Athena style)

void Hydro::AddGravityFluxWithGflx(void) {
  MeshBlock *pmb=pmy_block;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> &x1flux=flux[X1DIR];
  AthenaArray<Real> &x2flux=flux[X2DIR];
  AthenaArray<Real> &x3flux=flux[X3DIR];
  AthenaArray<Real> &x1gflx=gflx[X1DIR];
  AthenaArray<Real> &x2gflx=gflx[X2DIR];
  AthenaArray<Real> &x3gflx=gflx[X3DIR];

  for (int n=IM1; n<=IM3; ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          x1flux(n,k,j,i) += x1gflx(n,k,j,i);
          if (i==ie) x1flux(n,k,j,i+1) += x1gflx(n,k,j,i+1);
          if (pmb->block_size.nx2 > 1) {
            x2flux(n,k,j,i) += x2gflx(n,k,j,i);
            if (j==je) x2flux(n,k,j+1,i) += x2gflx(n,k,j+1,i);
          }
          if (pmb->block_size.nx3 > 1) {
            x3flux(n,k,j,i) += x3gflx(n,k,j,i);
            if (k==ke) x3flux(n,k+1,j,i) += x3gflx(n,k+1,j,i);
          }
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Hydro::CalculateGravityFlux(AthenaArray<Real> &phi_in)
//  \brief Calcuates gravity flux to hydro flux
//         Implemented for one poisson solver call. (old Athena style)

void Hydro::CalculateGravityFlux(AthenaArray<Real> &phi_in) {
  MeshBlock *pmb=pmy_block;
  Coordinates *pco=pmb->pcoord;
  Real four_pi_G=pmb->pgrav->four_pi_G, grav_mean_rho=pmb->pgrav->grav_mean_rho;

  AthenaArray<Real> &x1flux=flux[X1DIR];
  AthenaArray<Real> &x2flux=flux[X2DIR];
  AthenaArray<Real> &x3flux=flux[X3DIR];
  AthenaArray<Real> &x1gflx=gflx[X1DIR];
  AthenaArray<Real> &x2gflx=gflx[X2DIR];
  AthenaArray<Real> &x3gflx=gflx[X3DIR];
  AthenaArray<Real> &x1gflx_old=gflx_old[X1DIR];
  AthenaArray<Real> &x2gflx_old=gflx_old[X2DIR];
  AthenaArray<Real> &x3gflx_old=gflx_old[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  Real phic, phil;
  Real gxl,gyl,gzl;
//----------------------------------------------------------------------------------------
// i-direction
// set the loop limits
  jl=js, ju=je, kl=ks, ku=ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (pmb->block_size.nx2 > 1) {
      if (pmb->block_size.nx3 == 1) // 2D
        jl=js-1, ju=je+1, kl=ks, ku=ke;
      else // 3D
        jl=js-1, ju=je+1, kl=ks-1, ku=ke+1;
    }
  }

  gxl = 0.0, gyl = 0.0, gzl = 0.0;
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=is; i<=ie+1; ++i) {
        Real dx1 = pco->dx1v(i);
        Real dx2 = pco->dx2v(j);
        Real dx3 = pco->dx3v(k);
        phic = phi_in(k,j,i);
        phil = 0.5*(phi_in(k,j,i-1) + phi_in(k,j,i  ));
        // gx, gy, and gz centered at L and R x1-faces
        gxl =       (phi_in(k,j  ,i-1) - phi_in(k,j  ,i  ))/dx1;
        if (pmb->block_size.nx2 > 1) { // 2D or 3D
          gyl = 0.25*((phi_in(k,j-1,i-1) - phi_in(k,j+1,i-1)) +
                      (phi_in(k,j-1,i  ) - phi_in(k,j+1,i  )))/dx2;
          if (pmb->block_size.nx3 > 1) // 3D
            gzl = 0.25*((phi_in(k-1,j,i-1) - phi_in(k+1,j,i-1)) +
                        (phi_in(k-1,j,i  ) - phi_in(k+1,j,i  )))/dx3;
        }

        x1gflx_old(IM1,k,j,i) = x1gflx(IM1,k,j,i);
        x1gflx_old(IM2,k,j,i) = x1gflx(IM2,k,j,i);
        x1gflx_old(IM3,k,j,i) = x1gflx(IM3,k,j,i);
        // momentum fluxes in x1-dir.
        // 2nd term is needed only if Jean's swindle used
        x1gflx(IM1,k,j,i) = 0.5*(gxl*gxl-gyl*gyl-gzl*gzl)/four_pi_G
                          + grav_mean_rho*phil;
        x1gflx(IM2,k,j,i) = gxl*gyl/four_pi_G;
        x1gflx(IM3,k,j,i) = gxl*gzl/four_pi_G;
        // energy source term is included as a source term separately
        // see hydro/srctrm/gravity.cpp
      }
    }
  }
//----------------------------------------------------------------------------------------
// j-direction

  gxl = 0.0, gyl = 0.0, gzl = 0.0;
  if (pmb->block_size.nx2 > 1) {
    // set the loop limits
    il=is, iu=ie, kl=ks, ku=ke;
    if (MAGNETIC_FIELDS_ENABLED) {
      if (pmb->block_size.nx3 == 1) // 2D
        il=is-1, iu=ie+1, kl=ks, ku=ke;
      else // 3D
        il=is-1, iu=ie+1, kl=ks-1, ku=ke+1;
    }
    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
        for (int i=il; i<=iu; i++) {
          Real dx1 = pco->dx1v(i);
          Real dx2 = pco->dx2v(j);
          Real dx3 = pco->dx3v(k);
          phic = phi_in(k,j,i);
          phil = 0.5*(phi_in(k,j-1,i)+phi_in(k,j  ,i));

          // gx, gy, and gz centered at L and R x2-faces
          gxl = 0.25*((phi_in(k,j-1,i-1) - phi_in(k,j-1,i+1)) +
                      (phi_in(k,j  ,i-1) - phi_in(k,j  ,i+1)))/dx1;
          gyl =       (phi_in(k,j-1,i  ) - phi_in(k,j  ,i  ))/dx2;
          if (pmb->block_size.nx3 > 1)
            gzl = 0.25*((phi_in(k-1,j-1,i) - phi_in(k+1,j-1,i)) +
                        (phi_in(k-1,j  ,i) - phi_in(k+1,j  ,i)))/dx3;

          x2gflx_old(IM1,k,j,i) = x2gflx(IM1,k,j,i);
          x2gflx_old(IM2,k,j,i) = x2gflx(IM2,k,j,i);
          x2gflx_old(IM3,k,j,i) = x2gflx(IM3,k,j,i);
          // momentum fluxes in x2-dir.
          // 2nd term is needed only if Jean's swindle used
          x2gflx(IM1,k,j,i) = gyl*gxl/four_pi_G;
          x2gflx(IM2,k,j,i) = 0.5*(gyl*gyl-gxl*gxl-gzl*gzl)/four_pi_G
                            + grav_mean_rho*phil;
          x2gflx(IM3,k,j,i) = gyl*gzl/four_pi_G;
          // energy source term is included as a source term separately
          // see hydro/srctrm/gravity.cpp
        }
      }
    }
  } // 2D or 3D
//----------------------------------------------------------------------------------------
// k-direction

  gxl = 0.0, gyl = 0.0, gzl = 0.0;
  if (pmb->block_size.nx3 > 1) {
    // set the loop limits
    il=is, iu=ie, jl=js, ju=je;
    if (MAGNETIC_FIELDS_ENABLED)
      il=is-1, iu=ie+1, jl=js-1, ju=je+1;

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; i++) {
          Real dx1 = pco->dx1v(i);
          Real dx2 = pco->dx2v(j);
          Real dx3 = pco->dx3v(k);
          phic = phi_in(k,j,i);
          phil = 0.5*(phi_in(k-1,j,i)+phi_in(k  ,j,i));

          // gx, gy, and gz centered at L and R x3-faces
          gxl = 0.25*((phi_in(k-1,j,i-1) - phi_in(k-1,j,i+1)) +
                      (phi_in(k  ,j,i-1) - phi_in(k  ,j,i+1)))/dx1;
          gyl = 0.25*((phi_in(k-1,j-1,i) - phi_in(k-1,j+1,i)) +
                      (phi_in(k  ,j-1,i) - phi_in(k  ,j+1,i)))/dx2;
          gzl =       (phi_in(k-1,j  ,i) - phi_in(k  ,j  ,i))/dx3;

          x3gflx_old(IM1,k,j,i) = x3gflx(IM1,k,j,i);
          x3gflx_old(IM2,k,j,i) = x3gflx(IM2,k,j,i);
          x3gflx_old(IM3,k,j,i) = x3gflx(IM3,k,j,i);
          // momentum fluxes in x3-dir.
          // 2nd term is needed only if Jean's swindle used
          x3gflx(IM1,k,j,i) = gzl*gxl/four_pi_G;
          x3gflx(IM2,k,j,i) = gzl*gyl/four_pi_G;
          x3gflx(IM3,k,j,i) = 0.5*(gzl*gzl-gxl*gxl-gyl*gyl)/four_pi_G
                            + grav_mean_rho*phil;
          // energy source term is included as a source term separately
          // see hydro/srctrm/gravity.cpp
        }
      }
    }
  } //3D
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void Hydro::CorrectGravityFlux(void)
//  \brief Correct the fluxes using the updated gravitational fluxes
//         Implemented for one poisson solver call. (old Athena style)

void Hydro::CorrectGravityFlux(void) {
  MeshBlock *pmb=pmy_block;
  AthenaArray<Real> &x1gflx=gflx[X1DIR];
  AthenaArray<Real> &x2gflx=gflx[X2DIR];
  AthenaArray<Real> &x3gflx=gflx[X3DIR];
  AthenaArray<Real> &x1gflx_old=gflx_old[X1DIR];
  AthenaArray<Real> &x2gflx_old=gflx_old[X2DIR];
  AthenaArray<Real> &x3gflx_old=gflx_old[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> x1area, x2area, x2area_p1, x3area, x3area_p1, vol, dflx;
  x1area.InitWithShallowCopy(x1face_area_);
  x2area.InitWithShallowCopy(x2face_area_);
  x2area_p1.InitWithShallowCopy(x2face_area_p1_);
  x3area.InitWithShallowCopy(x3face_area_);
  x3area_p1.InitWithShallowCopy(x3face_area_p1_);
  vol.InitWithShallowCopy(cell_volume_);
  dflx.InitWithShallowCopy(dflx_);

#pragma omp for schedule(static)
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      // calculate x1-flux divergence
      pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);
      for (int n=IM1; n<=IM3; ++n) {
        for (int i=is; i<=ie; ++i) {
          dflx(n,i) = x1area(i+1)*(x1gflx(n,k,j,i+1)-x1gflx_old(n,k,j,i+1))
                      - x1area(i)*(x1gflx(n,k,j,i)-x1gflx_old(n,k,j,i));
        }
      }

      // calculate x2-flux divergence
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
        for (int n=IM1; n<=IM3; ++n) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) += x2area_p1(i)*(x2gflx(n,k,j+1,i)-x2gflx_old(n,k,j+1,i))
                       - x2area(i)*(x2gflx(n,k,j,i)-x2gflx_old(n,k,j,i));
          }
        }
      }

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);
        for (int n=IM1; n<=IM3; ++n) {
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            dflx(n,i) += x3area_p1(i)*(x3gflx(n,k+1,j,i)-x3gflx_old(n,k+1,j,i))
                       - x3area(i)*(x3gflx(n,k,j,i)-x3gflx_old(n,k,j,i));
          }
        }
      }

      // update conserved variables
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      for (int n=IM1; n<=IM3; ++n) {
        for (int i=is; i<=ie; ++i) {
          u(n,k,j,i) -= 0.5*(pmb->pmy_mesh->dt)*dflx(n,i)/vol(i);
        }
      }
    }
  }

  return;
}
