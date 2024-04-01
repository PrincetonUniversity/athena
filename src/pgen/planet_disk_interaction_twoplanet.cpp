#include <algorithm>  // min
#include <cmath>      // sqrt
#include <cstdlib>    // srand
#include <cstring>    // strcmp()
#include <fstream>
#include <iostream>   // endl
#include <limits>
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
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../parameter_input.hpp"

using namespace std;

namespace {
void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k);
Real DenProfileCyl(const Real rad, const Real phi, const Real z);
Real PoverR(const Real rad, const Real phi, const Real z);
Real VelProfileCyl(const Real rad, const Real phi, const Real z);
// problem parameters which are useful to make global to this file
Real gm0, r0, rho0, dslope, p0_over_r0, pslope, gamma_gas, gm_planet, alpha, nu_iso, scale, z, phi, r, rp, rp2, phip, d, dfloor, Omega0, cosine_term, sine_term, epsilon, R_H;
} // namespace

// User-defined boundary conditions for disk simulations
void Steady_State_Inner(MeshBlock *pmb, Coordinates *pco,
                  AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void Steady_State_Outer(MeshBlock *pmb, Coordinates *pco,
                  AthenaArray<Real> &prim, FaceField &b,
                  Real time, Real dt,
                  int il, int iu, int jl, int ju, int kl, int ku, int ngh);                  
void DiskInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);
void DiskOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//! \brief Function to initialize problem-specific data in mesh class.  Can also be used
//! to initialize variables which are global to (and therefore can be passed to) other
//! functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  Real x1, x2, x3;
  // Get parameters for gravitatonal potential of central star mass
  gm0 = pin->GetOrAddReal("problem","GM",0.0);
  r0 = pin->GetOrAddReal("problem","r0",1.0);

  // Get parameters for initial density and velocity
  rho0 = pin->GetReal("problem","rho0");
  dslope = pin->GetOrAddReal("problem","dslope",0.0);

  // Get parameters for gravitational potential of orbiting protoplanet
  gm_planet = pin -> GetOrAddReal("problem", "planetgm", 0.0);
  rp = pin -> GetOrAddReal("problem", "ptosr", 1.0);
  rp2 = pin -> GetOrAddReal("problem", "ptosr2", 1.1);

  // Get viscosity parameters and scale ratio
  alpha = pin -> GetOrAddReal("problem", "alpha", 0.0);
  nu_iso = pin -> GetOrAddReal("problem", "nu_iso", 0.0);
  scale = pin -> GetOrAddReal("hudro", "iso_sound_speed", 0.04);

  // Get parameters of initial pressure and cooling parameters
  if (NON_BAROTROPIC_EOS) {
    p0_over_r0 = pin->GetOrAddReal("problem","p0_over_r0",0.0025);
    pslope = pin->GetOrAddReal("problem","pslope",0.0);
    gamma_gas = pin->GetReal("hydro","gamma");
  } else {
    p0_over_r0=SQR(pin->GetReal("hydro","iso_sound_speed"));
  }
  Real float_min = std::numeric_limits<float>::min();
  dfloor=pin->GetOrAddReal("hydro","dfloor",(1024*(float_min)));

  Omega0 = pin->GetOrAddReal("orbital_advection","Omega0",0.0);

  // enroll user-defined boundary condition
  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, Steady_State_Inner);
  }
  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, Steady_State_Outer);
  }
  if (mesh_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x2, DiskInnerX2);
  }
  if (mesh_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x2, DiskOuterX2);
  }
  if (mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x3, DiskInnerX3);
  }
  if (mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x3, DiskOuterX3);
  }
  
  void Planet(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
              const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc,
              AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar);
  EnrollUserExplicitSourceFunction(Planet);

  void Viscosity(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, 
            int is, int ie, int js, int je, int ks, int ke);
  EnrollViscosityCoefficient(Viscosity);

  Real Torque(MeshBlock *pmb, int iout);
  Real Torque2(MeshBlock *pmb, int iout);
  AllocateUserHistoryOutput(2);
  EnrollUserHistoryOutput(0, Torque, "first planet torque");
  EnrollUserHistoryOutput(1, Torque2, "secondplanet torque");
  return;
}
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
    AllocateUserOutputVariables(2);
    return;
}
  
//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//! \brief Initializes Keplerian accretion disk.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real den, vel;
  Real x1,x2,x3;

  OrbitalVelocityFunc &vK = porb->OrbitalVelocity;
  for (int k=ks; k<=ke; ++k) {
    z = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      phi = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        r = pcoord->x1v(i);
        GetCylCoord(pcoord,rad,phi,z,i,j,k); // convert to cylindrical coordinates
        Real surface_density = rho0 / sqrt(r/rp);
        Real v_r = -3.0/2.0 * alpha * pow(scale,2) * sqrt((gm0+gm_planet)/r);
        Real v_phi = r * sqrt(1-0.5*pow(scale,2)) * sqrt(gm0+gm_planet)* sqrt(1 / pow(r,3));
        phydro->u(IDN,k,j,i) = surface_density;
        phydro->u(IM1,k,j,i) = surface_density * v_r;
        if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          phydro->u(IM2,k,j,i) = surface_density * v_phi;
          phydro->u(IM3,k,j,i) = 0.0;
        }
        if (NON_BAROTROPIC_EOS) {
          Real p_over_r = PoverR(rad,phi,z);
          phydro->u(IEN,k,j,i) = p_over_r*phydro->u(IDN,k,j,i)/(gamma_gas - 1.0);
          phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i))+SQR(phydro->u(IM2,k,j,i))
                                       + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  return;
}
//Center-of-Mass for planet-star system
void Planet(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
            const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bbc,
            AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar) {
  for (int planetnumber = 1; planetnumber <= 2; ++planetnumber) {
    for (int k = pmb->ks; k <= pmb->ke; ++k) {
      z = pmb->pcoord->x3v(k);
      for (int j = pmb->js; j <= pmb->je; ++j) {
        phi = pmb->pcoord->x2v(j);
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          r = pmb->pcoord->x1v(i);
          Real period;
          Real phip;
          Real com;
          Real rp_value;

          if (planetnumber == 1) {
            period = 2 * M_PI * sqrt(pow(rp, 3) / gm0);
            phip = 2 * (M_PI / period) * time;
            //com = (gm_planet * rp)/ (1 + gm_planet);
            rp_value = rp; //center of mass is calculated since the two planet system orbits around the barycenter rather than the central star
          } else {
            period = 2 * M_PI * sqrt(pow(rp2, 3) / gm0);
            phip = 2 * (M_PI / period) * time;
            //com = ((gm_planet * rp) + (gm_planet * rp2)) / (1 + (2.0 * gm_planet));
            rp_value = rp2;
          }

          Real d = sqrt(pow(rp_value, 2) + pow(r, 2) - 2 * rp_value * r * cos(phi - phip));
          Real dens = prim(IDN, k, j, i);
          Real velocity_x = prim(IVX, k, j, i);
          Real velocity_y = prim(IVY, k, j, i);
          epsilon = 0.3;
          R_H = cbrt(gm_planet / 3);
          Real F_g = -(dens) * ((gm_planet * d) / (sqrt(pow(pow(d, 2) + pow(epsilon, 2) * pow(R_H, 2), 3))));
          Real cosine_term = (pow(r, 2) * (pow(cos(phi), 2)) - r * rp_value * cos(phi) * cos(phip) + pow(r, 2) * (pow(sin(phi), 2)) - r * rp_value * sin(phi) * sin(phip)) / (r * d);
          Real sine_term = (r * rp_value * cos(phi) * sin(phip) - r * rp_value * sin(phi) * cos(phip)) / (r * d);
          Real Fg_x = F_g * cosine_term;
          Real Fg_y = -F_g * sine_term;
          Real delta_momentum_x = Fg_x * dt;
          Real delta_momentum_y = Fg_y * dt;
          cons(IM1, k, j, i) += delta_momentum_x;
          cons(IM2, k, j, i) += delta_momentum_y;
          if (NON_BAROTROPIC_EOS) cons(IEN, k, j, i) += (Fg_x * velocity_x + Fg_y * velocity_y) * dt;
          Real gamma = (rho0 * p0_over_r0) / (pow(r0, dslope));
          Real beta = rho0 / (pow(r0, dslope));
          Real pressure_0 = gamma * pow(r, pslope + dslope);
          Real surface_density_0 = beta * pow(r, dslope);
          Real pressure = dens * (pressure_0 / surface_density_0);
          if (NON_BAROTROPIC_EOS) cons(IEN, k, j, i) += 3.0 / 2.0 * (pressure - prim(IPR, k, j, i));
        }
      }
    }
  }
}

void Viscosity(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, 
               int is, int ie, int js, int je, int ks, int ke) {
    if (phdif->nu_iso > 0.0) {
      for (int k = ks; k <= ke; ++k) {
        z = pmb->pcoord->x3v(k);
        for (int j = js; j <= je; ++j) {
          phi = pmb->pcoord->x2v(j);
          for (int i = is; i <= ie; ++i) {
            r = pmb->pcoord->x1v(i);
            Real omega = sqrt(gm0/(pow(r,3)));
            Real sound_speed = scale * omega*r;
            Real kinematic_viscosity = alpha * sound_speed * (sound_speed/omega); 
            phdif->nu(HydroDiffusion::DiffProcess::iso,k,j,i) = kinematic_viscosity;
        }
      }
    }
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  Real time1 = pmy_mesh -> time;
  for (int k = ks; k <= ke; ++k) {
    z = pcoord->x3v(k);
    for (int j = js; j <= je; ++j) {
      phi = pcoord->x2v(j);
      for (int i = is; i <= ie; ++i) {
        r = pcoord->x1v(i);
        Real period = 2*M_PI*sqrt(pow(rp,3)/gm0);
        phip = 2*(M_PI / period)*time1;
        d = sqrt(pow(rp,2) + pow(r,2) - 2*rp*r*cos(phi - phip));
        epsilon = 0.3;
        R_H = cbrt(gm_planet/3);
        Real g_mag = -1*((gm_planet*d) / (sqrt(pow(pow(d,2) + pow(epsilon,2)*pow(R_H,2), 3))));
        cosine_term = (pow(r,2)*(pow(cos(phi),2)) - r*rp*cos(phi)*cos(phip) + pow(r,2)*(pow(sin(phi),2)) - r*rp*sin(phi)*sin(phip)) / (r*d);
        sine_term = (r*rp*cos(phi)*sin(phip) - r*rp*sin(phi)*cos(phip)) / (r*d);
        user_out_var(0,k,j,i) = g_mag*cosine_term;
        user_out_var(1,k,j,i) = -g_mag*sine_term;
      }
    }
  }
}

Real Torque(MeshBlock *pmb, int iout) { //This torque is only calculated for first planet
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  Real sum_torque = 0;
  for(int k=ks; k<=ke; k++) {
    z = pmb->pcoord->x3v(k);
    for(int j=js; j<=je; j++) {
      phi = pmb->pcoord->x2v(j);
      for(int i=is; i<=ie; i++) {
        r = pmb->pcoord->x1v(i);
        Real d = sqrt(pow(rp,2) + pow(r,2) - 2*rp*r*cos(phi - phip));
        Real g_mag = -1*((gm_planet*d) / (sqrt(pow(pow(d,2) + pow(epsilon,2)*pow(R_H,2), 3))));
        Real dens = pmb->phydro->u(IDN,k,j,i);
        Real volume = pmb ->pcoord->GetCellVolume(k,j,i);
        Real sine_term = (r*rp*cos(phi)*sin(phip) - r*rp*sin(phi)*cos(phip)) / (r*d);
        sum_torque +=  dens * volume *r * -g_mag * sine_term;
      }
    }
  }
  return sum_torque;
}

Real Torque2 (MeshBlock *pmb, int iout) { 
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  Real sum_torque2 = 0;
  for(int k=ks; k<=ke; k++) {
    z = pmb->pcoord->x3v(k); 
    for(int j=js; j<=je; j++) {
      phi = pmb->pcoord->x2v(j);
      for(int i=is; i<=ie; i++) {
        r = pmb->pcoord->x1v(i);
        Real d = sqrt(pow(rp2,2) + pow(r,2) - 2*rp2*r*cos(phi - phip));
        Real g_mag = -1*((gm_planet*d) / (sqrt(pow(pow(d,2) + pow(epsilon,2)*pow(R_H,2), 3))));
        Real dens = pmb->phydro->u(IDN,k,j,i);
        Real volume = pmb ->pcoord->GetCellVolume(k,j,i);
        Real sine_term = (r*rp2*cos(phi)*sin(phip) - r*rp2*sin(phi)*cos(phip)) / (r*d);
        sum_torque2 +=  dens * volume *r * g_mag * sine_term;       
      }
    }
  }
  return sum_torque2;
}

namespace {
//----------------------------------------------------------------------------------------
//! transform to cylindrical coordinate

void GetCylCoord(Coordinates *pco,Real &rad,Real &phi,Real &z,int i,int j,int k) {
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    rad=pco->x1v(i);
    phi=pco->x2v(j);
    z=pco->x3v(k);
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    rad=std::abs(pco->x1v(i)*std::sin(pco->x2v(j)));
    phi=pco->x3v(k);
    z=pco->x1v(i)*std::cos(pco->x2v(j));
  }
  return;
}

//----------------------------------------------------------------------------------------
//! computes density in cylindrical coordinates

Real DenProfileCyl(const Real rad, const Real phi, const Real z) {
  Real den;
  Real p_over_r = p0_over_r0;
  if (NON_BAROTROPIC_EOS) p_over_r = PoverR(rad, phi, z);
  Real denmid = rho0*std::pow(rad/r0,dslope);
  Real dentem = denmid*std::exp(gm0/p_over_r*(1./std::sqrt(SQR(rad)+SQR(z))-1./rad));
  den = dentem;
  return std::max(den,dfloor);
}

//----------------------------------------------------------------------------------------
//! computes pressure/density in cylindrical coordinates

Real PoverR(const Real rad, const Real phi, const Real z) {
  Real poverr;
  poverr = p0_over_r0*std::pow(rad/r0, pslope);
  return poverr;
}

//----------------------------------------------------------------------------------------
//! computes rotational velocity in cylindrical coordinates

Real VelProfileCyl(const Real rad, const Real phi, const Real z) {
  Real p_over_r = PoverR(rad, phi, z);
  Real vel = (dslope+pslope)*p_over_r/(gm0/rad) + (1.0+pslope)
             - pslope*rad/std::sqrt(rad*rad+z*z);
  vel = std::sqrt(gm0/rad)*std::sqrt(vel) - rad*Omega0;
  return vel;
}
} // namespace

//----------------------------------------------------------------------------------------
//! User-defined Boundary and Initial Conditions
void Steady_State_Inner(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    z = pmb->pcoord->x3v(k);
    for (int j=jl; j<=ju; ++j) {
      phi = pmb->pcoord->x2v(j);
      for (int i=1; i<=ngh; ++i) {
        r = pmb->pcoord->x1v(il-i);
        Real gamma = (rho0*p0_over_r0) / (pow(r0, dslope));
        Real beta = rho0/(pow(r0, dslope));
        Real pressure_0 = gamma * pow(r, pslope+dslope);
        Real surface_density_0 = beta * pow(r, dslope);
        Real surface_density = rho0 / sqrt(r/rp);
        Real pressure = surface_density * (pressure_0/surface_density_0);
        Real v_r = -3.0/2.0 * alpha * pow(scale,2) * sqrt((gm0+gm_planet)/r);
        Real v_phi = r * sqrt(1-0.5*pow(scale,2)) * sqrt(gm0+gm_planet)* sqrt(1 / pow(r,3));
        prim(IDN,k,j,il-i) = surface_density;
        prim(IPR,k,j,il-i) = pressure;
        prim(IVX,k,j,il-i) = v_r;
        prim(IVY,k,j,il-i) = v_phi;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined Boundary and Initial Conditions
void Steady_State_Outer(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &prim, FaceField &b,
                    Real time, Real dt,
                    int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    z = pmb->pcoord->x3v(k);
    for (int j=jl; j<=ju; ++j) {
      phi = pmb->pcoord->x2v(j);
      for (int i=1; i<=ngh; ++i) {
        r = pmb->pcoord->x1v(iu+1);
        Real gamma = (rho0*p0_over_r0) / (pow(r0, dslope));
        Real beta = rho0/(pow(r0, dslope));
        Real pressure_0 = gamma * pow(r, pslope+dslope);
        Real surface_density_0 = beta * pow(r, dslope);
        Real surface_density = rho0 / sqrt(r/rp);
        Real pressure = surface_density * (pressure_0/surface_density_0);
        Real v_r = -3.0/2.0 * alpha * pow(scale,2) * sqrt((gm0+gm_planet)/r);
        Real v_phi = r * sqrt(1-0.5*pow(scale,2)) * sqrt(gm0+gm_planet)* sqrt(1 / pow(r,3));
        prim(IDN,k,j,iu+1) = surface_density;
        prim(IPR,k,j,iu+1) = pressure;
        prim(IVX,k,j,iu+1) = v_r;
        prim(IVY,k,j,iu+1) = v_phi;
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

/*void DiskInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,il-i,j,k);
          prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(il-i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,il-i) = 0.0;
          prim(IM2,k,j,il-i) = vel;
          prim(IM3,k,j,il-i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,il-i,j,k);
          prim(IDN,k,j,il-i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(il-i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,il-i) = 0.0;
          prim(IM2,k,j,il-i) = 0.0;
          prim(IM3,k,j,il-i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,il-i) = PoverR(rad, phi, z)*prim(IDN,k,j,il-i);
        }
      }
    }
  }
}*/

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

/*void DiskOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
          prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(iu+i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,iu+i) = 0.0;
          prim(IM2,k,j,iu+i) = vel;
          prim(IM3,k,j,iu+i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=1; i<=ngh; ++i) {
          GetCylCoord(pco,rad,phi,z,iu+i,j,k);
          prim(IDN,k,j,iu+i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(iu+i), pco->x2v(j), pco->x3v(k));
          prim(IM1,k,j,iu+i) = 0.0;
          prim(IM2,k,j,iu+i) = 0.0;
          prim(IM3,k,j,iu+i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,j,iu+i) = PoverR(rad, phi, z)*prim(IDN,k,j,iu+i);
        }
      }
    }
  }
}*/

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskInnerX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,jl-j,k);
          prim(IDN,k,jl-j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(jl-j), pco->x3v(k));
          prim(IM1,k,jl-j,i) = 0.0;
          prim(IM2,k,jl-j,i) = vel;
          prim(IM3,k,jl-j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,jl-j,i) = PoverR(rad, phi, z)*prim(IDN,k,jl-j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,jl-j,k);
          prim(IDN,k,jl-j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(jl-j), pco->x3v(k));
          prim(IM1,k,jl-j,i) = 0.0;
          prim(IM2,k,jl-j,i) = 0.0;
          prim(IM3,k,jl-j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,jl-j,i) = PoverR(rad, phi, z)*prim(IDN,k,jl-j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskOuterX2(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,ju+j,k);
          prim(IDN,k,ju+j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(ju+j), pco->x3v(k));
          prim(IM1,k,ju+j,i) = 0.0;
          prim(IM2,k,ju+j,i) = vel;
          prim(IM3,k,ju+j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,ju+j,i) = PoverR(rad, phi, z)*prim(IDN,k,ju+j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=1; j<=ngh; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,ju+j,k);
          prim(IDN,k,ju+j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(ju+j), pco->x3v(k));
          prim(IM1,k,ju+j,i) = 0.0;
          prim(IM2,k,ju+j,i) = 0.0;
          prim(IM3,k,ju+j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,k,ju+j,i) = PoverR(rad, phi, z)*prim(IDN,k,ju+j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskInnerX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,kl-k);
          prim(IDN,kl-k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(kl-k));
          prim(IM1,kl-k,j,i) = 0.0;
          prim(IM2,kl-k,j,i) = vel;
          prim(IM3,kl-k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,kl-k,j,i) = PoverR(rad, phi, z)*prim(IDN,kl-k,j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,kl-k);
          prim(IDN,kl-k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(kl-k));
          prim(IM1,kl-k,j,i) = 0.0;
          prim(IM2,kl-k,j,i) = 0.0;
          prim(IM3,kl-k,j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,kl-k,j,i) = PoverR(rad, phi, z)*prim(IDN,kl-k,j,i);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------
//! User-defined boundary Conditions: sets solution in ghost zones to initial values

void DiskOuterX3(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, FaceField &b,
                 Real time, Real dt,
                 int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  Real rad(0.0), phi(0.0), z(0.0);
  Real vel;
  OrbitalVelocityFunc &vK = pmb->porb->OrbitalVelocity;
  if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,ku+k);
          prim(IDN,ku+k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(ku+k));
          prim(IM1,ku+k,j,i) = 0.0;
          prim(IM2,ku+k,j,i) = vel;
          prim(IM3,ku+k,j,i) = 0.0;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    for (int k=1; k<=ngh; ++k) {
      for (int j=jl; j<=ju; ++j) {
        for (int i=il; i<=iu; ++i) {
          GetCylCoord(pco,rad,phi,z,i,j,ku+k);
          prim(IDN,ku+k,j,i) = DenProfileCyl(rad,phi,z);
          vel = VelProfileCyl(rad,phi,z);
          if (pmb->porb->orbital_advection_defined)
            vel -= vK(pmb->porb, pco->x1v(i), pco->x2v(j), pco->x3v(ku+k));
          prim(IM1,ku+k,j,i) = 0.0;
          prim(IM2,ku+k,j,i) = 0.0;
          prim(IM3,ku+k,j,i) = vel;
          if (NON_BAROTROPIC_EOS)
            prim(IEN,ku+k,j,i) = PoverR(rad, phi, z)*prim(IDN,ku+k,j,i);
        }
      }
    }
  }
}