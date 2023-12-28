#ifndef SCALARS_SCALARS_HPP_
#define SCALARS_SCALARS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file scalars.hpp
//! \brief definitions for PassiveScalars class

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/cc/bvals_cc.hpp"
//chemistry headers
#include "../chemistry/network/network.hpp"
#include "../chemistry/ode_wrapper.hpp"
#include  CHEMNETWORK_HEADER //ChemNetwork class

class MeshBlock;
class ParameterInput;

//! \class PassiveScalars
//! \brief
//!
//! \todo (felker): consider renaming to Scalars
class PassiveScalars {
 public:
  // TODO(felker): pin is currently only used for checking ssprk5_4, otherwise unused.
  // Leaving as ctor parameter in case of run-time "nscalars" option
  PassiveScalars(MeshBlock *pmb, ParameterInput *pin);
  ~PassiveScalars();

  // public data:
  // "conserved vars" = passive scalar mass
  AthenaArray<Real> s, s1, s2;  // (no more than MAX_NREGISTER allowed)
  AthenaArray<Real> s0, s_fl_div;  // rkl2 STS memory registers;
  // "primitive vars" = (density-normalized) mass fraction/concentration of each species
  AthenaArray<Real> r;  // , r1;
  AthenaArray<Real> s_flux[3];  // face-averaged flux vector

  // fourth-order intermediate quantities
  AthenaArray<Real> mass_flux_fc[3];  // deep copy of Hydro intermediate flux quantities
  AthenaArray<Real> s_cc, r_cc;       // cell-centered approximations
  // (only needed for 4th order EOS evaluations that have explicit dependence on species
  // concentration)

  // storage for SMR/AMR
  // TODO(KGF): remove trailing underscore or revert to private:
  AthenaArray<Real> coarse_s_, coarse_r_;
  int refinement_idx{-1};

  CellCenteredBoundaryVariable sbvar;

  //chemistry variables
  //chemistry source term
  //s(ispec, k, j, i). read in s1(i, ispec), and loop over i,
  //maybe parallelize i with openmpi later.
  AthenaArray<Real> r_copy; //abundance of species copy at intermediate step
  AthenaArray<Real> h; //next stepsize in chemistry solver
  ChemNetwork chemnet; //pointer to chemical network
  ODEWrapper odew; //pointer to ode solver

  // public functions:
  // KGF: use inheritance for these functions / overall class?
  void AddFluxDivergence(const Real wght, AthenaArray<Real> &s_out);
  void AddFluxDivergence_STS(const Real wght, int stage,
                             AthenaArray<Real> &s_out, AthenaArray<Real> &s_flux_div_out);
  void CalculateFluxes(AthenaArray<Real> &s, const int order);
  void CalculateFluxes_STS();

  //! \note
  //! for now, not creating subfolder "scalars_diffusion/", nor class ScalarDiffusion
  //! that is would have an instance contained within PassiveScalars like HydroDiffusion
  //! approach. Consider creating an encapsulated class as these features are generalized.
  Real nu_scalar_iso; //, nu_scalar_aniso;          // diffusion coeff
  bool scalar_diffusion_defined;
  AthenaArray<Real> diffusion_flx[3];
  // AthenaArray<Real> nu_scalar;               // diffusion array

  //! \note
  //! No need for nu_scalar array, nor counterpart to HydroDiffusion::CalcDiffusionFlux
  //! wrapper function since, currently:
  //! - nu_scalar_iso must be constant across the mesh
  //!   (does not depend on local fluid or field variables),
  //! - there is only one type of
  //!   passive scalar diffusion process (nu_scalar_aniso disabled, no "eta"l, etc.)
  //! - nu_scalar_iso is identical for all NSCALARS
  void DiffusiveFluxIso(const AthenaArray<Real> &prim_r, const AthenaArray<Real> &w,
                        AthenaArray<Real> *flx_out);
  Real NewDiffusionDt();

 private:
  MeshBlock* pmy_block;
  // scratch space used to compute fluxes
  // 2D scratch arrays
  AthenaArray<Real> rl_, rr_, rlb_;
  // 1D scratch arrays
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_;
  AthenaArray<Real> dflx_;

  // fourth-order
  // 4D scratch arrays
  AthenaArray<Real> scr1_nkji_, scr2_nkji_;
  AthenaArray<Real> rl3d_, rr3d_;
  // 1D scratch arrays
  AthenaArray<Real> laplacian_l_fc_, laplacian_r_fc_;

  void ComputeUpwindFlux(const int k, const int j, const int il,
                         const int iu, // CoordinateDirection dir,
                         AthenaArray<Real> &rl, AthenaArray<Real> &rr,
                         AthenaArray<Real> &mass_flx,
                         AthenaArray<Real> &flx_out);
  void AddDiffusionFluxes();
  // TODO(felker): dedpulicate these arrays and the same named ones in HydroDiffusion
  AthenaArray<Real> dx1_, dx2_, dx3_;
};
#endif // SCALARS_SCALARS_HPP_
