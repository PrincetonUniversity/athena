#ifndef ORBITAL_ADVECTION_ORBITAL_ADVECTION_HPP_
#define ORBITAL_ADVECTION_ORBITAL_ADVECTION_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file  orbital_advection.hpp
//! \brief definitions of the OrbitalAdvection class and related functions

// C/C++ headers

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../bvals/orbital/bvals_orbital.hpp"

// Forward declarations
class PassiveScalars;

Real CartOrbitalVelocity(OrbitalAdvection *porb, Real x_, Real y_, Real z_);
Real CartOrbitalVelocity_x(OrbitalAdvection *porb, Real x_, Real y_, Real z_);

Real CylOrbitalVelocity2D(OrbitalAdvection *porb, Real x_, Real y_, Real z_);
Real CylOrbitalVelocity2D_r(OrbitalAdvection *porb, Real x_, Real y_, Real z_);

Real CylOrbitalVelocity3D(OrbitalAdvection *porb, Real x_, Real y_, Real z_);
Real CylOrbitalVelocity3D_r(OrbitalAdvection *porb, Real x_, Real y_, Real z_);
Real CylOrbitalVelocity3D_z(OrbitalAdvection *porb, Real x_, Real y_, Real z_);

Real SphOrbitalVelocity(OrbitalAdvection *porb, Real x_, Real y_, Real z_);
Real SphOrbitalVelocity_r(OrbitalAdvection *porb, Real x_, Real y_, Real z_);
Real SphOrbitalVelocity_t(OrbitalAdvection *porb, Real x_, Real y_, Real z_);

Real ZeroOrbitalVelocity(OrbitalAdvection *porb, Real x_, Real y_, Real z_);

// idntifiers for variable to be transformed in ConvertOrbitalSystem
enum class OrbitalTransform {none=0, prim=1, cons=2, all=3};

//! \class OrbitalAdvection
//! \brief data and functions for orbital advection
class OrbitalAdvection{
  friend class OrbitalBoundaryCommunication;
 public:
  // functions
  OrbitalAdvection(MeshBlock *pmb, ParameterInput *pin);
  ~OrbitalAdvection();

  void InitializeOrbitalAdvection();
  void SetOrbitalAdvectionCC(const AthenaArray<Real> &u, const AthenaArray<Real> &s);
  void SetOrbitalAdvectionFC(const FaceField &b);
  void CalculateOrbitalAdvectionCC(const Real dt, AthenaArray<Real> &u,
                                   AthenaArray<Real> &s);
  void CalculateOrbitalAdvectionFC(const Real dt, EdgeField &e);
  void ConvertOrbitalSystem(const AthenaArray<Real> &w0, const AthenaArray<Real> &u0,
                            const OrbitalTransform trans);
  void ResetOrbitalSystemConversionFlag();
  Real NewOrbitalAdvectionDt();
  void RemapFluxPlm(AthenaArray<Real> &pflux_, const AthenaArray<Real> &pbuf_,
                    const Real eps_, const int osgn_, const int k, const int j,
                    const int il, const int iu, const int shift_ = 0);
  void RemapFluxPlm(AthenaArray<Real> &pflux_, const AthenaArray<Real> &pbuf_,
                    const Real eps_, const int osgn_, const int k, const int j,
                    const int il, const int iu, const int nl, const int nu,
                    const int shift_ = 0);
  void RemapFluxPpm(AthenaArray<Real> &pflux_, AthenaArray<Real> &pbuf_,
                    const Real eps_, const int osgn_, const int k, const int j,
                    const int il, const int iu, const int shift_ = 0);
  void RemapFluxPpm(AthenaArray<Real> &pflux_, AthenaArray<Real> &pbuf_,
                    const Real eps_, const int osgn_, const int k, const int j,
                    const int il, const int iu, const int nl, const int nu,
                    const int shift_ = 0);

  // accessor to orbital velocity
  OrbitalVelocityFunc OrbitalVelocity, OrbitalVelocityDerivative[2];

  // flag
  int orbital_direction; //!> the direction of orbital motion x2(=1), x3 (=2)
  int orbital_splitting_order; //!> order of the orbital splitting method
  bool orbital_advection_defined; //!> flag for the orbital advection system
  bool orbital_advection_active; //!> flag for solving orbital advection
  bool orbital_refinement; //!> flag for refinement in the orbital direction
  bool orbital_uniform_mesh; //!> true: uniform grid, false: un-uniform grid

  AthenaArray<Real> vKc, vKf[2]; // Orbital Velocity (cell-centered and face)
  AthenaArray<Real> dvKc1, dvKc2; // Derivatives of vKc
  Real Omega0, qshear, shboxcoord; // parameters for shearing box in cartesian
  Real gm; // central gravity in cylindrical/spherical_polar

  // TODO(tomo-ono): Consider replace these buffers.
  AthenaArray<Real> w_orb, u_orb; // buffer for orbital advection system output

  OrbitalBoundaryCommunication *orb_bc;

 private:
  // Private Functions
  void SetVKc();
  void SetVKf();
  void SetVKcCoarse();
  void SetVKfCoarse();
  void SetDvKc();
  void SetOrbitalEdgeCC(const Real dt, int *ssize[2], int *rsize[2]);
  void SetOrbitalEdgeFC(const Real dt, int *ssize[2], int *rsize[2]);

  // ptr about this
  MeshBlock *pmb_;        // ptr to this meshblock
  Mesh *pm_;              // ptr to Mesh
  Hydro *ph_;             // ptr to Hydro
  Field *pf_;             // ptr to Field
  Coordinates *pco_;      // ptr to Coordinates
  BoundaryValues *pbval_; // ptr to Boundaryvalues
  PassiveScalars *ps_;    // ptr to PassiveScalars

  // max/min of orbital velocity
  Real vK_max, vK_min;
  // restriction on dt from orbital advection
  Real min_dt;

  // meshblock size
  int nc1, nc2, nc3;
  int onx;
  int xorder, xgh;

  // grids
  AthenaArray<Real> orc, orf[2]; // orbital residual of cell-centered values
  AthenaArray<int>  ofc, off[2]; // orbital offset of cell-centered values

  // uniform mesh
  Real dx;

  // orbital blocks with two meshblock orbital length for variables
  AthenaArray<Real> orbital_cons, orbital_scalar;
  AthenaArray<Real> orbital_b1, orbital_b2;

  // For Orbital Remapping
  // Orbital Velocity in coarse meshblock (cell-centered and face)
  AthenaArray<Real> vKc_coarse, vKf_coarse[2];
  // orbital offset of cell-centered values
  AthenaArray<int>  ofc_coarse, off_coarse[2];

  AthenaArray<Real> u_coarse_send, u_coarse_recv, u_temp;
  AthenaArray<Real> s_coarse_send, s_coarse_recv, s_temp;
  AthenaArray<Real> b1_coarse_send, b2_coarse_send;
  FaceField b_temp, b_coarse_recv;
  int max_ofc_coarse, min_ofc_coarse;
  int max_off_coarse, min_off_coarse;

  // pencil(1D) buffer
  AthenaArray<Real> hbuf;  // pencil buffer for shallow copy for hydro calculation
  AthenaArray<Real> pflux; // pencil buffer for flux

  // buffer for ppm
  AthenaArray<Real> s_src[5], d_src[13]; // s_src for deep copy, d_src for shallow copy

  // flag for orbital system conversion
  int orbital_system_conversion_done; //
};
#endif // ORBITAL_ADVECTION_ORBITAL_ADVECTION_HPP_
