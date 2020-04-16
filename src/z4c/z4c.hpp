#ifndef Z4c_HPP
#define Z4c_HPP
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file z4c.hpp
//  \brief definitions for the Z4c class
//
// Convention: tensor names are followed by tensor type suffixes:
//    _u --> contravariant component
//    _d --> covariant component
// For example g_dd is a tensor, or tensor-like object, with two covariant indices.

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../athena_tensor.hpp"
#include "../finite_differencing.hpp"

#include "../bvals/cc/bvals_cc.hpp"
#include "../bvals/vc/bvals_vc.hpp"

class MeshBlock;
class ParameterInput;

// Indexes for variables in AthenaArray
#define NDIM (3) // Manifold dimension

//! \class Z4c
//  \brief Z4c data and functions
class Z4c {

public:
  // Indexes of evolved variables
  enum {
    I_Z4c_chi,
    I_Z4c_gxx, I_Z4c_gxy, I_Z4c_gxz, I_Z4c_gyy, I_Z4c_gyz, I_Z4c_gzz,
    I_Z4c_Khat,
    I_Z4c_Axx, I_Z4c_Axy, I_Z4c_Axz, I_Z4c_Ayy, I_Z4c_Ayz, I_Z4c_Azz,
    I_Z4c_Gamx, I_Z4c_Gamy, I_Z4c_Gamz,
    I_Z4c_Theta,
    I_Z4c_alpha,
    I_Z4c_betax, I_Z4c_betay, I_Z4c_betaz,
    N_Z4c
  };
  // Names of Z4c variables
  static char const * const Z4c_names[N_Z4c];
  // Indexes of ADM variables
  enum {
    I_ADM_gxx, I_ADM_gxy, I_ADM_gxz, I_ADM_gyy, I_ADM_gyz, I_ADM_gzz,
    I_ADM_Kxx, I_ADM_Kxy, I_ADM_Kxz, I_ADM_Kyy, I_ADM_Kyz, I_ADM_Kzz,
    I_ADM_psi4,
    N_ADM
  };
  // Names of ADM variables
  static char const * const ADM_names[N_ADM];
  // Indexes of Constraint variables
  enum {
    I_CON_C,
    I_CON_H,
    I_CON_M,
    I_CON_Z,
    I_CON_Mx, I_CON_My, I_CON_Mz,
    N_CON,
  };
  // Names of costraint variables
  static char const * const Constraint_names[N_CON];
  // Indexes of matter fields
  enum {
    I_MAT_rho,
    I_MAT_Sx, I_MAT_Sy, I_MAT_Sz,
    I_MAT_Sxx, I_MAT_Sxy, I_MAT_Sxz, I_MAT_Syy, I_MAT_Syz, I_MAT_S_zz,
    N_MAT
  };
  // Names of matter variables
  static char const * const Matter_names[N_MAT];

public:
  Z4c(MeshBlock *pmb, ParameterInput *pin);
  ~Z4c();

  MeshBlock * pmy_block;     // pointer to MeshBlock containing this Z4c

  // public data storage
  struct {
    AthenaArray<Real> u;     // solution of Z4c evolution system
    AthenaArray<Real> u1;    // solution at intermediate steps
    AthenaArray<Real> u2;    // solution at intermediate steps
    AthenaArray<Real> rhs;   // Z4c rhs
    //DEBUG Needed for Traditional RK4
    AthenaArray<Real> rhs1;    // intermediate storage
    //ENDDEBUG
    AthenaArray<Real> adm;   // ADM variables
    AthenaArray<Real> con;   // constraints
    AthenaArray<Real> mat;   // matter variables
  } storage;

  // aliases for variables and RHS
  struct Z4c_vars {
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> chi;       // conf. factor
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> Khat;      // trace extr. curvature
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> Theta;     // Theta var in Z4c
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> alpha;     // lapse
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> Gam_u;     // Gamma functions (BSSN)
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> beta_u;    // shift
    AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> g_dd;      // conf. 3-metric
    AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> A_dd;      // conf. traceless extr. curvature
  };
  Z4c_vars z4c;
  Z4c_vars rhs;

  // aliases for the ADM variables
  struct ADM_vars {
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> psi4;      // conformal factor
    AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> g_dd;      // 3-metric
    AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> K_dd;      // curvature
  };
  ADM_vars adm;

  // aliases for the constraints
  struct Constraint_vars {
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> C;         // Z constraint monitor
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> H;         // hamiltonian constraint
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> M;         // norm squared of M_d
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> Z;         // Z constraint violation
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> M_d;       // momentum constraint
  };
  Constraint_vars con;

  // aliases for the matter variables
  struct Matter_vars {
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> rho;       // matter energy density
    AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> S_d;       // matter momentum density
    AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> S_dd;      // matter stress tensor
  };
  Matter_vars mat;

  // user settings and options
  struct {
    Real chi_psi_power;   // chi = psi^N, N = chi_psi_power
    Real chi_div_floor;   // puncture's floor value for chi, use max(chi, chi_div_floor) in non-differentiated chi
    Real diss;            // amount of numerical dissipation
    Real eps_floor;       // a small number O(10^-12)
    // Constraint damping parameters
    Real damp_kappa1;
    Real damp_kappa2;
    // Gauge conditions for the lapse
    Real lapse_oplog;
    Real lapse_harmonicf;
    Real lapse_harmonic;
    Real lapse_advect;
    // Gauge condition for the shift
    Real shift_advect;
    Real shift_eta;
    // Single puncture parameters
    Real punc_ADM_mass;
    // AwA parameters
    Real AwA_amplitude; // amplitude parameter
    Real AwA_d_x; // d_x (width) parameter
    Real AwA_d_y; // d_y (width) parameter
    int AwA_rho; // Resolution index
    int AwA_direction; // direction of the test
  } opt;


  // boundary and grid data
#if PREFER_VC
  VertexCenteredBoundaryVariable ubvar;
#else
  CellCenteredBoundaryVariable ubvar;
#endif
  AthenaArray<Real> empty_flux[3];

  // storage for SMR/AMR
  // BD: this should perhaps be combined with the above stuct.
  AthenaArray<Real> coarse_u_;
  int refinement_idx{-1};

  // for seamless CC/VC switching
  struct MB_info {
    int il, iu, jl, ju, kl, ku;        // local block iter.
    int nn1, nn2, nn3;                 // number of nodes (simplify switching)

    AthenaArray<Real> x1, x2, x3;      // for CC / VC grid switch
    AthenaArray<Real> cx1, cx2, cx3;   // for CC / VC grid switch (coarse)
  };

  MB_info mbi;

public:
  // scheduled functions
  //
  // compute new timestep on a MeshBlock
  Real NewBlockTimeStep(void);
  // compute the RHS given the Z4c and matter variables
  void Z4cRHS(AthenaArray<Real> & u, AthenaArray<Real> & mat, AthenaArray<Real> & rhs);
  // compute the boundary RHS given the Z4c and matter variables
  void Z4cBoundaryRHS(AthenaArray<Real> & u, AthenaArray<Real> & mat, AthenaArray<Real> & rhs);
  // compute linear combination of states
  void WeightedAve(AthenaArray<Real> &u_out, AthenaArray<Real> &u_in1,
                   AthenaArray<Real> &u_in2, const Real wght[3]);
  // add RHS to state
  void AddZ4cRHS(AthenaArray<Real> & rhs, Real const wght, AthenaArray<Real> &u_out);
  // compute Z4c variables from ADM variables
  void ADMToZ4c(AthenaArray<Real> & u_adm, AthenaArray<Real> & u);
  // compute ADM variables from Z4c variables
  void Z4cToADM(AthenaArray<Real> & u, AthenaArray<Real> & u_adm);
  // enforce algebraic constraints on the solution
  void AlgConstr(AthenaArray<Real> & u);
  // compute ADM constraints
  void ADMConstraints(AthenaArray<Real> & u_con, AthenaArray<Real> & u_adm,
                      AthenaArray<Real> & u_mat, AthenaArray<Real> & u_z4c);

  // utility functions
  //
  // set ADM aliases given u_adm
  void SetADMAliases(AthenaArray<Real> & u_adm, ADM_vars & adm);
  // set constraint aliases for a given u_con
  void SetConstraintAliases(AthenaArray<Real> & u_con, Constraint_vars & con);
  // set matter aliases given a state
  void SetMatterAliases(AthenaArray<Real> & u_mat, Matter_vars & mat);
  // set Z4c aliases given a state
  void SetZ4cAliases(AthenaArray<Real> & u, Z4c_vars & z4c);

  // compute spatial determinant of a 3x3  matrix
  Real SpatialDet(Real const gxx, Real const gxy, Real const gxz,
      Real const gyy, Real const gyz, Real const gzz);
  Real SpatialDet(AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> const & g,
                  int const k, int const j, int const i) {
    return SpatialDet(g(0,0,k,j,i), g(0,1,k,j,i), g(0,2,k,j,i),
                      g(1,1,k,j,i), g(1,2,k,j,i), g(2,2,k,j,i));
  }
  // compute inverse of a 3x3 matrix
  void SpatialInv(Real const detginv,
                  Real const gxx, Real const gxy, Real const gxz,
                  Real const gyy, Real const gyz, Real const gzz,
                  Real * uxx, Real * uxy, Real * uxz,
                  Real * uyy, Real * uyz, Real * uzz);
  // compute trace of a rank 2 covariant spatial tensor
  Real Trace(Real const detginv,
             Real const gxx, Real const gxy, Real const gxz,
             Real const gyy, Real const gyz, Real const gzz,
             Real const Axx, Real const Axy, Real const Axz,
             Real const Ayy, Real const Ayz, Real const Azz);

  // additional global functions

  // setup a Minkowski spacetime
  void ADMMinkowski(AthenaArray<Real> & u_adm);
  // set the gauge condition to geodesic slicing
  void GaugeGeodesic(AthenaArray<Real> & u);
  // set the matter variables to zero
  void MatterVacuum(AthenaArray<Real> & u_adm);

  // initial data for the AwA tests
  void ADMRobustStability(AthenaArray<Real> & u_adm);
  void GaugeRobStab(AthenaArray<Real> & u);
  void ADMLinearWave1(AthenaArray<Real> & u_adm);
  void ADMLinearWave2(AthenaArray<Real> & u_adm);
  void ADMGaugeWave1(AthenaArray<Real> & u_adm);
  void ADMGaugeWave1_shifted(AthenaArray<Real> & u_adm);
  void ADMGaugeWave2(AthenaArray<Real> & u_adm);
  void GaugeGaugeWave1(AthenaArray<Real> & u);
  void GaugeGaugeWave1_shifted(AthenaArray<Real> & u);
  void GaugeGaugeWave2(AthenaArray<Real> & u);
  void GaugeSimpleGaugeWave(AthenaArray<Real> & u);

  // initial data for a single BH
  void ADMOnePuncture(AthenaArray<Real> & u_adm);
  void GaugePreCollapsedLapse(AthenaArray<Real> & u);
  void ADMOnePunctureSpin(AthenaArray<Real> & u_adm);

  // initial data for binary BHs
  void ADMTwoPunctures(AthenaArray<Real> & u_adm);

private:
  AthenaArray<Real> dt1_,dt2_,dt3_;  // scratch arrays used in NewTimeStep

  // auxiliary tensors
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> r;           // radial coordinate
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> detg;        // det(g)
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> chi_guarded; // bounded version of chi
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> oopsi4;      // 1/psi4
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> A;           // trace of A
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> AA;          // trace of AA
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> R;           // Ricci scalar
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> Ht;          // tilde H
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> K;           // trace of extrinsic curvature
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> KK;          // K^a_b K^b_a
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> Ddalpha;     // Trace of Ddalpha_dd
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> S;           // Trace of S_ik
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> M_u;         // momentum constraint
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> Gamma_u;     // Gamma computed from the metric
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> DA_u;        // Covariant derivative of A
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> s_u;         // x^i/r where r is the coord. radius
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> g_uu;        // inverse of conf. metric
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> A_uu;        // inverse of A
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> AA_dd;       // g^cd A_ac A_db
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> R_dd;        // Ricci tensor
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> Rphi_dd;     // Ricci tensor, conformal contribution
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> Kt_dd;       // conformal extrinsic curvature
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> K_ud;        // extrinsic curvature
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> Ddalpha_dd;  // 2nd differential of the lapse
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> Ddphi_dd;    // 2nd differential of phi
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 3> Gamma_ddd;   // Christoffel symbols of 1st kind
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 3> Gamma_udd;   // Christoffel symbols of 2nd kind
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 3> DK_ddd;      // differential of K
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 3> DK_udd;      // differential of K

  // auxiliary derivatives
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> dbeta;       // d_a beta^a
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> dalpha_d;    // lapse 1st drvts
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> ddbeta_d;    // 2nd "divergence" of beta
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> dchi_d;      // chi 1st drvts
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> dphi_d;      // phi 1st drvts
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> dK_d;        // K 1st drvts
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> dKhat_d;     // Khat 1st drvts
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> dTheta_d;    // Theta 1st drvts
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> ddalpha_dd;  // lapse 2nd drvts
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 2> dbeta_du;    // shift 1st drvts
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> ddchi_dd;    // chi 2nd drvts
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 2> dGam_du;     // Gamma 1st drvts
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 3> dg_ddd;      // metric 1st drvts
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 3> dg_duu;      // inverse metric 1st drvts
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 3> dK_ddd;      // K 1st drvts
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 3> dA_ddd;      // A 1st drvts
  AthenaTensor<Real, TensorSymm::ISYM2, NDIM, 3> ddbeta_ddu; // shift 2nd drvts
  AthenaTensor<Real, TensorSymm::SYM22, NDIM, 4> ddg_dddd;   // metric 2nd drvts

  // auxialiry Lie derivatives along the shift vector
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> Lchi;        // Lie derivative of chi
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> LKhat;       // Lie derivative of Khat
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> LTheta;      // Lie derivative of Theta
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 0> Lalpha;      // Lie derivative of the lapse
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> LGam_u;      // Lie derivative of Gamma
  AthenaTensor<Real, TensorSymm::NONE, NDIM, 1> Lbeta_u;     // Lie derivative of the shift
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> Lg_dd;       // Lie derivative of conf. 3-metric
  AthenaTensor<Real, TensorSymm::SYM2, NDIM, 2> LA_dd;       // Lie derivative of A

private:
  void Z4cSommerfeld_(AthenaArray<Real> & u, AthenaArray<Real> & rhs,
      int const is, int const ie, int const js, int const je, int const ks, int const ke);

private:
  struct {
    // 1st derivative stecil
    typedef FDCenteredStencil<1, NGHOST-1> s1;
    // 2nd derivative stencil
    typedef FDCenteredStencil<2, NGHOST-1> s2;
    // dissipation operator
    typedef FDCenteredStencil<
      FDDissChoice<NGHOST-1>::degree,
      FDDissChoice<NGHOST-1>::nghost
      > sd;
    // left-biased derivative
    typedef FDLeftBiasedStencil<
        FDBiasedChoice<1, NGHOST-1>::degree,
        FDBiasedChoice<1, NGHOST-1>::nghost,
        FDBiasedChoice<1, NGHOST-1>::lopsize
      > sl;
    // right-biased derivative
    typedef FDRightBiasedStencil<
        FDBiasedChoice<1, NGHOST-1>::degree,
        FDBiasedChoice<1, NGHOST-1>::nghost,
        FDBiasedChoice<1, NGHOST-1>::lopsize
      > sr;

    int stride[3];
    Real idx[3];
    Real diss;

    // 1st derivative (high order centered)
    inline Real Dx(int dir, Real & u) {
      Real * pu = &u - s1::offset*stride[dir];

      Real out(0.);
      for(int n1 = 0; n1 < s1::nghost; ++n1) {
        int const n2  = s1::width - n1 - 1;
        Real const c1 = s1::coeff[n1] * pu[n1*stride[dir]];
        Real const c2 = s1::coeff[n2] * pu[n2*stride[dir]];
        out += (c1 + c2);
      }
      out += s1::coeff[s1::nghost] * pu[s1::nghost*stride[dir]];

      return out * idx[dir];
    }
    // 1st derivative 2nd order centered
    inline Real Ds(int dir, Real & u) {
      Real * pu = &u;
      return 0.5 * idx[dir] * (pu[stride[dir]] - pu[-stride[dir]]);
    }
    // Advective derivative
    // The advective derivative is for an equation in the form
    //    d_t u = vx d_x u
    // So negative vx means advection from the *left* to the *right*, so we use
    // *left* biased FD stencils
    inline Real Lx(int dir, Real & vx, Real & u) {
      Real * pu = &u;

      Real dl(0.);
      for(int n = 0; n < sl::width; ++n) {
        dl += sl::coeff[n] * pu[(n - sl::offset)*stride[dir]];
      }

      Real dr(0.);
      for(int n = sr::width-1; n >= 0; --n) {
        dr += sr::coeff[n] * pu[(n - sr::offset)*stride[dir]];
      }

      return ((vx < 0) ? (vx * dl) : (vx * dr)) * idx[dir];
    }
    // Homogeneous 2nd derivative
    inline Real Dxx(int dir, Real & u) {
      Real * pu = &u - s2::offset*stride[dir];

      Real out(0.);
      for(int n1 = 0; n1 < s2::nghost; ++n1) {
        int const n2  = s2::width - n1 - 1;
        Real const c1 = s2::coeff[n1] * pu[n1*stride[dir]];
        Real const c2 = s2::coeff[n2] * pu[n2*stride[dir]];
        out += (c1 + c2);
      }
      out += s2::coeff[s2::nghost] * pu[s2::nghost*stride[dir]];

      return out * SQR(idx[dir]);
    }
    // Mixed 2nd derivative
    inline Real Dxy(int dirx, int diry, Real & u) {
      Real * pu = &u - s1::offset*(stride[dirx] + stride[diry]);
      Real out(0.);

      for(int nx1 = 0; nx1 < s1::nghost; ++nx1) {
        int const nx2 = s1::width - nx1 - 1;
        for(int ny1 = 0; ny1 < s1::nghost; ++ny1) {
          int const ny2 = s1::width - ny1 - 1;

          Real const c11 = s1::coeff[nx1] * s1::coeff[ny1] * pu[nx1*stride[dirx] + ny1*stride[diry]];
          Real const c12 = s1::coeff[nx1] * s1::coeff[ny2] * pu[nx1*stride[dirx] + ny2*stride[diry]];
          Real const c21 = s1::coeff[nx2] * s1::coeff[ny1] * pu[nx2*stride[dirx] + ny1*stride[diry]];
          Real const c22 = s1::coeff[nx2] * s1::coeff[ny2] * pu[nx2*stride[dirx] + ny2*stride[diry]];

          Real const ca = (1./6.)*((c11 + c12) + (c21 + c22));
          Real const cb = (1./6.)*((c11 + c21) + (c12 + c22));
          Real const cc = (1./6.)*((c11 + c22) + (c12 + c21));

          out += ((ca + cb) + cc) + ((ca + cc) + cb);
        }
        int const ny = s1::nghost;

        Real const c1 = s1::coeff[nx1] * s1::coeff[ny] * pu[nx1*stride[dirx] + ny*stride[diry]];
        Real const c2 = s1::coeff[nx2] * s1::coeff[ny] * pu[nx2*stride[dirx] + ny*stride[diry]];

        out += (c1 + c2);
      }
      int const nx = s1::nghost;
      for(int ny1 = 0; ny1 < s1::nghost; ++ny1) {
        int const ny2 = s1::width - ny1 - 1;

        Real const c1 = s1::coeff[nx] * s1::coeff[ny1] * pu[nx*stride[dirx] + ny1*stride[diry]];
        Real const c2 = s1::coeff[nx] * s1::coeff[ny2] * pu[nx*stride[dirx] + ny2*stride[diry]];

        out += (c1 + c2);
      }
      int const ny = s1::nghost;
      out += s1::coeff[nx] * s1::coeff[ny] * pu[nx*stride[dirx] + ny*stride[diry]];

      return out * idx[dirx] * idx[diry];
    }
    // Kreiss-Oliger dissipation operator
    inline Real Diss(int dir, Real & u) {
      Real * pu = &u - sd::offset*stride[dir];

      Real out(0.);
      for(int n1 = 0; n1 < sd::nghost; ++n1) {
        int const n2  = sd::width - n1 - 1;
        Real const c1 = sd::coeff[n1] * pu[n1*stride[dir]];
        Real const c2 = sd::coeff[n2] * pu[n2*stride[dir]];
        out += (c1 + c2);
      }
      out += sd::coeff[sd::nghost] * pu[sd::nghost*stride[dir]];

      return out * idx[dir] * diss;
    }
  } FD;
};

#endif // Z4c_HPP
