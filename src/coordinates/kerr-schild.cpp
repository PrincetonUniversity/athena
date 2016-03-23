// Kerr spacetime, Kerr-Schild coordinates
// Notes:
//   coordinates: t, r, \theta, \phi
//   parameters:
//     M (mass), M > 0
//     a (spin), -M < a < M
//   metric:
//     ds^2 = -(1 - 2 M r/\Sigma) dt^2
//            + 4 M r/\Sigma dt dr
//            - 4 M a r/\Sigma dt d\phi
//            + (1 + 2 M r/\Sigma) dr^2
//            - 2 a (1 + 2 M r/\Sigma) \sin^2\theta dr d\phi
//            + \Sigma d\theta^2
//            + (r^2 + a^2 + (2 M a^2 r/\Sigma) \sin^2\theta) \sin^2\theta d\phi^2
//     \Delta = r^2 - 2 M r + a^2
//     \Sigma = r^2 + a^2 cos^2\theta
//     \Xi = r^2 - a^2 cos^2\theta
//   other "Kerr-Schild" coordinates exist

// Primary header
#include "coordinates.hpp"

// C++ headers
#include <cmath>  // abs(), cos(), log(), sin(), sqrt()

// Athena++ headers
#include "../athena.hpp"           // enums, macros
#include "../athena_arrays.hpp"    // AthenaArray
#include "../mesh.hpp"             // MeshBlock
#include "../parameter_input.hpp"  // ParameterInput
#include "../hydro/hydro.hpp"      // Hydro
#include "../hydro/eos/eos.hpp"    // HydroEqnOfState

//--------------------------------------------------------------------------------------

// Constructor
// Inputs:
//   pmb: pointer to MeshBlock containing this grid
//   pin: pointer to runtime inputs
//   flag: indicator that this is special refinement case
Coordinates::Coordinates(MeshBlock *pmb, ParameterInput *pin, int flag)
{
  // Set pointer to host MeshBlock and note active zone boundaries
  pmy_block = pmb;
  cflag = flag;
  int is, ie, js, je, ks, ke, ng;
  if(cflag == 0)
  {
    is = pmb->is;
    ie = pmb->ie;
    js = pmb->js;
    je = pmb->je;
    ks = pmb->ks;
    ke = pmb->ke;
    ng = NGHOST;
  }
  else
  {
    is = pmb->cis;
    ie = pmb->cie;
    js = pmb->cjs;
    je = pmb->cje;
    ks = pmb->cks;
    ke = pmb->cke;
    ng = pmb->cnghost;
  }

  // Set face-centered positions and distances
  AllocateAndSetBasicCoordinates();

  // Set parameters
  bh_mass_ = pin->GetReal("coord", "m");
  bh_spin_ = pin->GetReal("coord", "a");
  const Real &m = bh_mass_;
  const Real &a = bh_spin_;

  // Initialize volume-averaged positions and spacings: r-direction
  for (int i = is-ng; i <= ie+ng; ++i)
  {
    Real r_m = x1f(i);
    Real r_p = x1f(i+1);
    x1v(i) = 0.5 * (r_m + r_p);  // approximate
  }
  for (int i = is-ng; i <= ie+ng-1; ++i)
    dx1v(i) = x1v(i+1) - x1v(i);

  // Initialize volume-averaged positions and spacings: theta-direction
  if (pmb->block_size.nx2 == 1)  // no extent
  {
    Real theta_m = x2f(js);
    Real theta_p = x2f(js+1);
    x2v(js) = 0.5 * (theta_m + theta_p);  // approximate
    dx2v(js) = dx2f(js);
  }
  else  // extended
  {
    for (int j = js-ng; j <= je+ng; ++j)
    {
      Real theta_m = x2f(j);
      Real theta_p = x2f(j+1);
      x2v(j) = 0.5 * (theta_m + theta_p);  // approximate
    }
    for (int j = js-ng; j <= je+ng-1; ++j)
      dx2v(j) = x2v(j+1) - x2v(j);
  }

  // Initialize volume-averaged positions and spacings: phi-direction
  if (pmb->block_size.nx3 == 1)  // no extent
  {
    Real phi_m = x3f(ks);
    Real phi_p = x3f(ks+1);
    x3v(ks) = 0.5 * (phi_m + phi_p);
    dx3v(ks) = dx3f(ks);
  }
  else  // extended
  {
    for (int k = ks-ng; k <= ke+ng; ++k)
    {
      Real phi_m = x3f(k);
      Real phi_p = x3f(k+1);
      x3v(k) = 0.5 * (phi_m + phi_p);
    }
    for (int k = ks-ng; k <= ke+ng-1; ++k)
      dx3v(k) = x3v(k+1) - x3v(k);
  }

  // Prepare for MHD mesh refinement
  if (pmb->pmy_mesh->multilevel && MAGNETIC_FIELDS_ENABLED)
  {
    for (int i = is-NGHOST; i <= ie+NGHOST; ++i)
      x1s2(i) = x1s3(i) = x1v(i);
    if (pmb->block_size.nx2 == 1)
      x2s1(js) = x2s3(js) = x2v(js);
    else
      for (int j = js-NGHOST; j <= je+NGHOST; ++j)
        x2s1(j) = x2s3(j) = x2v(j);
    if (pmb->block_size.nx3 == 1)
      x3s1(ks) = x3s2(ks) = x3v(ks);
    else
      for (int k = ks-NGHOST; k <= ke+NGHOST; ++k)
        x3s1(k) = x3s2(k) = x3v(k);
  }

  // Allocate arrays for intermediate geometric quantities: r-direction
  int n_cells_1 = pmb->block_size.nx1 + 2*ng;
  coord_vol_i1_.NewAthenaArray(n_cells_1);
  coord_vol_i2_.NewAthenaArray(n_cells_1);
  coord_area1_i1_.NewAthenaArray(n_cells_1+1);
  coord_area2_i1_.NewAthenaArray(n_cells_1);
  coord_area2_i2_.NewAthenaArray(n_cells_1);
  coord_area3_i1_.NewAthenaArray(n_cells_1);
  coord_area3_i2_.NewAthenaArray(n_cells_1);
  coord_len1_i1_.NewAthenaArray(n_cells_1);
  coord_len1_i2_.NewAthenaArray(n_cells_1);
  coord_len2_i1_.NewAthenaArray(n_cells_1+1);
  coord_len3_i1_.NewAthenaArray(n_cells_1+1);
  coord_width1_i1_.NewAthenaArray(n_cells_1);
  coord_width2_i1_.NewAthenaArray(n_cells_1);
  coord_src_i1_.NewAthenaArray(n_cells_1);
  metric_cell_i1_.NewAthenaArray(n_cells_1);
  metric_face1_i1_.NewAthenaArray(n_cells_1+1);
  metric_face2_i1_.NewAthenaArray(n_cells_1);
  metric_face3_i1_.NewAthenaArray(n_cells_1);
  g_.NewAthenaArray(NMETRIC, n_cells_1+1);
  gi_.NewAthenaArray(NMETRIC, n_cells_1+1);

  // Allocate arrays for intermediate geometric quantities: theta-direction
  int n_cells_2 = (pmb->block_size.nx2 > 1) ? pmb->block_size.nx2 + 2*ng : 1;
  coord_vol_j1_.NewAthenaArray(n_cells_2);
  coord_vol_j2_.NewAthenaArray(n_cells_2);
  coord_area1_j1_.NewAthenaArray(n_cells_2);
  coord_area1_j2_.NewAthenaArray(n_cells_2);
  coord_area2_j1_.NewAthenaArray(n_cells_2+1);
  coord_area2_j2_.NewAthenaArray(n_cells_2+1);
  coord_area3_j1_.NewAthenaArray(n_cells_2);
  coord_area3_j2_.NewAthenaArray(n_cells_2);
  coord_len1_j1_.NewAthenaArray(n_cells_2+1);
  coord_len1_j2_.NewAthenaArray(n_cells_2+1);
  coord_len2_j1_.NewAthenaArray(n_cells_2);
  coord_len2_j2_.NewAthenaArray(n_cells_2);
  coord_len3_j1_.NewAthenaArray(n_cells_2+1);
  coord_len3_j2_.NewAthenaArray(n_cells_2+1);
  coord_width2_j1_.NewAthenaArray(n_cells_2);
  coord_width3_j1_.NewAthenaArray(n_cells_2);
  coord_width3_j2_.NewAthenaArray(n_cells_2);
  coord_width3_j3_.NewAthenaArray(n_cells_2);
  coord_src_j1_.NewAthenaArray(n_cells_2);
  coord_src_j2_.NewAthenaArray(n_cells_2);
  metric_cell_j1_.NewAthenaArray(n_cells_2);
  metric_cell_j2_.NewAthenaArray(n_cells_2);
  metric_face1_j1_.NewAthenaArray(n_cells_2);
  metric_face1_j2_.NewAthenaArray(n_cells_2);
  metric_face2_j1_.NewAthenaArray(n_cells_2+1);
  metric_face2_j2_.NewAthenaArray(n_cells_2+1);
  metric_face3_j1_.NewAthenaArray(n_cells_2);
  metric_face3_j2_.NewAthenaArray(n_cells_2);

  // Allocate arrays for intermediate geometric quantities: phi-direction
  int n_cells_3 = (pmb->block_size.nx3 > 1) ? pmb->block_size.nx3 + 2*ng : 1;
  coord_vol_k1_.NewAthenaArray(n_cells_3);
  coord_area1_k1_.NewAthenaArray(n_cells_3);
  coord_area2_k1_.NewAthenaArray(n_cells_3);
  coord_len3_k1_.NewAthenaArray(n_cells_3);
  coord_width3_k1_.NewAthenaArray(n_cells_3);

  // Allocate arrays for intermediate geometric quantities: r-theta-direction
  coord_width3_ji1_.NewAthenaArray(n_cells_2, n_cells_1);
  trans_face1_ji1_.NewAthenaArray(n_cells_2, n_cells_1+1);
  trans_face1_ji2_.NewAthenaArray(n_cells_2, n_cells_1+1);
  trans_face1_ji3_.NewAthenaArray(n_cells_2, n_cells_1+1);
  trans_face1_ji4_.NewAthenaArray(n_cells_2, n_cells_1+1);
  trans_face1_ji5_.NewAthenaArray(n_cells_2, n_cells_1+1);
  trans_face1_ji6_.NewAthenaArray(n_cells_2, n_cells_1+1);
  trans_face1_ji7_.NewAthenaArray(n_cells_2, n_cells_1+1);
  trans_face2_ji1_.NewAthenaArray(n_cells_2+1, n_cells_1);
  trans_face2_ji2_.NewAthenaArray(n_cells_2+1, n_cells_1);
  trans_face2_ji3_.NewAthenaArray(n_cells_2+1, n_cells_1);
  trans_face2_ji4_.NewAthenaArray(n_cells_2+1, n_cells_1);
  trans_face2_ji5_.NewAthenaArray(n_cells_2+1, n_cells_1);
  trans_face2_ji6_.NewAthenaArray(n_cells_2+1, n_cells_1);
  trans_face3_ji1_.NewAthenaArray(n_cells_2, n_cells_1);
  trans_face3_ji2_.NewAthenaArray(n_cells_2, n_cells_1);
  trans_face3_ji3_.NewAthenaArray(n_cells_2, n_cells_1);
  trans_face3_ji4_.NewAthenaArray(n_cells_2, n_cells_1);
  trans_face3_ji5_.NewAthenaArray(n_cells_2, n_cells_1);
  trans_face3_ji6_.NewAthenaArray(n_cells_2, n_cells_1);

  // Calculate intermediate geometric quantities: r-direction
  int il = is - ng;
  int iu = ie + ng;
  for (int i = il; i <= iu; ++i)
  {
    // Useful quantities
    Real r_c = x1v(i);
    Real r_m = x1f(i);
    Real r_p = x1f(i+1);
    Real r_c_sq = SQR(r_c);
    Real r_m_sq = SQR(r_m);
    Real r_p_sq = SQR(r_p);
    Real r_m_cu = r_m * r_m_sq;
    Real rm_m = std::sqrt(r_m_sq + SQR(m));
    Real rm_p = std::sqrt(r_p_sq + SQR(m));

    // Volumes, areas, lengths, and widths
    coord_vol_i1_(i) = dx1f(i);
    coord_vol_i2_(i) = r_m_sq + r_m * r_p + r_p_sq;
    coord_area1_i1_(i) = SQR(r_m);
    if (i == iu)
      coord_area1_i1_(i+1) = SQR(r_p);
    coord_area2_i1_(i) = coord_vol_i1_(i);
    coord_area2_i2_(i) = coord_vol_i2_(i);
    coord_area3_i1_(i) = coord_vol_i1_(i);
    coord_area3_i2_(i) = coord_vol_i2_(i);
    coord_len1_i1_(i) = coord_vol_i1_(i);
    coord_len1_i2_(i) = coord_vol_i2_(i);
    coord_len2_i1_(i) = coord_area1_i1_(i);
    if (i == iu)
      coord_len2_i1_(i+1) = coord_area1_i1_(i+1);
    coord_len3_i1_(i) = coord_area1_i1_(i);
    if (i == iu)
      coord_len3_i1_(i+1) = coord_area1_i1_(i+1);
    coord_width1_i1_(i) = rm_p - rm_m + m * std::log((rm_p + r_p) / (rm_m + r_m));
    coord_width2_i1_(i) = x1v(i);

    // Source terms
    coord_src_i1_(i) = r_c;

    // Metric coefficients
    metric_cell_i1_(i) = r_c;
    metric_face1_i1_(i) = r_m;
    if (i == iu)
      metric_face1_i1_(i+1) = r_p;
    metric_face2_i1_(i) = r_c;
    metric_face3_i1_(i) = r_c;
  }

  // Calculate intermediate geometric quantities: theta-direction
  int jl, ju;
  if (n_cells_2 > 1)  // extended
  {
    jl = js - ng;
    ju = je + ng;
  }
  else  // no extent
  {
    jl = js;
    ju = js;
  }
  for (int j = jl; j <= ju; ++j)
  {
    // Useful quantities
    Real theta_c = x2v(j);
    Real theta_m = x2f(j);
    Real theta_p = x2f(j+1);
    Real sin_c = std::sin(theta_c);
    Real sin_m = std::sin(theta_m);
    Real sin_p = std::sin(theta_p);
    Real cos_c = std::cos(theta_c);
    Real cos_m = std::cos(theta_m);
    Real cos_p = std::cos(theta_p);
    Real sin_c_sq = SQR(sin_c);
    Real sin_m_sq = SQR(sin_m);
    Real sin_p_sq = SQR(sin_p);
    Real cos_c_sq = SQR(cos_c);
    Real cos_m_sq = SQR(cos_m);
    Real cos_p_sq = SQR(cos_p);
    Real sin_m_cu = sin_m_sq*sin_m;
    Real sin_p_cu = SQR(sin_p)*sin_p;

    // Volumes, areas, lengths, and widths
    coord_vol_j1_(j) = std::abs(cos_m - cos_p);
    coord_vol_j2_(j) = SQR(a) * (cos_m_sq + cos_m * cos_p + cos_p_sq);
    coord_area1_j1_(j) = coord_vol_j1_(j);
    coord_area1_j2_(j) = coord_vol_j2_(j);
    coord_area2_j1_(j) = std::abs(sin_m);
    if (j == ju)
      coord_area2_j1_(j+1) = std::abs(sin_p);
    coord_area2_j2_(j) = SQR(a) * cos_m_sq;
    if (j == ju)
      coord_area2_j2_(j+1) = SQR(a) * cos_p_sq;
    coord_area3_j1_(j) = coord_vol_j1_(j);
    coord_area3_j2_(j) = coord_vol_j2_(j);
    coord_len1_j1_(j) = coord_area2_j1_(j);
    if (j == ju)
      coord_len1_j1_(j+1) = coord_area2_j1_(j+1);
    coord_len1_j2_(j) = coord_area2_j2_(j);
    if (j == ju)
      coord_len1_j2_(j+1) = coord_area2_j2_(j+1);
    coord_len2_j1_(j) = coord_vol_j1_(j);
    coord_len2_j2_(j) = coord_vol_j2_(j);
    coord_len3_j1_(j) = coord_area2_j1_(j);
    if (j == ju)
      coord_len3_j1_(j+1) = coord_area2_j1_(j+1);
    coord_len3_j2_(j) = coord_area2_j2_(j);
    if (j == ju)
      coord_len3_j2_(j+1) = coord_area2_j2_(j+1);
    coord_width2_j1_(j) = dx2f(j);
    coord_width3_j1_(j) = std::abs(sin_c);
    coord_width3_j2_(j) = SQR(a) * sin_c_sq;
    coord_width3_j3_(j) = SQR(a) * cos_c_sq;

    // Source terms
    coord_src_j1_(j) = sin_c;
    coord_src_j2_(j) = cos_c;

    // Metric coefficients
    metric_cell_j1_(j) = sin_c_sq;
    metric_cell_j2_(j) = cos_c_sq;
    metric_face1_j1_(j) = sin_c_sq;
    metric_face1_j2_(j) = cos_c_sq;
    metric_face2_j1_(j) = sin_m_sq;
    if (j == ju)
      metric_face2_j1_(j) = sin_p_sq;
    metric_face2_j2_(j) = cos_m_sq;
    if (j == ju)
      metric_face2_j2_(j) = cos_p_sq;
    metric_face3_j1_(j) = sin_c_sq;
    metric_face3_j2_(j) = cos_c_sq;
  }

  // Calculate intermediate geometric quantities: phi-direction
  int kl, ku;
  if (n_cells_3 > 1)  // extended
  {
    kl = ks - ng;
    ku = ke + ng;
  }
  else  // no extent
  {
    kl = ks;
    ku = ks;
  }
  for (int k = kl; k <= ku; ++k)
  {
    // Volumes, areas, lengths, and widths
    coord_vol_k1_(k) = dx3f(k);
    coord_area1_k1_(k) = coord_vol_k1_(k);
    coord_area2_k1_(k) = coord_vol_k1_(k);
    coord_len3_k1_(k) = coord_vol_k1_(k);
    coord_width3_k1_(k) = coord_vol_k1_(k);
  }

  // Calculate intermediate geometric quantities: r-theta-direction
  for (int j = jl; j <= ju; ++j)
    for (int i = il; i <= iu; ++i)
    {
      // Useful quantities
      Real a2 = SQR(a);
      Real r_c = x1v(i);
      Real r_m = x1f(i);
      Real r_p = x1f(i+1);
      Real r_c_sq = SQR(r_c);
      Real r_m_sq = SQR(r_m);
      Real r_p_sq = SQR(r_p);
      Real r_m_qu = SQR(r_m_sq);
      Real r_p_qu = SQR(r_p_sq);
      Real theta_c = x2v(j);
      Real theta_m = x2f(j);
      Real theta_p = x2f(j+1);
      Real sin_c = std::sin(theta_c);
      Real sin_m = std::sin(theta_m);
      Real sin_p = std::sin(theta_p);
      Real cos_c = std::cos(theta_c);
      Real cos_m = std::cos(theta_m);
      Real cos_p = std::cos(theta_p);
      Real sin_c_sq = SQR(sin_c);
      Real sin_m_sq = SQR(sin_m);
      Real sin_p_sq = SQR(sin_p);
      Real cos_c_sq = SQR(cos_c);
      Real cos_m_sq = SQR(cos_m);
      Real cos_p_sq = SQR(cos_p);
      Real delta_m = r_m_sq - 2.0*m*r_m + a2;
      Real delta_p = r_p_sq - 2.0*m*r_p + a2;
      Real sigma_cc = r_c_sq + a2 * cos_c_sq;
      Real sigma_cm = r_c_sq + a2 * cos_m_sq;
      Real sigma_cp = r_c_sq + a2 * cos_p_sq;
      Real sigma_mc = r_m_sq + a2 * cos_c_sq;
      Real sigma_pc = r_p_sq + a2 * cos_c_sq;

      // Volumes, areas, lengths, and widths
      coord_width3_ji1_(j,i) = std::sqrt(r_c_sq + a2
          + 2.0 * m * a2 * r_c * sin_c_sq / (r_c_sq + a2 * cos_c_sq));

      // Coordinate transformations
      trans_face1_ji1_(j,i) = 1.0 / std::sqrt(1.0 + 2.0*m*r_m/sigma_mc);
      if (i == iu)
        trans_face1_ji1_(j,i+1) = 1.0 / std::sqrt(1.0 + 2.0*m*r_p/sigma_pc);
      trans_face1_ji2_(j,i) = std::sqrt((sigma_mc + 2.0*m*r_m)
          / (r_m_sq + a2 + 2.0*m*a2*r_m/sigma_mc * sin_c_sq));
      if (i == iu)
        trans_face1_ji2_(j,i+1) = std::sqrt((sigma_pc + 2.0*m*r_p)
            / (r_p_sq + a2 + 2.0*m*a2*r_p/sigma_pc * sin_c_sq));
      trans_face1_ji3_(j,i) = std::sqrt(sigma_mc);
      if (i == iu)
        trans_face1_ji3_(j,i+1) = std::sqrt(sigma_pc);
      trans_face1_ji4_(j,i) = std::abs(sin_c)
          * std::sqrt(r_m_sq + a2 + 2.0*m*a2*r_m/sigma_mc * sin_c_sq);
      if (i == iu)
        trans_face1_ji4_(j,i+1) = std::abs(sin_c)
            * std::sqrt(r_p_sq + a2 + 2.0*m*a2*r_p/sigma_pc * sin_c_sq);
      trans_face1_ji5_(j,i) = 2.0*m*r_m * std::sqrt(sigma_mc
          / ((sigma_mc + 2.0*m*r_m)
          * (r_m_qu + a2*r_m_sq + 2.0*m*a2*r_m + delta_m*a2*cos_c_sq)));
      if (i == iu)
        trans_face1_ji5_(j,i+1) = 2.0*m*r_p * std::sqrt(sigma_pc
            / ((sigma_pc + 2.0*m*r_p)
            * (r_p_qu + a2*r_p_sq + 2.0*m*a2*r_p + delta_p*a2*cos_c_sq)));
      trans_face1_ji6_(j,i) = -2.0*m*a*r_m * std::abs(sin_c)
          * std::sqrt(r_m_sq + a2 + 2.0*m*a2*r_m/sigma_mc * sin_c_sq)
          / (r_m_qu + a2*r_m_sq + 2.0*m*a2*r_m + delta_m*a2*cos_c_sq);
      if (i == iu)
        trans_face1_ji6_(j,i+1) = -2.0*m*a*r_p * std::abs(sin_c)
            * std::sqrt(r_p_sq + a2 + 2.0*m*a2*r_p/sigma_pc * sin_c_sq)
            / (r_p_qu + a2*r_p_sq + 2.0*m*a2*r_p + delta_p*a2*cos_c_sq);
      trans_face1_ji7_(j,i) = -a * (sigma_mc + 2.0*m*r_m) * std::abs(sin_c)
          * std::sqrt(r_m_sq + a2 + 2.0*m*a2*r_m/sigma_mc * sin_c_sq)
          / (r_m_qu + a2*r_m_sq + 2.0*m*a2*r_m + delta_m*a2*cos_c_sq);
      if (i == iu)
        trans_face1_ji7_(j,i+1) = -a * (sigma_pc + 2.0*m*r_p) * std::abs(sin_c)
            * std::sqrt(r_p_sq + a2 + 2.0*m*a2*r_p/sigma_pc * sin_c_sq)
            / (r_p_qu + a2*r_p_sq + 2.0*m*a2*r_p + delta_p*a2*cos_c_sq);
      trans_face2_ji1_(j,i) = 1.0 / std::sqrt(1.0 + 2.0*m*r_c/sigma_cm);
      if (j == ju)
        trans_face2_ji1_(j+1,i) = 1.0 / std::sqrt(1.0 + 2.0*m*r_c/sigma_cp);
      trans_face2_ji2_(j,i) = std::sqrt(1.0 + 2.0*m*r_c/sigma_cm);
      if (j == ju)
        trans_face2_ji2_(j+1,i) = std::sqrt(1.0 + 2.0*m*r_c/sigma_cp);
      trans_face2_ji3_(j,i) = std::sqrt(sigma_cm);
      if (j == ju)
        trans_face2_ji3_(j+1,i) = std::sqrt(sigma_cp);
      trans_face2_ji4_(j,i) = std::sqrt(sigma_cm) * std::abs(sin_m);
      if (j == ju)
        trans_face2_ji4_(j+1,i) = std::sqrt(sigma_cp) * std::abs(sin_p);
      trans_face2_ji5_(j,i) = 2.0*m*r_c
          / (sigma_cm * std::sqrt(1.0 + 2.0*m*r_c/sigma_cm));
      if (j == ju)
        trans_face2_ji5_(j+1,i) = 2.0*m*r_c
            / (sigma_cp * std::sqrt(1.0 + 2.0*m*r_c/sigma_cp));
      trans_face2_ji6_(j,i) = -a * sin_m_sq * std::sqrt(1.0 + 2.0*m*r_c/sigma_cm);
      if (j == ju)
        trans_face2_ji6_(j+1,i) = -a * sin_p_sq * std::sqrt(1.0 + 2.0*m*r_c/sigma_cp);
      trans_face3_ji1_(j,i) = 1.0 / std::sqrt(1.0 + 2.0*m*r_c/sigma_cc);
      trans_face3_ji2_(j,i) = std::sqrt(1.0 + 2.0*m*r_c/sigma_cc);
      trans_face3_ji3_(j,i) = std::sqrt(sigma_cc);
      trans_face3_ji4_(j,i) = std::sqrt(sigma_cc) * std::abs(sin_c);
      trans_face3_ji5_(j,i) = 2.0*m*r_c
          / (sigma_cc * std::sqrt(1.0 + 2.0*m*r_c/sigma_cc));
      trans_face3_ji6_(j,i) = -a * sin_c_sq * std::sqrt(1.0 + 2.0*m*r_c/sigma_cc);
    }
}

//--------------------------------------------------------------------------------------

// Destructor
Coordinates::~Coordinates()
{
  DeleteBasicCoordinates();
  coord_vol_i1_.DeleteAthenaArray();
  coord_vol_i2_.DeleteAthenaArray();
  coord_area1_i1_.DeleteAthenaArray();
  coord_area2_i1_.DeleteAthenaArray();
  coord_area2_i2_.DeleteAthenaArray();
  coord_area3_i1_.DeleteAthenaArray();
  coord_area3_i2_.DeleteAthenaArray();
  coord_len1_i1_.DeleteAthenaArray();
  coord_len1_i2_.DeleteAthenaArray();
  coord_len2_i1_.DeleteAthenaArray();
  coord_len3_i1_.DeleteAthenaArray();
  coord_width1_i1_.DeleteAthenaArray();
  coord_width2_i1_.DeleteAthenaArray();
  coord_src_i1_.DeleteAthenaArray();
  metric_cell_i1_.DeleteAthenaArray();
  metric_face1_i1_.DeleteAthenaArray();
  metric_face2_i1_.DeleteAthenaArray();
  metric_face3_i1_.DeleteAthenaArray();
  g_.DeleteAthenaArray();
  gi_.DeleteAthenaArray();
  coord_vol_j1_.DeleteAthenaArray();
  coord_vol_j2_.DeleteAthenaArray();
  coord_area1_j1_.DeleteAthenaArray();
  coord_area1_j2_.DeleteAthenaArray();
  coord_area2_j1_.DeleteAthenaArray();
  coord_area2_j2_.DeleteAthenaArray();
  coord_area3_j1_.DeleteAthenaArray();
  coord_area3_j2_.DeleteAthenaArray();
  coord_len1_j1_.DeleteAthenaArray();
  coord_len1_j2_.DeleteAthenaArray();
  coord_len2_j1_.DeleteAthenaArray();
  coord_len2_j2_.DeleteAthenaArray();
  coord_len3_j1_.DeleteAthenaArray();
  coord_len3_j2_.DeleteAthenaArray();
  coord_width2_j1_.DeleteAthenaArray();
  coord_width3_j1_.DeleteAthenaArray();
  coord_width3_j2_.DeleteAthenaArray();
  coord_width3_j3_.DeleteAthenaArray();
  coord_src_j1_.DeleteAthenaArray();
  coord_src_j2_.DeleteAthenaArray();
  metric_cell_j1_.DeleteAthenaArray();
  metric_cell_j2_.DeleteAthenaArray();
  metric_face1_j1_.DeleteAthenaArray();
  metric_face1_j2_.DeleteAthenaArray();
  metric_face2_j1_.DeleteAthenaArray();
  metric_face2_j2_.DeleteAthenaArray();
  metric_face3_j1_.DeleteAthenaArray();
  metric_face3_j2_.DeleteAthenaArray();
  coord_vol_k1_.DeleteAthenaArray();
  coord_area1_k1_.DeleteAthenaArray();
  coord_area2_k1_.DeleteAthenaArray();
  coord_len3_k1_.DeleteAthenaArray();
  coord_width3_k1_.DeleteAthenaArray();
  coord_width3_ji1_.DeleteAthenaArray();
  trans_face1_ji1_.DeleteAthenaArray();
  trans_face1_ji2_.DeleteAthenaArray();
  trans_face1_ji3_.DeleteAthenaArray();
  trans_face1_ji4_.DeleteAthenaArray();
  trans_face1_ji5_.DeleteAthenaArray();
  trans_face1_ji6_.DeleteAthenaArray();
  trans_face1_ji7_.DeleteAthenaArray();
  trans_face2_ji1_.DeleteAthenaArray();
  trans_face2_ji2_.DeleteAthenaArray();
  trans_face2_ji3_.DeleteAthenaArray();
  trans_face2_ji4_.DeleteAthenaArray();
  trans_face2_ji5_.DeleteAthenaArray();
  trans_face2_ji6_.DeleteAthenaArray();
  trans_face3_ji1_.DeleteAthenaArray();
  trans_face3_ji2_.DeleteAthenaArray();
  trans_face3_ji3_.DeleteAthenaArray();
  trans_face3_ji4_.DeleteAthenaArray();
  trans_face3_ji5_.DeleteAthenaArray();
  trans_face3_ji6_.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------

// Function for computing cell volumes
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   volumes: 1D array of cell volumes
// Notes:
//   \Delta V = 1/3 * (r_+ - r_-) |\cos\theta_- - \cos\theta_+| (\phi_+ - \phi_-)
//       * (r_-^2 + r_- r_+ + r_+^2
//       + a^2 (\cos^2\theta_- + \cos\theta_- \cos\theta_+ + \cos^2\theta_+))
//   cf. GetCellVolume()
void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &volumes)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    volumes(i) = GetCellVolume(k, j, i);
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing single cell volume
// Inputs:
//   k,j,i: phi-, theta-, and r-indices
// Outputs:
//   returned value: cell volume
// Notes:
//   \Delta V = 1/3 * (r_+ - r_-) |\cos\theta_- - \cos\theta_+| (\phi_+ - \phi_-)
//       * (r_-^2 + r_- r_+ + r_+^2
//       + a^2 (\cos^2\theta_- + \cos\theta_- \cos\theta_+ + \cos^2\theta_+))
//   cf. CellVolume()
Real Coordinates::GetCellVolume(const int k, const int j, const int i)
{
  return 1.0/3.0 * coord_vol_i1_(i) * coord_vol_j1_(j) * coord_vol_k1_(k)
      * (coord_vol_i2_(i) + coord_vol_j2_(j));
}

//--------------------------------------------------------------------------------------

// Function for computing areas orthogonal to r
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to r
// Notes:
//   \Delta A = 1/3 * |\cos\theta_- - \cos\theta_+| (\phi_+ - \phi_-)
//       * (3 r_-^2 + a^2 (\cos^2\theta_- + \cos\theta_- \cos\theta_+ + \cos^2\theta_+))
//   cf. GetFace1Area()
void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    areas(i) = GetFace1Area(k, j, i);
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing single area orthogonal to r
// Inputs:
//   k,j,i: phi-, theta-, and r-indices
// Outputs:
//   returned value: interface area orthogonal to r
// Notes:
//   \Delta A = 1/3 * |\cos\theta_- - \cos\theta_+| (\phi_+ - \phi_-)
//       * (3 r_-^2 + a^2 (\cos^2\theta_- + \cos\theta_- \cos\theta_+ + \cos^2\theta_+))
//   cf. Face1Area()
Real Coordinates::GetFace1Area(const int k, const int j, const int i)
{
  return coord_area1_j1_(j) * coord_area1_k1_(k)
      * (coord_area1_i1_(i) + 1.0/3.0 * coord_area1_j2_(j));
}

//--------------------------------------------------------------------------------------

// Function for computing areas orthogonal to theta
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to theta
// Notes:
//   \Delta A = 1/3 (r_+ - r_-) |\sin\theta_-| (\phi_+ - \phi_-)
//       * (r_-^2 + r_- r_+ + r_+^2 + 3 a^2 \cos^2\theta_-)
void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    areas(i) = coord_area2_i1_(i) * coord_area2_j1_(j) * coord_area2_k1_(k)
               * (1.0/3.0 * coord_area2_i2_(i) + coord_area2_j2_(j));
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing areas orthogonal to phi
// Inputs:
//   k: phi-index (unused)
//   j: theta-index
//   il,iu: r-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to phi
// Notes:
//   \Delta A = 1/3 (r_+ - r_-) |\cos\theta_- - \cos\theta_+|
//       * (r_-^2 + r_- r_+ + r_+^2
//       + a^2 (\cos^2\theta_- + \cos\theta_- \cos\theta_+ + \cos^2\theta_+))
void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &areas)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    areas(i) = 1.0/3.0 * coord_area3_i1_(i) * coord_area3_j1_(j)
             * (coord_area3_i2_(i) + coord_area3_j2_(j));
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the r-direction
// Inputs:
//   k: phi-index (unused)
//   j: theta-index
//   il,iu: r-index bounds
// Outputs:
//   lengths: 1D array of edge lengths along r
// Notes:
//   \Delta L = 1/3 (r_+ - r_-) |\sin\theta_-|
//       * (r_-^2 + r_- r_+ + r_+^2 + 3 a^2 \cos^2\theta_-)
void Coordinates::Edge1Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &lengths)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    lengths(i) = coord_len1_i1_(i) * coord_len1_j1_(j)
        * (1.0/3.0 * coord_len1_i2_(i) + coord_len1_j2_(j));
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the theta-direction
// Inputs:
//   k: phi-index (unused)
//   j: theta-index
//   il,iu: r-index bounds
// Outputs:
//   lengths: 1D array of edge lengths along r
// Notes:
//   \Delta L = 1/3 * |\cos\theta_- - \cos\theta_+|
//       * (3 r_-^2 + a^2 (\cos^2\theta_- + \cos\theta_- \cos\theta_+ + \cos^2\theta_+))
//   cf. GetEdge2Length()
void Coordinates::Edge2Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &lengths)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    lengths(i) = GetEdge2Length(k, j, i);
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing single length of edge in the theta-direction
// Inputs:
//   k: phi-index (unused)
//   j,i: theta- and r-indices
// Outputs:
//   returned value: length of edge along theta
// Notes:
//   \Delta L = 1/3 * |\cos\theta_- - \cos\theta_+|
//       * (3 r_-^2 + a^2 (\cos^2\theta_- + \cos\theta_- \cos\theta_+ + \cos^2\theta_+))
//   cf. Edge2Length()
Real Coordinates::GetEdge2Length(const int k, const int j, const int i)
{
  return coord_len2_j1_(j) * (coord_len2_i1_(i) + 1.0/3.0 * coord_len2_j2_(j));
}

//--------------------------------------------------------------------------------------

// Function for computing lengths of edges in the phi-direction
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   lengths: 1D array of edge lengths along r
// Notes:
//   \Delta L = (\phi_+ - \phi_-) |\sin\theta_-| (r_-^2 + a^2 \cos^2\theta_-)
//   cf. GetEdge3Length()
void Coordinates::Edge3Length(const int k, const int j, const int il, const int iu,
  AthenaArray<Real> &lengths)
{
  #pragma simd
  for (int i = il; i <= iu; ++i)
    lengths(i) = GetEdge3Length(k, j, i);
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing single length of edge in the phi-direction
// Inputs:
//   k,j,i: phi-, theta-, and r-indices
// Outputs:
//   returned value: length of edge along phi
// Notes:
//   \Delta L = (\phi_+ - \phi_-) |\sin\theta_-| (r_-^2 + a^2 \cos^2\theta_-)
//   cf. Edge3Length()
Real Coordinates::GetEdge3Length(const int k, const int j, const int i)
{
  return coord_len3_k1_(k) * coord_len3_j1_(j)
      * (coord_len3_i1_(i) + coord_len3_j2_(j));
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the r-direction
// Inputs:
//   k,j: phi- and theta-indices (unused)
//   i: r-index
// Outputs:
//   returned value: r-width of cell (i,j,k)
// Notes:
//   \Delta W >= \sqrt{r_+^2 + M^2} - \sqrt{r_-^2 + M^2}
//       + M \log{(\sqrt{r_+^2 + M^2} + r_+) / (\sqrt{r_-^2 + M^2} + r_-)}
Real Coordinates::CenterWidth1(const int k, const int j, const int i)
{
  return coord_width1_i1_(i);
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the theta-direction
// Inputs:
//   k: phi-index (unused)
//   j,i: theta- and r-indices
// Outputs:
//   returned value: theta-width of cell (i,j,k)
// Notes:
//   \Delta W >= r (\theta_+ - \theta_-)
Real Coordinates::CenterWidth2(const int k, const int j, const int i)
{
  return coord_width2_i1_(i) * coord_width2_j1_(j);
}

//--------------------------------------------------------------------------------------

// Function for computing widths of cells in the phi-direction
// Inputs:
//   k,j,i: phi-, theta-, and r-indices
// Outputs:
//   returned value: phi-width of cell (i,j,k)
// Notes:
//   \Delta W = |\sin\theta| (\phi_+ - \phi_-)
//       * \sqrt{r^2 + a^2 + 2 M a^2 r \sin^2\theta / (r^2 + a^2 \cos^2\theta)}
Real Coordinates::CenterWidth3(const int k, const int j, const int i)
{
  return coord_width3_j1_(j) * coord_width3_k1_(k) * coord_width3_ji1_(j,i);
}

//--------------------------------------------------------------------------------------

// Function for computing source terms using fluxes
// Inputs:
//   dt: size of timestep
//   flux: 3D array of fluxes
//   prim: 3D array of primitive values at beginning of half timestep
//   bb_cc: 3D array of cell-centered magnetic fields
// Outputs:
//   cons: source terms added to 3D array of conserved variables
// Notes:
//   all source terms computed in this function
void Coordinates::CoordSrcTerms(const Real dt, const AthenaArray<Real> *flux,
    const AthenaArray<Real> &prim, const AthenaArray<Real> &bb_cc,
    AthenaArray<Real> &cons)
{
  // Extract ratio of specific heats
  const Real gamma_adi = pmy_block->phydro->peos->GetGamma();

  // Extract useful quantities that do not depend on location
  const Real &m = bh_mass_;
  Real m2 = SQR(m);
  const Real &a = bh_spin_;
  Real a2 = SQR(a);
  Real a3 = a * a2;

  // Go through cells
  for (int k = pmy_block->ks; k <= pmy_block->ke; ++k)
  {
    #pragma omp parallel for schedule(static)
    for (int j = pmy_block->js; j <= pmy_block->je; ++j)
    {
      // Extract useful quantities that do not depend on r
      const Real &sin = coord_src_j1_(j);
      Real sin2 = SQR(sin);
      Real sin3 = sin * sin2;
      const Real &cos = coord_src_j2_(j);
      Real cos2 = SQR(cos);
      Real cos4 = SQR(cos2);
      Real cot = cos/sin;

      // Calculate metric coefficients
      CellMetric(k, j, pmy_block->is, pmy_block->ie, g_, gi_);

      // Go through 1D slice
      #pragma simd
      for (int i = pmy_block->is; i <= pmy_block->ie; ++i)
      {
        // Extract metric coefficients
        const Real &g_00 = g_(I00,i);
        const Real &g_01 = g_(I01,i);
        const Real &g_02 = 0.0;
        const Real &g_03 = g_(I03,i);
        const Real &g_10 = g_(I01,i);
        const Real &g_11 = g_(I11,i);
        const Real &g_12 = 0.0;
        const Real &g_13 = g_(I13,i);
        const Real &g_20 = 0.0;
        const Real &g_21 = 0.0;
        const Real &g_22 = g_(I22,i);
        const Real &g_23 = 0.0;
        const Real &g_30 = g_(I03,i);
        const Real &g_31 = g_(I13,i);
        const Real &g_32 = 0.0;
        const Real &g_33 = g_(I33,i);
        const Real &g01 = gi_(I01,i);
        const Real &g02 = 0.0;
        const Real &g03 = 0.0;
        Real alpha = std::sqrt(-1.0/gi_(I00,i));

        // Extract primitives
        const Real &rho = prim(IDN,k,j,i);
        const Real &pgas = prim(IEN,k,j,i);
        const Real &uu1 = prim(IVX,k,j,i);
        const Real &uu2 = prim(IVY,k,j,i);
        const Real &uu3 = prim(IVZ,k,j,i);

        // Calculate 4-velocity
        Real tmp = g_11*uu1*uu1 + 2.0*g_12*uu1*uu2 + 2.0*g_13*uu1*uu3
                 + g_22*uu2*uu2 + 2.0*g_23*uu2*uu3
                 + g_33*uu3*uu3;
        Real gamma = std::sqrt(1.0 + tmp);
        Real u0 = gamma / alpha;
        Real u1 = uu1 - alpha * gamma * g01;
        Real u2 = uu2 - alpha * gamma * g02;
        Real u3 = uu3 - alpha * gamma * g03;
        Real u_0 = g_00*u0 + g_01*u1 + g_02*u2 + g_03*u3;
        Real u_1 = g_10*u0 + g_11*u1 + g_12*u2 + g_13*u3;
        Real u_2 = g_20*u0 + g_21*u1 + g_22*u2 + g_23*u3;
        Real u_3 = g_30*u0 + g_31*u1 + g_32*u2 + g_33*u3;

        // Extract and calculate magnetic field
        Real b0 = 0.0, b1 = 0.0, b2 = 0.0, b3 = 0.0;
        Real b_0 = 0.0, b_1 = 0.0, b_2 = 0.0, b_3 = 0.0;
        Real b_sq = 0.0;
        if (MAGNETIC_FIELDS_ENABLED)
        {
          const Real &bb1 = bb_cc(IB1,k,j,i);
          const Real &bb2 = bb_cc(IB2,k,j,i);
          const Real &bb3 = bb_cc(IB3,k,j,i);
          b0 =
                g_10*bb1*u0 + g_11*bb1*u1 + g_12*bb1*u2 + g_13*bb1*u3
              + g_20*bb2*u0 + g_21*bb2*u1 + g_22*bb2*u2 + g_23*bb2*u3
              + g_30*bb3*u0 + g_31*bb3*u1 + g_32*bb3*u2 + g_33*bb3*u3;
          b1 = (bb1 + b0 * u1) / u0;
          b2 = (bb2 + b0 * u2) / u0;
          b3 = (bb3 + b0 * u3) / u0;
          b_0 = g_00*b0 + g_01*b1 + g_02*b2 + g_03*b3;
          b_1 = g_10*b0 + g_11*b1 + g_12*b2 + g_13*b3;
          b_2 = g_20*b0 + g_21*b1 + g_22*b2 + g_23*b3;
          b_3 = g_30*b0 + g_31*b1 + g_32*b2 + g_33*b3;
          b_sq = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3;
        }

        // Calculate stress-energy tensor
        Real wtot = rho + gamma_adi/(gamma_adi-1.0) * pgas + b_sq;
        Real ptot = pgas + 0.5*b_sq;
        Real t0_0 = wtot*u0*u_0 - b0*b_0 + ptot;
        Real t0_1 = wtot*u0*u_1 - b0*b_1;
        Real t0_2 = wtot*u0*u_2 - b0*b_2;
        Real t0_3 = wtot*u0*u_3 - b0*b_3;
        Real t1_0 = wtot*u1*u_0 - b1*b_0;
        Real t1_1 = wtot*u1*u_1 - b1*b_1 + ptot;
        Real t1_2 = wtot*u1*u_2 - b1*b_2;
        Real t1_3 = wtot*u1*u_3 - b1*b_3;
        Real t2_0 = wtot*u2*u_0 - b2*b_0;
        Real t2_1 = wtot*u2*u_1 - b2*b_1;
        Real t2_2 = wtot*u2*u_2 - b2*b_2 + ptot;
        Real t2_3 = wtot*u2*u_3 - b2*b_3;
        Real t3_0 = wtot*u3*u_0 - b3*b_0;
        Real t3_1 = wtot*u3*u_1 - b3*b_1;
        Real t3_2 = wtot*u3*u_2 - b3*b_2;
        Real t3_3 = wtot*u3*u_3 - b3*b_3 + ptot;

        // Extract remaining useful quantities
        const Real &r = coord_src_i1_(i);
        Real r2 = SQR(r);
        Real r3 = r * r2;
        Real delta = r2 - 2.0*m*r + a2;
        Real sigma = r2 + a2 * cos2;
        Real sigma2 = SQR(sigma);
        Real sigma3 = sigma * sigma2;
        Real xi = r2 - a2 * cos2;

        // Calculate connection coefficients
        Real gamma0_00 = 2.0*m2*xi*r/sigma3;
        Real gamma0_01 = m*xi/sigma2 * (1.0 + 2.0*m*r/sigma);
        Real gamma0_02 = -2.0*m*a2*r/sigma2 * sin * cos;
        Real gamma0_03 = -2.0*m2*a*r*xi/sigma3 * sin2;
        Real gamma0_11 = 2.0*m*xi/sigma2 * (1.0 + m*r/sigma);
        Real gamma0_12 = -2.0*m*a2*r/sigma2 * sin * cos;
        Real gamma0_13 = -m*a*xi/sigma2 * (1.0 + 2.0*m*r/sigma) * sin2;
        Real gamma0_22 = -2.0*m*r2/sigma;
        Real gamma0_23 = 2.0*m*a3*r/sigma2 * sin3 * cos;
        Real gamma0_33 = -2.0*m*r/sigma * (r - m*a2*xi/sigma2 * sin2) * sin2;
        Real gamma1_00 = m*delta*xi/sigma3;
        Real gamma1_01 = m*xi/sigma3 * (delta-sigma);
        Real gamma1_02 = 0.0;
        Real gamma1_03 = -m*a*delta*xi/sigma3 * sin2;
        Real gamma1_11 = m*xi/sigma3 * (delta-2.0*sigma);
        Real gamma1_12 = -a2/sigma * sin * cos;
        Real gamma1_13 = a/sigma3 * (r3 * (r2+2.0*m2)
            + a2 * (r*a2*cos4 + (2.0*r3 + m*(delta-sigma)) * cos2 - m*r2*sin2));
        Real gamma1_22 = -r*delta/sigma;
        Real gamma1_23 = 0.0;
        Real gamma1_33 = -delta/sigma * (r - m*a2*xi/sigma2 * sin2) * sin2;
        Real gamma2_00 = -2.0*m*a2*r/sigma3 * sin * cos;
        Real gamma2_01 = -2.0*m*a2*r/sigma3 * sin * cos;
        Real gamma2_02 = 0.0;
        Real gamma2_03 = 2.0*m*a*r/sigma3 * (r2+a2) * sin * cos;
        Real gamma2_11 = -2.0*m*a2*r/sigma3 * sin * cos;
        Real gamma2_12 = r/sigma;
        Real gamma2_13 = a/sigma * (1.0 + 2.0*m*r/sigma2 * (r2+a2)) * sin * cos;
        Real gamma2_22 = -a2/sigma * sin * cos;
        Real gamma2_23 = 0.0;
        Real gamma2_33 = -1.0/sigma * (delta + 2.0*m*r/sigma2 * SQR(r2+a2)) * sin * cos;
        Real gamma3_00 = m*a*xi/sigma3;
        Real gamma3_01 = m*a*xi/sigma3;
        Real gamma3_02 = -2.0*m*a*r/sigma2 * cot;
        Real gamma3_03 = -m*a2*xi/sigma3 * sin2;
        Real gamma3_11 = m*a*xi/sigma3;
        Real gamma3_12 = -a/sigma2 * (2.0*m*r + sigma) * cot;
        Real gamma3_13 = 1.0/sigma3 * (r*sigma2 - m*a2*xi * sin2);
        Real gamma3_22 = -a*r/sigma;
        Real gamma3_23 = 2.0*m*a2*r/sigma2 * sin * cos + cot;
        Real gamma3_33 = a/sigma3 * (m*a2*xi * sin2 - r*sigma2) * sin2;
        const Real &gamma0_10 = gamma0_01;
        const Real &gamma0_20 = gamma0_02;
        const Real &gamma0_21 = gamma0_12;
        const Real &gamma0_30 = gamma0_03;
        const Real &gamma0_31 = gamma0_13;
        const Real &gamma0_32 = gamma0_23;
        const Real &gamma1_10 = gamma1_01;
        const Real &gamma1_20 = gamma1_02;
        const Real &gamma1_21 = gamma1_12;
        const Real &gamma1_30 = gamma1_03;
        const Real &gamma1_31 = gamma1_13;
        const Real &gamma1_32 = gamma1_23;
        const Real &gamma2_10 = gamma2_01;
        const Real &gamma2_20 = gamma2_02;
        const Real &gamma2_21 = gamma2_12;
        const Real &gamma2_30 = gamma2_03;
        const Real &gamma2_31 = gamma2_13;
        const Real &gamma2_32 = gamma2_23;
        const Real &gamma3_10 = gamma3_01;
        const Real &gamma3_20 = gamma3_02;
        const Real &gamma3_21 = gamma3_12;
        const Real &gamma3_30 = gamma3_03;
        const Real &gamma3_31 = gamma3_13;
        const Real &gamma3_32 = gamma3_23;

        // Calculate source terms
        Real s_0 = gamma0_00*t0_0 + gamma0_10*t1_0 + gamma0_20*t2_0 + gamma0_30*t3_0
                 + gamma1_00*t0_1 + gamma1_10*t1_1 + gamma1_20*t2_1 + gamma1_30*t3_1
                 + gamma2_00*t0_2 + gamma2_10*t1_2 + gamma2_20*t2_2 + gamma2_30*t3_2
                 + gamma3_00*t0_3 + gamma3_10*t1_3 + gamma3_20*t2_3 + gamma3_30*t3_3;
        Real s_1 = gamma0_01*t0_0 + gamma0_11*t1_0 + gamma0_21*t2_0 + gamma0_31*t3_0
                 + gamma1_01*t0_1 + gamma1_11*t1_1 + gamma1_21*t2_1 + gamma1_31*t3_1
                 + gamma2_01*t0_2 + gamma2_11*t1_2 + gamma2_21*t2_2 + gamma2_31*t3_2
                 + gamma3_01*t0_3 + gamma3_11*t1_3 + gamma3_21*t2_3 + gamma3_31*t3_3;
        Real s_2 = gamma0_02*t0_0 + gamma0_12*t1_0 + gamma0_22*t2_0 + gamma0_32*t3_0
                 + gamma1_02*t0_1 + gamma1_12*t1_1 + gamma1_22*t2_1 + gamma1_32*t3_1
                 + gamma2_02*t0_2 + gamma2_12*t1_2 + gamma2_22*t2_2 + gamma2_32*t3_2
                 + gamma3_02*t0_3 + gamma3_12*t1_3 + gamma3_22*t2_3 + gamma3_32*t3_3;
        Real s_3 = gamma0_03*t0_0 + gamma0_13*t1_0 + gamma0_23*t2_0 + gamma0_33*t3_0
                 + gamma1_03*t0_1 + gamma1_13*t1_1 + gamma1_23*t2_1 + gamma1_33*t3_1
                 + gamma2_03*t0_2 + gamma2_13*t1_2 + gamma2_23*t2_2 + gamma2_33*t3_2
                 + gamma3_03*t0_3 + gamma3_13*t1_3 + gamma3_23*t2_3 + gamma3_33*t3_3;

        // Extract conserved quantities
        Real &m_0 = cons(IEN,k,j,i);
        Real &m_1 = cons(IM1,k,j,i);
        Real &m_2 = cons(IM2,k,j,i);
        Real &m_3 = cons(IM3,k,j,i);

        // Add source terms to conserved quantities
        m_0 += dt * s_0;
        m_1 += dt * s_1;
        m_2 += dt * s_2;
        m_3 += dt * s_3;
      }
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing cell-centered metric coefficients
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::CellMetric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  // Extract useful quantities that do not depend on r
  const Real &m = bh_mass_;
  const Real &a = bh_spin_;
  Real a2 = SQR(a);
  const Real &sin2 = metric_cell_j1_(j);
  const Real &cos2 = metric_cell_j2_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract remaining useful quantities
    const Real &r = metric_cell_i1_(i);
    Real r2 = SQR(r);
    Real delta = r2 - 2.0*m*r + a2;
    Real sigma = r2 + a2 * cos2;

    // Set covariant metric coefficients
    g(I00,i) = -(1.0 - 2.0*m*r/sigma);
    g(I01,i) = 2.0*m*r/sigma;
    g(I03,i) = -2.0*m*a*r/sigma * sin2;
    g(I11,i) = 1.0 + 2.0*m*r/sigma;
    g(I13,i) = -(1.0 + 2.0*m*r/sigma) * a * sin2;
    g(I22,i) = sigma;
    g(I33,i) = (r2 + a2 + 2.0*m*a2*r/sigma * sin2) * sin2;

    // Set contravariant metric coefficients
    g_inv(I00,i) = -(1.0 + 2.0*m*r/sigma);
    g_inv(I01,i) = 2.0*m*r/sigma;
    g_inv(I11,i) = delta/sigma;
    g_inv(I13,i) = a/sigma;
    g_inv(I22,i) = 1.0/sigma;
    g_inv(I33,i) = 1.0 / (sigma * sin2);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: r-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face1Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  // Extract useful quantities that do not depend on r
  const Real &m = bh_mass_;
  const Real &a = bh_spin_;
  Real a2 = SQR(a);
  const Real &sin2 = metric_face1_j1_(j);
  const Real &cos2 = metric_face1_j2_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract remaining useful quantities
    const Real &r = metric_face1_i1_(i);
    Real r2 = SQR(r);
    Real delta = r2 - 2.0*m*r + a2;
    Real sigma = r2 + a2 * cos2;

    // Set covariant metric coefficients
    g(I00,i) = -(1.0 - 2.0*m*r/sigma);
    g(I01,i) = 2.0*m*r/sigma;
    g(I03,i) = -2.0*m*a*r/sigma * sin2;
    g(I11,i) = 1.0 + 2.0*m*r/sigma;
    g(I13,i) = -(1.0 + 2.0*m*r/sigma) * a * sin2;
    g(I22,i) = sigma;
    g(I33,i) = (r2 + a2 + 2.0*m*a2*r/sigma * sin2) * sin2;

    // Set contravariant metric coefficients
    g_inv(I00,i) = -(1.0 + 2.0*m*r/sigma);
    g_inv(I01,i) = 2.0*m*r/sigma;
    g_inv(I11,i) = delta/sigma;
    g_inv(I13,i) = a/sigma;
    g_inv(I22,i) = 1.0/sigma;
    g_inv(I33,i) = 1.0 / (sigma * sin2);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: theta-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face2Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  // Extract useful quantities that do not depend on r
  const Real &m = bh_mass_;
  const Real &a = bh_spin_;
  Real a2 = SQR(a);
  const Real &sin2 = metric_face2_j1_(j);
  const Real &cos2 = metric_face2_j2_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract remaining useful quantities
    const Real &r = metric_face2_i1_(i);
    Real r2 = SQR(r);
    Real delta = r2 - 2.0*m*r + a2;
    Real sigma = r2 + a2 * cos2;

    // Set covariant metric coefficients
    g(I00,i) = -(1.0 - 2.0*m*r/sigma);
    g(I01,i) = 2.0*m*r/sigma;
    g(I03,i) = -2.0*m*a*r/sigma * sin2;
    g(I11,i) = 1.0 + 2.0*m*r/sigma;
    g(I13,i) = -(1.0 + 2.0*m*r/sigma) * a * sin2;
    g(I22,i) = sigma;
    g(I33,i) = (r2 + a2 + 2.0*m*a2*r/sigma * sin2) * sin2;

    // Set contravariant metric coefficients
    g_inv(I00,i) = -(1.0 + 2.0*m*r/sigma);
    g_inv(I01,i) = 2.0*m*r/sigma;
    g_inv(I11,i) = delta/sigma;
    g_inv(I13,i) = a/sigma;
    g_inv(I22,i) = 1.0/sigma;
    g_inv(I33,i) = 1.0 / (sigma * sin2);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for computing face-centered metric coefficients: phi-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
void Coordinates::Face3Metric(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> &g, AthenaArray<Real> &g_inv)
{
  // Extract useful quantities that do not depend on r
  const Real &m = bh_mass_;
  const Real &a = bh_spin_;
  Real a2 = SQR(a);
  const Real &sin2 = metric_face3_j1_(j);
  const Real &cos2 = metric_face3_j2_(j);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract remaining useful quantities
    const Real &r = metric_face3_i1_(i);
    Real r2 = SQR(r);
    Real delta = r2 - 2.0*m*r + a2;
    Real sigma = r2 + a2 * cos2;

    // Set covariant metric coefficients
    g(I00,i) = -(1.0 - 2.0*m*r/sigma);
    g(I01,i) = 2.0*m*r/sigma;
    g(I03,i) = -2.0*m*a*r/sigma * sin2;
    g(I11,i) = 1.0 + 2.0*m*r/sigma;
    g(I13,i) = -(1.0 + 2.0*m*r/sigma) * a * sin2;
    g(I22,i) = sigma;
    g(I33,i) = (r2 + a2 + 2.0*m*a2*r/sigma * sin2) * sin2;

    // Set contravariant metric coefficients
    g_inv(I00,i) = -(1.0 + 2.0*m*r/sigma);
    g_inv(I01,i) = 2.0*m*r/sigma;
    g_inv(I11,i) = delta/sigma;
    g_inv(I13,i) = a/sigma;
    g_inv(I22,i) = 1.0/sigma;
    g_inv(I33,i) = 1.0 / (sigma * sin2);
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: r-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   bb1: 3D array of normal components B^1 of magnetic field, in global coordinates
//   prim_l: 1D array of left primitives, using global coordinates
//   prim_r: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
//   bbx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   expects \tilde{u}^1/\tilde{u}^2/\tilde{u}^3 in IVX/IVY/IVZ slots
//   expects B^1 in bb1
//   expects B^2/B^3 in IBY/IBZ slots
//   puts \tilde{u}^x/\tilde{u}^y/\tilde{u}^z in IVX/IVY/IVZ
//   puts B^x in bbx
//   puts B^y/B^z in IBY/IBZ slots
//   \tilde{u}^\hat{i} = u^\hat{i}
void Coordinates::PrimToLocal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb1, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &bbx)
{
  // Calculate metric coefficients
  if (MAGNETIC_FIELDS_ENABLED)
    Face1Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract transformation coefficients
    const Real &mt_0 = trans_face1_ji1_(j,i);
    const Real &mx_0 = trans_face1_ji5_(j,i);
    const Real &mx_1 = trans_face1_ji2_(j,i);
    const Real &my_2 = trans_face1_ji3_(j,i);
    const Real &mz_0 = trans_face1_ji6_(j,i);
    const Real &mz_1 = trans_face1_ji7_(j,i);
    const Real &mz_3 = trans_face1_ji4_(j,i);

    // Extract global projected 4-velocities
    Real uu0_l = 0.0;
    Real uu1_l = prim_l(IVX,i);
    Real uu2_l = prim_l(IVY,i);
    Real uu3_l = prim_l(IVZ,i);
    Real uu0_r = 0.0;
    Real uu1_r = prim_r(IVX,i);
    Real uu2_r = prim_r(IVY,i);
    Real uu3_r = prim_r(IVZ,i);

    // Transform projected 4-velocities
    Real uux_l = mx_0*uu0_l + mx_1*uu1_l;
    Real uuy_l = my_2*uu2_l;
    Real uuz_l = mz_0*uu0_l + mz_1*uu1_l + mz_3*uu3_l;
    Real uux_r = mx_0*uu0_r + mx_1*uu1_r;
    Real uuy_r = my_2*uu2_r;
    Real uuz_r = mz_0*uu0_r + mz_1*uu1_r + mz_3*uu3_r;

    // Set local projected 4-velocities
    prim_l(IVX,i) = uux_l;
    prim_l(IVY,i) = uuy_l;
    prim_l(IVZ,i) = uuz_l;
    prim_r(IVX,i) = uux_r;
    prim_r(IVY,i) = uuy_r;
    prim_r(IVZ,i) = uuz_r;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      // Extract metric coefficients
      const Real &g_00 = g_(I00,i);
      const Real &g_01 = g_(I01,i);
      const Real &g_02 = 0.0;
      const Real &g_03 = g_(I03,i);
      const Real &g_10 = g_(I01,i);
      const Real &g_11 = g_(I11,i);
      const Real &g_12 = 0.0;
      const Real &g_13 = g_(I13,i);
      const Real &g_20 = 0.0;
      const Real &g_21 = 0.0;
      const Real &g_22 = g_(I22,i);
      const Real &g_23 = 0.0;
      const Real &g_30 = g_(I03,i);
      const Real &g_31 = g_(I13,i);
      const Real &g_32 = 0.0;
      const Real &g_33 = g_(I33,i);
      const Real &g01 = gi_(I01,i);
      const Real &g02 = 0.0;
      const Real &g03 = 0.0;
      Real alpha = std::sqrt(-1.0/gi_(I00,i));

      // Calculate global 4-velocities
      Real tmp = g_11*uu1_l*uu1_l + 2.0*g_12*uu1_l*uu2_l + 2.0*g_13*uu1_l*uu3_l
               + g_22*uu2_l*uu2_l + 2.0*g_23*uu2_l*uu3_l
               + g_33*uu3_l*uu3_l;
      Real gamma_l = std::sqrt(1.0 + tmp);
      Real u0_l = gamma_l / alpha;
      Real u1_l = uu1_l - alpha * gamma_l * g01;
      Real u2_l = uu2_l - alpha * gamma_l * g02;
      Real u3_l = uu3_l - alpha * gamma_l * g03;
      tmp = g_11*uu1_r*uu1_r + 2.0*g_12*uu1_r*uu2_r + 2.0*g_13*uu1_r*uu3_r
          + g_22*uu2_r*uu2_r + 2.0*g_23*uu2_r*uu3_r
          + g_33*uu3_r*uu3_r;
      Real gamma_r = std::sqrt(1.0 + tmp);
      Real u0_r = gamma_r / alpha;
      Real u1_r = uu1_r - alpha * gamma_r * g01;
      Real u2_r = uu2_r - alpha * gamma_r * g02;
      Real u3_r = uu3_r - alpha * gamma_r * g03;

      // Extract global magnetic fields
      const Real &bb1_l = bb1(k,j,i);
      const Real &bb1_r = bb1(k,j,i);
      Real &bb2_l = prim_l(IBY,i);
      Real &bb3_l = prim_l(IBZ,i);
      Real &bb2_r = prim_r(IBY,i);
      Real &bb3_r = prim_r(IBZ,i);

      // Calculate global 4-magnetic fields
      Real b0_l = g_10*bb1_l*u0_l + g_11*bb1_l*u1_l + g_12*bb1_l*u2_l + g_13*bb1_l*u3_l
                + g_20*bb2_l*u0_l + g_21*bb2_l*u1_l + g_22*bb2_l*u2_l + g_23*bb2_l*u3_l
                + g_30*bb3_l*u0_l + g_31*bb3_l*u1_l + g_32*bb3_l*u2_l + g_33*bb3_l*u3_l;
      Real b1_l = (bb1_l + b0_l * u1_l) / u0_l;
      Real b2_l = (bb2_l + b0_l * u2_l) / u0_l;
      Real b3_l = (bb3_l + b0_l * u3_l) / u0_l;
      Real b0_r = g_10*bb1_r*u0_r + g_11*bb1_r*u1_r + g_12*bb1_r*u2_r + g_13*bb1_r*u3_r
                + g_20*bb2_r*u0_r + g_21*bb2_r*u1_r + g_22*bb2_r*u2_r + g_23*bb2_r*u3_r
                + g_30*bb3_r*u0_r + g_31*bb3_r*u1_r + g_32*bb3_r*u2_r + g_33*bb3_r*u3_r;
      Real b1_r = (bb1_r + b0_r * u1_r) / u0_r;
      Real b2_r = (bb2_r + b0_r * u2_r) / u0_r;
      Real b3_r = (bb3_r + b0_r * u3_r) / u0_r;

      // Transform 4-velocities
      Real ut_l = mt_0*u0_l;
      Real ux_l = mx_0*u0_l + mx_1*u1_l;
      Real uy_l = my_2*u2_l;
      Real uz_l = mz_0*u0_l + mz_1*u1_l + mz_3*u3_l;
      Real ut_r = mt_0*u0_r;
      Real ux_r = mx_0*u0_r + mx_1*u1_r;
      Real uy_r = my_2*u2_r;
      Real uz_r = mz_0*u0_r + mz_1*u1_r + mz_3*u3_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_0*b0_l + mx_1*b1_l;
      Real by_l = my_2*b2_l;
      Real bz_l = mz_0*b0_l + mz_1*b1_l + mz_3*b3_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_0*b0_r + mx_1*b1_r;
      Real by_r = my_2*b2_r;
      Real bz_r = mz_0*b0_r + mz_1*b1_r + mz_3*b3_r;

      // Set local magnetic fields
      Real bbx_l = ut_l * bx_l - ux_l * bt_l;
      Real bbx_r = ut_r * bx_r - ux_r * bt_r;
      bbx(i) = 0.5 * (bbx_l + bbx_r);
      bb2_l = ut_l * by_l - uy_l * bt_l;
      bb3_l = ut_l * bz_l - uz_l * bt_l;
      bb2_r = ut_r * by_r - uy_r * bt_r;
      bb3_r = ut_r * bz_r - uz_r * bt_r;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: theta-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   bb2: 3D array of normal components B^2 of magnetic field, in global coordinates
//   prim_l: 1D array of left primitives, using global coordinates
//   prim_r: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
//   bbx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   expects \tilde{u}^1/\tilde{u}^2/\tilde{u}^3 in IVX/IVY/IVZ slots
//   expects B^2 in bb2
//   expects B^3/B^1 in IBY/IBZ slots
//   puts \tilde{u}^x/\tilde{u}^y/\tilde{u}^z in IVY/IVZ/IVX slots
//   puts B^x in bbx
//   puts B^y/B^z in IBY/IBZ slots
//   \tilde{u}^\hat{i} = u^\hat{i}
void Coordinates::PrimToLocal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb2, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &bbx)
{
  // Calculate metric coefficients
  if (MAGNETIC_FIELDS_ENABLED)
    Face2Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract transformation coefficients
    const Real &mt_0 = trans_face2_ji1_(j,i);
    const Real &mx_2 = trans_face2_ji3_(j,i);
    const Real &my_3 = trans_face2_ji4_(j,i);
    const Real &mz_0 = trans_face2_ji5_(j,i);
    const Real &mz_1 = trans_face2_ji2_(j,i);
    const Real &mz_3 = trans_face2_ji6_(j,i);

    // Extract global projected 4-velocities
    Real uu0_l = 0.0;
    Real uu1_l = prim_l(IVX,i);
    Real uu2_l = prim_l(IVY,i);
    Real uu3_l = prim_l(IVZ,i);
    Real uu0_r = 0.0;
    Real uu1_r = prim_r(IVX,i);
    Real uu2_r = prim_r(IVY,i);
    Real uu3_r = prim_r(IVZ,i);

    // Transform projected 4-velocities
    Real uux_l = mx_2*uu2_l;
    Real uuy_l = my_3*uu3_l;
    Real uuz_l = mz_0*uu0_l + mz_1*uu1_l + mz_3*uu3_l;
    Real uux_r = mx_2*uu2_r;
    Real uuy_r = my_3*uu3_r;
    Real uuz_r = mz_0*uu0_r + mz_1*uu1_r + mz_3*uu3_r;

    // Set local projected 4-velocities
    prim_l(IVY,i) = uux_l;
    prim_l(IVZ,i) = uuy_l;
    prim_l(IVX,i) = uuz_l;
    prim_r(IVY,i) = uux_r;
    prim_r(IVZ,i) = uuy_r;
    prim_r(IVX,i) = uuz_r;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      // Extract metric coefficients
      const Real &g_00 = g_(I00,i);
      const Real &g_01 = g_(I01,i);
      const Real &g_02 = 0.0;
      const Real &g_03 = g_(I03,i);
      const Real &g_10 = g_(I01,i);
      const Real &g_11 = g_(I11,i);
      const Real &g_12 = 0.0;
      const Real &g_13 = g_(I13,i);
      const Real &g_20 = 0.0;
      const Real &g_21 = 0.0;
      const Real &g_22 = g_(I22,i);
      const Real &g_23 = 0.0;
      const Real &g_30 = g_(I03,i);
      const Real &g_31 = g_(I13,i);
      const Real &g_32 = 0.0;
      const Real &g_33 = g_(I33,i);
      const Real &g01 = gi_(I01,i);
      const Real &g02 = 0.0;
      const Real &g03 = 0.0;
      Real alpha = std::sqrt(-1.0/gi_(I00,i));

      // Calculate global 4-velocities
      Real tmp = g_11*uu1_l*uu1_l + 2.0*g_12*uu1_l*uu2_l
                   + 2.0*g_13*uu1_l*uu3_l
               + g_22*uu2_l*uu2_l + 2.0*g_23*uu2_l*uu3_l
               + g_33*uu3_l*uu3_l;
      Real gamma_l = std::sqrt(1.0 + tmp);
      Real u0_l = gamma_l / alpha;
      Real u1_l = uu1_l - alpha * gamma_l * g01;
      Real u2_l = uu2_l - alpha * gamma_l * g02;
      Real u3_l = uu3_l - alpha * gamma_l * g03;
      tmp = g_11*uu1_r*uu1_r + 2.0*g_12*uu1_r*uu2_r + 2.0*g_13*uu1_r*uu3_r
          + g_22*uu2_r*uu2_r + 2.0*g_23*uu2_r*uu3_r
          + g_33*uu3_r*uu3_r;
      Real gamma_r = std::sqrt(1.0 + tmp);
      Real u0_r = gamma_r / alpha;
      Real u1_r = uu1_r - alpha * gamma_r * g01;
      Real u2_r = uu2_r - alpha * gamma_r * g02;
      Real u3_r = uu3_r - alpha * gamma_r * g03;

      // Extract global magnetic fields
      const Real &bb2_l = bb2(k,j,i);
      const Real &bb2_r = bb2(k,j,i);
      Real &bb3_l = prim_l(IBY,i);
      Real &bb1_l = prim_l(IBZ,i);
      Real &bb3_r = prim_r(IBY,i);
      Real &bb1_r = prim_r(IBZ,i);

      // Calculate global 4-magnetic fields
      Real b0_l = g_10*bb1_l*u0_l + g_11*bb1_l*u1_l + g_12*bb1_l*u2_l + g_13*bb1_l*u3_l
                + g_20*bb2_l*u0_l + g_21*bb2_l*u1_l + g_22*bb2_l*u2_l + g_23*bb2_l*u3_l
                + g_30*bb3_l*u0_l + g_31*bb3_l*u1_l + g_32*bb3_l*u2_l + g_33*bb3_l*u3_l;
      Real b1_l = (bb1_l + b0_l * u1_l) / u0_l;
      Real b2_l = (bb2_l + b0_l * u2_l) / u0_l;
      Real b3_l = (bb3_l + b0_l * u3_l) / u0_l;
      Real b0_r = g_10*bb1_r*u0_r + g_11*bb1_r*u1_r + g_12*bb1_r*u2_r + g_13*bb1_r*u3_r
                + g_20*bb2_r*u0_r + g_21*bb2_r*u1_r + g_22*bb2_r*u2_r + g_23*bb2_r*u3_r
                + g_30*bb3_r*u0_r + g_31*bb3_r*u1_r + g_32*bb3_r*u2_r + g_33*bb3_r*u3_r;
      Real b1_r = (bb1_r + b0_r * u1_r) / u0_r;
      Real b2_r = (bb2_r + b0_r * u2_r) / u0_r;
      Real b3_r = (bb3_r + b0_r * u3_r) / u0_r;

      // Transform 4-velocities
      Real ut_l = mt_0*u0_l;
      Real ux_l = mx_2*u2_l;
      Real uy_l = my_3*u3_l;
      Real uz_l = mz_0*u0_l + mz_1*u1_l + mz_3*u3_l;
      Real ut_r = mt_0*u0_r;
      Real ux_r = mx_2*u2_r;
      Real uy_r = my_3*u3_r;
      Real uz_r = mz_0*u0_r + mz_1*u1_r + mz_3*u3_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_2*b2_l;
      Real by_l = my_3*b3_l;
      Real bz_l = mz_0*b0_l + mz_1*b1_l + mz_3*b3_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_2*b2_r;
      Real by_r = my_3*b3_r;
      Real bz_r = mz_0*b0_r + mz_1*b1_r + mz_3*b3_r;

      // Set local magnetic fields
      Real bbx_l = ut_l * bx_l - ux_l * bt_l;
      Real bbx_r = ut_r * bx_r - ux_r * bt_r;
      bbx(i) = 0.5 * (bbx_l + bbx_r);
      bb3_l = ut_l * by_l - uy_l * bt_l;
      bb1_l = ut_l * bz_l - uz_l * bt_l;
      bb3_r = ut_r * by_r - uy_r * bt_r;
      bb1_r = ut_r * bz_r - uz_r * bt_r;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming primitives to locally flat frame: phi-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   bb3: 3D array of normal components B^3 of magnetic field, in global coordinates
//   prim_l: 1D array of left primitives, using global coordinates
//   prim_r: 1D array of right primitives, using global coordinates
// Outputs:
//   prim_l: values overwritten in local coordinates
//   prim_r: values overwritten in local coordinates
//   bbx: 1D array of normal magnetic fields, in local coordinates
// Notes:
//   expects \tilde{u}^1/\tilde{u}^2/\tilde{u}^3 in IVX/IVY/IVZ slots
//   expects B^3 in bb3
//   expects B^1/B^2 in IBY/IBZ slots
//   puts \tilde{u}^x/\tilde{u}^y/\tilde{u}^z in IVZ/IVX/IVY slots
//   puts B^x in bbx
//   puts B^y/B^z in IBY/IBZ slots
//   \tilde{u}^\hat{i} = u^\hat{i}
void Coordinates::PrimToLocal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &bb3, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &bbx)
{
  // Calculate metric coefficients
  if (MAGNETIC_FIELDS_ENABLED)
    Face3Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract transformation coefficients
    const Real &mt_0 = trans_face3_ji1_(j,i);
    const Real &mx_3 = trans_face3_ji4_(j,i);
    const Real &my_0 = trans_face3_ji5_(j,i);
    const Real &my_1 = trans_face3_ji2_(j,i);
    const Real &my_3 = trans_face3_ji6_(j,i);
    const Real &mz_2 = trans_face3_ji3_(j,i);

    // Extract global projected 4-velocities
    Real uu0_l = 0.0;
    Real uu1_l = prim_l(IVX,i);
    Real uu2_l = prim_l(IVY,i);
    Real uu3_l = prim_l(IVZ,i);
    Real uu0_r = 0.0;
    Real uu1_r = prim_r(IVX,i);
    Real uu2_r = prim_r(IVY,i);
    Real uu3_r = prim_r(IVZ,i);

    // Transform projected 4-velocities
    Real uux_l = mx_3*uu3_l;
    Real uuy_l = my_0*uu0_l + my_1*uu1_l + my_3*uu3_l;
    Real uuz_l = mz_2*uu2_l;
    Real uux_r = mx_3*uu3_r;
    Real uuy_r = my_0*uu0_r + my_1*uu1_r + my_3*uu3_r;
    Real uuz_r = mz_2*uu2_r;

    // Set local projected 4-velocities
    prim_l(IVZ,i) = uux_l;
    prim_l(IVX,i) = uuy_l;
    prim_l(IVY,i) = uuz_l;
    prim_r(IVZ,i) = uux_r;
    prim_r(IVX,i) = uuy_r;
    prim_r(IVY,i) = uuz_r;

    // Transform magnetic field if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      // Extract metric coefficients
      const Real &g_00 = g_(I00,i);
      const Real &g_01 = g_(I01,i);
      const Real &g_02 = 0.0;
      const Real &g_03 = g_(I03,i);
      const Real &g_10 = g_(I01,i);
      const Real &g_11 = g_(I11,i);
      const Real &g_12 = 0.0;
      const Real &g_13 = g_(I13,i);
      const Real &g_20 = 0.0;
      const Real &g_21 = 0.0;
      const Real &g_22 = g_(I22,i);
      const Real &g_23 = 0.0;
      const Real &g_30 = g_(I03,i);
      const Real &g_31 = g_(I13,i);
      const Real &g_32 = 0.0;
      const Real &g_33 = g_(I33,i);
      const Real &g01 = gi_(I01,i);
      const Real &g02 = 0.0;
      const Real &g03 = 0.0;
      Real alpha = std::sqrt(-1.0/gi_(I00,i));

      // Calculate global 4-velocities
      Real tmp = g_11*uu1_l*uu1_l + 2.0*g_12*uu1_l*uu2_l
                   + 2.0*g_13*uu1_l*uu3_l
               + g_22*uu2_l*uu2_l + 2.0*g_23*uu2_l*uu3_l
               + g_33*uu3_l*uu3_l;
      Real gamma_l = std::sqrt(1.0 + tmp);
      Real u0_l = gamma_l / alpha;
      Real u1_l = uu1_l - alpha * gamma_l * g01;
      Real u2_l = uu2_l - alpha * gamma_l * g02;
      Real u3_l = uu3_l - alpha * gamma_l * g03;
      tmp = g_11*uu1_r*uu1_r + 2.0*g_12*uu1_r*uu2_r + 2.0*g_13*uu1_r*uu3_r
          + g_22*uu2_r*uu2_r + 2.0*g_23*uu2_r*uu3_r
          + g_33*uu3_r*uu3_r;
      Real gamma_r = std::sqrt(1.0 + tmp);
      Real u0_r = gamma_r / alpha;
      Real u1_r = uu1_r - alpha * gamma_r * g01;
      Real u2_r = uu2_r - alpha * gamma_r * g02;
      Real u3_r = uu3_r - alpha * gamma_r * g03;

      // Extract global magnetic fields
      const Real &bb3_l = bb3(k,j,i);
      const Real &bb3_r = bb3(k,j,i);
      Real &bb1_l = prim_l(IBY,i);
      Real &bb2_l = prim_l(IBZ,i);
      Real &bb1_r = prim_r(IBY,i);
      Real &bb2_r = prim_r(IBZ,i);

      // Calculate global 4-magnetic fields
      Real b0_l = g_10*bb1_l*u0_l + g_11*bb1_l*u1_l + g_12*bb1_l*u2_l + g_13*bb1_l*u3_l
                + g_20*bb2_l*u0_l + g_21*bb2_l*u1_l + g_22*bb2_l*u2_l + g_23*bb2_l*u3_l
                + g_30*bb3_l*u0_l + g_31*bb3_l*u1_l + g_32*bb3_l*u2_l + g_33*bb3_l*u3_l;
      Real b1_l = (bb1_l + b0_l * u1_l) / u0_l;
      Real b2_l = (bb2_l + b0_l * u2_l) / u0_l;
      Real b3_l = (bb3_l + b0_l * u3_l) / u0_l;
      Real b0_r = g_10*bb1_r*u0_r + g_11*bb1_r*u1_r + g_12*bb1_r*u2_r + g_13*bb1_r*u3_r
                + g_20*bb2_r*u0_r + g_21*bb2_r*u1_r + g_22*bb2_r*u2_r + g_23*bb2_r*u3_r
                + g_30*bb3_r*u0_r + g_31*bb3_r*u1_r + g_32*bb3_r*u2_r + g_33*bb3_r*u3_r;
      Real b1_r = (bb1_r + b0_r * u1_r) / u0_r;
      Real b2_r = (bb2_r + b0_r * u2_r) / u0_r;
      Real b3_r = (bb3_r + b0_r * u3_r) / u0_r;

      // Transform 4-velocities
      Real ut_l = mt_0*u0_l;
      Real ux_l = mx_3*u3_l;
      Real uy_l = my_0*u0_l + my_1*u1_l + my_3*u3_l;
      Real uz_l = mz_2*u2_l;
      Real ut_r = mt_0*u0_r;
      Real ux_r = mx_3*u3_r;
      Real uy_r = my_0*u0_r + my_1*u1_r + my_3*u3_r;
      Real uz_r = mz_2*u2_r;

      // Transform 4-magnetic fields
      Real bt_l = mt_0*b0_l;
      Real bx_l = mx_3*b3_l;
      Real by_l = my_0*b0_l + my_1*b1_l + my_3*b3_l;
      Real bz_l = mz_2*b2_l;
      Real bt_r = mt_0*b0_r;
      Real bx_r = mx_3*b3_r;
      Real by_r = my_0*b0_r + my_1*b1_r + my_3*b3_r;
      Real bz_r = mz_2*b2_r;

      // Set local magnetic fields
      Real bbx_l = ut_l * bx_l - ux_l * bt_l;
      Real bbx_r = ut_r * bx_r - ux_r * bt_r;
      bbx(i) = 0.5 * (bbx_l + bbx_r);
      bb1_l = ut_l * by_l - uy_l * bt_l;
      bb2_l = ut_l * bz_l - uz_l * bt_l;
      bb1_r = ut_r * by_r - uy_r * bt_r;
      bb2_r = ut_r * bz_r - uz_r * bt_r;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: r-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM1/IM2/IM3 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts r-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts r-fluxes of B2/B3 in IBY/IBZ slots
void Coordinates::FluxToGlobal1(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux)
{
  // Calculate metric coefficients
  Face1Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract transformation coefficients
    const Real &mt_0 = trans_face1_ji1_(j,i);
    const Real &mx_0 = trans_face1_ji5_(j,i);
    const Real &mx_1 = trans_face1_ji2_(j,i);
    const Real &my_2 = trans_face1_ji3_(j,i);
    const Real &mz_1 = trans_face1_ji7_(j,i);
    const Real &mz_3 = trans_face1_ji4_(j,i);
    Real m0_t = 1.0/mt_0;
    Real m1_t = -mx_0/(mt_0*mx_1);
    Real m1_x = 1.0/mx_1;
    Real m2_y = 1.0/my_2;
    Real m3_x = -mz_1/(mx_1*mz_3);
    Real m3_z = 1.0/mz_3;

    // Extract local conserved quantities and fluxes
    Real dt = cons(IDN,i);
    Real ttt = cons(IEN,i);
    Real ttx = cons(IM1,i);
    Real tty = cons(IM2,i);
    Real ttz = cons(IM3,i);
    Real dx = flux(IDN,i);
    Real txt = flux(IEN,i);
    Real txx = flux(IM1,i);
    Real txy = flux(IM2,i);
    Real txz = flux(IM3,i);

    // Transform stress-energy tensor
    Real t10 = m1_t*m0_t*ttt + m1_x*m0_t*txt;
    Real t11 = m1_t*m1_t*ttt + m1_t*m1_x*ttx + m1_x*m1_t*txt + m1_x*m1_x*txx;
    Real t12 = m1_t*m2_y*tty + m1_x*m2_y*txy;
    Real t13 = m1_t*m3_x*ttx + m1_t*m3_z*ttz + m1_x*m3_x*txx + m1_x*m3_z*txz;

    // Extract metric coefficients
    const Real &g_00 = g_(I00,i);
    const Real &g_01 = g_(I01,i);
    const Real &g_03 = g_(I03,i);
    const Real &g_10 = g_(I01,i);
    const Real &g_11 = g_(I11,i);
    const Real &g_13 = g_(I13,i);
    const Real &g_22 = g_(I22,i);
    const Real &g_30 = g_(I03,i);
    const Real &g_31 = g_(I13,i);
    const Real &g_33 = g_(I33,i);

    // Extract global fluxes
    Real &d1 = flux(IDN,i);
    Real &t1_0 = flux(IEN,i);
    Real &t1_1 = flux(IM1,i);
    Real &t1_2 = flux(IM2,i);
    Real &t1_3 = flux(IM3,i);

    // Set fluxes
    d1 = m1_t*dt + m1_x*dx;
    t1_0 = g_00*t10 + g_01*t11 + g_03*t13;
    t1_1 = g_10*t10 + g_11*t11 + g_13*t13;
    t1_2 = g_22*t12;
    t1_3 = g_30*t10 + g_31*t11 + g_33*t13;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      const Real &fxt = bbx(i);
      const Real &fyt = cons(IBY,i);
      const Real &fzt = cons(IBZ,i);
      Real fyx = flux(IBY,i);
      Real fzx = flux(IBZ,i);
      Real &f21 = flux(IBY,i);
      Real &f31 = flux(IBZ,i);
      f21 = m2_y*m1_t*fyt + m2_y*m1_x*fyx;
      f31 = m3_x*m1_t*fxt + m3_z*m1_t*fzt + m3_z*m1_x*fzx;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: theta-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM2/IM3/IM1 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts theta-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts theta-fluxes of B3/B1 in IBY/IBZ slots
void Coordinates::FluxToGlobal2(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux)
{
  // Calculate metric coefficients
  Face2Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract transformation coefficients
    const Real &mt_0 = trans_face2_ji1_(j,i);
    const Real &mx_2 = trans_face2_ji3_(j,i);
    const Real &my_3 = trans_face2_ji4_(j,i);
    const Real &mz_0 = trans_face2_ji5_(j,i);
    const Real &mz_1 = trans_face2_ji2_(j,i);
    const Real &mz_3 = trans_face2_ji6_(j,i);
    Real m0_t = 1.0/mt_0;
    Real m1_t = -mz_0/(mt_0*mz_1);
    Real m1_y = -mz_3/(mz_1*my_3);
    Real m1_z = 1.0/mz_1;
    Real m2_x = 1.0/mx_2;
    Real m3_y = 1.0/my_3;

    // Extract local conserved quantities and fluxes
    Real dt = cons(IDN,i);
    Real ttt = cons(IEN,i);
    Real ttx = cons(IM1,i);
    Real tty = cons(IM2,i);
    Real ttz = cons(IM3,i);
    Real dx = flux(IDN,i);
    Real txt = flux(IEN,i);
    Real txx = flux(IM2,i);
    Real txy = flux(IM3,i);
    Real txz = flux(IM1,i);

    // Transform stress-energy tensor
    Real t20 = m2_x*m0_t*txt;
    Real t21 = m2_x*m1_t*txt + m2_x*m1_y*txy + m2_x*m1_z*txz;
    Real t22 = m2_x*m2_x*txx;
    Real t23 = m2_x*m3_y*txy;

    // Extract metric coefficients
    const Real &g_00 = g_(I00,i);
    const Real &g_01 = g_(I01,i);
    const Real &g_03 = g_(I03,i);
    const Real &g_10 = g_(I01,i);
    const Real &g_11 = g_(I11,i);
    const Real &g_13 = g_(I13,i);
    const Real &g_22 = g_(I22,i);
    const Real &g_30 = g_(I03,i);
    const Real &g_31 = g_(I13,i);
    const Real &g_33 = g_(I33,i);

    // Extract global fluxes
    Real &d2 = flux(IDN,i);
    Real &t2_0 = flux(IEN,i);
    Real &t2_1 = flux(IM1,i);
    Real &t2_2 = flux(IM2,i);
    Real &t2_3 = flux(IM3,i);

    // Set fluxes
    d2 = m2_x*dx;
    t2_0 = g_00*t20 + g_01*t21 + g_03*t23;
    t2_1 = g_10*t20 + g_11*t21 + g_13*t23;
    t2_2 = g_22*t22;
    t2_3 = g_30*t20 + g_31*t21 + g_33*t23;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      Real ftx = -bbx(i);
      Real fyx = flux(IBY,i);
      Real fzx = flux(IBZ,i);
      Real &f32 = flux(IBY,i);
      Real &f12 = flux(IBZ,i);
      f32 = m3_y*m2_x*fyx;
      f12 = m1_t*m2_x*ftx + m1_y*m2_x*fyx + m1_z*m2_x*fzx;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming fluxes to global frame: phi-interface
// Inputs:
//   k,j: phi- and theta-indices
//   il,iu: r-index bounds
//   cons: array of conserved quantities in 1D, using local coordinates
//   bbx: 1D array of longitudinal magnetic fields, in local coordinates
//   flux: array of fluxes in 1D, using local coordinates
// Outputs:
//   flux: values overwritten in global coordinates
// Notes:
//   expects values and x-fluxes of Mx/My/Mz in IM3/IM1/IM2 slots
//   expects values and x-fluxes of By/Bz in IBY/IBZ slots
//   puts phi-fluxes of M1/M2/M3 in IM1/IM2/IM3 slots
//   puts phi-fluxes of B1/B2 in IBY/IBZ slots
void Coordinates::FluxToGlobal3(const int k, const int j, const int il, const int iu,
    const AthenaArray<Real> &cons, const AthenaArray<Real> &bbx,
    AthenaArray<Real> &flux)
{
  // Calculate metric coefficients
  Face3Metric(k, j, il, iu, g_, gi_);

  // Go through 1D block of cells
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract transformation coefficients
    const Real &mt_0 = trans_face3_ji1_(j,i);
    const Real &mx_3 = trans_face3_ji4_(j,i);
    const Real &my_0 = trans_face3_ji5_(j,i);
    const Real &my_1 = trans_face3_ji2_(j,i);
    const Real &my_3 = trans_face3_ji6_(j,i);
    const Real &mz_2 = trans_face3_ji3_(j,i);
    Real m0_t = 1.0/mt_0;
    Real m1_t = -my_0/(mt_0*my_1);
    Real m1_x = -my_3/(my_1*mx_3);
    Real m1_y = 1.0/my_1;
    Real m2_z = 1.0/mz_2;
    Real m3_x = 1.0/mx_3;

    // Extract local conserved quantities and fluxes
    Real dt = cons(IDN,i);
    Real ttt = cons(IEN,i);
    Real ttx = cons(IM1,i);
    Real tty = cons(IM2,i);
    Real ttz = cons(IM3,i);
    Real dx = flux(IDN,i);
    Real txt = flux(IEN,i);
    Real txx = flux(IM3,i);
    Real txy = flux(IM1,i);
    Real txz = flux(IM2,i);

    // Transform stress-energy tensor
    Real t30 = m3_x*m0_t*txt;
    Real t31 = m3_x*m1_t*txt + m3_x*m1_x*txx + m3_x*m1_y*txy;
    Real t32 = m3_x*m2_z*txz;
    Real t33 = m3_x*m3_x*txx;

    // Extract metric coefficients
    const Real &g_00 = g_(I00,i);
    const Real &g_01 = g_(I01,i);
    const Real &g_03 = g_(I03,i);
    const Real &g_10 = g_(I01,i);
    const Real &g_11 = g_(I11,i);
    const Real &g_13 = g_(I13,i);
    const Real &g_22 = g_(I22,i);
    const Real &g_30 = g_(I03,i);
    const Real &g_31 = g_(I13,i);
    const Real &g_33 = g_(I33,i);

    // Extract global fluxes
    Real &d3 = flux(IDN,i);
    Real &t3_0 = flux(IEN,i);
    Real &t3_1 = flux(IM1,i);
    Real &t3_2 = flux(IM2,i);
    Real &t3_3 = flux(IM3,i);

    // Set fluxes
    d3 = m3_x*dx;
    t3_0 = g_00*t30 + g_01*t31 + g_03*t33;
    t3_1 = g_10*t30 + g_11*t31 + g_13*t33;
    t3_2 = g_22*t32;
    t3_3 = g_30*t30 + g_31*t31 + g_33*t33;

    // Transform magnetic fluxes if necessary
    if (MAGNETIC_FIELDS_ENABLED)
    {
      Real ftx = -bbx(i);
      Real fyx = flux(IBY,i);
      Real fzx = flux(IBZ,i);
      Real &f13 = flux(IBY,i);
      Real &f23 = flux(IBZ,i);
      f13 = m1_t*m3_x*ftx + m1_y*m3_x*fyx;
      f23 = m2_z*m3_x*fzx;
    }
  }
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Boyer-Lindquist to Kerr-Schild
// Inputs:
//   a0_bl,a1_bl,a2_bl,a3_bl: upper 4-vector components in Boyer-Lindquist coordinates
//   k,j,i: indices of cell in which transformation is desired
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in Kerr-Schild coordinates
void Coordinates::TransformVectorCell(
    Real a0_bl, Real a1_bl, Real a2_bl, Real a3_bl, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  const Real &m = bh_mass_;
  const Real &a = bh_spin_;
  const Real &r = x1v(i);
  Real delta = SQR(r) - 2.0*m*r + SQR(a);
  *pa0 = a0_bl + 2.0*m*r/delta * a1_bl;
  *pa1 = a1_bl;
  *pa2 = a2_bl;
  *pa3 = a3_bl + a/delta * a1_bl;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Boyer-Lindquist to Kerr-Schild
// Inputs:
//   a0_bl,a1_bl,a2_bl,a3_bl: upper 4-vector components in Boyer-Lindquist coordinates
//   k,j,i: indices of x1-face in which transformation is desired
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in Kerr-Schild coordinates
void Coordinates::TransformVectorFace1(
    Real a0_bl, Real a1_bl, Real a2_bl, Real a3_bl, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  const Real &m = bh_mass_;
  const Real &a = bh_spin_;
  const Real &r = x1f(i);
  Real delta = SQR(r) - 2.0*m*r + SQR(a);
  *pa0 = a0_bl + 2.0*m*r/delta * a1_bl;
  *pa1 = a1_bl;
  *pa2 = a2_bl;
  *pa3 = a3_bl + a/delta * a1_bl;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Boyer-Lindquist to Kerr-Schild
// Inputs:
//   a0_bl,a1_bl,a2_bl,a3_bl: upper 4-vector components in Boyer-Lindquist coordinates
//   k,j,i: indices of x2-face in which transformation is desired
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in Kerr-Schild coordinates
void Coordinates::TransformVectorFace2(
    Real a0_bl, Real a1_bl, Real a2_bl, Real a3_bl, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  const Real &m = bh_mass_;
  const Real &a = bh_spin_;
  const Real &r = x1v(i);
  Real delta = SQR(r) - 2.0*m*r + SQR(a);
  *pa0 = a0_bl + 2.0*m*r/delta * a1_bl;
  *pa1 = a1_bl;
  *pa2 = a2_bl;
  *pa3 = a3_bl + a/delta * a1_bl;
  return;
}

//--------------------------------------------------------------------------------------

// Function for transforming 4-vector from Boyer-Lindquist to Kerr-Schild
// Inputs:
//   a0_bl,a1_bl,a2_bl,a3_bl: upper 4-vector components in Boyer-Lindquist coordinates
//   k,j,i: indices of x3-face in which transformation is desired
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to upper 4-vector components in Kerr-Schild coordinates
void Coordinates::TransformVectorFace3(
    Real a0_bl, Real a1_bl, Real a2_bl, Real a3_bl, int k, int j, int i,
    Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  const Real &m = bh_mass_;
  const Real &a = bh_spin_;
  const Real &r = x1v(i);
  Real delta = SQR(r) - 2.0*m*r + SQR(a);
  *pa0 = a0_bl + 2.0*m*r/delta * a1_bl;
  *pa1 = a1_bl;
  *pa2 = a2_bl;
  *pa3 = a3_bl + a/delta * a1_bl;
  return;
}

//--------------------------------------------------------------------------------------

// Function for raising covariant components of a vector
// Inputs:
//   a_0,a_1,a_2,a_3: covariant components of vector
//   k,j,i: indices of cell in which transformation is desired
// Outputs:
//   pa0,pa1,pa2,pa3: pointers to contravariant 4-vector components
void Coordinates::RaiseVectorCell(Real a_0, Real a_1, Real a_2, Real a_3, int k, int j,
    int i, Real *pa0, Real *pa1, Real *pa2, Real *pa3)
{
  // Extract geometric quantities
  const Real &m = bh_mass_;
  const Real &a = bh_spin_;
  Real asq = SQR(a);
  const Real &r = metric_cell_i1_(i);
  Real rsq = SQR(r);
  const Real &sin2 = metric_cell_j1_(j);
  const Real &cos2 = metric_cell_j2_(j);
  Real delta = rsq - 2.0*m*r + asq;
  Real sigma = rsq + asq * cos2;

  // Calculate metric coefficients
  Real g00 = -(1.0 + 2.0*m*r/sigma);
  Real g01 = 2.0*m*r/sigma;
  Real g02 = 0.0;
  Real g03 = 0.0;
  Real g11 = delta/sigma;
  Real g12 = 0.0;
  Real g13 = a/sigma;
  Real g22 = 1.0/sigma;
  Real g23 = 0.0;
  Real g33 = 1.0/(sigma*sin2);
  const Real &g10 = g01;
  const Real &g20 = g02;
  const Real &g21 = g12;
  const Real &g30 = g03;
  const Real &g31 = g13;
  const Real &g32 = g23;

  // Set raised components
  *pa0 = g00*a_0 + g01*a_1 + g02*a_2 + g03*a_3;
  *pa1 = g10*a_0 + g11*a_1 + g12*a_2 + g13*a_3;
  *pa2 = g20*a_0 + g21*a_1 + g22*a_2 + g23*a_3;
  *pa3 = g30*a_0 + g31*a_1 + g32*a_2 + g33*a_3;
  return;
}

//--------------------------------------------------------------------------------------

// Function for lowering contravariant components of a vector
// Inputs:
//   a0,a1,a2,a3: contravariant components of vector
//   k,j,i: indices of cell in which transformation is desired
// Outputs:
//   pa_0,pa_1,pa_2,pa_3: pointers to covariant 4-vector components
void Coordinates::LowerVectorCell(Real a0, Real a1, Real a2, Real a3, int k, int j,
    int i, Real *pa_0, Real *pa_1, Real *pa_2, Real *pa_3)
{
  // Extract geometric quantities
  const Real &m = bh_mass_;
  const Real &a = bh_spin_;
  Real asq = SQR(a);
  const Real &r = metric_cell_i1_(i);
  Real rsq = SQR(r);
  const Real &sin2 = metric_cell_j1_(j);
  const Real &cos2 = metric_cell_j2_(j);
  Real sigma = rsq + asq * cos2;

  // Calculate metric coefficients
  Real g_00 = -(1.0 - 2.0*m*r/sigma);
  Real g_01 = 2.0*m*r/sigma;
  Real g_02 = 0.0;
  Real g_03 = -2.0*m*a*r/sigma * sin2;
  Real g_11 = 1.0 + 2.0*m*r/sigma;
  Real g_12 = 0.0;
  Real g_13 = -(1.0 + 2.0*m*r/sigma) * a * sin2;
  Real g_22 = sigma;
  Real g_23 = 0.0;
  Real g_33 = (rsq + asq + 2.0*m*asq*r/sigma * sin2) * sin2;
  const Real &g_10 = g_01;
  const Real &g_20 = g_02;
  const Real &g_21 = g_12;
  const Real &g_30 = g_03;
  const Real &g_31 = g_13;
  const Real &g_32 = g_23;

  // Set lowered components
  *pa_0 = g_00*a0 + g_01*a1 + g_02*a2 + g_03*a3;
  *pa_1 = g_10*a0 + g_11*a1 + g_12*a2 + g_13*a3;
  *pa_2 = g_20*a0 + g_21*a1 + g_22*a2 + g_23*a3;
  *pa_3 = g_30*a0 + g_31*a1 + g_32*a2 + g_33*a3;
  return;
}

//--------------------------------------------------------------------------------------

// Function for returning Boyer-Lindquist coordinates of given cell
// Inputs:
//   x1,x2,x3: Kerr-Schild coordinates to be converted
// Outputs:
//   pr: pointer to stored value of r
//   ptheta: pointer to stored value of theta
//   pphi: pointer to stored value of phi
// Notes:
//   Kerr-Schild (x1,x2) match Boyer-Lindquist (r,theta)
//   returns x3 = phi
//       not correct
//       not a problem if used for initializing azimuthally symmetric problem
//   ignores x0/t
//       surfaces of constant x0 are not surfaces of constant t
//       not a problem if used for initializing stationary problem
void Coordinates::GetBoyerLindquistCoordinates(Real x1, Real x2, Real x3,
    Real *pr, Real *ptheta, Real *pphi)
{
  *pr = x1;
  *ptheta = x2;
  *pphi = x3;
  return;
}
