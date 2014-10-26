// Minkowski spacetime, Cartesian coordinates
// Notes:
//   coordinates: t, x, y, z
//   metric: ds^2 = -dt^2 + dx^2 + dy^2 + dz^2

// Primary header
#include "coordinates.hpp"

// Athena headers
#include "../athena.hpp"         // enums, macros, Real
#include "../athena_arrays.hpp"  // AthenaArray
#include "../mesh.hpp"           // MeshBlock

// Constructor
// Inputs:
//   pb: pointer to block containing this grid
//   pin: pointer to runtime inputs
Coordinates::Coordinates(MeshBlock *pb, ParameterInput *pin)
{
  // Set pointer to host MeshBlock
  pmy_block = pb;

  // Initialize volume-averated positions and spacings: x-direction
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST; i++)
    pb->x1v(i) = 0.5 * (pb->x1f(i) + pb->x1f(i+1));
  for (int i = pb->is-NGHOST; i <= pb->ie+NGHOST-1; i++)
    pb->dx1v(i) = pb->x1v(i+1) - pb->x1v(i);

  // Initialize volume-averated positions and spacings: y-direction
  if (pb->block_size.nx2 == 1)  // no extent
  {
    pb->x2v(pb->js) = 0.5 * (pb->x2f(pb->js) + pb->x2f(pb->js+1));
    pb->dx2v(pb->js) = pb->dx2f(pb->js);
  }
  else  // extended
  {
    for (int j = pb->js-NGHOST; j <= pb->je+NGHOST; j++)
      pb->x2v(j) = 0.5 * (pb->x2f(j) + pb->x2f(j+1));
    for (int j = pb->js-NGHOST; j <= pb->je+NGHOST-1; j++)
      pb->dx2v(j) = pb->x2v(j+1) - pb->x2v(j);
  }

  // Initialize volume-averated positions and spacings: z-direction
  if (pb->block_size.nx3 == 1)  // no extent
  {
    pb->x3v(pb->ks) = 0.5 * (pb->x3f(pb->ks) + pb->x3f(pb->ks+1));
    pb->dx3v(pb->ks) = pb->dx3f(pb->ks);
  }
  else  // extended
  {
    for (int k = pb->ks-NGHOST; k <= pb->ke+NGHOST; k++)
      pb->x3v(k) = 0.5 * (pb->x3f(k) + pb->x3f(k+1));
    for (int k = pb->ks-NGHOST; k <= pb->ke+NGHOST-1; k++)
      pb->dx3v(k) = pb->x3v(k+1) - pb->x3v(k);
  }
}

// Destructor
Coordinates::~Coordinates()
{
}

// Function for computing areas orthogonal to x
// Inputs:
//   k: z-index
//   j: y-index
//   il, iu: x-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to x
// Notes:
//   \Delta A = \Delta y * \Delta z
void Coordinates::Face1Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> *pareas)
{
  Real &delta_y = pmy_block->dx2f(j);
  Real &delta_z = pmy_block->dx3f(k);
#pragma simd
  for (int i = il; i <= iu; i++)
  {
    Real &area = (*pareas)(i);
    area = delta_y * delta_z;
  }
  return;
}

// Function for computing areas orthogonal to y
// Inputs:
//   k: z-index
//   j: y-index
//   il, iu: x-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to y
// Notes:
//   \Delta A = \Delta x * \Delta z
void Coordinates::Face2Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> *pareas)
{
  Real &delta_z = pmy_block->dx3f(k);
#pragma simd
  for (int i = il; i <= iu; i++)
  {
    Real &area = (*pareas)(i);
    Real &delta_x = pmy_block->dx1f(i);
    area = delta_x * delta_z;
  }
  return;
}

// Function for computing areas orthogonal to z
// Inputs:
//   k: z-index
//   j: y-index
//   il, iu: x-index bounds
// Outputs:
//   areas: 1D array of interface areas orthogonal to z
// Notes:
//   \Delta A = \Delta x * \Delta y
void Coordinates::Face3Area(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> *pareas)
{
  Real &delta_y = pmy_block->dx2f(j);
#pragma simd
  for (int i = il; i <= iu; i++)
  {
    Real &area = (*pareas)(i);
    Real &delta_x = pmy_block->dx1f(i);
    area = delta_x * delta_y;
  }
  return;
}

// Function for computing cell volumes
// Inputs:
//   k: z-index
//   j: y-index
//   il, iu: x-index bounds
// Outputs:
//   volumes: 1D array of cell volumes
// Notes:
//   \Delta V = \Delta x * \Delta y * \Delta z
void Coordinates::CellVolume(const int k, const int j, const int il, const int iu,
    AthenaArray<Real> *pvolumes)
{
  Real &delta_y = pmy_block->dx2f(j);
  Real &delta_z = pmy_block->dx3f(k);
#pragma simd
  for (int i = il; i <= iu; i++)
  {
    Real &volume = (*pvolumes)(i);
    Real &delta_x = pmy_block->dx1f(i);
    volume = delta_x * delta_y * delta_z;
  }
  return;
}

// Function for computing displacement vector
// Inputs:
//   pt1, pt2: endpoints
// Outputs:
//   returned value: vector from pt1 to pt2
// Notes:
//   currently not implemented correctly
//   TODO: implement?
ThreeVector VectorBetweenPoints(const ThreeVector pt1, const ThreeVector pt2)
{
  ThreeVector r;
  r.x1 = p1.x1 - p2.x1;
  r.x2 = p1.x2 - p2.x2;
  r.x3 = p1.x3 - p2.x3;
  return r;
}

// Function for computing source terms
// Inputs:
//   dt: size of timestep
//   prim: full grid of primitive values at beginning of half timestep
//   cons: full grid of conserved variables at end of half timestep
// Outputs:
//   cons: source terms added
// Notes:
//   source terms all vanish identically
void Coordinates::CoordinateSourceTerms(Real dt, const AthenaArray<Real> &prim,
    AthenaArray<Real> &cons)
{
  return;
}

// Function for computing metric terms
// Inputs:
//   k: z-index
//   j: y-index
// Outputs:
//   g: array of metric components in 1D
//   g_inv: array of inverse metric components in 1D
// Notes:
//   both arrays assumed to be 0-initialized
void Coordinates::CellMetric(const int k, const int j, AthenaArray<Real> &g,
    AthenaArray<Real> &g_inv)
{
#pragma simd
  for (int i = pmy_block->is-NGHOST; i <= pmy_block->ie+NGHOST; i++)
  {
    // TODO: should 0's be set explicitly?
    g(I00,i) = -1.0;
    g(I11,i) = 1.0;
    g(I22,i) = 1.0;
    g(I33,i) = 1.0;
    g_inv(I00,i) = -1.0;
    g_inv(I11,i) = 1.0;
    g_inv(I22,i) = 1.0;
    g_inv(I33,i) = 1.0;
  }
  return;
}

// Function for transforming primitives to locally flat frame: x-interface
// Inputs:
//   k: z-index
//   j: y-index
//   pprim: pointer to array of primitives in 1D, using global coordinates
// Outputs:
//   pprim: pointer to values in local coordinates
// Notes:
//   transformation is trivial
void Coordinates::PrimToLocal1(const int k, const int j, AthenaArray<Real> *pprim)
{
  return;
}

// Function for transforming primitives to locally flat frame: y-interface
// Inputs:
//   k: z-index
//   j: y-index
//   pprim: pointer to array of primitives in 1D, using global coordinates
// Outputs:
//   pprim: pointer to values in local coordinates
// Notes:
//   transformation is trivial
void Coordinates::PrimToLocal2(const int k, const int j, AthenaArray<Real> *pprim)
{
  return;
}

// Function for transforming primitives to locally flat frame: z-interface
// Inputs:
//   k: z-index
//   j: y-index
//   pprim: pointer to array of primitives in 1D, using global coordinates
// Outputs:
//   pprim: pointer to values in local coordinates
// Notes:
//   transformation is trivial
void Coordinates::PrimToLocal3(const int k, const int j, AthenaArray<Real> *pprim)
{
  return;
}

// Function for transforming fluxes to global frame: x-interface
// Inputs:
//   k: z-index
//   j: y-index
//   pflux: pointer to array of fluxes in 1D, using local coordinates
// Outputs:
//   pflux: pointer to values overwritten in global coordinates
// Notes:
//   transformation is trivial except for sign change from lowering time index
void Coordinates::FluxToGlobal1(const int k, const int j, AthenaArray<Real> *pflux)
{
  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie+1; i++)
  {
    // Extract fluxes for reading and writing
    Real &txt = (*pflux)(IEN,i);
    Real &t10 = (*pflux)(IEN,i);

    // Set new fluxes
    t10 = -txt;
  }
  return;
}

// Function for transforming fluxes to global frame: y-interface
// Inputs:
//   k: z-index
//   j: y-index
//   pflux: pointer to array of fluxes in 1D, using local coordinates
// Outputs:
//   pflux: pointer to values overwritten in global coordinates
// Notes:
//   transformation is trivial except for sign change from lowering time index
void Coordinates::FluxToGlobal2(const int k, const int j, AthenaArray<Real> *pflux)
{
  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
  {
    // Extract fluxes for reading and writing
    Real &tyt = (*pflux)(IEN,i);
    Real &t20 = (*pflux)(IEN,i);

    // Set new fluxes
    t20 = -tyt;
  }
  return;
}

// Function for transforming fluxes to global frame: z-interface
// Inputs:
//   k: z-index
//   j: y-index
//   pflux: pointer to array of fluxes in 1D, using local coordinates
// Outputs:
//   pflux: pointer to values overwritten in global coordinates
// Notes:
//   transformation is trivial except for sign change from lowering time index
void Coordinates::FluxToGlobal3(const int k, const int j, AthenaArray<Real> *pflux)
{
  // Go through 1D block of cells
#pragma simd
  for (int i = pmy_block->is; i <= pmy_block->ie; i++)
  {
    // Extract fluxes for reading and writing
    Real &tzt = (*pflux)(IEN,i);
    Real &t30 = (*pflux)(IEN,i);

    // Set new fluxes
    t30 = -tzt;
  }
  return;
}
