// Tests for adiabatic hydro variable inversion in Schwarzschild

// Primary header
#include "test_inversion_schwarzschild_adiabatic_hydro.hpp"

// gtest headers
#include "gtest/gtest.h"

// Athena headers
#include "../../../src/athena.hpp"                   // enums, macros
#include "../../../src/athena_arrays.hpp"            // array access
#include "../../../src/mesh.hpp"                     // Mesh, MeshDomain, MeshBlock
#include "../../../src/coordinates/coordinates.hpp"  // Coordinates
#include "../../../src/fluid/fluid.hpp"              // Fluid
#include "../../../src/fluid/eos/eos.hpp"            // FluidEqnOfState

// Test inversion with no motion, perfect initial guess
TEST_F(AdiabaticHydroInversionTest1, TestStaticPerfectGuess)
{
  // Set primitives and corresponding conserved quantities
  set_primitives(1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
  mesh->pdomain->pblock->pcoord->PrimToCons(prim_original, cons);

  // Check conserved-to-primitive inversion
  mesh->pdomain->pblock->pfluid->pf_eos->ConservedToPrimitive(cons, prim_close,
      prim_returned);
  for (int n = 0; n < NFLUID; n++)
    for (int k = ks; k <= ke; k++)
      for (int j = js; j <= je; j++)
        for (int i = is; i <= ie; i++)
          EXPECT_DOUBLE_TOL(prim_original(n,k,j,i), prim_returned(n,k,j,i));
}

// Test inversion with no motion
TEST_F(AdiabaticHydroInversionTest1, TestStatic)
{
  // Set primitives and corresponding conserved quantities
  set_primitives(1.0, 1.0, 0.0, 0.0, 0.0, 1.0e-6, 1.01);
  mesh->pdomain->pblock->pcoord->PrimToCons(prim_original, cons);

  // Check conserved-to-primitive inversion
  mesh->pdomain->pblock->pfluid->pf_eos->ConservedToPrimitive(cons, prim_close,
      prim_returned);
  for (int n = 0; n < NFLUID; n++)
    for (int k = ks; k <= ke; k++)
      for (int j = js; j <= je; j++)
        for (int i = is; i <= ie; i++)
          EXPECT_DOUBLE_TOL(prim_original(n,k,j,i), prim_returned(n,k,j,i));
}

// Test inversion with no motion and very low pressure
TEST_F(AdiabaticHydroInversionTest1, TestStaticLowPressure)
{
  // Set primitives and corresponding conserved quantities
  set_primitives(1.0, 1.0e-6, 0.0, 0.0, 0.0, 1.0e-8, 1.01);
  mesh->pdomain->pblock->pcoord->PrimToCons(prim_original, cons);

  // Check conserved-to-primitive inversion
  mesh->pdomain->pblock->pfluid->pf_eos->ConservedToPrimitive(cons, prim_close,
      prim_returned);
  for (int n = 0; n < NFLUID; n++)
    for (int k = ks; k <= ke; k++)
      for (int j = js; j <= je; j++)
        for (int i = is; i <= ie; i++)
          EXPECT_DOUBLE_TOL(prim_original(n,k,j,i), prim_returned(n,k,j,i));
}

// Test inversion with motion, perfect initial guess
TEST_F(AdiabaticHydroInversionTest1, TestMovingPerfectGuess)
{
  // Set primitives and corresponding conserved quantities
  set_primitives(1.0, 1.0, 0.05, 0.02, -0.03, 0.0, 1.0);
  mesh->pdomain->pblock->pcoord->PrimToCons(prim_original, cons);

  // Check conserved-to-primitive inversion
  mesh->pdomain->pblock->pfluid->pf_eos->ConservedToPrimitive(cons, prim_close,
      prim_returned);
  for (int n = 0; n < NFLUID; n++)
    for (int k = ks; k <= ke; k++)
      for (int j = js; j <= je; j++)
        for (int i = is; i <= ie; i++)
          EXPECT_DOUBLE_TOL(prim_original(n,k,j,i), prim_returned(n,k,j,i));
}

// Test inversion with motion
TEST_F(AdiabaticHydroInversionTest1, TestMoving)
{
  // Set primitives and corresponding conserved quantities
  set_primitives(1.0, 1.0, 0.05, 0.02, -0.03, 1.0e-6, 1.01);
  mesh->pdomain->pblock->pcoord->PrimToCons(prim_original, cons);

  // Check conserved-to-primitive inversion
  mesh->pdomain->pblock->pfluid->pf_eos->ConservedToPrimitive(cons, prim_close,
      prim_returned);
  for (int n = 0; n < NFLUID; n++)
    for (int k = ks; k <= ke; k++)
      for (int j = js; j <= je; j++)
        for (int i = is; i <= ie; i++)
          EXPECT_DOUBLE_TOL(prim_original(n,k,j,i), prim_returned(n,k,j,i));
}

// Test inversion with motion and very low pressure
TEST_F(AdiabaticHydroInversionTest1, TestMovingLowPressure)
{
  // Set primitives and corresponding conserved quantities
  set_primitives(1.0, 1.0e-6, 0.05, 0.02, -0.03, 1.0e-8, 1.01);
  mesh->pdomain->pblock->pcoord->PrimToCons(prim_original, cons);

  // Check conserved-to-primitive inversion
  mesh->pdomain->pblock->pfluid->pf_eos->ConservedToPrimitive(cons, prim_close,
      prim_returned);
  for (int n = 0; n < NFLUID; n++)
    for (int k = ks; k <= ke; k++)
      for (int j = js; j <= je; j++)
        for (int i = is; i <= ie; i++)
          EXPECT_DOUBLE_TOL(prim_original(n,k,j,i), prim_returned(n,k,j,i));
}
