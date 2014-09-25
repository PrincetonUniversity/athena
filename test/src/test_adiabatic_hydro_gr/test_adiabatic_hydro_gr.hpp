// Class for testing variable inversion in adiabatic hydro (Minkowski Cartesian metric)

#ifndef TEST_ADIABATIC_HYDRO_GR_HPP
#define TEST_ADIABATIC_HYDRO_GR_HPP

// Primary headers
#include "../../../src/fluid/eos/eos.hpp"
#include "../test.hpp"

// C++ headers
#include <cmath>   // sqrt()
#include <string>  // string

// Athena headers
#include "../../../src/athena.hpp"           // enums, macros, Real
#include "../../../src/athena_arrays.hpp"    // AthenaArray
#include "../../../src/mesh.hpp"             // Mesh, MeshDomain, MeshBlock
#include "../../../src/fluid/fluid.hpp"      // Fluid
#include "../../../src/parameter_input.hpp"  // ParameterInput

namespace {

// Test class
class AdiabaticHydroGRTest : public GeneralTest
{
 protected:

  // Workspace
  ParameterInput *inputs;
  std::string input_file;
  Mesh *mesh;
  FluidEqnOfState *pf_eos;
  AthenaArray<Real> cons, prim;
  Real prim_expected[NFLUID];

  // Constructor
  AdiabaticHydroGRTest(std::string input) : GeneralTest(1.0e-5, 1.0e-5)
  {
    input_file = input;
  }

  // Destructor
  virtual ~AdiabaticHydroGRTest() {};

  // Function invoked before each test
  virtual void SetUp()
  {
    inputs = new ParameterInput;
    inputs->LoadFromFile(input_file);
    const char *argv[] = {"athena", "-i", input_file.c_str()};
    inputs->ModifyFromCmdline(3, const_cast<char **>(argv));
    mesh = new Mesh(inputs);
    mesh->ForAllDomains(init_fluid, inputs);
    MeshBlock *pblock = mesh->pdomain->pblock;
    pblock->is = NGHOST;
    pblock->ie = -NGHOST;
    pblock->js = pblock->je = 0;
    pblock->ks = pblock->ke = 0;
    pf_eos = pblock->pfluid->pf_eos;
    cons.NewAthenaArray(NFLUID,1,1,1);
    prim.NewAthenaArray(NFLUID,1,1,1);
  }

  // Function invoked after each test
  virtual void TearDown()
  {
    cons.DeleteAthenaArray();
    prim.DeleteAthenaArray();
    delete mesh;
    delete inputs;
  }

  // Function for copying expected primitives into malleable array
  void set_primitives()
  {
    prim(IDN,0,0,0) = prim_expected[IDN];
    prim(IEN,0,0,0) = prim_expected[IEN];
    prim(IM1,0,0,0) = prim_expected[IM1];
    prim(IM2,0,0,0) = prim_expected[IM2];
    prim(IM3,0,0,0) = prim_expected[IM3];
    return;
  }
};

// Gamma = 4/3
class AdiabaticHydroGRTest1 : public AdiabaticHydroGRTest
{
 protected:

  // Constructor
  AdiabaticHydroGRTest1()
    : AdiabaticHydroGRTest("inputs/athinput.adiabatic_hydro_gr_a") {};

  // Destructor
  virtual ~AdiabaticHydroGRTest1() {};

  // Function for calculating conserved quantities from primitives
  void primitive_to_conserved()
  {
    // Extract primitives
    Real rho = prim_expected[IDN];
    Real pgas = prim_expected[IEN];
    Real v1 = prim_expected[IM1];
    Real v2 = prim_expected[IM2];
    Real v3 = prim_expected[IM3];

    // Set conserved quantites
    Real gamma_prime = 4.0;
    Real gamma_lorentz_sq = 1.0 / (1.0 - v1*v1 - v2*v2 - v3*v3);
    Real gamma_lorentz = std::sqrt(gamma_lorentz_sq);
    Real rho_h = rho + gamma_prime * pgas;
    cons(IDN,0,0,0) = gamma_lorentz * rho;
    cons(IEN,0,0,0) = -gamma_lorentz_sq * rho_h + pgas;
    cons(IM1,0,0,0) = gamma_lorentz_sq * rho_h * v1;
    cons(IM2,0,0,0) = gamma_lorentz_sq * rho_h * v2;
    cons(IM3,0,0,0) = gamma_lorentz_sq * rho_h * v3;
    return;
  }
};

// Gamma = 5/3
class AdiabaticHydroGRTest2 : public AdiabaticHydroGRTest,
                              public ::testing::WithParamInterface<Real>
{
 protected:

  // Constructor
  AdiabaticHydroGRTest2()
    : AdiabaticHydroGRTest("inputs/athinput.adiabatic_hydro_gr_b") {};

  // Destructor
  virtual ~AdiabaticHydroGRTest2() {};

  // Function for calculating conserved quantities from primitives
  void primitive_to_conserved()
  {
    // Extract primitives
    Real rho = prim_expected[IDN];
    Real pgas = prim_expected[IEN];
    Real v1 = prim_expected[IM1];
    Real v2 = prim_expected[IM2];
    Real v3 = prim_expected[IM3];

    // Set conserved quantites
    Real gamma_prime = 2.5;
    Real gamma_lorentz_sq = 1.0 / (1.0 - v1*v1 - v2*v2 - v3*v3);
    Real gamma_lorentz = std::sqrt(gamma_lorentz_sq);
    Real rho_h = rho + gamma_prime * pgas;
    cons(IDN,0,0,0) = gamma_lorentz * rho;
    cons(IEN,0,0,0) = -gamma_lorentz_sq * rho_h + pgas;
    cons(IM1,0,0,0) = gamma_lorentz_sq * rho_h * v1;
    cons(IM2,0,0,0) = gamma_lorentz_sq * rho_h * v2;
    cons(IM3,0,0,0) = gamma_lorentz_sq * rho_h * v3;
    return;
  }
};

}

#endif  // TEST_ADIABATIC_HYDRO_GR_HPP
