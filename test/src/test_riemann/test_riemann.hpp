// Class for testing Riemann solvers

#ifndef TEST_RIEMANN_HPP
#define TEST_RIEMANN_HPP

// Primary headers
#include "../../../src/fluid/integrators/integrators.hpp"
#include "../test.hpp"

// C++ headers
#include <string>  // string

// Athena headers
#include "../../../src/athena.hpp"           // macros, Real
#include "../../../src/athena_arrays.hpp"    // AthenaArray
#include "../../../src/mesh.hpp"             // Mesh, MeshDomain, MeshBlock
#include "../../../src/parameter_input.hpp"  // ParameterInput

namespace {

// Test class
class RiemannTest : public GeneralTest
{
 protected:

  // Workspace
  ParameterInput *inputs;
  std::string input_file;
  Mesh *mesh;
  FluidIntegrator *pfluid_integrator;
  AthenaArray<Real> prim_left, prim_right, flux;
  Real flux_expected[NVAR];

  // Constructor
  RiemannTest(std::string input) : GeneralTest(1.0e-5, 1.0e-5)
  {
    input_file = input;
  }

  // Destructor
  virtual ~RiemannTest() {};

  // Function invoked before each test
  virtual void SetUp()
  {
    inputs = new ParameterInput;
    inputs->LoadFromFile(input_file);
    const char *argv[] = {"athena", "-i", input_file.c_str()};
    inputs->ModifyFromCmdline(3, const_cast<char **>(argv));
    mesh = new Mesh(inputs);
    mesh->ForAllDomains(init_fluid, inputs);
    pfluid_integrator = new FluidIntegrator(mesh->pdomain->pblock->pfluid);
    prim_left.NewAthenaArray(NVAR,1,1,1);
    prim_right.NewAthenaArray(NVAR,1,1,1);
    flux.NewAthenaArray(NVAR,1,1,1);
  }

  // Function invoked after each test
  virtual void TearDown()
  {
    //prim_left.DeleteAthenaArray();
    //prim_right.DeleteAthenaArray();
    //flux.DeleteAthenaArray();
    //delete pfluid_integrator;
    //delete mesh;
    delete inputs;
  }
};

// HLLC SR, Gamma = 5/3
class HLLCSRTest : public RiemannTest
{
 protected:
  HLLCSRTest() : RiemannTest("inputs/athinput.riemann_a") {};
  virtual ~HLLCSRTest() {};
};

}

#endif  // TEST_RIEMANN_HPP
