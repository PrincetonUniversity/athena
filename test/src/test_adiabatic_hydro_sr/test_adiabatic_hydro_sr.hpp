// Class for testing variable inversion in adiabatic hydro

#ifndef TEST_ADIABATIC_HYDRO_SR_HPP
#define TEST_ADIABATIC_HYDRO_SR_HPP

// Primary headers
#include "../../../src/fluid.hpp"
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
class AdiabaticHydroSRTest : public GeneralTest
{
 protected:

  // Workspace
  ParameterInput *inputs;
  std::string input_file;
  Mesh *mesh;
  Fluid *pfluid;
  AthenaArray<Real> cons, prim;
  Real prim_expected[NVAR];

  // Constructor
  AdiabaticHydroSRTest(std::string input) : GeneralTest(1.0e-5, 1.0e-5)
  {
    input_file = input;
  }

  // Destructor
  virtual ~AdiabaticHydroSRTest() {};

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
    pfluid = pblock->pfluid;
    cons.NewAthenaArray(NVAR,1,1,1);
    prim.NewAthenaArray(NVAR,1,1,1);
  }

  // Function invoked after each test
  virtual void TearDown()
  {
    //prim.DeleteAthenaArray();
    //cons.DeleteAthenaArray();
    //delete mesh;
    delete inputs;
  }
};

// Gamma = 4/3
class AdiabaticHydroSRTest1 : public AdiabaticHydroSRTest
{
 protected:
  AdiabaticHydroSRTest1() : AdiabaticHydroSRTest("inputs/athinput.adiabatic_hydro_sr_a") {};
  virtual ~AdiabaticHydroSRTest1() {};
};

}

#endif  // TEST_ADIABATIC_HYDRO_SR_HPP
