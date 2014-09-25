// Class for testing adiabatic hydro variable inversion in Schwarzschild

#ifndef TEST_INVERSION_SCHWARZSCHILD_ADIABATIC_HYDRO_HPP
#define TEST_INVERSION_SCHWARZSCHILD_ADIABATIC_HYDRO_HPP

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
class AdiabaticHydroInversionTest : public GeneralTest
{
 protected:

  // Workspace
  ParameterInput *inputs;
  std::string input_file;
  Mesh *mesh;
  int is, ie, js, je, ks, ke;
  AthenaArray<Real> prim_original, prim_close, prim_returned, cons;

  // Constructor
  AdiabaticHydroInversionTest(std::string input) : GeneralTest(1.0e-5, 1.0e-5)
  {
    input_file = input;
  }

  // Destructor
  virtual ~AdiabaticHydroInversionTest() {};

  // Function invoked before each test
  virtual void SetUp()
  {
    inputs = new ParameterInput;
    inputs->LoadFromFile(input_file);
    const char *argv[] = {"athena", "-i", input_file.c_str()};
    inputs->ModifyFromCmdline(3, const_cast<char **>(argv));
    mesh = new Mesh(inputs);
    mesh->ForAllDomains(init_fluid, inputs);
    is = mesh->pdomain->pblock->is;
    ie = mesh->pdomain->pblock->ie;
    js = mesh->pdomain->pblock->js;
    je = mesh->pdomain->pblock->je;
    ks = mesh->pdomain->pblock->ks;
    ke = mesh->pdomain->pblock->ke;
    int n1 = mesh->pdomain->pblock->pfluid->w.GetDim1();
    int n2 = mesh->pdomain->pblock->pfluid->w.GetDim2();
    int n3 = mesh->pdomain->pblock->pfluid->w.GetDim3();
    int n4 = mesh->pdomain->pblock->pfluid->w.GetDim4();
    prim_original.NewAthenaArray(n4, n3, n2, n1);
    prim_close.NewAthenaArray(n4, n3, n2, n1);
    prim_returned.NewAthenaArray(n4, n3, n2, n1);
    cons.NewAthenaArray(n4, n3, n2, n1);
  }

  // Function invoked after each test
  virtual void TearDown()
  {
    prim_original.DeleteAthenaArray();
    prim_close.DeleteAthenaArray();
    prim_returned.DeleteAthenaArray();
    cons.DeleteAthenaArray();
    delete mesh;
    delete inputs;
  }

  // Function for setting primitives everywhere
  void set_primitives(Real rho, Real pgas, Real v1, Real v2, Real v3, Real offset,
      Real scale)
  {
    for (int k = ks; k <= ke; k++)
      for (int j = js; j <= je; j++)
        for (int i = is; i <= ie; i++)
        {
          prim_original(IDN,k,j,i) = rho;
          prim_original(IEN,k,j,i) = pgas;
          prim_original(IM1,k,j,i) = v1;
          prim_original(IM2,k,j,i) = v2;
          prim_original(IM3,k,j,i) = v3;
          prim_close(IDN,k,j,i) = (rho + offset) * scale;
          prim_close(IEN,k,j,i) = (pgas + offset) * scale;
          prim_close(IM1,k,j,i) = (v1 + offset) * scale;
          prim_close(IM2,k,j,i) = (v2 + offset) * scale;
          prim_close(IM3,k,j,i) = (v3 + offset) * scale;
        }
    return;
  }
};

// 2D, Gamma = 5/3
class AdiabaticHydroInversionTest1 : public AdiabaticHydroInversionTest
{
 protected:

  // Constructor
  AdiabaticHydroInversionTest1()
    : AdiabaticHydroInversionTest("inputs/athinput.geodesic_2d_odd") {};

  // Destructor
  virtual ~AdiabaticHydroInversionTest1() {};
};

}

#endif  // TEST_INVERSION_SCHWARZSCHILD_ADIABATIC_HYDRO_HPP
