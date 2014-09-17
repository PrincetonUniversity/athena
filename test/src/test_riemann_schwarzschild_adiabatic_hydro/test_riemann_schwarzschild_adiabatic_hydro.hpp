// Class for testing HLLE GR Riemann solver

#ifndef TEST_RIEMANN_GR_HPP
#define TEST_RIEMANN_GR_HPP

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
  int is, ie, js, je, jm, ks, ncells1;
  AthenaArray<Real> prim_left, prim_right, flux, flux_expected;

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
    is = mesh->pdomain->pblock->is;
    ie = mesh->pdomain->pblock->ie;
    js = mesh->pdomain->pblock->js;
    je = mesh->pdomain->pblock->je;
    jm = (js + je + 1) / 2;
    ks = mesh->pdomain->pblock->ks;
    ncells1 = mesh->pdomain->pblock->block_size.nx1 + 2*NGHOST;
    prim_left.NewAthenaArray(NVAR,ncells1);
    prim_right.NewAthenaArray(NVAR,ncells1);
    flux.NewAthenaArray(NVAR,ncells1);
    flux_expected.NewAthenaArray(NVAR,ncells1);
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

// 2D geodesic infall, even number of cells in theta-direction
class Geodesic2DTest : public RiemannTest
{
 protected:
  Geodesic2DTest() : RiemannTest("inputs/athinput.geodesic_2d_even") {};
  virtual ~Geodesic2DTest() {};
};

}

#endif  // TEST_RIEMANN_GR_HPP
