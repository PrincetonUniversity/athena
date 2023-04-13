#ifndef PARTICLES_PARTICLES_HPP_
#define PARTICLES_PARTICLES_HPP_
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//======================================================================================
//! \file particles.hpp
//! \brief defines classes for particle dynamics.
//======================================================================================

// C/C++ Standard Libraries
#include <string>
#include <vector>

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../outputs/outputs.hpp"
#include "../parameter_input.hpp"
#include "particle_buffer.hpp"
#include "particle_gravity.hpp"
#include "particle_mesh.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// Forward definitions
class ParticleGravity;

//--------------------------------------------------------------------------------------
//! \struct Neighbor
//  \brief defines a structure for links to neighbors

struct Neighbor {
  NeighborBlock *pnb;
  MeshBlock *pmb;
  Neighbor *next, *prev;

  Neighbor() : pnb(NULL), pmb(NULL), next(NULL), prev(NULL) {}
};

//----------------------------------------------------------------------------------------
//! \struct ParticleParameters
//! \brief container for parameters read from `<particle?>` block in the input file

struct ParticleParameters {
  int block_number, ipar;
  bool table_output, gravity;
  std::string block_name;
  std::string partype;

  ParticleParameters() : block_number(0), ipar(-1), table_output(false), gravity(false) {}
};

//--------------------------------------------------------------------------------------
//! \class Particles
//! \brief defines the base class for all implementations of particles.

class Particles {
friend class MeshBlock;  // Make writing initial conditions possible.
friend class OutputType;
friend class ParticleGravity;
friend class ParticleMesh;

 public:
  // Class methods
  static void AMRCoarseToFine(Particles *pparc, Particles *pparf, MeshBlock* pmbf);
  static void AMRFineToCoarse(Particles *pparc, Particles *pparf);
  static void Initialize(Mesh *pm, ParameterInput *pin);
  static void PostInitialize(Mesh *pm, ParameterInput *pin);
  static void FindDensityOnMesh(Mesh *pm, bool include_momentum);
  static void FormattedTableOutput(Mesh *pm, OutputParameters op);
  static void GetHistoryOutputNames(std::string output_names[], int ipar);
  static std::int64_t GetTotalNumber(Mesh *pm);

  // Class constant
  static const int NHISTORY = 8;  //!> number of variables in history output
  // number of particle containers
  static int num_particles, num_particles_grav, num_particles_output;

  // Constructor
  Particles(MeshBlock *pmb, ParameterInput *pin, ParticleParameters *pp);

  // Destructor
  virtual ~Particles();

  // virtual public methods that can be overriden
  virtual AthenaArray<Real> GetMassDensity() const;
  virtual void Integrate(int step);

  // Accessor
  Real GetMaximumWeight() const;
  AthenaArray<Real> GetVelocityField() const;
  bool IsGravity() { return isgravity_; }

  // Instance methods
  bool CheckInMeshBlock(Real x1, Real x2, Real x3);
  void AddHistoryOutput(Real data_sum[], int pos);
  void ClearBoundary();
  void ClearNeighbors();
  void LinkNeighbors(MeshBlockTree &tree, int64_t nrbx1, int64_t nrbx2, int64_t nrbx3,
                     int root_level);
  void AddOneParticle(Real x1, Real x2, Real x3, Real v1, Real v2, Real v3);
  void RemoveOneParticle(int k);
  void LoadParticleBuffer(ParticleBuffer *ppb, int k);
#ifdef MPI_PARALLEL
  void SendParticleBuffer(ParticleBuffer& send, int dst);
  void ReceiveParticleBuffer(int nb_rank, ParticleBuffer& recv,
                             enum BoundaryStatus& bstatus);
#endif
  void SendToNeighbors();
  void SetPositionIndices();
  bool ReceiveFromNeighbors();

  void StartReceivingParticlesShear();
  void SendParticlesShear();
  int FindTargetGidAlongX2(Real x2);
  void ClearBoundaryShear();
  bool ReceiveFromNeighborsShear();

  void DepositPMtoMesh(int stage);

  Real NewBlockTimeStep();
  virtual void FindLocalDensityOnMesh(bool include_momentum);

  // output individual particle history
  void OutputParticles(bool header);
  void OutputParticles(bool header, int kid);
  void OutputOneParticle(std::ostream &os, int k, bool header);
  void ToggleParHstOutFlag();

  std::size_t GetSizeInBytes();
  void UnpackParticlesForRestart(char *mbdata, std::size_t &os);
  void PackParticlesForRestart(char *&pdata);

  ParticleMesh *ppm;  //!> ptr to particle-mesh

 protected:
  // Class variables
  static bool initialized;  //!> whether or not the class is initialized
  static ParameterInput *pinput;

  // Instance variables
  std::string input_block_name, partype;

  int nint;          //!> numbers of integer particle properties
  int nreal;         //!> numbers of real particle properties
  int naux;          //!> number of auxiliary particle properties
  int nwork;         //!> number of working arrays for particles
  int nint_buf, nreal_buf; //!> number of properties for buffer

  int ipid;                 //!> index for the particle ID
  int ixp, iyp, izp;        // indices for the position components
  int ivpx, ivpy, ivpz;     // indices for the velocity components

  int ixp0, iyp0, izp0;     // indices for beginning position components
  int ivpx0, ivpy0, ivpz0;  // indices for beginning velocity components

  int ixi1, ixi2, ixi3;     // indices for position indices

  int imom1, imom2, imom3;  // indices for momentum components on mesh

  int imass, ish;
  int igx, igy, igz; // indices for gravity force

  // Instance methods
  virtual void AssignShorthands();  //!> Needs to be called everytime
                                    //!> intprop, realprop, & auxprop are resized
                                    //!> Be sure to call back when derived.
  virtual void AllocateMemory();    //!> Needs to be called in the derived class init

  int AddIntProperty();
  int AddRealProperty();
  int AddAuxProperty();
  int AddWorkingArray();

  void UpdateCapacity(int new_nparmax);  //!> Change the capacity of particle arrays
  void ConvertToDensity(bool include_momentum);
  void SaveStatus(); // x->x0, v->v0

  // Instance variables
  // std::uint64_t npar;     //!> number of particles
  // std::uint64_t nparmax;  //!> maximum number of particles per meshblock
  int npar;     //!> number of particles
  int nparmax;  //!> maximum number of particles per meshblock
  int my_ipar_;
  bool isgravity_; //!> flag for gravity
  bool parhstout_; //!> flag for individual particle history output
  Real mass;   //!> common mass of particle
  Real cfl_par;  //!> CFL number for particles

                               // Data attached to the particles:
  AthenaArray<int> intprop;    //!>   integer properties
  AthenaArray<Real> realprop;  //!>   real properties
  AthenaArray<Real> auxprop;   //!>   auxiliary properties (communicated when
                               //!>     particles moving to another meshblock)
  AthenaArray<Real> work;      //!>   working arrays (not communicated)

  std::vector<std::string> intfieldname, realfieldname, auxfieldname;

  ParticleGravity *ppgrav; //!> ptr to particle-gravity
                                       // Shorthands:
  AthenaArray<int> pid;                //!>   particle ID
  AthenaArray<Real> xp, yp, zp;        //   position
  AthenaArray<Real> vpx, vpy, vpz;     //   velocity
  AthenaArray<Real> xi1, xi2, xi3;     //   position indices in local meshblock
  AthenaArray<Real> xp0, yp0, zp0;     //   beginning position
  AthenaArray<Real> vpx0, vpy0, vpz0;  //   beginning velocity

  MeshBlock* pmy_block;  //!> MeshBlock pointer
  Mesh* pmy_mesh;        //!> Mesh pointer

  // shearing box parameters
  Real Omega_0_, qshear_, qomL;
  int ShBoxCoord_;
  bool orbital_advection_defined_;

 private:
  // Class method
  static void ProcessNewParticles(Mesh *pmesh, int ipar);

  // Instance methods
  // pure virtual methods
  virtual void SourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc)=0;
  virtual void UserSourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc)=0;
  virtual void ReactToMeshAux(Real t, Real dt, const AthenaArray<Real>& meshsrc)=0;
  virtual void DepositToMesh(Real t, Real dt, const AthenaArray<Real>& meshsrc,
                             AthenaArray<Real>& meshdst)=0;

  int CountNewParticles() const;
  void ApplyBoundaryConditions(int k, Real &x1, Real &x2, Real &x3);
  void EulerStep(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void FlushReceiveBuffer(ParticleBuffer& recv);
  void GetPositionIndices(int npar,
                          const AthenaArray<Real>& xp,
                          const AthenaArray<Real>& yp,
                          const AthenaArray<Real>& zp,
                          AthenaArray<Real>& xi1,
                          AthenaArray<Real>& xi2,
                          AthenaArray<Real>& xi3);
  void SetNewParticleID(int id);
  struct Neighbor* FindTargetNeighbor(
      int ox1, int ox2, int ox3, int xi1, int xi2, int xi3);
  void ApplyBoundaryConditionsShear(int k, Real &x1, Real &x2, Real &x3);

  // Class variable
  static std::vector<int> idmax;

  // Instance variables
  bool active1_, active2_, active3_;  // active dimensions
  int my_particle_num_;

  // MeshBlock-to-MeshBlock communication:
  BoundaryValues *pbval_;                            //!> ptr to my BoundaryValues
  Neighbor neighbor_[3][3][3];                       //!> links to neighbors
  ParticleBuffer recv_[56], recv_sh_[8];             //!> particle receive buffers
  enum BoundaryStatus bstatus_[56], bstatus_recv_sh_[8];  //!> boundary status
#ifdef MPI_PARALLEL
  static MPI_Comm my_comm;   //!> my MPI communicator
  ParticleBuffer send_[56],send_sh_[8];  //!> particle send buffers
  enum BoundaryStatus bstatus_send_sh_[8];  //!> comm. flags
#endif
};

//--------------------------------------------------------------------------------------
//! \fn Real Particles::GetMaximumWeight()
//! \brief returns the maximum weight on the mesh.

inline Real Particles::GetMaximumWeight() const {
  return ppm->FindMaximumWeight();
}

//--------------------------------------------------------------------------------------
//! \fn Real Particles::GetMaximumWeight()
//! \brief returns the maximum weight on the mesh.

inline AthenaArray<Real> Particles::GetMassDensity() const {
  return ppm->weight;
}

//--------------------------------------------------------------------------------------
//! \class DustParticles
//! \brief defines the class for dust particles that interact with the gas via drag
//!        force.

class DustParticles : public Particles {
friend class MeshBlock;

 public:
  void SetOneParticleMass(Real new_mass);
  bool GetBackReaction() { return backreaction; }
  bool GetDragForce() { return dragforce; }
  bool GetVariableTaus() { return variable_taus; }
  Real GetOneParticleMass() { return mass; }
  Real GetStoppingTime() { return taus0; }

  //!Constructor
  DustParticles(MeshBlock *pmb, ParameterInput *pin, ParticleParameters *pp);

  // Destructor
  ~DustParticles();

  // Instance method
  Real NewBlockTimeStep();

 private:
  bool backreaction;   //!> turn on/off back reaction
  bool dragforce;      //!> turn on/off drag force
  bool variable_taus;  //!> whether or not the stopping time is variable

  int iwx, iwy, iwz;         // indices for working arrays
  int idpx1, idpx2, idpx3;   // indices for momentum change
  int itaus;                 //!> index for stopping time

  Real taus0;  //!> constant/default stopping time (in code units)

  // Instance methods.
  void AssignShorthands() override;
  void SourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) override;
  void UserSourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) override;
  void ReactToMeshAux(Real t, Real dt, const AthenaArray<Real>& meshsrc) override;
  void DepositToMesh(Real t, Real dt, const AthenaArray<Real>& meshsrc,
                     AthenaArray<Real>& meshdst) override;
  void UserStoppingTime(Real t, Real dt, const AthenaArray<Real>& meshsrc);

  // Instance variables
  AthenaArray<Real> wx, wy, wz;        // shorthand for working arrays
  AthenaArray<Real> dpx1, dpx2, dpx3;  // shorthand for momentum change
  AthenaArray<Real> taus;              // shorthand for stopping time
};

//--------------------------------------------------------------------------------------
//! \class TracerParticles
//! \brief defines the class for velocity Tracer particles

class TracerParticles : public Particles {
friend class MeshBlock;

 public:
  //!Constructor
  TracerParticles(MeshBlock *pmb, ParameterInput *pin, ParticleParameters *pp);

  // Destructor
  ~TracerParticles();

  // Instance method
  void SetOneParticleMass(Real new_mass);
  Real GetOneParticleMass() { return mass; }

 private:
  int iwx, iwy, iwz;         // indices for working arrays

  // Instance methods.
  void AssignShorthands() override;
  void SourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) override;
  void UserSourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) override;
  void ReactToMeshAux(Real t, Real dt, const AthenaArray<Real>& meshsrc) override;
  void DepositToMesh(Real t, Real dt, const AthenaArray<Real>& meshsrc,
                     AthenaArray<Real>& meshdst) override;

  // Instance variables
  AthenaArray<Real> wx, wy, wz;        // shorthand for working arrays
};

//--------------------------------------------------------------------------------------
//! \class StarParticles
//! \brief defines the class for Star particles

class StarParticles : public Particles {
friend class MeshBlock;

 public:
  //!Constructor
  StarParticles(MeshBlock *pmb, ParameterInput *pin, ParticleParameters *pp);

  // Destructor
  ~StarParticles();

  // override integrator
  void Integrate(int step) override;

  void AddOneParticle(Real mass, Real x1, Real x2, Real x3,
                      Real v1, Real v2, Real v3);
  AthenaArray<Real> GetMassDensity() const override;
  void FindLocalDensityOnMesh(bool include_momentum) override;

 private:
  Real dt_old;
  int imetal, iage; // indices for additional Real properties
  int igas;                // indices for additional Aux properties
  // Instance methods.
  void AssignShorthands() override;
  void SourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) override;
  void UserSourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) override;
  void ReactToMeshAux(Real t, Real dt, const AthenaArray<Real>& meshsrc) override;
  void DepositToMesh(Real t, Real dt, const AthenaArray<Real>& meshsrc,
                     AthenaArray<Real>& meshdst) override;

  void Kick(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void Drift(Real t, Real dt);
  void BorisKick(Real t, Real dt);
  void Age(Real t, Real dt);

  void ExertTidalForce(Real t, Real dt);
  void PointMass(Real t, Real dt, Real gm);
  void ConstantAcceleration(Real t, Real dt, Real g1, Real g2, Real g3);

  // Instance variables
  AthenaArray<Real> mp, mzp, tage;        // shorthand for real properties
  AthenaArray<Real> fgas;                     // shorthand for aux properties
};

//--------------------------------------------------------------------------------------
//! \fn Real Particles::GetMaximumWeight()
//! \brief returns the maximum weight on the mesh.

inline AthenaArray<Real> StarParticles::GetMassDensity() const {
  return ppm->density;
}

#endif  // PARTICLES_PARTICLES_HPP_
