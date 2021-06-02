#ifndef PARTICLES_HPP
#define PARTICLES_HPP

// C++ headers
#include <vector>
#include <string>

// Athena++ classes headers
#include "../athena.hpp"

class MeshBlock;
class MaterialPoint;
class ParticleBuffer;

class Particles {
public:
// data
  MeshBlock* pmy_block;
  std::string myname;
  Particles *prev, *next;
  ParticleBuffer *ppb;
  AthenaArray<Real> c, dc;
  std::vector<MaterialPoint> mp;

// functions
  Particles(MeshBlock *pmb, ParameterInput *pin);
  Particles(MeshBlock *pmb, ParameterInput *pin, std::string name);
  virtual ~Particles();
  Particles(Particles const& other);
  Particles& operator=(Particles const& other);

  template<typename T> Particles* AddParticles(T const& other) {
    T* pt = new T(other);
    Particles *p = this;
    while (p->next != nullptr) p = p->next;
    p->next = pt;
    p->next->prev = p;
    p->next->next = nullptr;
    return p->next;
  }

  Particles* FindParticle(std::string name);
  void Translate(Real dt);
  void AddHydroVelocities(AthenaArray<Real> const &w);
  void SetMeshLocation();
  void AggregateMass(AthenaArray<Real> &c_sum);
  void Particulate(AthenaArray<Real> &c_new, int seeds_per_cell);

  virtual void TimeIntegrate(Real time, Real dt) {
    Translate(dt);
  }

protected:
  AthenaArray<Real> vol_;
  std::vector<Real> coordinates_;
  std::vector<int> lengths_;
};

class TwoPhaseCloudParticles : public Particles {
public:
  TwoPhaseCloudParticles(MeshBlock *pmb, ParameterInput *pin);
  ~TwoPhaseCloudParticles() {}

protected:
  int max_number_;
  int seeds_per_cell_;
};

#endif
