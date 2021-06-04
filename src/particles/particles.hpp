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
  std::vector<MaterialPoint> mp, mp1;

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
  void TranslateEuler(std::vector<MaterialPoint> &mp, Real dt);
  void ExchangeHydro(AthenaArray<Real> &du, AthenaArray<Real> const &w);
  void AggregateMass(AthenaArray<Real> &c_sum);
  int Categories() const {
    return c.GetDim4();
  }
  std::string CategoryName(int i) const {
    return cnames_.at(i);
  }

  virtual void Particulate(std::vector<MaterialPoint> &mp, AthenaArray<Real> &c_dif);
  virtual void TimeIntegrate(std::vector<MaterialPoint> &mp, Real time, Real dt);
  virtual void WeightedAverage(std::vector<MaterialPoint> &mp_out,
    std::vector<MaterialPoint> const& mp_in, Real ave_wghts[]);

protected:
  AthenaArray<Real> vol_;
  std::vector<Real> coordinates_;
  std::vector<int> lengths_;
  std::vector<std::string> cnames_;
};

class TwoPhaseCloudParticles : public Particles {
public:
  TwoPhaseCloudParticles(MeshBlock *pmb, ParameterInput *pin, std::string name);
  ~TwoPhaseCloudParticles() {}
  //void TimeIntegrate(std::vector<MaterialPoint> &mp, Real time, Real dt);
  //void Particulate(AthenaArray<Real> &c_dif);

protected:
  int max_number_;
  int seeds_per_cell_;
};

#endif
