#ifndef PARTICLES_HPP
#define PARTICLES_HPP

// C++ headers
#include <vector>
#include <string>

// Athena++ classes headers
#include "../athena.hpp"
#include "material_point.hpp"

class MeshBlock;
class ParticleBuffer;

class Particles {
public:
// data
  MeshBlock* pmy_block;
  std::string myname;
  Particles *prev, *next;
  ParticleBuffer *ppb;
  AthenaArray<Real> c, c1;
  std::vector<MaterialPoint> mp, mp1;

// functions
  Particles(MeshBlock *pmb, ParameterInput *pin);
  Particles(MeshBlock *pmb, ParameterInput *pin, std::string name, int nct);
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

  std::string CategoryName(int i) const {
    return cnames_.at(i);
  }

  int GetNextId() {
    int id;
    if (available_ids_.size() > 0) {
      id = available_ids_.back();
      available_ids_.pop_back();
    } else
      id = mp.size()+1;
    return id;
  }

  Real GetMassRatio(int i, Real mu) const {
    return mu_[i]/mu;
  }

  int CountParticlesInCell(int t, int k, int j, int i) {
    int num = 0;
    MaterialPoint *pc = pcell_(t,k,j,i);
    while (pc != nullptr) {
      pc = pc->next;
      num++;
    }
    return num;
  }

  Real GetTotalCv(int k, int j, int i) {
    Real cvt = 0.;
    for (int t = 0; t < c.GetDim4(); ++t)
      cvt += c(t,k,j,i)*cc_[t];
    return cvt;
  }

  Real GetMolarDensity(int k, int j ,int i) {
    Real mol = 0.;
    for (int t = 0; t < c.GetDim4(); ++t)
      mol += c(t,k,j,i)/mu_[t];
    return mol;
  }

  Particles* FindParticle(std::string name);
  void Initialize();
  void AggregateDensity(AthenaArray<Real> &c1, std::vector<MaterialPoint> const& mp);
  void Particulate(std::vector<MaterialPoint> &mp, AthenaArray<Real> const& c);

  virtual void ExchangeHydro(std::vector<MaterialPoint> &mp, AthenaArray<Real> &du,
    AthenaArray<Real> const &w);
  virtual void TimeIntegrate(std::vector<MaterialPoint> &mp, Real time, Real dt);
  virtual void WeightedAverage(std::vector<MaterialPoint> &mp_out,
    std::vector<MaterialPoint> const& mp_in, Real ave_wghts[]);

protected:
  std::vector<Real> coordinates_;
  std::vector<int> dims_;
  std::vector<std::string> cnames_;
  std::vector<int> available_ids_;
  //! heat capacity
  std::vector<Real> cc_;
  //! mean molecular weight
  std::vector<Real> mu_;
  AthenaArray<MaterialPoint*> pcell_;

  int seeds_per_cell_;
  int nmax_per_cell_;
  Real density_floor_;
};

class SimpleCloudParticles : public Particles {
public:
  SimpleCloudParticles(MeshBlock *pmb, ParameterInput *pin, std::string name);
  ~SimpleCloudParticles() {}
  void ExchangeHydro(std::vector<MaterialPoint> &mp, AthenaArray<Real> &du,
    AthenaArray<Real> const &w);
};

class TwoPhaseCloudParticles : public Particles {
public:
  TwoPhaseCloudParticles(MeshBlock *pmb, ParameterInput *pin, std::string name);
  ~TwoPhaseCloudParticles() {}
  void ExchangeHydro(std::vector<MaterialPoint> &mp, AthenaArray<Real> &du,
    AthenaArray<Real> const &w);
  //void TimeIntegrate(std::vector<MaterialPoint> &mp, Real time, Real dt);
};

#endif
