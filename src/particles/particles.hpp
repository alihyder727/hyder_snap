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
  AthenaArray<Real> u;
  std::vector<MaterialPoint> mp, mp1;

// functions
  Particles(MeshBlock *pmb, ParameterInput *pin);
  Particles(MeshBlock *pmb, ParameterInput *pin, std::string name, int nct);
  virtual ~Particles();
  Particles(Particles const& other);
  //Particles& operator=(Particles const& other);

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

  Real GetMolecularWeight(int n) const {
    return mu_[n];
  }

  Real GetCv(int n) const {
    return cc_[n];
  }

  Real GetDensityFloor() const {
    return density_floor_;
  }

  int GetNextId() {
    int id;
    if (available_ids_.size() > 0) {
      id = available_ids_.back();
      available_ids_.pop_back();
    } else
      id = mp.size() + mp1.size() +1;
    return id;
  }

  Particles* FindParticle(std::string name);
  void Initialize();
  void AggregateDensity(AthenaArray<Real> &c, std::vector<MaterialPoint> &mp);
  void Particulate(std::vector<MaterialPoint> &mp, AthenaArray<Real> const& c);

  virtual void ExchangeHydro(std::vector<MaterialPoint> &mp, AthenaArray<Real> &du,
    AthenaArray<Real> const &w, Real dt);
  virtual void TimeIntegrate(std::vector<MaterialPoint> &mp, Real time, Real dt);
  virtual void WeightedAverage(std::vector<MaterialPoint> &mp_out,
    std::vector<MaterialPoint> const& mp_in, Real ave_wghts[]);

protected:
  std::vector<Real> xface_;
  std::vector<Real> xcenter_;
  std::vector<int> dims_;
  std::vector<std::string> cnames_;
  std::vector<int> available_ids_;
  //! heat capacity
  std::vector<Real> cc_;
  //! mean molecular weight
  std::vector<Real> mu_;
  //! copy of c before doing chemistry
  AthenaArray<Real> u1_;
  //! linked list of particles in cell
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
    AthenaArray<Real> const &w, Real dt);
};

class TwoPhaseCloudParticles : public Particles {
public:
  TwoPhaseCloudParticles(MeshBlock *pmb, ParameterInput *pin, std::string name);
  ~TwoPhaseCloudParticles() {}
  void ExchangeHydro(std::vector<MaterialPoint> &mp, AthenaArray<Real> &du,
    AthenaArray<Real> const &w, Real dt);
  //void TimeIntegrate(std::vector<MaterialPoint> &mp, Real time, Real dt);
};

#endif
