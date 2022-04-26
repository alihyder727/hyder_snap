#ifndef RADIATION_HPP
#define RADIATION_HPP

// C/C++ headers
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../astronomy/celestrial_body.hpp"

class MeshBlock;
class ParameterInput;
class Absorber;
class Radiation;

#ifdef RT_DISORT
#undef SQR
extern "C" {
  #include "rtsolver/cdisort213/cdisort.h"
}
#endif

struct Spectrum {
  Real rad, wav1, wav2, wght;
};

struct Direction {
  Real mu, phi;
};

namespace RadiationFlags {
  const uint64_t None = 0LL;
  const uint64_t Dynamic = 1LL << 0;
  const uint64_t LineByLine = 1LL << 1;
  const uint64_t CorrelatedK = 1LL << 2;
  const uint64_t Planck = 1LL << 3;
  const uint64_t Star = 1LL << 4;
  const uint64_t Sphere = 1LL << 5;
  const uint64_t FluxOnly = 1LL << 6;
}

class RadiationBand {
public:
  // data
  std::string myname;
  int num_bins;
  uint64_t bflags;
  Radiation *pmy_rad;
  RadiationBand *prev, *next;
  Absorber *pabs;

  // spectra
  std::vector<Spectrum> spec;
  Real wmin, wmax;

  // band radiation results
  AthenaArray<Real> btau, bssa, bpmom;
  AthenaArray<Real> bflxup, bflxdn;
  //! btoa is a reference to radiance in Radiation
  AthenaArray<Real> btoa;

  // outgoing rays
  std::vector<Direction> rayOutput;

  // functions
  RadiationBand(Radiation *prad); // delayed initialization
  RadiationBand(Radiation *prad, std::string name, ParameterInput *pin);
  ~RadiationBand();

  void addAbsorber(std::string name, std::string file, ParameterInput *pin);
  void addAbsorber(Absorber *pab);
  void setSpectralProperties(AthenaArray<Real> const& w, int k, int j, int il, int iu);
  void calculateRadiativeFlux(Direction rayInput, Real dist_au, int k, int j, int il, int iu);
  void calculateRadiance(Direction rayInput, Real dist_au, int k, int j, int il, int iu);

#ifdef RT_DISORT
  void init_disort(ParameterInput *pin);
  void free_disort();
  disort_state ds;
  disort_output ds_out;
#endif

protected:
  Real **tau_, **ssa_, ***pmom_, *tem_, *temf_;
  Real **flxup_, **flxdn_;
  Real **toa_;
  Real alpha_;  // T ~ Ts*(\tau/\tau_s)^\alpha at lower boundary
};

class Radiation {
public:
  // constants
  static Real const hPlanck;
  static Real const hPlanck_cgs;
  static Real const cLight;
  static Real const cLight_cgs;
  static Real const stefanBoltzmann;

  // data
  MeshBlock *pmy_block;
  RadiationBand *pband;
  uint64_t rflags;
  Real cooldown, current;
  CelestrialBody *planet;

  // incomming rays and outgoing radiance
  std::vector<Direction> rayInput;
  AthenaArray<Real> radiance;

  // functions
  Radiation(MeshBlock *pmb); // delayed initialization
  Radiation(MeshBlock *pmb, ParameterInput *pin);
  ~Radiation();
  RadiationBand* getBand(int n);
  int getNumBands();
  void calculateRadiativeFluxes(AthenaArray<Real> const& w, Real time,
    int k, int j, int il, int iu);
  void calculateRadiances(AthenaArray<Real> const& w, Real time,
    int k, int j, int il, int iu);
  void addRadiativeFluxes(AthenaArray<Real>& x1flux,
    int k, int j, int il, int iu);
  void readRadiationBands(ParameterInput *pin, int &b);
  int getTotalNumberOutgoingRays();
  //std::vector<Direction> getIncomingRays();

  // restart functions
  size_t getRestartDataSizeInBytes();
  size_t dumpRestartData(char *pdst);
  size_t loadRestartData(char *psrc);

protected:
  Real stellarDistance_au_;
};

#endif
