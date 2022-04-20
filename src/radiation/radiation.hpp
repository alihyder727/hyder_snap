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
  Real rad, wav, wgt;
};

struct Direction {
  Real mu, phi;
};

enum class RadiationType {fixed = 0, dynamic = 1, band = 2};

class RadiationBand {
public:
  // data
  std::string myname;
  Radiation *pmy_rad;
  RadiationBand *prev, *next;
  Absorber *pabs;

  // spectra
  Spectrum *spec;
  int nspec, npmom;   // number of spectra and Legendre moments

  // band radiation results
  AthenaArray<Real> btau, bssa, bpmom;
  AthenaArray<Real> bflxup, bflxdn;
  AthenaArray<Real> btoa;

  // functions
  RadiationBand(Radiation *prad); // delayed initialization
  RadiationBand(Radiation *prad, std::string name, ParameterInput *pin);
  ~RadiationBand();
  void AddAbsorber(std::string name, std::string file, ParameterInput *pin);
  void AddAbsorber(Absorber *pab);
  void SetSpectralProperties(AthenaArray<Real> const& w, int k, int j, int il, int iu);
  void RadtranFlux(Direction const rin, Real dist_au,
    int k, int j, int il, int iu);
  void RadtranRadiance(Direction const rin, Direction const *rout, int nrout, Real dist_au,
    int k, int j, int il, int iu);

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
  Real cooldown, current;
  CelestrialBody *planet;
  RadiationType rtype;

  // functions
  Radiation(MeshBlock *pmb); // delayed initialization
  Radiation(MeshBlock *pmb, ParameterInput *pin);
  ~Radiation();
  RadiationBand* GetBand(int n);
  int GetNumBands();
  std::vector<Direction> GetOutgoingRays();
  std::vector<Direction> GetIncomingRays();
  void CalculateFluxes(AthenaArray<Real> const& w, Real time,
    int k, int j, int il, int iu);
  void CalculateRadiances(AthenaArray<Real> const& w, Real time,
    int k, int j, int il, int iu);
  void AddRadiativeFluxes(AthenaArray<Real>& x1flux,
    int k, int j, int il, int iu);

  // restart functions
  size_t RestartDataSizeInBytes();
  size_t DumpRestartData(char *pdst);
  size_t LoadRestartData(char *psrc);

protected:
  // reserved incoming and outgoing rays
  Direction *rin_, *rout_;
  int nrin_, nrout_;
  Real stellar_dist_au_;
};

#endif
