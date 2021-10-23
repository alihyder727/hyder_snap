#ifndef JUNO_MWR_DATA_HPP_
#define JUNO_MWR_DATA_HPP_

// C/C++ headers
#include <vector>

// Athena++ headers
#include "../athena.hpp"

class MeshBlock;
class ParameterInput;
class Coordinates;
class Hydro;
class Thermodynamics;
class Radiation;

struct TPData {
  Real P, T, ERR;
};

class RadioData {
public:
  // data
  Real **target;
  Real ***icov;
  Real *z1, *p1, *t1;

  Real Tstd, Tlen, NH3std, NH3len;
  Real grav;
  int is, ie, js, je, ks, ke, nx1;
  //int nwave, nangle;
  int iH2O, iNH3;
  std::vector<Real> pdiv, zdiv, zfrac;
  std::vector<Real> TpSample, NH3pSample;
  std::vector<TPData> tpdata;

  Thermodynamics *pthermo;
  Hydro *phydro;
  Radiation *prad;
  Coordinates *pcoord;

  // functions
  RadioData(MeshBlock *pmb, ParameterInput *pin,
    Real grav_, int iH2O_, int iNH3_,
    Real *z1_, Real *p1_, Real *t1_, int nx1_, bool testrun_ = false);
  ~RadioData();
  void ReadObservationFile(char const *fname);
  void WriteObservationFile(char const *fname);

  bool IsTestRun() {
    return testrun_;
  }

private:
  bool testrun_;
};

#endif
