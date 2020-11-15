#ifndef CELESTRIAL_BODY_HPP
#define CELESTRIAL_BODY_HPP

// C/C++ headers
#include <string>

// Athena++ headers
#include "../athena.hpp"

class ParameterInput;
struct float_triplet;

class CelestrialBody {
public:
// data
  CelestrialBody *parent;
  std::string name;
  Real re;        // equatorial radius [km]
  Real rp;        // polar radius [km]
  Real omega;     // rotation angular velocity [rad/s]
  Real obliq;     // obliquity [deg]
  Real orbit_a;   // orbital semi-major axis [au]
  Real orbit_e;   // orbital eccentricity [deg]
  Real orbit_i;   // orbital inclination to ecliptic [deg]
  Real orbit_p;   // orbital period [day]
  Real equinox;   // vernal equinox
  Real grav_eq;   // equatorial gravity at surface [m/s]

// functions
  CelestrialBody(ParameterInput *pin);
  CelestrialBody(ParameterInput *pin, std::string myname);
  ~CelestrialBody();

  void ReadSpectraFile(std::string sfile);
  void ParentZenithAngle(Real *mu, Real *phi, Real time, Real colat, Real lon);
  Real ParentInsolationFlux(Real wav, Real dist = 1.);
  Real ParentDistanceInAu(Real time);

protected:
  void ReadCelestrialData_(ParameterInput *pin, std::string myname);

// emission spectra data
  float_triplet *spec_;
  int nspec_;
  int il_;  // search pointer
};

Real GetGravity(char const *name, Real pclat);

#endif
