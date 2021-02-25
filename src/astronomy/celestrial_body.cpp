// C/C++ headers
#include <stdexcept>
#include <sstream>

// Athena++ headers
#include "../math/core.h"
#include "../math/interpolation.h"
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../radiation/radiation.hpp"
#include "../utils/utils.hpp"
#include "celestrial_body.hpp"

void CelestrialBody::ReadCelestrialData_(ParameterInput *pin, std::string myname)
{
  char entry[80];
  sprintf(entry, "%s.re", name.c_str());
  re = km2m(pin->GetOrAddReal("astronomy", entry, 0.));

  sprintf(entry, "%s.rp", name.c_str());
  rp = km2m(pin->GetOrAddReal("astronomy", entry, re));

  sprintf(entry, "%s.obliq", name.c_str());
  obliq = deg2rad(pin->GetOrAddReal("astronomy", entry, 0.));

  sprintf(entry, "%s.omega", name.c_str());
  omega = pin->GetOrAddReal("astronomy", entry, 0.);

  sprintf(entry, "%s.orbit_a", name.c_str());
  orbit_a = au2m(pin->GetOrAddReal("astronomy", entry, 0.));

  sprintf(entry, "%s.orbit_e", name.c_str());
  orbit_e = pin->GetOrAddReal("astronomy", entry, 0.);

  sprintf(entry, "%s.orbit_i", name.c_str());
  orbit_i = deg2rad(pin->GetOrAddReal("astronomy", entry, 0.));

  sprintf(entry, "%s.orbit_p", name.c_str());
  orbit_p = day2sec(pin->GetOrAddReal("astronomy", entry, 0.));

  sprintf(entry, "%s.equinox", name.c_str());
  equinox = pin->GetOrAddReal("astronomy", entry, 0.);

  sprintf(entry, "%s.grav_eq", name.c_str());
  grav_eq = pin->GetOrAddReal("astronomy", entry, 0.);
}

CelestrialBody::CelestrialBody(ParameterInput *pin):
  parent(NULL), spec_(NULL), il_(-1)
{
  std::stringstream msg;
  name = pin->GetOrAddString("astronomy", "planet", "unknown");
  ReadCelestrialData_(pin, name);

  std::string parent_name = pin->GetOrAddString("astronomy", name + ".parent", "");
  if (!parent_name.empty())
    parent = new CelestrialBody(pin, parent_name);
  else
    parent = NULL;

  std::string sfile = pin->GetOrAddString("astronomy", name + ".spec_file", "");
  if (!sfile.empty()) {
    if (!FileExists(sfile)) {
      msg << "### FATAL ERROR in initializing CelestrialBody"
          << std::endl << "Cannot open spectral file " << sfile;
      ATHENA_ERROR(msg);
    } else 
      ReadSpectraFile(sfile);
  } else {
    spec_ = new float_triplet [1];
    nspec_ = 1;
  }
}

CelestrialBody::CelestrialBody(ParameterInput *pin, std::string myname):
  parent(NULL), name(myname), spec_(NULL), il_(-1)
{
  std::stringstream msg;
  ReadCelestrialData_(pin, name);

  std::string parent_name = pin->GetOrAddString("astronomy", name + ".parent", "");
  if (!parent_name.empty())
    parent = new CelestrialBody(pin, parent_name);
  else
    parent = NULL;

  std::string sfile = pin->GetOrAddString("astronomy", name + ".spec_file", "");
  if (!sfile.empty()) {
    if (!FileExists(sfile)) {
      msg << "### FATAL ERROR in initializing CelestrialBody"
          << std::endl << "Cannot open spectral file " << sfile;
      ATHENA_ERROR(msg);
    } else 
      ReadSpectraFile(sfile);
  } else {
    spec_ = new float_triplet [1];
    nspec_ = 1;
  }
}

CelestrialBody::~CelestrialBody()
{
  if (parent != NULL)
    delete parent;
  delete [] spec_;
}

void CelestrialBody::ReadSpectraFile(std::string sfile)
{
  AthenaArray<Real> spectrum;
  ReadDataTable(spectrum, sfile);

  nspec_ = spectrum.GetDim2();

  if (spec_ != NULL) delete [] spec_;
  spec_ = new float_triplet [nspec_];
  for (int i = 0; i < nspec_; ++i) {
    spec_[i].x = spectrum(i,0);
    spec_[i].y = spectrum(i,1);
  }

  spline(nspec_, spec_, 0., 0.);
}

void CelestrialBody::ParentZenithAngle(Real *mu, Real *phi, Real time, Real colat, Real lon)
{
  Real lat = M_PI/2. - colat;
  *mu = cos((time*(omega - 2.*M_PI/orbit_p)) - lon + M_PI)*cos(lat);
  *phi = 0.;
}

Real CelestrialBody::ParentInsolationFlux(Real wav, Real dist)
{
  Real dx;
  // assuming ascending in wav
  if ((il_ >= 0) && (wav < parent->spec_[il_].x)) il_ = -1;
  il_ = find_place_in_table(parent->nspec_, parent->spec_, wav, &dx, il_);
  Real flux = splint(wav, parent->spec_ + il_, dx);
  return flux/(dist*dist);
}

Real CelestrialBody::ParentDistanceInAu(Real time)
{
  Real orbit_b = sqrt(1. - orbit_e*orbit_e)*orbit_a;
  Real orbit_c = orbit_b*orbit_b/orbit_a;
  return m2au(orbit_c/(1. + orbit_e*cos(2.*M_PI/orbit_p*time - equinox - M_PI/2.)));
}
