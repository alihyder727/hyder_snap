// C/C++ headers
#include <sstream>
#include <stdexcept>

// Athena++ header
//#include "../math_funcs.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "../math/core.h"
#include "../debugger/debugger.hpp"
#include "radiation.hpp"
#include "radiation_utils.hpp"  // setRadiationFlags

Real const Radiation::hPlanck = 6.63E-34;
Real const Radiation::hPlanck_cgs = 6.63E-27;
Real const Radiation::cLight = 3.E8;
Real const Radiation::cLight_cgs = 3.E10;
Real const Radiation::stefanBoltzmann = 5.670374419E-8;

Radiation::Radiation(MeshBlock *pmb):
  pmy_block(pmb), pband(nullptr), rflags(0LL),
  cooldown(0.), current(0.), planet(nullptr), stellarDistance_au_(1.)
{}

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin):
  pmy_block(pmb), pband(nullptr), rflags(0LL)
{
  pmb->pdebug->Enter("Radiation");
  std::stringstream &msg = pmb->pdebug->msg;

  // radiation flags
  setRadiationFlags(&rflags, pin->GetOrAddString("radiation", "flags", ""));

  // distance to parent star
  stellarDistance_au_ = pin->GetOrAddReal("radiation", "distance_au", 1.);
  msg << "- stellar distance = " << stellarDistance_au_ << " au" << std::endl;

  // radiation bands
  int b = 1;
  readRadiationBands(pin, b);

  // incoming radiation direction (mu,phi) in degree
  std::string str = pin->GetOrAddString("radiation", "indir", "(0.,0.)");
  readRadiationDirections(rayInput, str);

  // output radiance
  radiance_units = pin->GetOrAddString("radiation", "radiance_units", "w/(m^2.cm^{-1}.sr)");
  int nout = 0;
  RadiationBand *p = pband;
  while (p != nullptr) {
    nout += p->rayOutput.size();
    p = p->next;
  }
  if (nout > 0) {
    radiance.NewAthenaArray(nout, pmb->ncells3, pmb->ncells2);
    // set band toa
    p = pband;
    nout = 0;
    while (p != nullptr) {
      p->btoa.InitWithShallowSlice(radiance,3,nout,p->rayOutput.size());
      nout += p->rayOutput.size();
      p = p->next;
    }
  }

  // time control
  cooldown = pin->GetOrAddReal("radiation", "dt", 0.);
  current = 0.;

  planet = new CelestrialBody(pin);
  pmb->pdebug->Leave();
}

Radiation::~Radiation()
{
  if (pband != nullptr) {
    while (pband->prev != nullptr) // should not be true
      delete pband->prev;
    while (pband->next != nullptr)
      delete pband->next;
    delete pband;
  }
  delete planet;
}

RadiationBand* Radiation::getBand(int n) {
  std::stringstream msg;
  RadiationBand* p = pband;
  int b = 0;
  while (p != nullptr) {
    if (b++ == n) break;
    p = p->next;
  }
  return p;
}

int Radiation::getNumBands() {
  int n = 0;
  RadiationBand* p = pband;
  while (p != nullptr) {
    p = p->next;
    n++;
  }
  return n;
}

void Radiation::calculateRadiativeFluxes(AthenaArray<Real> const& w, Real time,
  int k, int j, int il, int iu)
{
  pmy_block->pdebug->Call("Radiation::calculateRadiativeFluxes");
  Coordinates *pcoord = pmy_block->pcoord;
  Real dist = stellarDistance_au_;

  RadiationBand *p = pband;
  if (pband == nullptr) return;

  while (p != nullptr) {
    Direction ray;
    if (p->bflags & RadiationFlags::Dynamic) {
      planet->ParentZenithAngle(&ray.mu, &ray.phi, time, pcoord->x2v(j), pcoord->x3v(k));
      dist = planet->ParentDistanceInAu(time);
    } else {
      ray = rayInput[0];
    }

    // iu ~= ie + 1
    p->setSpectralProperties(w, k, j, il - NGHOST, iu + NGHOST - 1);
    p->calculateRadiativeFlux(ray, dist, k, j, il, iu);
    p = p->next;
  }
  pmy_block->pdebug->Leave();
}

void Radiation::calculateRadiances(AthenaArray<Real> const& w, Real time,
  int k, int j, int il, int iu)
{
  pmy_block->pdebug->Call("Radiation::calculateRadiances");
  Coordinates *pcoord = pmy_block->pcoord;
  Real dist = stellarDistance_au_;

  RadiationBand *p = pband;
  if (pband == nullptr) return;

  Direction ray;
  if (rflags & RadiationFlags::Dynamic) {
    planet->ParentZenithAngle(&ray.mu, &ray.phi, time, pcoord->x2v(j), pcoord->x3v(k));
    dist = planet->ParentDistanceInAu(time);
  }

  while (p != nullptr) {
    // iu ~= ie + 1
    p->setSpectralProperties(w, k, j, il - NGHOST, iu + NGHOST - 1);
    p->calculateRadiance(ray, dist, k, j, il, iu);
    p = p->next;
  }
  pmy_block->pdebug->Leave();
}

void Radiation::addRadiativeFluxes(AthenaArray<Real>& x1flux, 
  int k, int j, int il, int iu)
{
  RadiationBand *p = pband;
  if (pband == nullptr) return;

  MeshBlock *pmb = pmy_block;
  pmb->pdebug->Call("Radiation::addRadiativeFluxes");

  // x1-flux divergence
  p = pband;
  while (p != nullptr) {
#pragma omp simd
    for (int i = il; i <= iu; ++i)
      x1flux(IEN,k,j,i) += p->bflxup(k,j,i) - p->bflxdn(k,j,i);
    p = p->next;
  }
  pmb->pdebug->Leave();
}

void Radiation::readRadiationBands(ParameterInput *pin, int &b)
{
  char name[80];
  RadiationBand *plast = pband;
  while (plast != nullptr)
    plast = plast->next;

  while (true) {
    sprintf(name, "b%d", b);
    if (!pin->DoesParameterExist("radiation", name))
      break;
    RadiationBand* p = new RadiationBand(this, name, pin);
    if (plast == nullptr) {
      plast = p;
      pband = p;
    } else {
      plast->next = p;
      plast->next->prev = plast;
      plast->next->next = nullptr;
      plast = plast->next;
    }
    b++;
  }

  if (pin->DoesParameterExist("radiation", "bands")) {
    ParameterInput* pin_next = new ParameterInput;
    IOWrapper infile;
    infile.Open(pin->GetString("radiation", "bands").c_str(), IOWrapper::FileMode::read);
    pin_next->LoadFromFile(infile);
    infile.Close();
    readRadiationBands(pin_next, b);
    delete pin_next;
  }
}

int Radiation::getTotalNumberOutgoingRays() {
  int num = 0;
  RadiationBand *p = pband;
  while (p != nullptr) {
    num += p->rayOutput.size();
    p = p->next;
  }
  return num;
}

size_t Radiation::getRestartDataSizeInBytes()
{
  size_t size = 0;

  RadiationBand *p = pband;
  while (p != nullptr) {
    size += p->bflxup.GetSizeInBytes() + p->bflxdn.GetSizeInBytes();
    p = p->next;
  }

  return size;
}

size_t Radiation::dumpRestartData(char *pdst)
{
  RadiationBand *p = pband;
  int offset = 0;
  while (p != nullptr) {
    std::memcpy(pdst + offset, p->bflxup.data(), p->bflxup.GetSizeInBytes());
    offset += p->bflxup.GetSizeInBytes();
    std::memcpy(pdst + offset, p->bflxdn.data(), p->bflxdn.GetSizeInBytes());
    offset += p->bflxdn.GetSizeInBytes();
    p = p->next;
  }

  return getRestartDataSizeInBytes();
}

size_t Radiation::loadRestartData(char *psrc)
{
  RadiationBand *p = pband;
  int offset = 0;
  while (p != nullptr) {
    std::memcpy(p->bflxup.data(), psrc + offset, p->bflxup.GetSizeInBytes());
    offset += p->bflxup.GetSizeInBytes();
    std::memcpy(p->bflxdn.data(), psrc + offset, p->bflxdn.GetSizeInBytes());
    offset += p->bflxdn.GetSizeInBytes();
    p = p->next;
  }

  return getRestartDataSizeInBytes();
}
