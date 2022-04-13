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

Real const Radiation::hPlanck = 6.63E-34;
Real const Radiation::hPlanck_cgs = 6.63E-27;
Real const Radiation::c = 3.E8;
Real const Radiation::c_cgs = 3.E10;

Radiation::Radiation(MeshBlock *pmb):
  pmy_block(pmb), pband(nullptr), dynamic_(false), beam_(-1.),
  cooldown(0.), current(0.), nrin_(1), nrout_(1), dist_(1.), planet(nullptr)
{
  rin_ = new Direction [1];
  rout_ = new Direction [1];
}

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin)
{
  pmb->pdebug->Enter("Radiation");
  //ATHENA_LOG("Radiation");
  pmy_block = pmb;
  pband = nullptr;
  RadiationBand *plast = pband;

  // incoming radiation direction (mu,phi) in degree
  std::string str = pin->GetOrAddString("radiation", "indir", "(0.,0.)");
  std::vector<std::string> dstr = Vectorize<std::string>(str.c_str());
  nrin_ = dstr.size();
  rin_ = new Direction [nrin_];
  for (int i = 0; i < nrin_; ++i) {
    rin_[i].phi = 0.;
    sscanf(dstr[i].c_str(), "(%lf,%lf)", &rin_[i].mu, &rin_[i].phi);
    rin_[i].mu = cos(deg2rad(rin_[i].mu));
    rin_[i].phi = deg2rad(rin_[i].phi);
    //std::cout << rout[i].mu << " " << rout[i].phi << std::endl;
  }

  // outgoing radiation direction (mu,phi) in degree
  str = pin->GetOrAddString("radiation", "outdir", "(0.,0.)");
  dstr = Vectorize<std::string>(str.c_str());
  nrout_ = dstr.size();
  rout_ = new Direction [nrout_];
  for (int i = 0; i < nrout_; ++i) {
    rout_[i].phi = 0.;
    sscanf(dstr[i].c_str(), "(%lf,%lf)", &rout_[i].mu, &rout_[i].phi);
    rout_[i].mu = cos(deg2rad(rout_[i].mu));
    rout_[i].phi = deg2rad(rout_[i].phi);
    //std::cout << rout[i].mu << " " << rout[i].phi << std::endl;
  }

  // distance to parent star
  dist_ = pin->GetOrAddReal("radiation", "distance", 1.);

  int b = 1;
  char name[80];
  while (true) {
    sprintf(name, "b%d", b);
    try {
      pin->GetString("radiation", name);
    } catch (const std::runtime_error& e) {
      break;
    }
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

  dynamic_ = pin->GetOrAddBoolean("radiation", "dynamic", false);
  beam_ = pin->GetOrAddReal("radiation", "beam", -1.);
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

  delete[] rin_;
  delete[] rout_;
  delete planet;
}

RadiationBand* Radiation::GetBand(int n) {
  std::stringstream msg;
  RadiationBand* p = pband;
  int b = 0;
  while (p != nullptr) {
    if (b++ == n) break;
    p = p->next;
  }
  return p;
}

int Radiation::GetNumBands() {
  int n = 0;
  RadiationBand* p = pband;
  while (p != nullptr) {
    p = p->next;
    n++;
  }
  return n;
}

std::vector<Direction> Radiation::GetOutgoingRays() {
  std::vector<Direction> dir;
  for (int i = 0; i < nrout_; ++i)
    dir.push_back(rout_[i]);
  return dir;
}

std::vector<Direction> Radiation::GetIncomingRays() {
  std::vector<Direction> dir;
  for (int i = 0; i < nrin_; ++i)
    dir.push_back(rin_[i]);
  return dir;
}

void Radiation::CalculateFluxes(AthenaArray<Real> const& w, Real time,
  int k, int j, int il, int iu)
{
  pmy_block->pdebug->Call("Radiation::CalculateFluxes");
  Coordinates *pcoord = pmy_block->pcoord;

  RadiationBand *p = pband;
  if (pband == nullptr) return;

  if (dynamic_) {
    planet->ParentZenithAngle(&rin_->mu, &rin_->phi, time, pcoord->x2v(j), pcoord->x3v(k));
    dist_ = planet->ParentDistanceInAu(time);
  }

  while (p != nullptr) {
    // iu ~= ie + 1
    p->SetSpectralProperties(w, k, j, il - NGHOST, iu + NGHOST - 1);
    p->RadtranFlux(*rin_, dist_, k, j, il, iu);
    p = p->next;
  }
  pmy_block->pdebug->Leave();
}

void Radiation::CalculateRadiances(AthenaArray<Real> const& w, Real time,
  int k, int j, int il, int iu)
{
  pmy_block->pdebug->Call("Radiation::CalculateRadiances");
  Coordinates *pcoord = pmy_block->pcoord;

  RadiationBand *p = pband;
  if (pband == nullptr) return;

  if (dynamic_) {
    planet->ParentZenithAngle(&rin_->mu, &rin_->phi, time, pcoord->x2v(j), pcoord->x3v(k));
    dist_ = planet->ParentDistanceInAu(time);
  }

  while (p != nullptr) {
    // iu ~= ie + 1
    p->SetSpectralProperties(w, k, j, il - NGHOST, iu + NGHOST - 1);
    p->RadtranRadiance(*rin_, rout_, nrout_, dist_, k, j, il, iu);
    p = p->next;
  }
  pmy_block->pdebug->Leave();
}

void Radiation::AddRadiativeFluxes(AthenaArray<Real>& x1flux, 
  int k, int j, int il, int iu)
{
  RadiationBand *p = pband;
  if (pband == nullptr) return;

  MeshBlock *pmb = pmy_block;
  pmb->pdebug->Call("Radiation::AddRadiativeFluxes");

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
