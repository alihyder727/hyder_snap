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
#include "radiation.hpp"

Radiation::Radiation(MeshBlock *pmb):
  pmy_block(pmb), pband(NULL), dynamic_(false), beam_(-1.),
  cooldown(0.), current(0.), nrin_(1), nrout_(1), dist_(1.), planet(NULL)
{
  rin_ = new Direction [1];
  rout_ = new Direction [1];
}

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin)
{
  pmy_block = pmb;
  pband = NULL;
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
    if (plast == NULL) {
      plast = p;
      pband = p;
    } else {
      plast->next = p;
      plast->next->prev = plast;
      plast->next->next = NULL;
      plast = plast->next;
    }
    b++;
  }

  dynamic_ = pin->GetOrAddBoolean("radiation", "dynamic", false);
  beam_ = pin->GetOrAddReal("radiation", "beam", -1.);
  cooldown = pin->GetOrAddReal("radiation", "dt", 0.);
  current = 0.;

  planet = new CelestrialBody(pin);
}

Radiation::~Radiation()
{
  if (pband != NULL) {
    while (pband->prev != NULL) // should not be true
      delete pband->prev;
    while (pband->next != NULL)
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
  while (p != NULL) {
    if (b++ == n) break;
    p = p->next;
  }
  return p;
}

int Radiation::GetNumBands() {
  int n = 0;
  RadiationBand* p = pband;
  while (p != NULL) {
    p = p->next;
    n++;
  }
  return n;
}

void Radiation::CalculateFluxes(AthenaArray<Real> const& w, Real time,
  int k, int j, int il, int iu)
{
  Coordinates *pcoord = pmy_block->pcoord;

  RadiationBand *p = pband;
  if (pband == NULL) return;

  if (dynamic_) {
    planet->ParentZenithAngle(&rin_->mu, &rin_->phi, time, pcoord->x2v(j), pcoord->x3v(k));
    dist_ = planet->ParentDistanceInAu(time);
  }

  while (p != NULL) {
    p->SetSpectralProperties(w, k, j, il - NGHOST, iu + NGHOST - 1);
    p->RadtranFlux(*rin_, dist_, k, j, il, iu);
    p = p->next;
  }
}

void Radiation::CalculateRadiances(AthenaArray<Real> const& w, Real time,
  int k, int j, int il, int iu)
{
  Coordinates *pcoord = pmy_block->pcoord;

  RadiationBand *p = pband;
  if (pband == NULL) return;

  if (dynamic_) {
    planet->ParentZenithAngle(&rin_->mu, &rin_->phi, time, pcoord->x2v(j), pcoord->x3v(k));
    dist_ = planet->ParentDistanceInAu(time);
  }

  while (p != NULL) {
    p->SetSpectralProperties(w, k, j, il - NGHOST, iu + NGHOST - 1);
    p->RadtranRadiance(*rin_, rout_, nrout_, dist_, k, j, il, iu);
    p = p->next;
  }
}

void Radiation::AddRadiativeFluxes(AthenaArray<Real>& x1flux, 
  int k, int j, int il, int iu)
{
  RadiationBand *p = pband;
  if (pband == NULL) return;

  MeshBlock *pmb = pmy_block;

  // x1-flux divergence
  p = pband;
  while (p != NULL) {
#pragma omp simd
    for (int i = il; i <= iu; ++i)
      x1flux(IEN,k,j,i) += p->bflxup(k,j,i) - p->bflxdn(k,j,i);
    p = p->next;
  }
}
