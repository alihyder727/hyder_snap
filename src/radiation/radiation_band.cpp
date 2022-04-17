// C/C++ headers
#include <vector>
#include <stdexcept>
#include <type_traits>

// Athena++ header
#include "../parameter_input.hpp"
#include "absorber.hpp"
#include "radiation.hpp"
#include "../mesh/mesh.hpp"
#include "../debugger/debugger.hpp"
#include "../utils/utils.hpp" // Vectorize, ReadTabular, ReplaceChar

RadiationBand::RadiationBand(Radiation *prad):
  myname(""), npmom(0), nspec(1),
  pmy_rad(prad), prev(nullptr), next(nullptr), pabs(nullptr)
{
  spec = new Spectrum [1];
  tem_ = new Real [1];
}

RadiationBand::RadiationBand(Radiation *prad, std::string name, ParameterInput *pin)
{
  Debugger *pdbg = prad->pmy_block->pdebug;
  pdbg->Enter("RadiationBand " + name);
  std::stringstream &msg = pdbg->msg;

  myname = name;
  prev = nullptr;
  next = nullptr;
  pmy_rad = prad;

  // number of Legendre moments
  npmom = pin->GetOrAddInteger("radiation", "npmom", 0);

  // name radiation band in the format of "min_wave max_wave nwave"
  std::string str = pin->GetString("radiation", name);
  char default_file[80];
  sprintf(default_file, "kcoeff.%s.nc", str.c_str());
  ReplaceChar(default_file, ' ', '-');

  std::vector<Real> v = Vectorize<Real>(str.c_str());
  if (v.size() != 3) {
    msg << "### FATAL ERROR in construction function RadiationBand"
        << std::endl << "Length of '" << name << "' "
        << "must be 3.";
    ATHENA_ERROR(msg);
  }

  // set default wave number and weights
  nspec = (int)v[2];
  if (nspec < 1) {
    msg << "### FATAL ERROR in construction function RadiationBand"
        << "Length of some spectral band is not a positive number";
    ATHENA_ERROR(msg);
  }

  spec = new Spectrum [nspec];
  if (nspec == 1) {
    if (v[0] != v[1]) {
      msg << "### FATAL ERROR in construction function RadiationBand"
          << std::endl << "The first spectrum must equal the last spectrum "
          << "if the length of the spectral band is 1.";
      ATHENA_ERROR(msg);
    }
    spec[0].wav = v[0];
    spec[0].wgt = 1.;
  } else {
    Real dwave = (v[1] - v[0])/(nspec - 1);
    for (int i = 0; i < nspec; ++i) {
      spec[i].wav = v[0] + dwave*i;
      spec[i].wgt = (i == 0) || (i == nspec - 1) ? 0.5*dwave : dwave;
    }
  }

  // outgoing radiation direction (mu,phi) in degree
  str = pin->GetOrAddString("radiation", "outdir", "(0.,0.)");
  std::vector<std::string> dstr = Vectorize<std::string>(str.c_str());
  Real nrout = dstr.size();

  // allocate memory
  MeshBlock *pmb = prad->pmy_block;
  int ncells1 = pmb->ncells1;
  int ncells2 = pmb->ncells2;
  int ncells3 = pmb->ncells3;

  // spectral properties
  tem_ = new Real [ncells1];
  temf_ = new Real [ncells1+1];
  NewCArray(tau_, nspec, ncells1);
  std::fill(tau_[0], tau_[0] + nspec*ncells1, 0.);
  NewCArray(ssa_, nspec, ncells1);
  std::fill(ssa_[0], ssa_[0] + nspec*ncells1, 0.);
  NewCArray(pmom_, nspec, ncells1, npmom+1);
  std::fill(pmom_[0][0], pmom_[0][0] + nspec*ncells1*(npmom+1), 0.);
  NewCArray(flxup_, nspec, ncells1+1);
  NewCArray(flxdn_, nspec, ncells1+1);
  NewCArray(toa_, nspec, nrout);

  // band properties
  btau.NewAthenaArray(ncells3, ncells2, ncells1);
  bssa.NewAthenaArray(ncells3, ncells2, ncells1);
  bpmom.NewAthenaArray(npmom+1, ncells3, ncells2, ncells1);
  bflxup.NewAthenaArray(ncells3, ncells2, ncells1+1);
  bflxdn.NewAthenaArray(ncells3, ncells2, ncells1+1);
  btoa.NewAthenaArray(nrout, ncells3, ncells2);

  // absorbers
  char astr[1024];
  sprintf(astr, "%s.absorbers", name.c_str());
  str = pin->GetOrAddString("radiation", astr, "");
  std::vector<std::string> aname = Vectorize<std::string>(str.c_str());

  pabs = new Absorber(this);  // first one is empty

  for (int i = 0; i < aname.size(); ++i) {
    sprintf(astr, "%s.%s", name.c_str(), aname[i].c_str());
    std::string afile = pin->GetOrAddString("radiation", astr, default_file);
    AddAbsorber(aname[i], afile, pin);
  }

  if (pabs->next != nullptr) {
    pabs = pabs->next;
    delete pabs->prev;  // remove first one
  }

  // band parameters
  sprintf(astr, "%s.alpha", name.c_str());
  alpha_ = pin->GetOrAddReal("radiation", astr, 0.);

  msg << "- spectral range = " << spec[0].wav << " - " << spec[nspec-1].wav << std::endl
      << "- number of lines = " << nspec << std::endl;

  // initialize radiative transfer solver
#ifdef RT_DISORT
  init_disort(pin);
#endif

  pdbg->Leave();
}

RadiationBand::~RadiationBand()
{
  if (prev != nullptr) prev->next = next;
  if (next != nullptr) next->prev = prev;
  if (pabs != nullptr) {
    while (pabs->prev != nullptr)  // should not be true
      delete pabs->prev;
    while (pabs->next != nullptr)
      delete pabs->next;
    delete pabs;
  }

  delete[] spec;
  delete[] tem_;
  delete[] temf_;
  FreeCArray(tau_);
  FreeCArray(ssa_);
  FreeCArray(pmom_);
  FreeCArray(flxup_);
  FreeCArray(flxdn_);
  FreeCArray(toa_);

  // destroy radiative transfer solver
#ifdef RT_DISORT
  free_disort();
#endif
}

void RadiationBand::AddAbsorber(Absorber *pab) {
  // detach the current one
  if (pab->prev != nullptr) {
    pab->prev->next = nullptr;
    pab->prev = nullptr;
  }
  
  if (pabs == nullptr) { // new absorber
    pabs = pab;
  } else {  // attach to tail
    Absorber *p = pabs;
    while (p->next != nullptr) p = p->next;
    p->next = pab;
    p->next->prev = p;
  }

  pab->pmy_band = this;
}

// overide in the pgen file
void __attribute__((weak)) RadiationBand::AddAbsorber(
  std::string name, std::string file, ParameterInput *pin)
{}

// overide in rtsolver folder
void __attribute__((weak)) RadiationBand::RadtranFlux(
  Direction const rin, Real dist, int k, int j, int il, int iu)
{}

// overide in rtsolver folder
void __attribute__((weak)) RadiationBand::RadtranRadiance(
  Direction const rin, Direction const *rout, int nrout, Real dist,
  int k, int j, int il, int iu)
{}
