// C/C++ headers
#include <vector>
#include <stdexcept>
#include <type_traits>

// Athena++ header
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../debugger/debugger.hpp"
#include "../utils/utils.hpp" // Vectorize, ReadTabular, ReplaceChar
#include "absorber.hpp"
#include "radiation.hpp"
#include "radiation_utils.hpp" // readRadiationDirections

RadiationBand::RadiationBand(Radiation *prad):
  myname(""), num_bins(1), bflags(prad->rflags),
  pmy_rad(prad), prev(nullptr), next(nullptr), pabs(nullptr)
{}

RadiationBand::RadiationBand(Radiation *prad, std::string name, ParameterInput *pin):
  myname(name), prev(nullptr), next(nullptr), pmy_rad(prad)
{
  Debugger *pdbg = prad->pmy_block->pdebug;
  pdbg->Enter("RadiationBand " + myname);
  std::stringstream &msg = pdbg->msg;

  // band flags
  setRadiationFlags(&bflags, pin->GetOrAddString("radiation", myname + ".flags", ""));

  // number of Legendre moments
  int npmom = pin->GetOrAddInteger("radiation", "npmom", 0);

  // name radiation band in the format of "min_wave max_wave nbins"
  std::string str = pin->GetString("radiation", myname);
  char default_file[80];
  sprintf(default_file, "kcoeff.%s.nc", str.c_str());
  ReplaceChar(default_file, ' ', '-');

  std::vector<Real> val = Vectorize<Real>(str.c_str());
  if (val.size() != 3) {
    msg << "### FATAL ERROR in construction function RadiationBand"
        << std::endl << "Length of '" << myname << "' "
        << "must be 3.";
    ATHENA_ERROR(msg);
  }

  // set wavenumber and weights
  wmin = val[0];
  wmax = val[1];
  num_bins = (int)val[2];
  if (num_bins < 1) {
    msg << "### FATAL ERROR in function RadiationBand::RadiationBand"
        << "Length of some spectral band is not a positive number";
    ATHENA_ERROR(msg);
  }

  spec.resize(num_bins);
  Real dwave = (val[1] - val[0])/num_bins;
  for (int i = 0; i < num_bins; ++i) {
    spec[i].wav1 = val[0] + dwave*i;
    spec[i].wav2 = val[0] + dwave*(i+1);
    spec[i].wght = val[0] == val[1] ? 1. : dwave;
  }

  if (bflags & RadiationFlags::CorrelatedK) {
    val = Vectorize<Real>(pin->GetString("radiatin", myname + ".gpoints").c_str());
    if (val.size() != num_bins) {
      msg << "### FATAL ERROR in function RadiationBand::RadiationBand"
          << std::endl << "Number of gpoints does not equal " << num_bins;
      ATHENA_ERROR(msg);
    }

    for (int i = 0; i < num_bins; ++i) {
      spec[i].wav1 = val[i];
      spec[i].wav2 = val[i];
    }

    val = Vectorize<Real>(pin->GetString("radiatin", myname + ".weights").c_str());
    if (val.size() != num_bins) {
      msg << "### FATAL ERROR in function RadiationBand::RadiationBand"
          << std::endl << "Number of weights does not equal " << num_bins;
      ATHENA_ERROR(msg);
    }

    for (int i = 0; i < num_bins; ++i)
      spec[i].wght = val[i];
  }

  // outgoing radiation direction (mu,phi) in degree
  if (pin->DoesParameterExist("radiatin", myname + ".outdir")) {
    str = pin->GetString("radiation", myname + ".outdir");
  } else
    str = pin->GetOrAddString("radiation", "outdir", "(0.,0.)");
  readRadiationDirections(rayOutput, str);

  // allocate memory
  MeshBlock *pmb = prad->pmy_block;
  int ncells1 = pmb->ncells1;
  int ncells2 = pmb->ncells2;
  int ncells3 = pmb->ncells3;

  // spectral properties
  tem_ = new Real [ncells1];
  temf_ = new Real [ncells1+1];
  NewCArray(tau_, num_bins, ncells1);
  std::fill(tau_[0], tau_[0] + num_bins*ncells1, 0.);
  NewCArray(ssa_, num_bins, ncells1);
  std::fill(ssa_[0], ssa_[0] + num_bins*ncells1, 0.);
  NewCArray(pmom_, num_bins, ncells1, npmom+1);
  std::fill(pmom_[0][0], pmom_[0][0] + num_bins*ncells1*(npmom+1), 0.);
  NewCArray(flxup_, num_bins, ncells1+1);
  NewCArray(flxdn_, num_bins, ncells1+1);
  NewCArray(toa_, num_bins, rayOutput.size());

  // band properties
  btau.NewAthenaArray(ncells3, ncells2, ncells1);
  bssa.NewAthenaArray(ncells3, ncells2, ncells1);
  bpmom.NewAthenaArray(npmom+1, ncells3, ncells2, ncells1);
  bflxup.NewAthenaArray(ncells3, ncells2, ncells1+1);
  bflxdn.NewAthenaArray(ncells3, ncells2, ncells1+1);
  //! \todo btoa needs to be changed
  btoa.NewAthenaArray(rayOutput.size(), ncells3, ncells2);

  // absorbers
  str = pin->GetOrAddString("radiation", name + ".absorbers", "");
  std::vector<std::string> aname = Vectorize<std::string>(str.c_str());

  pabs = new Absorber(this);  // first one is empty

  char astr[1024];
  for (int i = 0; i < aname.size(); ++i) {
    sprintf(astr, "%s.%s", name.c_str(), aname[i].c_str());
    std::string afile = pin->GetOrAddString("radiation", astr, default_file);
    addAbsorber(aname[i], afile, pin);
  }

  if (pabs->next != nullptr) {
    pabs = pabs->next;
    delete pabs->prev;  // remove first one
  }

  // band parameters
  alpha_ = pin->GetOrAddReal("radiation", name + ".alpha", 0.);

  msg << "- spectral range = " << wmin << " - " << wmax << std::endl
      << "- number of spectral bins = " << num_bins << std::endl;

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

void RadiationBand::addAbsorber(Absorber *pab) {
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

/*std::vector<Direction> RadiationBand::getOutgoingRays() {
  std::vector<Direction> dir;
  for (int i = 0; i < nrayOutput_; ++i)
    dir.push_back(rayOutput_[i]);
  return dir;
}

std::vector<Direction> RadiationBand::getIncomingRays() {
  std::vector<Direction> dir;
  for (int i = 0; i < nrayInput_; ++i)
    dir.push_back(rayInput_[i]);
  return dir;
}*/

// overide in the pgen file
void __attribute__((weak)) RadiationBand::addAbsorber(
  std::string name, std::string file, ParameterInput *pin)
{}

// overide in rtsolver folder
void __attribute__((weak)) RadiationBand::calculateRadiativeFlux(
  Direction rayInput, Real dist, int k, int j, int il, int iu)
{}

// overide in rtsolver folder
void __attribute__((weak)) RadiationBand::calculateRadiance(
  Direction rayInput, Real dist, int k, int j, int il, int iu)
{}
