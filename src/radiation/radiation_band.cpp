// C/C++ headers
#include <vector>
#include <sstream>
#include <stdexcept>
#include <type_traits>

// Athena++ header
#include "../parameter_input.hpp"
#include "absorber.hpp"
#include "radiation.hpp"
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../particles/particles.hpp"
#include "../scalars/scalars.hpp"
#include "../utils/utils.hpp" // Vectorize, ReadTabular, ReplaceChar

RadiationBand::RadiationBand(Radiation *prad):
  myname(""), npmom(0), nspec(1),
  pmy_rad(prad), prev(NULL), next(NULL), pabs(NULL)
{
  spec = new Spectrum [1];
  tem_ = new Real [1];
}

RadiationBand::RadiationBand(Radiation *prad, std::string name, ParameterInput *pin)
{
  ATHENA_LOG("RadiationBand");
  std::stringstream msg;

  myname = name;
  prev = NULL;
  next = NULL;
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
  NewCArray(tau_, ncells1, nspec);
  std::fill(tau_[0], tau_[0] + ncells1*nspec, 0.);
  NewCArray(ssa_, ncells1, nspec);
  NewCArray(pmom_, ncells1, nspec, npmom+1);
  NewCArray(flxup_, ncells1+1, nspec);
  NewCArray(flxdn_, ncells1+1, nspec);
  NewCArray(toa_, nrout, nspec);

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

  if (pabs->next != NULL) {
    pabs = pabs->next;
    delete pabs->prev;  // remove first one
  }

  // band parameters
  sprintf(astr, "%s.alpha", name.c_str());
  alpha_ = pin->GetOrAddReal("radiation", astr, 0.);

  // initialize radiative transfer solver
#ifdef RT_DISORT
  init_disort(pin);
#endif

  std::cout << "- Finish setting spectral BAND " << name
            << " with RANGE " << spec[0].wav << " - " << spec[nspec-1].wav
            << " and LENGTH " << nspec << std::endl;
}

RadiationBand::~RadiationBand()
{
  if (prev != NULL) prev->next = next;
  if (next != NULL) next->prev = prev;
  if (pabs != NULL) {
    while (pabs->prev != NULL)  // should not be true
      delete pabs->prev;
    while (pabs->next != NULL)
      delete pabs->next;
    delete pabs;
  }

  delete[] spec;
  delete[] tem_;
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
  if (pab->prev != NULL) {
    pab->prev->next = NULL;
    pab->prev = NULL;
  }
  
  if (pabs == NULL) { // new absorber
    pabs = pab;
  } else {  // attach to tail
    Absorber *p = pabs;
    while (p->next != NULL) p = p->next;
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

// setting optical properties
void RadiationBand::SetSpectralProperties(AthenaArray<Real> const& w,
  int k, int j, int il, int iu)
{
  MeshBlock *pmb = pmy_rad->pmy_block;
  int is = pmb->is, ie = pmb->ie;

  // set tau, ssalb, pmom, etc...
  int ncells1 = pmy_rad->pmy_block->ncells1;
  std::fill(*tau_, *tau_ + ncells1*nspec, 0.);
  std::fill(*ssa_, *ssa_ + ncells1*nspec, 0.);
  std::fill(**pmom_, **pmom_ + ncells1*nspec*(npmom+1), 0.);

  Absorber *a = pabs;
  Thermodynamics *pthermo = pmy_rad->pmy_block->pthermo;
  Coordinates *pcoord = pmy_rad->pmy_block->pcoord;
  Hydro *phydro = pmy_rad->pmy_block->phydro;
  PassiveScalars *pscalars = pmy_rad->pmy_block->pscalars;
  Real *mypmom = new Real[1+npmom];

  int num_clouds = 0.;
  Particles *ppart = pmb->ppart;
  while (ppart != nullptr) {
    num_clouds += ppart->u.GetDim4();
    ppart = ppart->next;
  }
  Real *q = new Real [NHYDRO + num_clouds];
  Real c[1];
  Real s[NSCALARS];

  while (a != NULL) {
    for (int i = il; i <= iu; ++i) {
      for (int n = 0; n < NSCALARS; ++n) s[n] = pscalars->s(n,k,j,i);
      pthermo->PrimitiveToChemical(q, w.at(k,j,i));
      // molar concentration to molar mixing ratio
      // \todo TODO: do we need it?
      Real nmols = q[IPR]/(Thermodynamics::Rgas*q[IDN]);
      for (int n = 1; n <= NVAPOR; ++n) q[n] /= nmols;

      // molar density of clouds
      ppart = pmb->ppart;
      int ip  = 0;
      while (ppart != nullptr) {
        for (int n = 0; n < ppart->u.GetDim4(); ++n) {
          q[NHYDRO + ip] = ppart->u(n,k,j,i)/ppart->GetMolecularWeight(n);
          ip++;
        }
        ppart = ppart->next;
      }

      tem_[i] = q[IDN];
      //std::cout << i << " " << tem_[i] << std::endl;
      for (int m = 0; m < nspec; ++m) {
        Real kcoeff = a->Attenuation(spec[m].wav, q, c, s);  // 1/m
        Real dssalb = a->SingleScatteringAlbedo(spec[m].wav, q)*kcoeff;
        // tau 
        tau_[i][m] += kcoeff;
        // ssalb
        ssa_[i][m] += dssalb;
        // pmom
        a->PhaseMomentum(spec[m].wav, q, mypmom, npmom);
        for (int p = 0; p <= npmom; ++p)
          pmom_[i][m][p] += mypmom[p]*dssalb;
      }
    }
    a = a->next;
  }

  delete [] mypmom; 
  delete [] q;

  // absorption coefficiunts -> optical thickness
  for (int i = il; i <= iu; ++i) {
    for (int m = 0; m < nspec; ++m) {
      if (tau_[i][m] > 1e-6 && ssa_[i][m] > 1e-6) {  // has scattering
        for (int p = 0; p <= npmom; ++p)
          pmom_[i][m][p] /= ssa_[i][m];
        ssa_[i][m] /= tau_[i][m];
      } else {
        ssa_[i][m] = 0.;
        pmom_[i][m][0] = 1.;
        for (int p = 1; p <= npmom; ++p)
          pmom_[i][m][p] = 0.;
      }
      if (HYDROSTATIC) {
        Real grav = -phydro->hsrc.GetG1();
        // \delta z = \delt Z * P/(\rho g H), cli, TODO
        tau_[i][m] *= pcoord->dx1f(i)*w(IPR,k,j,i)/(w(IDN,k,j,i)*grav*phydro->scale_height);
      } else {
        tau_[i][m] *= pcoord->dx1f(i);
      }
    }

    // aggregated band properties
    btau(k,j,i) = 0; bssa(k,j,i) = 0;
    for (int p = 0; p <= npmom; ++p) bpmom(p,k,j,i) = 0.;

    for (int m = 0; m < nspec; ++m) {
      btau(k,j,i) += tau_[i][m];
      bssa(k,j,i) += ssa_[i][m]*tau_[i][m];
      for (int p = 0; p <= npmom; ++p)
        bpmom(p,k,j,i) += pmom_[i][m][p]*ssa_[i][m];
    }

    for (int p = 0; p <= npmom; ++p)
      bpmom(p,k,j,i) /= bssa(k,j,i);
    bssa(k,j,i) /= btau(k,j,i);
    btau(k,j,i) /= nspec;
  }
}
