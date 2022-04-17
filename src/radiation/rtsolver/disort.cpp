//! \file disort.cpp
//  \brief DISORT radiative transfer solver
//=====================================================

// C/C++ header
#include <cstdlib>
#include <sstream>
#include <algorithm>

// Athena++ headers
#include "../../reconstruct/interpolation.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../../utils/utils.hpp" // StringToArray, ReadTabular
#include "../../coordinates/coordinates.hpp"
#include "../../communicator/communicator.hpp"
#include "../../globals.hpp"
#include "../../debugger/debugger.hpp"
#include "../radiation_utils.hpp"
#include "../radiation.hpp"

// DISORT headers

#ifdef RT_DISORT

void RadiationBand::init_disort(ParameterInput *pin)
{
  Debugger *pdbg = pmy_rad->pmy_block->pdebug;
  pdbg->Enter("RadiationBand " + myname + "-disort");
  std::stringstream &msg = pdbg->msg;

  ds.nlyr = pmy_rad->pmy_block->pmy_mesh->mesh_size.nx1;
  ds.nmom = pin->GetInteger("radiation", "npmom");

  ds.nstr   = pin->GetInteger("disort", "nstr");
  ds.nphase = pin->GetInteger("disort", "nphase");
  ds.accur  = pin->GetOrAddReal("disort", "accur", 0.);

  // bottom boundary isotropic radiation
  if (pin->DoesParameterExist("radiation", myname + ".fluor_K")) {
    Real Tbot = pin->GetReal("radiation", myname + ".fluor_K");
    ds.bc.fluor = Radiation::stefanBoltzmann*pow(Tbot, 4.);
  } else
    ds.bc.fluor = pin->GetOrAddReal("disort", "fluor", 0.);

  msg << "- fluor = " << ds.bc.fluor << " w/m^2" << std::endl;

  // bottom boundary albedo
  if (pin->DoesParameterExist("radiation", myname + ".albedo"))
    ds.bc.albedo = pin->GetReal("disort", "albedo");
  else
    ds.bc.albedo = pin->GetOrAddReal("disort", "albedo", 0.);

  msg << "- albedo = " << ds.bc.albedo << std::endl;

  // top boundary isotropic radiation
  if (pin->DoesParameterExist("radiation", myname + ".fisot_K")) {
    Real Ttop = pin->GetReal("radiation", myname + ".fisot_K");
    ds.bc.fisot = Radiation::stefanBoltzmann*pow(Ttop, 4.);
  } else
    ds.bc.fisot = pin->GetOrAddReal("disort", "fisot", 0.);

  msg << "- fisot = " << ds.bc.fisot << " w/m^2" << std::endl;

  // top boundary direct beam
  if (pin->DoesParameterExist("radiation", myname + ".fbeam_K")) {
    Real Ttop = pin->GetReal("radiation", myname + ".fbeam_K");
    ds.bc.fbeam = Radiation::stefanBoltzmann*pow(Ttop, 4.);
  } else
    ds.bc.fbeam = pin->GetOrAddReal("disort", "fbeam", 0.);

  msg << "- fbeam = " << ds.bc.fbeam << " w/m^2" << std::endl;

  // top boundary emissivity
  if (pin->DoesParameterExist("radiation", myname + ".temis"))
    ds.bc.temis = pin->GetReal("radiation", myname + ".temis");
  else
    ds.bc.temis = pin->GetOrAddReal("disort", "temis", 0.);

  msg << "- temis = " << ds.bc.temis << std::endl;

  ds.flag.ibcnd = pin->GetOrAddBoolean("disort", "ibcnd", false);
  ds.flag.usrtau = pin->GetOrAddBoolean("disort", "usrtau", false);
  ds.flag.usrang = pin->GetOrAddBoolean("disort", "usrang",false);
  ds.flag.lamber = pin->GetOrAddBoolean("disort", "lamber", true);

  // calculate planck function
  if (pin->DoesParameterExist("radiation", myname + ".planck"))
    ds.flag.planck = pin->GetBoolean("radiation", myname + ".planck");
  else
    ds.flag.planck = pin->GetOrAddBoolean("disort", "planck", false);

  if (ds.flag.planck)
    msg << "- planck = true" << std::endl;
  else
    msg << "- planck = false" << std::endl;

  ds.flag.spher = pin->GetOrAddBoolean("disort", "spher", false);
  ds.flag.onlyfl = pin->GetOrAddBoolean("disort", "onlyfl", true);
  ds.flag.quiet = pin->GetOrAddBoolean("disort", "quiet", true);
  ds.flag.intensity_correction = 
    pin->GetOrAddBoolean("disort", "intensity_correction", true);
  ds.flag.old_intensity_correction = 
    pin->GetOrAddBoolean("disort", "old_intensity_correction", false);
  ds.flag.general_source = pin->GetOrAddBoolean("disort", "general_source", false);
  ds.flag.output_uum = pin->GetOrAddBoolean("disort", "output_uum", false);

  for (int i = 0; i < 5; ++i) ds.flag.prnt[i] = 0;
  std::string str;
  if (ds.flag.usrtau) {
    std::vector<Real> utau;
    str = pin->GetString("disort", "utau"); 
    utau = Vectorize<Real>(str.c_str());
    ds.ntau = utau.size();
    for (int i = 0; i < ds.ntau; ++i)
      ds.utau[i] = utau[i];
  }

  std::vector<Real> umu, phi;
  if (ds.flag.usrang) {
    str = pin->GetString("disort", "umu"); 
    umu = Vectorize<Real>(str.c_str());
    str = pin->GetString("disort", "phi"); 
    phi = Vectorize<Real>(str.c_str());
    ds.numu = umu.size();
    ds.nphi = phi.size();
  } else {
    ds.numu = 0;
    ds.nphi = 0;
  }

  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds, &ds_out);

  for (int i = 0; i < ds.numu; ++i)
    ds.umu[i] = umu[i];
  for (int i = 0; i < ds.nphi; ++i)
    ds.phi[i] = phi[i];

  pdbg->Leave();
}

void RadiationBand::free_disort()
{
  c_disort_state_free(&ds);
  c_disort_out_free(&ds, &ds_out);
}

void RadiationBand::RadtranRadiance(Direction const rin, Direction const *rout, int nrout,
  Real dist_au, int k, int j, int il, int iu)
{
  /* place holder for calculating radiance
  if (!ds.flag.onlyfl) {
    int count = 0;
    for (int j = 0; j < ds.nphi; ++j)
      for (int k = 0; k < ds.ntau; ++k)
        for (int l = 0; l < ds.numu; ++l, ++count)
          uu(n,j,k,l) = ds_out.uu[count];
  }*/
}

void RadiationBand::RadtranFlux(Direction const rin, Real dist_au, int k, int j, int il, int iu)
{
  MeshBlock *pmb = pmy_rad->pmy_block;
  std::stringstream msg;
  if (ds.flag.ibcnd != 0) {
    msg << "### FATAL ERROR in RadiationBand::RadtranFlux (disort): ";
    msg << "ibcnd != 0" << std::endl;
    ATHENA_ERROR(msg);
  }
  pmb->pcomm->setColor(X1DIR);

  int nblocks = pmb->pmy_mesh->mesh_size.nx1/pmb->block_size.nx1;
  Real *buf = new Real [(iu-il)*nblocks*(npmom+3)];
  if (ds.flag.planck) {
    pmb->pcomm->gatherData(temf_+il, buf, iu-il+1, X1DIR);
    for (int i = 0; i < (iu-il+1)*nblocks; ++i) {
      int m = i/(iu-il+1);
      ds.temper[m*(iu-il) + i%(iu-il+1)] = buf[i];
    }
    std::reverse(ds.temper, ds.temper + ds.nlyr+1);
  }

  ds.bc.umu0 = rin.mu > 1.E-3 ? rin.mu : 1.E-3;
  ds.bc.phi0 = rin.phi;
  if (ds.flag.planck) {
    ds.bc.btemp = ds.temper[ds.nlyr];
    ds.bc.ttemp = ds.temper[0];
  }

  // reset flx of this column 
  for (int i = il; i <= iu; ++i)
    bflxup(k,j,i) = bflxdn(k,j,i) = 0.;

  int r = pmb->pcomm->getRank(X1DIR);
  // loop over lines in the band
  for (int n = 0; n < nspec; ++n) {
    // stellar source function
    if (pmy_rad->rtype != RadiationType::band)
      ds.bc.fbeam = pmy_rad->planet->ParentInsolationFlux(spec[n].wav, dist_au);
    else {
      ds.bc.fbeam /= dist_au*dist_au;
      ds.bc.fisot /= dist_au*dist_au;
    }

    // planck source function
    ds.wvnmlo = spec[n].wav;
    ds.wvnmhi = spec[n].wav;

    // pack data
    int dsize = (npmom+3)*(iu-il);
    packSpectralProperties(buf+r*dsize, tau_[n]+il, ssa_[n]+il, pmom_[n][il], iu-il, npmom+1);
    pmb->pcomm->gatherDataInPlace(buf, dsize, X1DIR);
    unpackSpectralProperties(ds.dtauc, ds.ssalb, ds.pmom, buf, iu-il, npmom+1, nblocks, ds.nmom_nstr+1);

    // absorption
    std::reverse(ds.dtauc, ds.dtauc + ds.nlyr);

    // single scatering albedo
    std::reverse(ds.ssalb, ds.ssalb + ds.nlyr);

    // Legendre coefficients
    std::reverse(ds.pmom, ds.pmom + ds.nlyr*(ds.nmom_nstr+1));
    for (int i = 0; i < ds.nlyr; ++i)
      std::reverse(ds.pmom + i*(ds.nmom_nstr+1), ds.pmom + (i+1)*(ds.nmom_nstr+1));

    // run disort
    c_disort(&ds, &ds_out);

    // Counting index
    // Example, il = 0, iu = 2, ds.nlyr = 6, partition in to 3 blocks
    // face id   -> 0 - 1 - 2 - 3 - 4 - 5 - 6
    // cell id   -> | 0 | 1 | 2 | 3 | 4 | 5 |
    // disort id -> 6 - 5 - 4 - 3 - 2 - 1 - 0
    // blocks    -> ---------       *       *
    //           ->  r = 0  *       *       *
    //           ->         ---------       *
    //           ->           r = 1 *       *
    //           ->                 ---------
    //           ->                   r = 2
    // block r = 0 gets, 6 - 5 - 4
    // block r = 1 gets, 4 - 3 - 2
    // block r = 2 gets, 2 - 1 - 0
    // accumulate flux from lines
    for (int i = il; i <= iu; ++i) {
      int m = ds.nlyr - (r*(iu-il) + i - il);
      /*! \bug does not work for spherical geometry, need to scale area using
       * farea(il)/farea(i)
       */
      // flux up
      flxup_[n][i] = ds_out.rad[m].flup;

      /*! \bug does not work for spherical geomtry, need to scale area using
       * farea(il)/farea(i)
       */
      // flux down
      flxdn_[n][i] = ds_out.rad[m].rfldir + ds_out.rad[m].rfldn;
      /*if (Globals::my_rank == 0) {
        for (int i = 0; i <= ds.nlyr; ++i)
          std::cout << i << " " << flxup_[n][i] << " " << flxdn_[n][i] << std::endl;
      }*/
      bflxup(k,j,i) += spec[n].wgt*flxup_[n][i];
      bflxdn(k,j,i) += spec[n].wgt*flxdn_[n][i];
    }
  }
  delete[] buf;
}

#endif // RT_DISORT
