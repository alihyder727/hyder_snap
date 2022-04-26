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

  // number of legengre moments
  if (pin->DoesParameterExist("radiation", myname + ".npmom"))
    ds.nmom = pin->GetInteger("radiation", myname + ".npmom");
  else
    ds.nmom = pin->GetInteger("radiation", "npmom");

  // number of streams
  if (pin->DoesParameterExist("radiation", myname + ".nstr"))
    ds.nstr = pin->GetInteger("radiation", myname + ".nstr");
  else
    ds.nstr = pin->GetInteger("radiation", "nstr");

  // number of phases
  if (pin->DoesParameterExist("radiation", myname + ".nphase"))
    ds.nphase = pin->GetInteger("radiation", myname + ".nphase");
  else
    ds.nphase = pin->GetInteger("radiation", "nphase");

  // accuracy
  if (pin->DoesParameterExist("radiation", myname + ".accur"))
    ds.accur  = pin->GetReal("radiation", myname + ".accur");
  else
    ds.accur  = pin->GetOrAddReal("radiation", "accur", 0.);

  // bottom boundary isotropic radiation
  if (pin->DoesParameterExist("radiation", myname + ".fluor_K")) {
    Real Tbot = pin->GetReal("radiation", myname + ".fluor_K");
    ds.bc.fluor = Radiation::stefanBoltzmann*pow(Tbot, 4.);
  } else
    ds.bc.fluor = pin->GetOrAddReal("radiation", "fluor", 0.);

  msg << "- fluor = " << ds.bc.fluor << " w/m^2" << std::endl;

  // bottom boundary albedo
  if (pin->DoesParameterExist("radiation", myname + ".albedo"))
    ds.bc.albedo = pin->GetReal("radiation", "albedo");
  else
    ds.bc.albedo = pin->GetOrAddReal("radiation", "albedo", 0.);

  msg << "- albedo = " << ds.bc.albedo << std::endl;

  // top boundary isotropic radiation
  if (pin->DoesParameterExist("radiation", myname + ".fisot_K")) {
    Real Ttop = pin->GetReal("radiation", myname + ".fisot_K");
    ds.bc.fisot = Radiation::stefanBoltzmann*pow(Ttop, 4.);
  } else
    ds.bc.fisot = pin->GetOrAddReal("radiation", "fisot", 0.);

  msg << "- fisot = " << ds.bc.fisot << " w/m^2" << std::endl;

  // top boundary direct beam
  if (pin->DoesParameterExist("radiation", myname + ".fbeam_K")) {
    Real Ttop = pin->GetReal("radiation", myname + ".fbeam_K");
    ds.bc.fbeam = Radiation::stefanBoltzmann*pow(Ttop, 4.);
  } else
    ds.bc.fbeam = pin->GetOrAddReal("radiation", "fbeam", 0.);

  msg << "- fbeam = " << ds.bc.fbeam << " w/m^2" << std::endl;

  // top boundary emissivity
  if (pin->DoesParameterExist("radiation", myname + ".temis"))
    ds.bc.temis = pin->GetReal("radiation", myname + ".temis");
  else
    ds.bc.temis = pin->GetOrAddReal("radiation", "temis", 0.);

  msg << "- temis = " << ds.bc.temis << std::endl;

  ds.flag.ibcnd = pin->GetOrAddBoolean("disort", "ibcnd", false);
  ds.flag.usrtau = pin->GetOrAddBoolean("disort", "usrtau", false);
  ds.flag.usrang = pin->GetOrAddBoolean("disort", "usrang",false);
  ds.flag.lamber = pin->GetOrAddBoolean("disort", "lamber", true);

  // calculate planck function
  ds.flag.planck = (bflags & RadiationFlags::Planck) > 0 ? true : false;
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

void RadiationBand::calculateRadiance(Direction ray, Real dist_au,
    int k, int j, int il, int iu)
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

void RadiationBand::calculateRadiativeFlux(Direction const ray, Real dist_au,
    int k, int j, int il, int iu)
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
  Real *buf = new Real [(iu-il)*nblocks*(ds.nmom_nstr+4)];
  if (ds.flag.planck) {
    pmb->pcomm->gatherData(temf_+il, buf, iu-il+1);
    for (int i = 0; i < (iu-il+1)*nblocks; ++i) {
      int m = i/(iu-il+1);
      ds.temper[m*(iu-il) + i%(iu-il+1)] = buf[i];
    }
    std::reverse(ds.temper, ds.temper + ds.nlyr+1);
  }
  //for (int i = 0; i <= ds.nlyr; ++i)
  //  std::cout << ds.temper[i] << std::endl;

  ds.bc.umu0 = ray.mu > 1.E-3 ? ray.mu : 1.E-3;
  ds.bc.phi0 = ray.phi;
  if (ds.flag.planck) {
    ds.bc.btemp = ds.temper[ds.nlyr];
    ds.bc.ttemp = ds.temper[0];
  }

  // reset flx of this column 
  for (int i = il; i <= iu; ++i)
    bflxup(k,j,i) = bflxdn(k,j,i) = 0.;

  Coordinates *pcoord = pmy_rad->pmy_block->pcoord;
  AthenaArray<Real> farea(iu+1), vol(iu+1);
  pcoord->Face1Area(k, j, il, iu, farea);
  pcoord->CellVolume(k, j, il, iu, vol);

  if (bflags & RadiationFlags::CorrelatedK) {
    // stellar source function
    if (bflags & RadiationFlags::Star)
      ds.bc.fbeam = pmy_rad->planet->ParentInsolationFlux(wmin, wmax, dist_au);
    // planck source function
    ds.wvnmlo = wmin;
    ds.wvnmhi = wmax;
  }

  int r = pmb->pcomm->getRank(X1DIR);
  // loop over bins in the band
  for (int n = 0; n < num_bins; ++n) {
    if (!(bflags & RadiationFlags::CorrelatedK)) {
      // stellar source function
      if (bflags & RadiationFlags::Star)
        ds.bc.fbeam = pmy_rad->planet->ParentInsolationFlux(spec[n].wav1, spec[n].wav2, dist_au);
      // planck source function
      ds.wvnmlo = spec[n].wav1;
      ds.wvnmhi = spec[n].wav2;
    }

    // pack data
    int npmom = bpmom.GetDim4() - 1;
    int dsize = (npmom+3)*(iu-il);
    std::fill(buf + r*dsize, buf + (r+1)*dsize, 0.);
    packSpectralProperties(buf+r*dsize, tau_[n]+il, ssa_[n]+il, pmom_[n][il], iu-il, npmom+1);
    pmb->pcomm->gatherDataInPlace(buf, dsize);
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
      bflxup(k,j,i) += spec[n].wght*flxup_[n][i];
      bflxdn(k,j,i) += spec[n].wght*flxdn_[n][i];
    }

    // spherical correction by XIZ
    // xiz 2022 flux scaling so that the heating rate is the same as the plane-parallel scheme
    // volheating scaling: first calculate flux divergence from DISORT using Plane-parallel in a cell
    // then mulitpled by the cell volume divided by dx1f
    // then solve for F using F1*S1-F2*S2 = volheating
    // the top fluxes are the still the same as the plane-paralell values
    Real volh, bflxup1 = bflxup(k,j,iu), bflxdn1 = bflxdn(k,j,iu);
    for (int i = iu-1; i >= il; --i) {
      // upward
      volh = (bflxup1 - bflxup(k,j,i))/pcoord->dx1f(i)*vol(i);
      bflxup1 = bflxup(k,j,i);
      bflxup(k,j,i) = (bflxup(k,j,i+1)*farea(i+1) - volh)/farea(i);

      // downward
      volh = (bflxdn1 - bflxdn(k,j,i))/pcoord->dx1f(i)*vol(i);
      bflxdn1 = bflxdn(k,j,i);
      bflxdn(k,j,i) = (bflxdn(k,j,i+1)*farea(i+1) - volh)/farea(i);
    }
  }
  delete[] buf;
}

#endif // RT_DISORT
