//! \file disort.cpp
//  \brief DISORT radiative transfer solver
//=====================================================

// C/C++ header
#include <cstdlib>
#include <sstream>

// Athena++ headers
#include "../../reconstruct/interpolation.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../../utils/utils.hpp" // StringToArray, ReadTabular
#include "../radiation.hpp"
#include "../../coordinates/coordinates.hpp"

// DISORT headers

#ifdef RT_DISORT

#undef SQR
extern "C" {
  #include "cdisort213/cdisort.h"
}

void RadiationBand::init_disort(ParameterInput *pin)
{
  ds = (disort_state*)malloc(sizeof(disort_state));
  ds_out = (disort_output*)malloc(sizeof(disort_output));

  ds->nlyr = pmy_rad->pmy_block->block_size.nx1;
  ds->nmom = pin->GetInteger("radiation", "npmom");

  ds->nstr   = pin->GetInteger("disort", "nstr");
  ds->nphase = pin->GetInteger("disort", "nphase");
  ds->accur  = pin->GetOrAddReal("disort", "accur", 0.);

  ds->bc.fisot = pin->GetOrAddReal("disort", "fisot", 0.);
  ds->bc.albedo = pin->GetReal("disort", "albedo");
  ds->bc.temis = pin->GetReal("disort", "temis");

  ds->flag.ibcnd = pin->GetOrAddBoolean("disort", "ibcnd", false);
  ds->flag.usrtau = pin->GetOrAddBoolean("disort", "usrtau", false);
  ds->flag.usrang = pin->GetOrAddBoolean("disort", "usrang",false);
  ds->flag.lamber = pin->GetOrAddBoolean("disort", "lamber", true);
  ds->flag.planck = pin->GetOrAddBoolean("disort", "planck", false);
  ds->flag.spher = pin->GetOrAddBoolean("disort", "spher", false);
  ds->flag.onlyfl = pin->GetOrAddBoolean("disort", "onlyfl", true);
  ds->flag.quiet = pin->GetOrAddBoolean("disort", "quiet", true);
  ds->flag.intensity_correction = 
    pin->GetOrAddBoolean("disort", "intensity_correction", true);
  ds->flag.old_intensity_correction = 
    pin->GetOrAddBoolean("disort", "old_intensity_correction", false);
  ds->flag.general_source = pin->GetOrAddBoolean("disort", "general_source", false);
  ds->flag.output_uum = pin->GetOrAddBoolean("disort", "output_uum", false);

  for (int i = 0; i < 5; ++i) ds->flag.prnt[i] = 0;
  std::string str;
  if (ds->flag.usrtau) {
    std::vector<Real> utau;
    str = pin->GetString("disort", "utau"); 
    utau = Vectorize<Real>(str.c_str());
    ds->ntau = utau.size();
    for (int i = 0; i < ds->ntau; ++i)
      ds->utau[i] = utau[i];
  }

  std::vector<Real> umu, phi;
  if (ds->flag.usrang) {
    str = pin->GetString("disort", "umu"); 
    umu = Vectorize<Real>(str.c_str());
    str = pin->GetString("disort", "phi"); 
    phi = Vectorize<Real>(str.c_str());
    ds->numu = umu.size();
    ds->nphi = phi.size();
  } else {
    ds->numu = 0;
    ds->nphi = 0;
  }

  c_disort_state_alloc(ds);
  c_disort_out_alloc(ds, ds_out);

  for (int i = 0; i < ds->numu; ++i)
    ds->umu[i] = umu[i];
  for (int i = 0; i < ds->nphi; ++i)
    ds->phi[i] = phi[i];
}

void RadiationBand::free_disort()
{
  c_disort_state_free(ds);
  c_disort_out_free(ds, ds_out);
  free(ds);
  free(ds_out);
}


void RadiationBand::RadtranRadiance(Direction const rin, Direction const *rout, int nrout,
  Real dist, int k, int j, int il, int iu)
{
  /* place holder for calculating radiance
  if (!ds->flag.onlyfl) {
    int count = 0;
    for (int j = 0; j < ds->nphi; ++j)
      for (int k = 0; k < ds->ntau; ++k)
        for (int l = 0; l < ds->numu; ++l, ++count)
          uu(n,j,k,l) = ds_out->uu[count];
  }*/
}

void RadiationBand::RadtranFlux(Direction const rin, Real dist, int k, int j, int il, int iu)
{
  std::stringstream msg;
  if (ds->flag.ibcnd != 0) {
    msg << "### FATAL ERROR in RadiationBand::RadiativeTransfer (disort): ";
    msg << "ibcnd != 0" << std::endl;
    ATHENA_ERROR(msg);
  }

  AthenaArray<Real> farea(iu+1);
  pmy_rad->pmy_block->pcoord->Face1Area(k, j, il, iu, farea);
   
  if (ds->flag.planck) {
    ds->temper[iu-il] = interp_weno5(tem_[il+1], tem_[il], tem_[il-1], tem_[il-2], tem_[il-3]);
    for (int i = il+1; i < iu; ++i)
      ds->temper[iu-i] = interp_cp6(tem_[i+2], tem_[i+1], tem_[i], tem_[i-1], tem_[i-2], tem_[i-3]);
    ds->temper[0] = interp_weno5(tem_[iu-2], tem_[iu-1], tem_[iu], tem_[iu+1], tem_[iu+2]);
  }

  ds->bc.umu0 = rin.mu > 1.E-3 ? rin.mu : 1.E-3;
  ds->bc.phi0 = rin.phi;
  if (ds->flag.planck) {
    ds->bc.btemp = ds->temper[iu-il];
    ds->bc.ttemp = ds->temper[0];
  }

  // reset flx of this column 
  for (int i = il; i <= iu; ++i)
    bflxup(k,j,i) = bflxdn(k,j,i) = 0.;

  // loop over lines in the band
  for (int n = 0; n < nspec; ++n) {
    // stellar source function
    if (pmy_rad->GetBeam() < 0.)
      ds->bc.fbeam = pmy_rad->planet->ParentInsolationFlux(spec[n].wav, dist);
    else
      ds->bc.fbeam = pmy_rad->GetBeam();

    // planck source function
    ds->wvnmlo = spec[n].wav;
    ds->wvnmhi = spec[n].wav;

    // absorption
#pragma omp simd
    for (int i = il; i < iu; ++i)
      ds->dtauc[iu-1-i] = tau_[i][n];

    // single scatering albedo
#pragma omp simd
    for (int i = il; i < iu; ++i)
      ds->ssalb[iu-1-i] = ssa_[i][n];

    // Legendre coefficients
#pragma omp simd
    for (int i = il; i < iu; ++i)
      for (int p = 0; p <= npmom; ++p)
        ds->pmom[(iu-1-i)*(ds->nmom_nstr + 1) + p] = pmom_[i][n][p];

    // run disort
    c_disort(ds, ds_out);

    // accumulate flux from lines
    for (int i = il; i <= iu; ++i) {
      flxup_[i][n] = ds_out->rad[iu-i].flup*farea(il)/farea(i);   // flux up
      flxdn_[i][n] = (ds_out->rad[iu-i].rfldir + ds_out->rad[iu-i].rfldn)
                     *farea(il)/farea(i);  // flux down
      bflxup(k,j,i) += spec[n].wgt*flxup_[i][n];
      bflxdn(k,j,i) += spec[n].wgt*flxdn_[i][n];
    }
  }
}

#endif // RT_DISORT
