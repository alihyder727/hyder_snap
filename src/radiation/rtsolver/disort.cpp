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
#include "../radiation.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../communicator/communicator.hpp"
#include "../../globals.hpp"

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

  ds->nlyr = pmy_rad->pmy_block->pmy_mesh->mesh_size.nx1;
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
  MeshBlock *pmb = pmy_rad->pmy_block;
  std::stringstream msg;
  if (ds->flag.ibcnd != 0) {
    msg << "### FATAL ERROR in RadiationBand::RadtranFlux (disort): ";
    msg << "ibcnd != 0" << std::endl;
    ATHENA_ERROR(msg);
  }

  //AthenaArray<Real> farea(iu+1);
  //pmy_rad->pmy_block->pcoord->Face1Area(k, j, il, iu, farea);
  
  int nblocks = pmb->pmy_mesh->mesh_size.nx1/pmb->block_size.nx1;
  Real *buf = new Real [(iu-il+1)*nblocks*(ds->nmom_nstr+1)];
  if (ds->flag.planck) {
    pmb->pcomm->gatherData(temf_+il, buf, iu-il+1, X1DIR);
    for (int i = 0; i < (iu-il+1)*nblocks; ++i) {
      int m = i/(iu-il+1);
      ds->temper[m*(iu-il) + i%(iu-il+1)] = buf[i];
    }
    std::reverse(ds->temper, ds->temper + ds->nlyr+1);
  }

  ds->bc.umu0 = rin.mu > 1.E-3 ? rin.mu : 1.E-3;
  ds->bc.phi0 = rin.phi;
  if (ds->flag.planck) {
    ds->bc.btemp = ds->temper[ds->nlyr];
    ds->bc.ttemp = ds->temper[0];
  }

  // reset flx of this column 
  for (int i = il; i <= iu; ++i)
    bflxup(k,j,i) = bflxdn(k,j,i) = 0.;

  pmb->pcomm->setColor(X1DIR);
  int r = pmb->pcomm->getRank(X1DIR);
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
      buf[i-il] = tau_[i][n];
    pmb->pcomm->gatherData(buf, ds->dtauc, iu-il, X1DIR);
    std::reverse(ds->dtauc, ds->dtauc + ds->nlyr);
    /*if (Globals::my_rank == 0)
      for (int i = 0; i < ds->nlyr; ++i)
        std::cout << i << " " << ds->dtauc[i] << std::endl;*/

    // single scatering albedo
#pragma omp simd
    for (int i = il; i < iu; ++i)
      buf[i-il] = ssa_[i][n];
    pmb->pcomm->gatherData(buf, ds->ssalb, iu-il, X1DIR);
    std::reverse(ds->ssalb, ds->ssalb + ds->nlyr);

    //! \bug npmom and nmom_nstr may not be consistent
    // Legendre coefficients
#pragma omp simd
    for (int i = il; i < iu; ++i)
      for (int p = 0; p <= npmom; ++p)
        buf[(i-il)*(ds->nmom_nstr+1) + p] = pmom_[i][n][p];
    pmb->pcomm->gatherData(buf, ds->pmom, (iu-il)*(ds->nmom_nstr+1), X1DIR);
    std::reverse(ds->pmom, ds->pmom + ds->nlyr*(ds->nmom_nstr+1));
    for (int i = 0; i < ds->nlyr; ++i)
      std::reverse(ds->pmom + i*(ds->nmom_nstr+1), ds->pmom + (i+1)*(ds->nmom_nstr+1));

    // run disort
    c_disort(ds, ds_out);

    // Counting index
    // Example, il = 0, iu = 2, ds->nlyr = 6, partition in to 3 blocks
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
      int m = ds->nlyr - (r*(iu-il) + i - il);
      /*! \bug does not work for spherical geometry, need to scale area using
       * farea(il)/farea(i)
       */
      // flux up
      flxup_[i][n] = ds_out->rad[m].flup;

      /*! \bug does not work for spherical geomtry, need to scale area using
       * farea(il)/farea(i)
       */
      // flux down
      flxdn_[i][n] = ds_out->rad[m].rfldir + ds_out->rad[m].rfldn;
      bflxup(k,j,i) += spec[n].wgt*flxup_[i][n];
      bflxdn(k,j,i) += spec[n].wgt*flxdn_[i][n];
    }
  }

  delete[] buf;
}

#endif // RT_DISORT
