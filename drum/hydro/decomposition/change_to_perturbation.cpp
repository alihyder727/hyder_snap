// C/C++ headers
#include <sstream>
#include <iomanip>

// Athena++ headers
#include "decomposition.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../thermodynamics/thermodynamics.hpp"
#include "../../globals.hpp"
#include "../../utils/utils.hpp"

inline void IntegrateUpwards(AthenaArray<Real>& psf, AthenaArray<Real> const& w, Coordinates *pco,
  Real grav, int kl, int ku, int jl, int ju, int il, int iu)
{
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i)
        psf(k,j,i+1) = psf(k,j,i) - grav*w(IDN,k,j,i)*pco->dx1f(i);

}

inline void IntegrateDownwards(AthenaArray<Real>& psf, AthenaArray<Real> const& w, Coordinates *pco,
  Real grav, int kl, int ku, int jl, int ju, int il, int iu)
{
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = iu; i >= il; --i)
        psf(k,j,i) = psf(k,j,i+1) + grav*w(IDN,k,j,i)*pco->dx1f(i);
}

void Decomposition::ChangeToPerturbation(AthenaArray<Real> &w, int kl, int ku, int jl, int ju)
{
  // Need to integrate upwards
  MeshBlock *pmb = pmy_hydro->pmy_block;
  Coordinates *pco = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pthermo;

  Real grav = -pmy_hydro->hsrc.GetG1();  // positive downward pointing
  Real gamma = pmb->peos->GetGamma();

  int is = pmb->is, ie = pmb->ie;
  
  if (grav == 0.) return;

  std::stringstream msg;
  if (NGHOST < 3) {
    msg << "### FATAL ERROR in function [Hydro::DecomposePressure]"
        << std::endl << "number of ghost cells (NGHOST) must be at least 3" <<std::endl;
    ATHENA_ERROR(msg);
  }

  Real **w1;
  NewCArray(w1, 2, NHYDRO);
  FindNeighbors();

  if (has_top_neighbor) {
    RecvFromTop(psf_, entropy_, kl, ku, jl, ju);
    //std::cout << "Rank " << Globals::my_rank << " receive from top" << std::endl;
  } else {
    //std::cout << "Rank " << Globals::my_rank << " start from top" << std::endl;
    // top layer polytropic index and entropy
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        entropy_(k,j) = log(w(IPR,k,j,ie)) - gamma*log(w(IDN,k,j,ie));

    // adiabatic extrapolation
    for (int k=kl; k<=ku; ++k)
      for (int j=jl; j<=ju; ++j) {
        Real P1 = w(IPR,k,j,ie);
        Real T1 = pthermo->Temp(w.at(k,j,ie));
        Real dz = pco->dx1f(ie);
        for (int n = 0; n < NHYDRO; ++n)
          w1[0][n] = w(n,k,j,ie);

        // adiabatic extrapolation for half a grid
        pthermo->ConstructAdiabat(w1, T1, P1, grav, dz/2., 2, Adiabat::reversible);
        psf_(k,j,ie+1) = w1[1][IPR];

        // outflow boundary condition
        if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::outflow)
          for (int n = 0; n < NHYDRO; ++n)
            for (int i = 1; i <= NGHOST; ++i)
              w(n,k,j,ie+i) = w1[1][n];
      }
  }
  IntegrateDownwards(psf_, w, pco, grav, kl, ku, jl, ju, is, ie);
  
  //if (has_bot_neighbor)
  //  SendToBottom(psf_, entropy_, kl, ku, jl, ju);

  // integrate ghost cells
  if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::reflect)
    IntegrateDownwards(psf_, w, pco, -grav, kl, ku, jl, ju, is - NGHOST, is - 1);
  else if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::outflow)
    IntegrateDownwards(psf_, w, pco,  0., kl, ku, jl, ju, is - NGHOST, is - 1);
  else { // block boundary
    IntegrateDownwards(psf_, w, pco,  grav, kl, ku, jl, ju, is - NGHOST, is - 1);
    SendToBottom(psf_, entropy_, kl, ku, jl, ju);
  }

  if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::reflect) {
    IntegrateUpwards(psf_, w, pco, -grav, kl, ku, jl, ju, ie + 1, ie + NGHOST);
  } else if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::outflow)
    IntegrateUpwards(psf_, w, pco,  0., kl, ku, jl, ju, ie + 1, ie + NGHOST);
  //else  // block boundary
  //  IntegrateUpwards(psf_, w, pco,  grav, kl, ku, jl, ju, ie + 1, ie + NGHOST);

  // decompose pressure
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      // 1. change density and pressure (including ghost cells)
      for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
        // save pressure and density
        pres_(k,j,i) = w(IPR,k,j,i);
        dens_(k,j,i) = w(IDN,k,j,i);

        // interpolate hydrostatic pressure, prevent divided by zero
        Real psv, dsv;
        if (fabs(psf_(k,j,i) - psf_(k,j,i+1)) < 1.E-6)
          psv = (psf_(k,j,i) + psf_(k,j,i+1))/2.;
        else
          psv = (psf_(k,j,i) - psf_(k,j,i+1))/log(psf_(k,j,i)/psf_(k,j,i+1));
        dsv = pow(psv, 1./gamma)*exp(-entropy_(k,j)/gamma);

        // change pressure/density to pertubation quantities
        w(IPR,k,j,i) -= psv;
        w(IDN,k,j,i) -= dsv;
      }
    }

  /* debug
  if (Globals::my_rank == 1) {
    //int km = (kl + ku)/2;
    //int jm = (jl + ju)/2;
    int km = kl;
    int jm = jl;
    std::cout << "my.gid = " << pmb->gid << std::endl;
    std::cout << "bblock.gid = " << bblock.snb.gid << std::endl;
    std::cout << "===== k = " << km << " j = " << jm << std::endl;
    for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
      if (i == is)
        std::cout << "-------- ";
      if (i == 0)
        std::cout << "i = " << "-1/2 ";
      else if (i == 1)
        std::cout << "i = " << "+1/2 ";
      else
        std::cout << "i = " << i-1 << "+1/2 ";
      std::cout << "psf = " << std::setprecision(8) << psf_(km,jm,i) << std::endl;
      std::cout << "i = " << i  << "    ";
      //std::cout << "pre = " << w(IPR,kl,jl,i) << std::endl;
      std::cout << " pre = " << w(IPR,km,jm,i) << " den = " << w(IDN,km,jm,i) << std::endl;
      if (i == ie)
        std::cout << "-------- ";
      if (i == ie + NGHOST) {
        std::cout << "i = " << i+1 << "+1/2 ";
        std::cout << "psf = " << psf_(km,jm,i+1) << std::endl;
      }
    }
    std::cout << "==========" << std::endl;
  }*/

  FreeCArray(w1);

  // finish send top pressure
  WaitToFinishSend();
}

void Decomposition::RestoreFromPerturbation(AthenaArray<Real> &w,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr, int k, int j, int il, int iu)
{
  MeshBlock *pmb = pmy_hydro->pmy_block;
  Real gamma = pmb->peos->GetGamma();
  int is = pmb->is, ie = pmb->ie;
  if (pmy_hydro->hsrc.GetG1() == 0.) return;
  
  for (int i = is - NGHOST; i <= ie + NGHOST; ++i) {
    w(IPR,k,j,i) = pres_(k,j,i);
    w(IDN,k,j,i) = dens_(k,j,i);
  }

  for (int i = il; i <= iu; ++i) {
    wr(IPR,i) += psf_(k,j,i);
    if (wr(IPR,i) < 0.) wr(IPR,i) = psf_(k,j,i);

    wl(IPR,i) += psf_(k,j,i);
    if (wl(IPR,i) < 0.) wl(IPR,i) = psf_(k,j,i);

    wr(IDN,i) += pow(psf_(k,j,i), 1./gamma)*exp(-entropy_(k,j)/gamma);
    if (wr(IDN,i) < 0.)
      wr(IDN,i) = pow(psf_(k,j,i), 1./gamma)*exp(-entropy_(k,j)/gamma);

    wl(IDN,i) += pow(psf_(k,j,i), 1./gamma)*exp(-entropy_(k,j)/gamma);
    if (wl(IDN,i) < 0.)
      wl(IDN,i) = pow(psf_(k,j,i), 1./gamma)*exp(-entropy_(k,j)/gamma);
  }

  /* debug
  //int km = (pmb->ks + pmb->ke)/2, jm = (pmb->js + pmb->je)/2;
  int km = pmb->ks-1, jm = pmb->js-1;
  if (Globals::my_rank == 1 && km == k & jm == j) {
    std::cout << "my.gid = " << pmb->gid << std::endl;
    std::cout << "bblock.gid = " << bblock.snb.gid << std::endl;
    std::cout << "===== k = " << km << " j = " << jm << std::endl;
    for (int i = il; i <= iu; ++i) {
      if (i == is)
        std::cout << "-------- ";
      if (i == 0)
        std::cout << "i = " << "-1/2 ";
      else if (i == 1)
        std::cout << "i = " << "+1/2 ";
      else
        std::cout << "i = " << i-1 << "+1/2 ";
      std::cout << "psf = " << psf_(km,jm,i) 
                << " wl(IPR) = " << wl(IPR,i)
                << " wr(IPR) = " << wr(IPR,i)
                << " wl(IDN) = " << wl(IDN,i)
                << " wr(IDN) = " << wr(IDN,i)
                << std::endl;
      std::cout << "i = " << i  << "    ";
      std::cout << " pre = " << w(IPR,km,jm,i) << " den = " << w(IDN,km,jm,i) << std::endl;
      if (i == ie)
        std::cout << "-------- ";
      if (i == ie + NGHOST) {
        std::cout << "i = " << i+1 << "+1/2 ";
        std::cout << "psf = " << psf_(km,jm,i+1) << std::endl;
      }
    }
    std::cout << "==========" << std::endl;
  }*/
}
