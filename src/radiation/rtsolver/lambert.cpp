//! \file lambert.cpp
//  \brief Beer-Lambert Radiative Transfer
//=====================================================

// C/C++ headers
#include <cmath>
#include <iostream>

// Athena++ headers
#include "../../mesh/mesh.hpp"
#include "../../reconstruct/interpolation.hpp"
#include "../../math/special.h"
#include "../radiation.hpp"

#ifdef RT_LAMBERT

void RadiationBand::RadtranRadiance(Direction const rin, Direction const *rout, 
  int nrout, Real dist, int k, int j, int il, int iu)
{
  MeshBlock *pmb = pmy_rad->pmy_block;
  Coordinates *pcoord = pmb->pcoord;

  int nlevels = pmb->ncells1;

  Real *taut = new Real [nlevels];
  /*Real *temf = new Real [nlevels];  // temperature at cell faces

  temf[il] = interp_weno5(tem_[il-2], tem_[il-1], tem_[il], tem_[il+1], tem_[il+2]);
  for (int i = il+1; i < iu; ++i)
    temf[i] = interp_cp4(tem_[i-2], tem_[i-1], tem_[i], tem_[i+1]);
  temf[iu] = interp_weno5(tem_[iu+1], tem_[iu], tem_[iu-1], tem_[iu-2], tem_[iu-3]);*/
  //for (int i = il; i <= iu; ++i)
  //  std::cout << i << " " << temf[i] << " " << tem_[i] << " " << temf[i+1] << std::endl;

  // integrate from top to bottom
  for (int m = 0; m < nrout; ++m) {
    btoa(m,k,j) = 0.;
    for (int n = 0; n < nspec; ++n) {
      taut[iu] = 0.;
      toa_[n][m] = 0.;
      for (int i = iu-1; i >= il; --i) {
        taut[i] = taut[i+1] + tau_[n][i]/rout[m].mu;
        toa_[n][m] += 0.5*(temf_[i+1]*exp(-taut[i+1]) + temf_[i]*exp(-taut[i]))*tau_[n][i]/rout[m].mu;
        //if (m == 0)
        //  std::cout << taut[i] << " " << temf_[i] << " " << temf_[i]*exp(-taut[i]) << " " << tau_[n][i] << std::endl;
      }
      //  std::cout << taut[il] << " " << gammq(alpha_, taut[il]) << std::endl;
      toa_[n][m] += temf_[il]*exp(-taut[il]);
      if ((alpha_ > 0) && (taut[il] < 1000.)) // correction for small optical opacity
        toa_[n][m] += temf_[il]*alpha_*gammq(alpha_, taut[il])*pow(taut[il], -alpha_)*tgamma(alpha_);
      btoa(m,k,j) += spec[n].wgt*toa_[n][m];
    }
  }

  delete[] taut;
  //delete[] temf;
}

#endif  // RT_LAMBERT
