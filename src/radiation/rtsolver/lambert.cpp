//! \file lambert.cpp
//  \brief Beer-Lambert Radiative Transfer
//=====================================================

// C/C++ headers
#include <cmath>

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
  Real *temf = new Real [nlevels];  // temperature at cell faces

  temf[il] = interp_weno5(tem_[il+1], tem_[il], tem_[il-1], tem_[il-2], tem_[il-3]);
  for (int i = il+1; i < iu; ++i)
    temf[i] = interp_cp6(tem_[i-3], tem_[i-2], tem_[i-1], tem_[i], tem_[i+1], tem_[i+2]);
  temf[iu] = interp_weno5(tem_[iu-2], tem_[iu-1], tem_[iu], tem_[iu+1], tem_[iu+2]);
  //for (int i = il; i <= iu; ++i)
    //std::cout << i << " " << temf[i] << std::endl;

  // integrate from top to bottom
  Real dtau;
  for (int m = 0; m < nrout; ++m) {
    btoa(m,k,j) = 0.;
    for (int n = 0; n < nspec; ++n) {
      taut[iu] = 0.;
      dtau = interp_weno5(tau_[iu-2][n], tau_[iu-1][n], tau_[iu][n], tau_[iu+1][n], tau_[iu+2][n]);
      toa_[m][n] = temf[iu]*dtau/rout[m].mu;
      for (int i = iu-1; i >= il; --i) {
        taut[i] = taut[i+1] + tau_[i][n]/rout[m].mu;
        dtau = interp_cp6(tau_[i-3][n], tau_[i-2][n], tau_[i-1][n], 
                          tau_[i][n], tau_[i+1][n], tau_[i+2][n]);
        toa_[m][n] += temf[i]*exp(-taut[i])*dtau/rout[m].mu;
        //if (m == 0)
        //  std::cout << taut[i] << " " << temf[i] << " " << temf[i]*exp(-taut[i]) << " " << tau_[i][n] << std::endl;
      }
      //  std::cout << taut[il] << " " << gammq(alpha_, taut[il]) << std::endl;
      toa_[m][n] += temf[il]*exp(-taut[il]);
      if ((alpha_ > 0) && (taut[il] < 1000.)) // correction for small optical opacity
        toa_[m][n] += temf[il]*alpha_*gammq(alpha_, taut[il])*pow(taut[il], -alpha_)*tgamma(alpha_);
      btoa(m,k,j) += spec[n].wgt*toa_[m][n];
    }
  }

  delete[] taut;
  delete[] temf;
}

#endif  // RT_LAMBERT
