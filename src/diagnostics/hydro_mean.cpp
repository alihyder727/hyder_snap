#include "diagnostics.hpp"
#include "../reconstruct/interpolation.hpp"
#include "../coordinates/coordinates.hpp"
#include "../thermodynamics/thermodynamics.hpp"

HydroMean::HydroMean(MeshBlock *pmb):Diagnostics(pmb, "mean")
{
  type = "VECTORS";
  data.NewAthenaArray(NHYDRO,ncells3_,ncells2_,ncells1_);
}

HydroMean::~HydroMean() {
  data.DeleteAthenaArray();
}

void HydroMean::Progress(AthenaArray<Real> const& w)
{
  MeshBlock *pmb = pmy_block_;

  if (ncycle == 0)
    std::fill(data.data(), data.data() + data.GetSize(), 0.);

  for (int n = 0; n < IPR; ++n)
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i)
          data(n,k,j,i) += w(n,k,j,i);

  // Temperature field is save in IPR
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i)
        data(IPR,k,j,i) += pmb->pthermo->GetTemp(w.at(k,j,i));

  ncycle++;
}

void HydroMean::Finalize(AthenaArray<Real> const& w)
{
  MeshBlock *pmb = pmy_block_;

  // hydro mean
  if (ncycle > 0) {
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = pmb->ks; k <= pmb->ke; ++k)
        for (int j = pmb->js; j <= pmb->je; ++j)
          for (int i = pmb->is; i <= pmb->ie; ++i)
            data(n,k,j,i) /= ncycle;
  } else {
    for (int n = 0; n < IPR; ++n)
      for (int k = pmb->ks; k <= pmb->ke; ++k)
        for (int j = pmb->js; j <= pmb->je; ++j)
          for (int i = pmb->is; i <= pmb->ie; ++i)
            data(n,k,j,i) = w(n,k,j,i);

    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i)
          data(IPR,k,j,i) = pmb->pthermo->GetTemp(w.at(k,j,i));
  }

  // clear cycle
  ncycle = 0;
}
