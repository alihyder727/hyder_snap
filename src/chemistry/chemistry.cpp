// Athena header files
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../reconstruct/interpolation.hpp"
#include "../utils/utils.hpp"
#include "chemistry.hpp"

Chemistry::Chemistry(MeshBlock *pmb, ParameterInput *pin)
{
  /*pmy_block_ = pmb;
  NewCArray(identity_, NMASS, NMASS);
  NewCArray(r1_, NMASS, NMASS);
  r0_ = new Real [NMASS];

  std::fill(*identity_, *identity_ + NMASS*NMASS, 0.);
  std::fill(*r1_, *r1_ + NMASS*NMASS, 0.);
  std::fill(r0_, r0_ + NMASS, 0.);
  for (int n = 0; n < NMASS; ++n)
    identity_[n][n] = 1.;

  min_dt = pin->GetOrAddReal("chemistry", "min_dt", 1.);
  last_time = pmb->pmy_mesh->time;

  char buf[80];
  int numk = pin->GetOrAddInteger("chemistry", "num_coeffs", 0) + 1;
  kc_.resize(numk);
  for (int i = 0; i < numk; ++i) {
    sprintf(buf, "k%d", i);
    kc_[i] = pin->GetOrAddReal("chemistry", buf, 0.);
  }

  order_ = pin->GetOrAddInteger("chemistry", "order", 1);
  max_iter_ = pin->GetOrAddInteger("chemistry", "max_iter", 5);
  ftol_ = pin->GetOrAddReal("chemistry", "ftol", 1.0E-10);
  gamma_ = 2. - sqrt(2.);

  // sedimentation velocity
  for (int n = 0; n <= NVAPOR; ++n) vsed_default_[n] = 0.;

  for (int n = 1+NVAPOR; n < NMASS; ++n) {
    sprintf(buf, "vsed%d", n);
    vsed_default_[n] = pin->GetOrAddReal("chemistry", buf, 0.);
  }*/
}

Chemistry::~Chemistry()
{
  //FreeCArray(identity_);
  //FreeCArray(r1_);
  //delete[] r0_;
}

void Chemistry::AddSedimentationFlux(AthenaArray<Real>& x1flux,
  AthenaArray<Real> const& wr, int k, int j, int il, int iu)
{
  /*MeshBlock *pmb = pmy_block_;
  Thermodynamics *pthermo = pmb->pthermo;

  Real vsed[NMASS], w1[NHYDRO];
  for (int i = il; i <= iu; ++i) {
    Real temp = pthermo->Temp(wr.at(i));
    for (int n = 0; n < NHYDRO; ++n)
      w1[n] = wr(n,i);
    SedimentationVelocity(vsed, w1, temp);

    for (int n = 1+NVAPOR; n < NMASS; ++n) {
      Real rho = w1[IDN]*w1[n];
      Real v1 = w1[IVX];
      Real v2 = w1[IVY];
      Real v3 = w1[IVZ];
      Real en = pthermo->GetCv(n)*temp + 0.5*(v1*v1 + v2*v2 + v3*v3)
                - pthermo->GetLatent(n);
      x1flux(n,k,j,i) += vsed[n]*rho;
      x1flux(IVX,k,j,i) += vsed[n]*rho*v1;
      x1flux(IVY,k,j,i) += vsed[n]*rho*v2;
      x1flux(IVZ,k,j,i) += vsed[n]*rho*v3;
      x1flux(IEN,k,j,i) += vsed[n]*rho*en;
    }
  }*/
}

void Chemistry::AddFrictionalHeating(AthenaArray<Real> &u, 
  AthenaArray<Real> const& w, Real dt)
{
  /*MeshBlock *pmb = pmy_block_;
  Real grav = pmb->phydro->hsrc.GetG1();

  Real vsed[NMASS], w1[NHYDRO];
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        for (int n = 0; n < NHYDRO; ++n)
          w1[n] = w(n,k,j,i);
        SedimentationVelocity(vsed, w1);
        for (int n = 1+NVAPOR; n < NMASS; ++n)
          u(IEN,k,j,i) += dt*w1[IDN]*w1[n]*vsed[n]*grav;
      }*/
}

void Chemistry::AssembleReactionMatrix(Real *r0, Real **r1, Real const q[], Real time)
{
  //std::fill(*r1_, *r1_ + NMASS*NMASS, 0.);
  //std::fill(r0_, r0_ + NMASS, 0.);
}

void Chemistry::SedimentationVelocity(Real vsed[], Real const w[], Real temp)
{
  //for (int n = 1+NVAPOR; n < NMASS; ++n)
  //  vsed[n] = vsed_default_[n];
}
