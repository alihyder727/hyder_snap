// Athena++ headers
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "../math/core.h"
#include "../debugger/debugger.hpp"
#include "../globals.hpp"
#include "physics.hpp"

TaskStatus Physics::TopSpongeLayer(AthenaArray<Real> &du,
  AthenaArray<Real> const& w, Real time, Real dt)
{
  MeshBlock *pmb = pmy_block;
  pmb->pdebug->Call("Physics::TopSpongerLayer");

  Coordinates *pcoord = pmb->pcoord;

  int iu = pmb->ie;
  Real x1min = pmb->pmy_mesh->mesh_size.x1min;
  Real x1max = pmb->pmy_mesh->mesh_size.x1max;

  while (iu >= pmb->is && (x1max - pcoord->x1v(iu) < width_top_)) iu--;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = iu+1; i <= pmb->ie; ++i) {
        Real eta = (pcoord->x1v(i) - x1min)/(x1max - x1min);
        Real eta_sp = 1. - width_top_/(x1max - x1min);
        Real scale = eta > eta_sp ? sqr(sin(M_PI/2.*(eta - eta_sp)/(1. - eta_sp))) : 0.;
        du(IVX,k,j,i) -= w(IDN,k,j,i)*w(IVX,k,j,i)/tau_top_*scale*dt;
        du(IVY,k,j,i) -= w(IDN,k,j,i)*w(IVY,k,j,i)/tau_top_*scale*dt;
        du(IVZ,k,j,i) -= w(IDN,k,j,i)*w(IVZ,k,j,i)/tau_top_*scale*dt;
      }
    }

#if DEBUG_LEVEL > 2
  pmb->pdebug->CheckConservation("du", du, pmb->is, pmb->ie, pmb->js, pmb->je, pmb->ks, pmb->ke);
#endif
  pmb->pdebug->Leave();
  return TaskStatus::success;
}

TaskStatus Physics::BotSpongeLayer(AthenaArray<Real> &du,
  AthenaArray<Real> const& w, Real time, Real dt)
{
  MeshBlock *pmb = pmy_block;
  pmb->pdebug->Call("Physics::BotSpongeLayer");

  Coordinates *pcoord = pmb->pcoord;

  int il = pmb->is;
  Real x1min = pmb->pmy_mesh->mesh_size.x1min;
  Real x1max = pmb->pmy_mesh->mesh_size.x1max;

  while (il <= pmb->ie + 1 && (pcoord->x1v(il) - x1min < width_bot_)) il++;

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i < il; ++i) {
        Real eta = (pcoord->x1v(i) - x1min)/(x1max - x1min);
        Real eta_sp = width_bot_/(x1max - x1min);
        if (eta > eta_sp) break;
        Real scale = sqr(sin(M_PI/2.*(eta_sp - eta)/eta_sp));
        du(IVX,k,j,i) -= w(IDN,k,j,i)*w(IVX,k,j,i)/tau_bot_*scale*dt;
        du(IVY,k,j,i) -= w(IDN,k,j,i)*w(IVY,k,j,i)/tau_bot_*scale*dt;
        du(IVZ,k,j,i) -= w(IDN,k,j,i)*w(IVZ,k,j,i)/tau_bot_*scale*dt;
      }
    }

#if (DEBUG_LEVEL > 1)
  pmb->pdebug->CheckConservation("du", du, pmb->is, pmb->ie, pmb->js, pmb->je, pmb->ks, pmb->ke);
#endif
  pmb->pdebug->Leave();
  return TaskStatus::success;
}

/*void __attribute__((weak)) HydroSourceTerms::SpongeLayer2(
  AthenaArray<Real> &u, AthenaArray<Real> const& w,
  Real dt, int il, int iu, int &jl, int &ju)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  Coordinates *pcoord = pmb->pcoord;

  int js, je, ks, ke;

  jl = js = pmb->js; ks = pmb->ks;
  ju = je = pmb->je; ke = pmb->ke;

  Real x2min = pmb->pmy_mesh->mesh_size.x2min;
  Real x2max = pmb->pmy_mesh->mesh_size.x2max;

  while (jl <= je + 1 && pcoord->x2v(jl) - x2min < width_i2_) jl++;
  while (ju >= js && x2max - pcoord->x2v(ju) < width_o2_) ju--;

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j < jl; ++j) {
      Real dx = pcoord->x2v(j) - pcoord->x2v(js);
      Real scale = 1. - pow(dx/width_i2_, 0.5);
      for (int i = il; i <= iu; ++i)
        for (int n = 0; n < 3; ++n)
          u(IVX+n,k,j,i) /= (1. + dt/tau_i2_*scale);
    }

    for (int j = ju+1; j <= je; ++j) {
      Real dx = pcoord->x2v(je) - pcoord->x2v(j);
      Real scale = 1. - pow(dx/width_o2_, 0.5);
      for (int i = il; i <= iu; ++i)
        for (int n = 0; n < 3; ++n)
          u(IVX+n,k,j,i) /= (1. + dt/tau_o2_*scale);
    }
  }
}

void __attribute__((weak)) HydroSourceTerms::SpongeLayer3(
  AthenaArray<Real> &u, AthenaArray<Real> const& w,
  Real dt, int il, int iu, int jl, int ju, int &kl, int &ku)
{
  MeshBlock *pmb = pmy_hydro_->pmy_block;
  Coordinates *pcoord = pmb->pcoord;

  int ks, ke;

  kl = ks = pmb->ks;
  ku = ke = pmb->ke;

  Real x3min = pmb->pmy_mesh->mesh_size.x3min;
  Real x3max = pmb->pmy_mesh->mesh_size.x3max;

  while (kl <= ke + 1 && pcoord->x3v(kl) - x3min < width_i3_) kl++;
  while (ku >= ks && x3max - pcoord->x3v(ku) < width_o3_) ku--;

  for (int k = ks; k < kl; ++k) {
    Real dx = pcoord->x3v(k) - pcoord->x3v(ks);
    Real scale = 1. - pow(dx/width_i3_, 0.5);
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i)
        for (int n = 0; n < 3; ++n)
          u(IVX+n,k,j,i) /= (1. + dt/tau_i3_*scale);
  }

  for (int k = ku+1; k <= ke; ++k) {
    Real dx = pcoord->x3v(ke) - pcoord->x3v(k);
    Real scale = 1. - pow(dx/width_i3_, 0.5);
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i)
        for (int n = 0; n < 3; ++n)
          u(IVX+n,k,j,i) /= (1. + dt/tau_i3_*scale);
  }
}*/
