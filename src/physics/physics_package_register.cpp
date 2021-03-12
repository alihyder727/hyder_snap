// C/C++ headers
#include <sstream>

// Athena++ headers
#include "../athena.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../task_list/task_manager.hpp"
#include "physics.hpp"

using namespace PhysicsPackageNames;

Physics::Physics(MeshBlock *pmb, ParameterInput *pin):
  pmy_block(pmb)
{
  std::stringstream msg;
  char package_names[1024], *p;
  std::string str = pin->GetOrAddString("physics", "packages", "");
  std::strcpy(package_names, str.c_str());
  p = std::strtok(package_names, " ,");

  ptm = new TaskManager<Physics>(this);

  PhysicsPackage pkg;
  while (p != NULL) {
    if (std::strcmp(p, "fix_bot_temperature") == 0) {
      pkg.id = FIX_BOT_TEMPERATURE;
      pkg.dependency = 0LL;
      pkg.incompatible = 0LL;
      pkg.Function = &Physics::RelaxBotTemperature;
      ptm->AddPackage(pkg, "fix_bot_temperature");

      tau_Tbot_ = pin->GetReal("physics", "fix_bot_temperature.tau");
      if (!hydro_bot_.IsAllocated())
        hydro_bot_.NewAthenaArray(NHYDRO, pmb->ncells3, pmb->ncells2);
      tem_bot_.InitWithShallowSlice(hydro_bot_, 3, IDN, 1);
    } else if (std::strcmp(p, "fix_bot_velocity") == 0) {
      pkg.id = FIX_BOT_VELOCITY;
      pkg.dependency = 0LL;
      pkg.incompatible = 0LL;
      // pkg.Function = &Physics::RelaxBotVelocity;
      ptm->AddPackage(pkg, "fix_bot_velocity");

      tau_Ubot_ = pin->GetReal("physics", "fix_bot_velocity.tau");
      if (!hydro_bot_.IsAllocated())
        hydro_bot_.NewAthenaArray(NHYDRO, pmb->ncells3, pmb->ncells2);
      vel_bot_.InitWithShallowSlice(hydro_bot_, 3, IVX, 3);
    } else if (std::strcmp(p, "fix_bot_composition") == 0) {
      pkg.id = FIX_BOT_COMPOSITION;
      pkg.dependency = 0LL;
      pkg.incompatible = 0LL;
      // pkg.Function = &Physics::RelaxBotComposition;
      ptm->AddPackage(pkg, "fix_bot_composition");

      tau_Ubot_ = pin->GetReal("physics", "fix_bot_composition.tau");
      if (!hydro_bot_.IsAllocated())
        hydro_bot_.NewAthenaArray(NHYDRO, pmb->ncells3, pmb->ncells2);
      com_bot_.InitWithShallowSlice(hydro_bot_, 3, IDN, NMASS);
    } else if (std::strcmp(p, "top_sponge_layer") == 0) {
      pkg.id = TOP_SPONGE_LAYER;
      pkg.dependency = 0LL;
      pkg.incompatible = 0LL;
      // pkg.Function = &Physics::TopSpongeLayer;
      ptm->AddPackage(pkg, "top_sponge_layer");

      tau_top_ = pin->GetReal("physics", "top_sponge_layer.tau");
      width_top_ = pin->GetReal("physics", "top_sponge_layer.width");
    } else if (std::strcmp(p, "bot_sponge_layer") == 0) {
      pkg.id = BOT_SPONGE_LAYER;
      pkg.dependency = 0LL;
      pkg.incompatible = 0LL;
      // pkg.Function = &Physics::BotSpongeLayer;
      ptm->AddPackage(pkg, "bot_sponge_layer");

      tau_bot_ = pin->GetReal("physics", "bot_sponge_layer.tau");
      width_bot_ = pin->GetReal("physics", "bot_sponge_layer.width");
    } else {
      msg << "### FATAL ERROR in function Physics::Physics"
          << std::endl << "Package '" << p << "' "
          << "unrecognized." << std::endl;
      ATHENA_ERROR(msg);
    }
    packages_.push_back(pkg);
    p = std::strtok(NULL, " ,");
  }
}

void Physics::Initialize(AthenaArray<Real> const& w)
{
  MeshBlock *pmb = pmy_block;

  // find top and bot neighbor
  NeighborBlock ntop, nbot;
  bool has_top_neighbor = false;
  bool has_bot_neighbor = false;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock& nb = pmb->pbval->neighbor[n];
    if ((nb.ni.ox1 == -1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
      nbot = nb;
      has_bot_neighbor = true;
    } if ((nb.ni.ox1 == 1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
      ntop = nb;
      has_top_neighbor = true;
    }
  }

  if (has_bot_neighbor)
    ptm->RemoveTask(FIX_BOT_TEMPERATURE | FIX_BOT_VELOCITY | FIX_BOT_COMPOSITION);

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j) {
      if (ptm->HasTask(FIX_BOT_TEMPERATURE))
        tem_bot_(k,j) = pmb->pthermo->Temp(w.at(k,j,pmb->is));
      if (ptm->HasTask(FIX_BOT_VELOCITY)) {
        vel_bot_(0,k,j) = w(IVX,k,j,pmb->is);
        vel_bot_(1,k,j) = w(IVY,k,j,pmb->is);
        vel_bot_(2,k,j) = w(IVZ,k,j,pmb->is);
      }
      if (ptm->HasTask(FIX_BOT_COMPOSITION))
        for (int n = 1; n < NMASS; ++n)
          com_bot_(n,k,j) = w(n,k,j,pmb->is);
    }
}
