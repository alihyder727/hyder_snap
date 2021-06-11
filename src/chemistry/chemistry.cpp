/** @file chemistry.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Jun 10, 2021 11:02:34 PDT
 * @bug No known bugs.
 */
// C/C++ headers
#include <string>
#include <cstring>

// Athena++ headers
#include "chemistry.hpp"
#include "kessler94.hpp"

Chemistry::Chemistry(MeshBlock *pmb, ParameterInput *pin) :
  pmy_block(pmb), pkessler94_(nullptr)
{
  char chem_names[1024], *p;
  std::string str = pin->GetOrAddString("chemistry", "chemistry", "");
  std::strcpy(chem_names, str.c_str());
  p = std::strtok(chem_names, " ,");

  while (p != NULL) {
    std::stringstream msg;
    std::string name;
    char *c = std::strchr(p, '.');
    if (c != NULL) name = c+1;
    else name = p;
    if (std::strncmp(p, "kessler94", 9) == 0) {
      AddToChemistry(pkessler94_, pin, name);
    } else {
      msg << "### FATAL ERROR in function Particles::Particles"
          << std::endl << "Particles '" << p << "' "
          << "unrecognized." << std::endl;
      ATHENA_ERROR(msg);
    }
    p = std::strtok(NULL, " ,");
  }
}

Chemistry::~Chemistry()
{
  if (pkessler94_ != nullptr) {
    while (pkessler94_->prev != nullptr)
      delete pkessler94_->prev;
    while (pkessler94_->next != nullptr)
      delete pkessler94_->next;
    delete pkessler94_;
  }
}

void Chemistry::TimeIntegrate(AthenaArray<Real> &u, Real time, Real dt)
{
  Particles *ppart = pkessler94_->pmy_part;
  pkessler94_->IntegrateDense(u, ppart->dc, ppart->c, time, dt);
}
