/** @file chemistry.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Jun 10, 2021 09:58:18 PDT
 * @bug No known bugs.
 */

#ifndef CHEMISTRY_HPP
#define CHEMISTRY_HPP

#include "../athena.hpp"
#include "chemistry_base.hpp"

class MeshBlock;
class ParameterInput;
class Kessler94;

class Chemistry {
public:
// data
  MeshBlock *pmy_block;

  Chemistry(MeshBlock *pmb, ParameterInput *pin);
  ~Chemistry();

  template<typename T>
  void AddToChemistry(ChemistryBase<T>* &pchem,
    ParameterInput *pin, std::string name) {
    T* pnew = new T(this, pin, name);
    if (pchem == nullptr) {
      pchem = static_cast<ChemistryBase<T>*>(pnew);
      pchem->prev = nullptr;
      pchem->next = nullptr;
    } else {
      ChemistryBase<T>* p = pchem;
      while (p->next != nullptr) p = p->next;
      p->next = static_cast<ChemistryBase<T>*>(pnew);
      p->next->prev = p;
      p->next->next = nullptr;
    }
  }

  void TimeIntegrate(Real time, Real dt) const;

protected:
  ChemistryBase<Kessler94> *pkessler94_;
};

#endif /* end of include guard CHEMISTRY_HPP */
