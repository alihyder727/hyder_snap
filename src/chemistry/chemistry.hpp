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
  ChemistryBase<T>* AddToChemistry(ChemistryBase<T> *pchem,
    ParameterInput *pin, std::string name) {
    T* pnew = new T(pmy_block, pin, name);
    ChemistryBase<T>* p = pchem;
    if (p == nullptr) {
      p = static_cast<ChemistryBase<T>*>(pnew);
      p->prev = nullptr;
      p->next = nullptr;
    } else {
      while (p->next != nullptr) p = p->next;
      p->next = static_cast<ChemistryBase<T>*>(pnew);
      p->next->prev = p;
      p->next->next = nullptr;
    }
    return p;
  }

  void TimeIntegrate(AthenaArray<Real> &u, Real time, Real dt);

protected:
  ChemistryBase<Kessler94> *pkessler94_;
};


#endif /* end of include guard CHEMISTRY_HPP */
