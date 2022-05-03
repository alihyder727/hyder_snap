#ifndef ABSORBER_HPP_
#define ABSORBER_HPP_

// C++ header
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../debugger/debugger.hpp"
#include "../mesh/mesh.hpp"
#include "../mesh/grid_data.hpp"
#include "radiation.hpp"

/**@file
 * @brief This file contains declaration of Absorber
*
* **Author** : Cheng Li, California Institute of Technology <br>
* **Contact** : cli@gps.caltech.edu <br>
* **Revision history** :
* - June 21 2016, start documenting this file
* - July 28 2016, merge scatter into absorber
* - June 24 2017, adapt to Athena++ framework
* - April 03 2019, merge to snap
* - July 27 2019, add multiple dependent molecules
*/

class RadiationBand;

class Absorber {
public:
  // data
  RadiationBand *pmy_band;
  std::string myname;
  Absorber *prev, *next;
  
  // functions
  Absorber(RadiationBand *pband):
    pmy_band(pband), myname("NULL"), prev(nullptr), next(nullptr), imol_(-1), mixr_(0) {}

  Absorber(RadiationBand *pband, std::string name):
    pmy_band(pband), myname(name), prev(nullptr), next(nullptr), imol_(-1), mixr_(0) {}

  Absorber(RadiationBand *pband, std::string name, int imol, Real mixr = 1.): 
    pmy_band(pband), myname(name), prev(nullptr), next(nullptr), imol_(imol), mixr_(mixr)
  {
    Debugger *pdebug = pband->pmy_rad->pmy_block->pdebug;
    pdebug->Enter("Absorber " + name);
    std::stringstream &msg = pdebug->msg;
    msg << "- molar mixing ratio = " << mixr
        << " and id = " << imol << std::endl;
    pdebug->Leave();
  }

  Absorber(RadiationBand *pband, std::string name, std::vector<int> imols, Real mixr = 1): 
    pmy_band(pband), myname(name), prev(nullptr), next(nullptr), imols_(imols), mixr_(mixr)
  {
    Debugger *pdebug = pband->pmy_rad->pmy_block->pdebug;
    pdebug->Enter("Absorber " + name);
    std::stringstream &msg = pdebug->msg;
    msg << "- molar mixing ratio = " << mixr
        << " and dependent id = ";
    for (int i = 0; i < imols.size(); ++i)
      msg << imols_[i] << " ";
    msg << std::endl;
    pdebug->Leave();
  }

  virtual ~Absorber() {
    if (prev != nullptr) prev->next = next;
    if (next != nullptr) next->prev = prev;
  }

  template<typename Ab> Absorber* addAbsorber(Ab const& a) {
    Ab* pa = new Ab(a);
    Absorber *p = this;
    while (p->next != nullptr) p = p->next;
    p->next = pa;
    p->next->prev = p;
    p->next->next = nullptr;
    return p->next;
  }

  //virtual void SaveCoefficient(std::string fname) const {}
  virtual void loadCoefficient(std::string fname, int bid = -1) {}
  virtual Real getAttenuation(Real wave1, Real wave2, GridData const& gdata) const { return 0.; }
  virtual Real getSingleScatteringAlbedo(Real wave1, Real wave2, GridData const& gdata) const { return 0.; }
  virtual void getPhaseMomentum(Real *pp, Real wave1, Real wave2, GridData const& gdata, int np) const {}

  //! \deprecated use setAttenuation instead
  //virtual Real AbsorptionCoefficient(Real wave, Real const prim[]) { return 0.;}
  //virtual Real SingleScatteringAlbedo(Real wave, Real const prim[]) const { return 0.; }
  //virtual void PhaseMomentum(Real wave, Real const prim[], Real *pp, int np) const {}

protected:
  int  imol_;       /**< id of dependent molecule */
  Real mixr_;       /**< mixing ratio for dependent molecule */
  std::vector<int>  imols_;     /**< id of dependent molecules */
};


/*class RosselandMean: public Absorber {
public:
  RosselandMean() : Absorber("") {}
  virtual ~RosselandMean() {}
  Real Attenuation(Real wave, Real const prim[]) const;
};*/


#endif
