#ifndef ABSORBER_HPP_
#define ABSORBER_HPP_

// C++ header
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>

// Athena++ headers
#include "../athena.hpp"

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
    pmy_band(pband), myname("NULL"), prev(NULL), next(NULL), imol_(-1), mixr_(0) {}

  Absorber(RadiationBand *pband, std::string name, int imol, Real mixr = 1.): 
    pmy_band(pband), myname(name), prev(NULL), next(NULL), imol_(imol), mixr_(mixr) {}

  Absorber(RadiationBand *pband, std::string name, std::vector<int> imols, Real mixr = 1): 
    pmy_band(pband), myname(name), prev(NULL), next(NULL), imols_(imols), mixr_(mixr) {}

  virtual ~Absorber() {
    if (prev != NULL) prev->next = next;
    if (next != NULL) next->prev = prev;
  }

  template<typename Ab> Absorber* AddAbsorber(Ab const& a) {
    Ab* pa = new Ab(a);
    Absorber *p = this;
    while (p->next != NULL) p = p->next;
    p->next = pa;
    p->next->prev = p;
    p->next->next = NULL;
    return p->next;
  }

  virtual void SaveCoefficient(std::string fname) const {}
  virtual void LoadCoefficient(std::string fname) {}
  virtual Real AbsorptionCoefficient(Real wave, Real const prim[]) const { return 0.; }
  virtual Real SingleScatteringAlbedo(Real wave, Real const prim[]) const { return 0.; }
  virtual void PhaseMomentum(Real wave, Real const prim[], Real *pp, int np) const {}

protected:
  int  imol_;       /**< id of dependent molecule */
  Real mixr_;       /**< mixing ratio for dependent molecule */
  std::vector<int>  imols_;     /**< id of dependent molecules */
};


/*class RosselandMean: public Absorber {
public:
  RosselandMean() : Absorber("") {}
  virtual ~RosselandMean() {}
  Real AbsorptionCoefficient(Real wave, Real const prim[]) const;
};*/


#endif
