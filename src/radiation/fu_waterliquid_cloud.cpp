// C/C++ headers
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cstring>
#include <cassert>  // assert

// Athena++ headers
#include "../math/interpolation.h"
#include "absorber.hpp"
#include "water_cloud.hpp"
#include "radiation_utils.hpp"  // getPhaseMomentum


// coefficient from fu and liou

// wavenumber range

Real liquid_wband[19] = { 0, 280, 400, 540, 670, 800, 980, 1100, 1250, 1400, 1700, 1900, 2200, 2857, 4000, 5263, 7692, 14286, 50000}; // cm-1

// radiuse of water cloud

Real re[8] = { 4.18,  5.36,  5.89,  6.16, 9.27,  9.84, 12.10, 31.23};

// fl for tau
Real fl[8] = { 0.05, 0.14, 0.22, 0.28, 0.50, 0.47, 1.00, 2.50};

// bz for tau
Real bz[18][8] = {{15.11,  40.25,  59.81,  72.43, 83.69,  73.99, 128.17, 120.91},
                   {15.74,  41.70,  61.52,  74.47, 85.78,  75.59, 130.46, 121.84},
                   {16.38,  43.52,  64.84,  77.97, 87.31,  77.36, 134.30, 124.06}, 
                   {17.57,  45.78,  66.44,  80.15, 90.49,  79.90, 137.56, 125.92},
                   {18.19,  46.63,  69.39,  82.20, 91.46,  79.99, 138.21, 126.08},
                   {21.30,  51.88,  77.77,  87.02, 94.91,  83.55, 143.46, 128.45},
                   {22.44,  57.35,  84.41, 103.50, 103.49,  84.17, 152.77, 132.07},
                   {18.32,  52.69,  76.67, 100.31, 105.46,  92.86, 157.82, 133.03},
                   {17.27,  50.44,  74.18,  96.76, 105.32,  95.25, 158.07, 134.48},
                   {13.73,  44.90,  67.70,  90.85, 109.16, 105.48, 163.11, 136.21},
                   {10.30,  36.28,  57.23,  76.43, 106.45, 104.90, 161.73, 136.62},
                   {7.16,  26.40,  43.51,  57.24, 92.55,  90.55, 149.10, 135.13},
                   {6.39,  21.00,  33.81,  43.36, 66.90,  63.58, 113.83, 125.65},
                   {10.33,  30.87,  47.63,  60.33, 79.54,  73.92, 127.46, 128.21},
                   {11.86,  35.64,  54.81,  69.85, 90.39,  84.16, 142.49, 135.25},
                   {10.27,  33.08,  51.81,  67.26, 93.24,  88.60, 148.71, 140.42},
                   {6.72,  24.09,  39.42,  51.68, 83.34,  80.72, 140.14, 143.57},
                   {3.92,  14.76,  25.32,  32.63, 60.85,  58.81, 112.30, 145.62}};

// wz for ssalb

Real wz[18][8] = {{.999999, .999999, .999999, .999999, .999998, .999999, .999998, .999997},
                  {.999753, .999700, .999667, .999646, .999492, .999470, .999344, .998667},
                  {.995914, .994967, .994379, .993842, .991385, .990753, .988908, .974831},
                  {.983761, .978981, .976568, .974700, .963466, .959934, .953865, .897690},
                  {.702949, .683241, .679723, .669045, .642616, .632996, .629776, .588820},
                  {.947343, .929619, .924806, .914557, .877169, .867047, .853661, .737426},
                  {.919356, .896274, .885924, .881097, .812772, .781637, .775418, .637341},
                  {.874717, .861122, .847850, .851677, .787171, .772952, .753143, .618656},
                  {.764750, .752410, .736529, .743435, .671272, .659392, .639492, .549941},
                  {.807536, .808700, .795994, .805489, .750577, .755524, .709472, .571989},
                  {.753346, .772026, .767273, .777079, .751264, .760973, .712536, .568286},
                  {.632722, .676332, .684631, .693552, .707986, .717724, .682430, .552867},
                  {.288885, .348489, .371653, .380367, .454540, .465769, .475409, .493881},
                  {.261827, .306283, .321340, .333051, .392917, .406876, .417450, .484593},
                  {.295804, .339929, .352494, .365502, .416229, .430369, .435267, .491356},
                  {.301214, .354746, .369346, .381906, .433602, .447397, .447406, .486968},
                  {.243714, .318761, .344642, .352770, .427906, .438979, .445972, .477264},
                  {.109012, .187230, .226849, .224976, .331382, .335917, .374882, .457067}};



// gz for gg


Real gz[18][8] = {{.838, .839, .844, .847, .849, .860, .853, .859},
                  {.809, .810, .819, .823, .823, .849, .833, .843},
                  {.774, .787, .781, .792, .812, .836, .815, .833},
                  {.801, .802, .793, .793, .814, .829, .818, .832},
                  {.877, .873, .879, .880, .885, .899, .891, .908},
                  {.783, .769, .777, .756, .764, .776, .770, .797},
                  {.818, .805, .824, .830, .815, .801, .820, .845},
                  {.810, .802, .826, .840, .829, .853, .840, .868},
                  {.774, .766, .799, .818, .815, .869, .834, .869},
                  {.734, .728, .767, .797, .796, .871, .818, .854},
                  {.693, .688, .736, .772, .780, .880, .808, .846},
                  {.643, .646, .698, .741, .759, .882, .793, .839},
                  {.564, .582, .637, .690, .719, .871, .764, .819},
                  {.466, .494, .546, .609, .651, .823, .701, .766},
                  {.375, .410, .455, .525, .583, .773, .637, .710},
                  {.262, .301, .334, .406, .485, .695, .545, .631},
                  {.144, .181, .200, .256, .352, .562, .413, .517},
                  {.060, .077, .088, .112, .181, .310, .222, .327}};



Real FuWaterLiquidCloud::getAttenuation(Real wave1, Real wave2,
    Real const q[], Real const c[], Real const s[]) const
{  
  static const Real Rgas = 8.314462;
  const Real mu = 29.E-3;  // FIXME: hard-wired to Earth

  Real pre = 10.; // choose in 4.18 ~ 31.23 micron 
  //Real pre = 10.0; // choose in 4.18 ~ 31.23 micron 
                   // effective radius of water cloud ( um )
  
  int j = locate(re, pre, 8);
  assert(j >= 0 && j <= 7);

  int iband = locate(liquid_wband, wave1, 19)+1;
  assert(iband >= 1 && iband <= 18);
  iband = 18 - iband;

  Real result=0.;
  Real k = (1.0/re[j+1]-1.0/re[j])/(1.0/pre-1.0/re[j]);
  Real dens   = q[IPR]*mu/(Rgas*q[IDN]);

  result = bz[iband-1][j]/fl[j] + (bz[iband-1][j+1]/fl[j+1]-bz[iband-1][j]/fl[j]) / k;
  
  return dens*q[imol_]*result;     // -> 1/m
}

Real FuWaterLiquidCloud::getSingleScateringAlbedo(Real wave1, Real wave2,
    Real const q[], Real const c[], Real const s[]) const
{
// from fu & liou code
  Real pre = 10.;
    
  int j = locate(re, pre, 8);
  assert(j >= 0 && j <= 7);
    
  int iband = locate(liquid_wband, wave1, 19)+1;
  assert(iband >= 1 && iband <= 18);
  iband = 18 - iband;
    
  Real ww=0.;
  ww = wz[iband-1][j]+(wz[iband-1][j+1]-wz[iband-1][j])/(re[j+1]-re[j])*(pre-re[j]);

  if (q[imol_]<2e-19) {
    return 0.0;
  } else {
    return ww;
  }
}

void FuWaterLiquidCloud::getPhaseMomentum(Real *pp, Real wave1, Real wave2,
    Real const q[], Real const c[], Real const s[], int np) const
{
// from fu & liou code
    
  Real pre = 10.;
    
  int j = locate(re, pre, 8);
  assert(j >= 0 && j <= 7);
    
  int iband = locate(liquid_wband, wave1, 19)+1;
  assert(iband >= 1 && iband <= 18);
  iband = 18 - iband;
    
  Real gg=0.;
    
  gg = gz[iband-1][j]+(gz[iband-1][j+1]-gz[iband-1][j])/(re[j+1]-re[j])*(pre-re[j]);
  
  if (q[imol_]<2e-19){
    getPhaseHenyeyGreenstein(pp, 0, 0.0, np); // 0 for HENYEY_GREENSTEIN
  } else {
    getPhaseHenyeyGreenstein(pp, 0, gg, np);  // 0 for HENYEY_GREENSTEIN
  }
}
